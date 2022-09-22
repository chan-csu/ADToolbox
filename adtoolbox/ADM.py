
from typing import Callable
import plotly
from sympy import plot
import numpy as np
import scipy.optimize
import scipy.integrate
import pandas as pd
import json
import os
import plotly.express as px
from dash import Dash, html, dcc
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output,dash_table
from ADToolBox import Database as Database
from ADToolBox import Reaction_Toolkit as Reaction_Toolkit
from dash.dash_table.Format import Format, Scheme, Sign, Symbol
import pandas as pd
from ADToolBox import Reaction
from collections import OrderedDict
import rich
from rich.console import Console
from rich.table import Table
from __init__ import Main_Dir,PKG_DATA
from rich.style import Style
import dash_escher
### Note ###
# The following code is a modified version of the code from the PyADM1 package
# It is extensively based on PyADM1, and I would like to thank the author of PyADM1 for this work
# As far as implementing the original ADM1 in python goes, I still consider this as a modification from
# the PyADM1 code.

# ----------


RT = Reaction_Toolkit(Reaction_DB=os.path.join(
    Main_Dir, "..", "Database", "reactions.json"))

class _Fake_Sol:
    def __init__(self, y,t):
        self.y = y
        self.t=t



class Model:

    """General Class For a Model
    Args:
        Model Parameters (dict): a dictionary which contains model parameters
        Base_Parameters (dict): a dictionary which contains base paramters
        Initial Conditions (dict): a dictionary containing inlet conditions for all species
        Inlet Conditions (dict): a dictionary containing inlet conditions for all species
        Reactions (list): a list containing all types of reactions
        Species (list): a list containing all species
        ODE_System (function): a function which solves inital value ODEs
        Build Stoiciometric Matrix (function): a function which outputs a matrix for of all reactions
        
    Returns:
        adtoolbox.ADM.Model: returns a model instance for downstream purposes.
    """
    


    def __init__(self, Model_Parameters: dict, Base_Parameters: dict, Initial_Conditions: dict, Inlet_Conditions:dict, Reactions: list, Species: list, ODE_System:Callable, Build_Stoiciometric_Matrix:Callable,Control_States:dict={},Metagenome_Report=None, Name="ADM1", Switch="DAE",simulation_time=30):
        self.Model_Parameters = Model_Parameters
        self.Base_Parameters = Base_Parameters
        
        # try:
        #     for items in Control_States.keys():
        #         Initial_Conditions[items]
        #         Initial_Conditions[items]=Control_States[items]
        
        # except TypeError:

        #     rich.print("[red]No acceptable Controled State passed!\nEverything will change according to the system's dynamics!") 
        #     self.Control_States={}
        # except KeyError:
            
        #     rich.print("[red]Control states must follow the ADM state variable notations!\nCheck the Control_States json file.")
        #     self.Control_States={}
        # else:
        #     rich.print("[green]The following states were set as control:\n")
        #     table = Table(title="Control States\n",header_style=Style(color="cyan", blink=True, bold=True))
        #     for i in Control_States.keys():
        #         table.add_column(i)
            
        #     table.add_row(*[str(Control_States[i]) for i in Control_States.keys()],style=Style(color="white", blink=True, bold=True))
        #     rich.print(table)
        for items in Control_States.keys():
            Initial_Conditions[items]
            Initial_Conditions[items]=Control_States[items]
        self.Control_States=Control_States


        self.Inlet_Conditions = np.array(
            [Inlet_Conditions[i+"_in"] for i in Species])[:, np.newaxis]
        self.Reactions = Reactions
        self.Species = Species
        self.Initial_Conditions = np.array(
            [Initial_Conditions[i] for i in Species])[:, np.newaxis]
        self._IC=Initial_Conditions
        self._InC=Inlet_Conditions
        self.Switch = Switch
        self.Name = Name
        self.Metagenome_Report = Metagenome_Report
        self.S = Build_Stoiciometric_Matrix(
            Base_Parameters, Model_Parameters, Reactions, Species)
        self.ODE_System = ODE_System
        self.Simtime=simulation_time

    # def Build_COBRA_Model(self, save_model=True):
    #     model = cobra.Model(self.Name)

    #     for i in range(len(self.Reactions)):
    #         Temp_Rxn = cobra.Reaction(
    #             self.Reactions[i].replace(" ", "_"), name=self.Reactions[i].replace(" ", "_"))
    #         Temp_Mets = np.where(self.S[:, i] != 0)
    #         Met_Dict = {}
    #         for met in Temp_Mets[0]:
    #             Metabolite = cobra.Metabolite(self.Species[met].replace(" ", "_"),
    #                                           name=self.Species[met].replace(" ", "_"), compartment="Model")
    #             Met_Dict[Metabolite] = self.S[met, i]
    #         Temp_Rxn.add_metabolites(Met_Dict)
    #         # print(Temp_Rxn.reaction)
    #         model.add_reaction(Temp_Rxn)

        # if save_model:
        #     cobra.io.save_json_model(model, self.Name+'.json')

    def Solve_Model(self, y0: np.ndarray, T_eval: np.ndarray, method="LSODA", switch="DAF")->scipy.integrate._ivp.ivp.OdeResult:
        """
        Function to solve the model. 
        Examples:
            >>> from ADM import Model
            >>> import numpy as np
            >>> reactions=['rxn1','rxn2']
            >>> species=['a','b','c']
            >>> Initial_Conditions={'a':.001,'b':.002,'c':.003}
            >>> Inlet_Conditions={'a_in':.001,'b_in':.002,'c_in':.003}
            >>> Model_Parameters={'k1':0.001,'k2':0.002}
            >>> Base_Parameters={'T':0.1}
            >>> def Build_Stoiciometric_Matrix(Base_Parameters,Model_Parameters,Reactions,Species):
            ...    S = np.zeros((len(Species), len(Reactions)))
            ...    S[[0,1],0]=[-1,0.001]
            ...    S[[1,2],1]=[-5,1]
            ...    return S
            >>> def ODE_System(t,c,Model1):
            ...    v = np.zeros((len(Model1.Reactions), 1))
            ...    v[0]=Model1.Model_Parameters['k1']*c[0]*Model1.Base_Parameters['T']/1000
            ...    v[1]=Model1.Model_Parameters['k2']*c[1]/1000
            ...    dCdt=np.matmul(Model1.S,v)
            ...    return dCdt[:, 0]
            >>> m= Model(Model_Parameters,Base_Parameters,Initial_Conditions,Inlet_Conditions,reactions,species,ODE_System,Build_Stoiciometric_Matrix)
            >>> m.Solve_Model((0,.1),m.Initial_Conditions[:, 0],np.linspace(0,0.1,10),method='RK45')['status']==0
            True
        
        Args:
            t_span (tuple): The range the time will include
            y0 (np.ndarray): The initial value at time = 0
            T_eval (np.ndarray): The time steps that will be evaluated
        
        Returns:
            scipy.integrate._ivp.ivp.OdeResult: Returns the results of the simulation being run and gives optimized paramters.
        """
        self.Info={"Fluxes":[]}
        try:
            C = scipy.integrate.solve_ivp(
                self.ODE_System, (0,self.Simtime), y0, t_eval=T_eval, method=method, args=[self])
        
        except Exception as e:
            print("Could not solve model, setting C to a very large value")
            C=_Fake_Sol(np.ones((y0.shape[0],T_eval.__len__()))*1e10,T_eval)
       
        return C

    
       #C = scipy.integrate.solve_ivp(
       #        self.ODE_System, t_span, y0, t_eval=T_eval, method=method, args=[self])
       #
       #return C

    def Plot(self, Sol: scipy.integrate._ivp.ivp.OdeResult, Type: str = "Line")-> None:
        """ A function which returns a plot of the solution from the ODE
        """
        Solution = {
            't': Sol.t,
        }
        for i in range(len(self.Species)):
            Solution[self.Species[i]] = Sol.y[i, :]
        Sol_DF = pd.DataFrame(Solution)

        if Type == "Line":
            fig = px.line(Sol_DF, x="t", y=Sol_DF.columns,
                          title="Concentration of species",text=Sol_DF.columns,textposition="bottom center")
            fig.update_layout(
                title={
                    'y': 0.95,
                    'x': 0.5,

                    "font_size": 30,
                    'xanchor': 'center',
                    'yanchor': 'top'}

            )
            fig.update_xaxes(
                title={
                    "text": "Time (Days)",
                    "font_size": 25,
                }
            )
            fig.update_yaxes(
                title={
                    "text": "Concentrations (kg COD/m^3)",
                    "font_size": 25,
                }
            )
            fig.show()

        elif Type == "Sankey":
            pass

    def Dash_App(self, Sol: scipy.integrate._ivp.ivp.OdeResult, Type: str = "Line")->None:
        """A fuction which opens a web browser that displays the results"""
        with open(os.path.join(PKG_DATA,"Modified_ADM_Map.json"),'rb') as f:
            MD=json.load(f)
        with open(os.path.join(PKG_DATA,"Modified_ADM_Model.json"),'rb') as f:
            MOD=json.load(f)


        app = Dash(__name__)
        colors = {
            'background': '#659dbd',
            'text': '#3e4444'
        }
        

        Solution = {
            't': Sol.t,
        }
        for i in range(len(self.Species)):
            Solution[self.Species[i]] = Sol.y[i, :]
        Sol_DF = pd.DataFrame(Solution)

        if Type == "Line":
            fig = px.line(Sol_DF, x="t", y=Sol_DF.columns,
                          title="Concentration of species")
            fig.update_layout(
            title={
            'y': 0.95,
            'x': 0.5,
            "font_size": 30,
            'xanchor': 'center',
            'yanchor': 'top'},
            legend=dict(font=dict(size= 20))

                )
            fig.update_xaxes(
            title={
            "text": "Time (Days)",
            "font_size": 25,
                },
                 tickfont_size=20
                )
            fig.update_yaxes(
            title={
            "text": "Concentrations (kg COD/m^3)",
            "font_size": 25,
             },
             tickfont_size=20
            
                )
            fig.update_traces(line=dict(width=5))

        with_Report=html.Div([html.Div(style={'backgroundColor': 'rgb(200, 200, 200)', 'height':'300px', "margin-top":'-50px'}, children=[html.H1(
            children='ADToolbox',
            style={
                'textAlign': 'center',
                'color': colors['text'],
                'font-size': '60px',
                # 'font-family': 'sans-serif',
                'padding-top': '30px',
                'color': "rgb(30, 30, 30)",

            }),html.H2(children="A toolbox for modeleing and optimization of anaerobic digestion process",
            style={
                'textAlign': 'center',
                'color': "rgb(30, 30, 30)",
                'font-family': 'sans-serif'

            })
            ]),

            html.H2(f"{self.Name} Concentration Plot", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }),
            

            dcc.Graph(figure=fig, id='Concentrations_Line_Plot',
            style={
                "backgroundColor": "#171010",
                "height":"600px"
            }
            ),



            html.Br(),

            html.H3("Base_Parameters", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }),
            
            dash_table.DataTable(
            id='Base_Parameters',
            columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self.Base_Parameters.keys())],
            data=pd.DataFrame(self.Base_Parameters,index=[0]).to_dict('records'),
            editable=True,
            style_table={'overflowX': 'scroll'},
            style_header={
            'backgroundColor': 'rgb(200, 200, 200)',
            'color': 'black',
            'font-size': '20px',
                },
            style_data={
            'backgroundColor': 'rgb(250, 250, 250)',
            'color': 'black',
            'font-size': '20px',
            'font-family': 'Trebuchet MS',
            }),


            html.Br(),
            html.H3("Model_Parameters", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }),



            dash_table.DataTable(
            id='Model_Parameters',
            columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self.Model_Parameters.keys())],
            data=pd.DataFrame(self.Model_Parameters,index=[0]).to_dict('records'),
            editable=True,
            style_table={'overflowX': 'scroll'},
            style_header={
            'backgroundColor': 'rgb(200, 200, 200)',
            'color': 'black',
            'font-size': '20px',
                },
            style_data={
            'backgroundColor': 'rgb(250, 250, 250)',
            'color': 'black',
            'font-size': '20px',
            'font-family': 'Trebuchet MS',
            }
            ,

            
                ),            

            html.Br(),

            html.H3("Initial_Conditions", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }),
            
            dash_table.DataTable(
            id='Initial_Conditions',
            columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self._IC.keys())],
            data=pd.DataFrame(self._IC,index=[0]).to_dict('records'),
            editable=True,
            style_table={'overflowX': 'scroll'},
            style_header={
            'backgroundColor': 'rgb(200, 200, 200)',
            'color': 'black',
            'font-size': '20px',
                },
            style_data={
            'backgroundColor': 'rgb(250, 250, 250)',
            'color': 'black',
            'font-size': '20px',
            'font-family': 'Trebuchet MS',
            }
            
            ),

            html.Br(),
            html.H3("Inlet_Conditions", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }),
            dash_table.DataTable(
            id='Inlet_Conditions',
            columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self._InC.keys())],
            data=pd.DataFrame(self._InC,index=[0]).to_dict('records'),
            editable=True,
            style_table={'overflowX': 'scroll'},
            style_header={
            'backgroundColor': 'rgb(200, 200, 200)',
            'color': 'black',
            'font-size': '20px',
                },
            style_data={
            'backgroundColor': 'rgb(250, 250, 250)',
            'color': 'black',
            'font-size': '20px',
            'font-family': 'Trebuchet MS',
            },
            export_format='csv',
            export_headers='display',
            merge_duplicate_headers=True),

            html.Br(),

            html.H2("Microbial Association", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }) ,
            

            dcc.Dropdown(self.Reactions,
                         self.Reactions[0], style={"width": "300px","font-size":25}, id="Drop_Down") ,

            dcc.Graph(figure=None, id="Annotation_Graph", style={
            "height": "650px"}),

            html.H2("Escher Map", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }) ,
            
            dcc.Dropdown(["Show Map","Hide Map"],
                         self.Reactions[0], style={"width": "300px","font-size":25}, id="Drop_Down_Escher"),
            html.Br(),
            html.Div(children=None,id="Escher_"),
            html.Div(children=None,id="Escher"),
            
        #     html.Div(children=[dash_escher.DashEscher(mapData=MD,modelData=MOD,
        #      options={
        #          'reaction_data':rxn_data
        #      }
        #      ,height='1000px',
        #  width='100%')])
         ])
##


        without_Report=html.Div([html.Div(style={'backgroundColor': colors['background'], 'height':'200px', "margin-top":'-30px'}, children=[html.H1(
            children='ADToolbox',
            style={
                'textAlign': 'center',
                'color': colors['text'],
                'font-size': '50px',
                'padding-top': '20px'
            })]),

            html.H2("Original ADM1 Line Plot", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }),
            

            dcc.Graph(figure=fig, id='Concentrations_Line_Plot'),

            html.Br(),

            html.H3("Base_Parameters", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }),
            
            dash_table.DataTable(
            id='Base_Parameters',
            columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self.Base_Parameters.keys())],
            data=pd.DataFrame(self.Base_Parameters,index=[0]).to_dict('records'),
            editable=True,
            style_table={'overflowX': 'scroll'}),


            html.Br(),
            html.H3("Model_Parameters", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }),



            dash_table.DataTable(
            id='Model_Parameters',
            columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self.Model_Parameters.keys())],
            data=pd.DataFrame(self.Model_Parameters,index=[0]).to_dict('records'),
            editable=True,
            style_table={'overflowX': 'scroll'}),            

            html.Br(),

            html.H3("Initial_Conditions", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }),
            
            dash_table.DataTable(
            id='Initial_Conditions',
            columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self._IC.keys())],
            data=pd.DataFrame(self._IC,index=[0]).to_dict('records'),
            editable=True,
            style_table={'overflowX': 'scroll'}),

            html.Br(),
            html.H3("Inlet_Conditions", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }),
            dash_table.DataTable(
            id='Inlet_Conditions',
            columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self._InC.keys())],
            data=pd.DataFrame(self._InC,index=[0]).to_dict('records'),
            editable=True,
            style_table={'overflowX': 'scroll'},
            style_header={
            'backgroundColor': 'rgb(30, 30, 30)',
            'color': 'white'
                },
            style_data={
            'backgroundColor': 'rgb(50, 50, 50)',
            'color': 'white'
                })])
        

        
        
        
        
        
        
        
        app.layout = with_Report if self.Metagenome_Report else without_Report





        @app.callback(Output(component_id="Escher_", component_property='children'), Input(component_id="Drop_Down_Escher", component_property='value'))
        def escher_wrapper(Drop_Down_Escher):
            if Drop_Down_Escher=="Show Map":
                Labels={}
                for i in range(0,len(Sol.t),int(len(Sol.t)/20)):
                    Labels[i]={'label':str(int(Sol.t[i])),'style':{'color': '#77b0b1'}}
                Labels[len(Sol.t)-1]=self.Simtime
                return [html.Br(),html.H2("Time (Day)",style={'textAlign': 'center'}),dcc.Slider(0, len(Sol.t),len(Sol.t)/20,value=0,id="Escher_Slider",marks=Labels),html.Br()]

        @app.callback(Output(component_id="Escher", component_property='children'), Input(component_id="Drop_Down_Escher", component_property='value'),
        Input(component_id="Escher_Slider", component_property='value'))        
        def draw_escher(Drop_Down_Escher,Escher_Slider):
            rxn_data={}
    
            for ind,i in enumerate(self.Reactions):
                rxn_data[i.replace(" ","_")]=self.Info["Fluxes"][int(Escher_Slider)][ind]
            
            if Drop_Down_Escher=="Show Map":
                return [dash_escher.DashEscher(mapData=MD,modelData=MOD,
            options={
             'reaction_data':rxn_data
            }
            ,height='1000px',
        width='100%')
             ]

        @app.callback(Output(component_id="Annotation_Graph", component_property='figure'), Input(component_id='Drop_Down', component_property='value'))
        def update_graph_fig(input_value):
            Reactions = self.Reactions
            id_list = []
            Label_List = []
            Parents_List = []
            Label_List.append(input_value)
            id_list.append("")
            Parents_List.append(id_list[0])

            Genome_List = self.Metagenome_Report[input_value]
            for item in Genome_List.keys():
                Label_List.append(item)
                id_list.append(item)
                Parents_List.append(input_value)

            for item in Genome_List.keys():
                C = 0
                if self.Metagenome_Report[input_value][item].__len__() > 0:
                    for ECs in self.Metagenome_Report[input_value][item].keys():
                        Label_List.append(ECs)
                        id_list.append(item+" : "+ECs)
                        Parents_List.append(item)
                        Label_List.append("<br>".join(
                            self.Metagenome_Report[input_value][item][ECs]))
                        id_list.append(item+" : "+ECs+" : " + str(C))
                        Parents_List.append(item+" : "+ECs)
                        C += 1

                else:
                    Label_List.append(
                        f"No Relevant EC Numbers were found for {item.split(';')[-1]} !")
                    Parents_List.append(item)

            fig = go.Figure(go.Treemap(
                ids=id_list,
                labels=Label_List,
                parents=Parents_List,
                values=list([1 for i in range(id_list.__len__())]),
                root_color="lightgray",
                maxdepth=2,
                marker_colorscale = 'RdBu'


            ))
            fig.update_layout(
            font=dict(size=25,color='white')


            )
            return fig

        @app.callback(Output(component_id='Concentrations_Line_Plot', component_property='figure'),
                    Input(component_id='Base_Parameters', component_property='data'),
                    Input(component_id='Model_Parameters', component_property='data'),
                    Input(component_id='Initial_Conditions', component_property='data'),
                    Input(component_id='Inlet_Conditions', component_property='data'))
        def update_graph_fig(Base_Parameters: dict, Model_Parameters:dict, Initial_Conditions: dict, Inlet_Conditions: dict,prevent_initial_call=True)->plotly.graph_objects.Figure:
            """A functions that imports Base_Parameters, Model_Parameters, Initial_Conditions, Inlet_Conditions into the plot.
            """
            if len(self.Control_States.keys()):
                for i in self.Control_States.keys():
                    self.Control_States[i]=Initial_Conditions[0][i]

            self.Base_Parameters = Base_Parameters[0]
            self.Model_Parameters = Model_Parameters[0]
            self.Initial_Conditions = np.array(
            [Initial_Conditions[0][i] for i in self.Species])[:, np.newaxis]
            self.Inlet_Conditions = np.array(
            [Inlet_Conditions[0][i+"_in"] for i in self.Species])[:, np.newaxis]
            Update_Sol = self.Solve_Model(
                         self.Initial_Conditions[:, 0], np.linspace(0, 30, 10000))


            Solution = {
                    't': Update_Sol.t,
                        }
            for i in range(len(self.Species)):
                Solution[self.Species[i]] = Update_Sol.y[i, :]
            Sol_DF = pd.DataFrame(Solution)

            fig = px.line(Sol_DF, x="t", y=Sol_DF.columns,
                  title="Concentration of species")
            fig.update_layout(
            title={
            'y': 0.95,
            'x': 0.5,
            "font_size": 30,
            'xanchor': 'center',
            'yanchor': 'top'},
            legend=dict(font=dict(size= 20))

                )
            fig.update_xaxes(
            title={
            "text": "Time (Days)",
            "font_size": 25,
                },
                 tickfont_size=20
                )
            fig.update_yaxes(
            title={
            "text": "Concentrations (kg COD/m^3)",
            "font_size": 25,
             },
             tickfont_size=20
            
                )
            fig.update_traces(line=dict(width=5))

            return fig
            


        app.run_server(port=8000, host='127.0.0.1')

    def S_to_DF(self,Address: str)->None:
        """Converts the matrix to a pandas data frame then to a csv"""

        DF = pd.DataFrame(self.S, columns=self.Reactions, index=self.Species)
        DF.to_csv(Address, header=True,
                  index=True)

    def CSV_Report(self,Sol: scipy.integrate._ivp.ivp.OdeResult ,Address: str)->None:
        """Converts the results to a pandas data frame then to a csv"""
        DF = pd.DataFrame(Sol.y, columns=Sol.t, index=self.Species)
        DF.to_csv(os.path.join(Address,self.Name+"_Report.csv"), header=True,
                  index=True)


def Build_ADM1_Stoiciometric_Matrix(Base_Parameters: dict, Model_Parameters: dict, Reactons: list, Species:list)-> np.ndarray:
    """This function builds the stoichiometric matrix for the original ADM Model.
    No testing is done.
    """
    if type(Base_Parameters)!= dict and type(Model_Parameters)!= dict:
        raise TypeError("Base_Parameters and Model_Parameters need to be dictionary")
    if type(Reactons)!= list and type(Species)!= list:
        raise TypeError("Reactions and Species must be list")


    S = np.zeros((len(Species), len(Reactons)))
    S[0, [1, 3, 4]] = [1, (1-Model_Parameters["f_fa_li"]), - 1]
    S[1, [2, 5]] = [1, -1]
    S[2, [3, 6]] = [(Model_Parameters["f_fa_li"]), - 1]
    S[3, [5, 7]] = [(1-Model_Parameters['Y_aa']) *
                    Model_Parameters['f_va_aa'], - 1]
    S[4, [4, 5, 8]] = [(1-Model_Parameters['Y_su'])*Model_Parameters['f_bu_su'],
                       (1-Model_Parameters['Y_aa'])*Model_Parameters["f_bu_aa"], - 1]
    S[5, [4, 5, 7, 9]] = [(1-Model_Parameters["Y_su"])*Model_Parameters['f_pro_su'],
                          (1-Model_Parameters['Y_aa'])*Model_Parameters["f_pro_aa"], (1 - Model_Parameters['Y_c4'])*0.54, -1]
    S[6, [4, 5, 6, 7, 8, 9, 10]] = [(1-Model_Parameters['Y_su'])*Model_Parameters['f_ac_su'],
                                    (1-Model_Parameters['Y_aa']) *
                                    Model_Parameters['f_ac_aa'],
                                    (1-Model_Parameters['Y_fa'])*0.7,
                                    (1-Model_Parameters['Y_c4'])*0.31,
                                    (1-Model_Parameters['Y_c4'])*0.8,
                                    (1-Model_Parameters['Y_pro'])*0.57,
                                    -1]
    S[7, [4, 5, 6, 7, 8, 9, 11, 25]] = [(1-Model_Parameters['Y_su'])*Model_Parameters['f_h2_su'],
                                        (1-Model_Parameters['Y_aa']) *
                                        Model_Parameters['f_h2_aa'],
                                        (1-Model_Parameters['Y_fa'])*0.3,
                                        (1-Model_Parameters['Y_c4'])*0.15,
                                        (1-Model_Parameters['Y_c4'])*0.2,
                                        (1-Model_Parameters['Y_pro'])*0.43,
                                        -1,
                                        -1]
    S[8, [10, 11, 26]] = [(1-Model_Parameters['Y_ac']),
                          (1-Model_Parameters['Y_h2']),
                          -1]
    s_1 = (-1 * Model_Parameters['C_xc'] + Model_Parameters['f_sI_xc'] * Model_Parameters['C_sI'] + Model_Parameters['f_ch_xc'] * Model_Parameters['C_ch'] +
           Model_Parameters['f_pr_xc'] * Model_Parameters['C_pr'] + Model_Parameters['f_li_xc'] * Model_Parameters['C_li'] + Model_Parameters['f_xI_xc'] * Model_Parameters['C_xI'])
    s_2 = (-1 * Model_Parameters['C_ch'] + Model_Parameters['C_su'])
    s_3 = (-1 * Model_Parameters['C_pr'] + Model_Parameters['C_aa'])
    s_4 = (-1 * Model_Parameters['C_li'] + (1 - Model_Parameters['f_fa_li']) *
           Model_Parameters['C_su'] + Model_Parameters['f_fa_li'] * Model_Parameters['C_fa'])
    s_5 = (-1 * Model_Parameters['C_su'] + (1 - Model_Parameters['Y_su']) * (Model_Parameters['f_bu_su'] * Model_Parameters['C_bu'] + Model_Parameters['f_pro_su']
                                                                             * Model_Parameters['C_pro'] + Model_Parameters['f_ac_su'] * Model_Parameters['C_ac']) + Model_Parameters['Y_su'] * Model_Parameters['C_bac'])
    s_6 = (-1 * Model_Parameters['C_aa'] + (1 - Model_Parameters['Y_aa']) * (Model_Parameters['f_va_aa'] * Model_Parameters['C_va'] + Model_Parameters['f_bu_aa'] * Model_Parameters['C_bu'] +
                                                                             Model_Parameters['f_pro_aa'] * Model_Parameters['C_pro'] + Model_Parameters['f_ac_aa'] * Model_Parameters['C_ac']) + Model_Parameters['Y_aa'] * Model_Parameters['C_bac'])
    s_7 = (-1 * Model_Parameters['C_fa'] + (1 - Model_Parameters['Y_fa']) * 0.7 *
           Model_Parameters['C_ac'] + Model_Parameters['Y_fa'] * Model_Parameters['C_bac'])
    s_8 = (-1 * Model_Parameters['C_va'] + (1 - Model_Parameters['Y_c4']) * 0.54 * Model_Parameters['C_pro'] + (
        1 - Model_Parameters['Y_c4']) * 0.31 * Model_Parameters['C_ac'] + Model_Parameters['Y_c4'] * Model_Parameters['C_bac'])
    s_9 = (-1 * Model_Parameters['C_bu'] + (1 - Model_Parameters['Y_c4']) * 0.8 *
           Model_Parameters['C_ac'] + Model_Parameters['Y_c4'] * Model_Parameters['C_bac'])
    s_10 = (-1 * Model_Parameters['C_pro'] + (1 - Model_Parameters['Y_pro']) * 0.57 *
            Model_Parameters['C_ac'] + Model_Parameters['Y_pro'] * Model_Parameters['C_bac'])
    s_11 = (-1 * Model_Parameters['C_ac'] + (1 - Model_Parameters['Y_ac']) *
            Model_Parameters['C_ch4'] + Model_Parameters['Y_ac'] * Model_Parameters['C_bac'])
    s_12 = ((1 - Model_Parameters['Y_h2']) * Model_Parameters['C_ch4'] +
            Model_Parameters['Y_h2'] * Model_Parameters['C_bac'])
    s_13 = (-1 * Model_Parameters['C_bac'] + Model_Parameters['C_xc'])
    S[9, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 27]] = [-s_1, -s_2, -s_3, -s_4, -
                                                                                    s_5, -s_6, -s_7, -s_8, -s_9, -s_10, -s_11, -s_12, -s_13, -s_13, -s_13, -s_13, -s_13, -s_13, -s_13, -1]
    S[10, [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]] = [Model_Parameters['N_xc']-Model_Parameters['f_xI_xc']*Model_Parameters['N_I']-Model_Parameters['f_sI_xc']*Model_Parameters['N_I']-Model_Parameters['f_pr_xc']*Model_Parameters['N_aa'],
                                                                        -Model_Parameters['Y_su']*Model_Parameters['N_bac'],
                                                                        Model_Parameters['N_aa']-Model_Parameters['Y_aa'] *
                                                                        Model_Parameters['N_bac'],
                                                                        -Model_Parameters['Y_fa']*Model_Parameters['N_bac'],
                                                                        -Model_Parameters['Y_c4']*Model_Parameters['N_bac'],
                                                                        -Model_Parameters['Y_c4']*Model_Parameters['N_bac'],
                                                                        -Model_Parameters['Y_pro']*Model_Parameters['N_bac'],
                                                                        -Model_Parameters['Y_ac']*Model_Parameters['N_bac'],
                                                                        -Model_Parameters['Y_h2']*Model_Parameters['N_bac'],
                                                                        Model_Parameters['N_bac'] -
                                                                        Model_Parameters['N_xc'],
                                                                        Model_Parameters['N_bac'] -
                                                                        Model_Parameters['N_xc'],
                                                                        Model_Parameters['N_bac'] -
                                                                        Model_Parameters['N_xc'],
                                                                        Model_Parameters['N_bac'] -
                                                                        Model_Parameters['N_xc'],
                                                                        Model_Parameters['N_bac'] -
                                                                        Model_Parameters['N_xc'],
                                                                        Model_Parameters['N_bac'] -
                                                                        Model_Parameters['N_xc'],
                                                                        Model_Parameters['N_bac']-Model_Parameters['N_xc']]
    S[11, 0] = Model_Parameters['f_sI_xc']
    S[12, [0, 12, 13, 14, 15, 16, 17, 18]] = [-1, 1, 1, 1, 1, 1, 1, 1]
    S[13, [0, 1]] = [Model_Parameters['f_ch_xc'], -1]
    S[14, [0, 2]] = [Model_Parameters['f_pr_xc'], -1]
    S[15, [0, 3]] = [Model_Parameters['f_li_xc'], -1]
    S[16, [4, 12]] = [Model_Parameters['Y_su'], -1]
    S[17, [5, 13]] = [Model_Parameters['Y_aa'], -1]
    S[18, [6, 14]] = [Model_Parameters['Y_fa'], -1]
    S[19, [7, 8, 15]] = [Model_Parameters['Y_c4'], Model_Parameters['Y_c4'], -1]
    S[20, [9, 16]] = [Model_Parameters['Y_pro'], -1]
    S[21, [10, 17]] = [Model_Parameters['Y_ac'], -1]
    S[22, [11, 18]] = [Model_Parameters['Y_h2'], -1]
    S[23, 0] = Model_Parameters['f_xI_xc']
    S[24, :] = 0
    S[25, :] = 0
    S[26, :] = 0
    S[27, 19] = -1
    S[28, 20] = -1
    S[29, 21] = -1
    S[30, 22] = -1
    S[31, 23] = -1
    S[32, :] = 0
    S[33, 24] = -1
    S[34, :] = 0
    S[35, 25] = Base_Parameters['V_liq']/Base_Parameters['V_gas']
    S[36, 26] = Base_Parameters['V_liq']/Base_Parameters['V_gas']
    S[37, 27] = Base_Parameters['V_liq']/Base_Parameters['V_gas']
    return S


def ADM1_ODE_Sys(t: float, c: np.ndarray, ADM1_Instance:Model)-> np.ndarray:
    """ The ODE system for the original ADM.
        No testing is done.

        Args:
            t (float):a matrix of zeros to be filled
            c (np.ndarray): an array of concentrations to be filled
            Model (Model): The model to calculate ODE with

        Returns:
            np.ndarray: The output is dCdt, the change of concentration with respect to time.
    """
    c[34] = c[10] - c[33]
    c[32] = c[9] - c[31]
    I_pH_aa = (ADM1_Instance.Model_Parameters["K_pH_aa"] ** ADM1_Instance.Model_Parameters['nn_aa'])/(np.power(
        c[26], ADM1_Instance.Model_Parameters['nn_aa']) + np.power(ADM1_Instance.Model_Parameters["K_pH_aa"], ADM1_Instance.Model_Parameters['nn_aa']))
    I_pH_ac = (ADM1_Instance.Model_Parameters['K_pH_ac'] ** ADM1_Instance.Model_Parameters["n_ac"])/(
        c[26] ** ADM1_Instance.Model_Parameters['n_ac'] + ADM1_Instance.Model_Parameters['K_pH_ac'] ** ADM1_Instance.Model_Parameters['n_ac'])
    I_pH_h2 = (ADM1_Instance.Model_Parameters['K_pH_h2']**ADM1_Instance.Model_Parameters['n_h2'])/(
        c[26] ** ADM1_Instance.Model_Parameters['n_h2'] + ADM1_Instance.Model_Parameters['K_pH_h2']**ADM1_Instance.Model_Parameters['n_h2'])
    I_IN_lim = 1 / (1+(ADM1_Instance.Model_Parameters['K_S_IN'] / c[10]))
    I_h2_fa = 1 / (1+(c[7] / ADM1_Instance.Model_Parameters['K_I_h2_fa']))
    I_h2_c4 = 1 / (1+(c[7]/ADM1_Instance.Model_Parameters['K_I_h2_c4']))
    I_h2_pro = (1/(1+(c[7]/ADM1_Instance.Model_Parameters['K_I_h2_pro'])))
    I_nh3 = 1/(1+(c[33]/ADM1_Instance.Model_Parameters['K_I_nh3']))
    I5 = (I_pH_aa * I_IN_lim)
    I6 = np.copy(I5)
    I7 = (I_pH_aa * I_IN_lim * I_h2_fa)
    I8 = (I_pH_aa * I_IN_lim * I_h2_c4)
    I9 = np.copy(I8)
    I10 = (I_pH_aa * I_IN_lim * I_h2_pro)
    I11 = (I_pH_ac * I_IN_lim * I_nh3)
    I12 = (I_pH_h2 * I_IN_lim)
    v = np.zeros((len(ADM1_Instance.Reactions), 1))
    v[0] = ADM1_Instance.Model_Parameters["k_dis"]*c[12]
    v[1] = ADM1_Instance.Model_Parameters['k_hyd_ch']*c[13]
    v[2] = ADM1_Instance.Model_Parameters['k_hyd_pr']*c[14]
    v[3] = ADM1_Instance.Model_Parameters['k_hyd_li']*c[15]
    v[4] = ADM1_Instance.Model_Parameters['k_m_su']*c[0] / \
        (ADM1_Instance.Model_Parameters['K_S_su']+c[0])*c[16]*I5
    v[5] = ADM1_Instance.Model_Parameters['k_m_aa']*c[1] / \
        (ADM1_Instance.Model_Parameters['K_S_aa']+c[1])*c[17]*I6
    v[6] = ADM1_Instance.Model_Parameters['k_m_fa']*c[2] / \
        (ADM1_Instance.Model_Parameters['K_S_fa']+c[2])*c[18]*I7
    v[7] = ADM1_Instance.Model_Parameters['k_m_c4']*c[3] / \
        (ADM1_Instance.Model_Parameters['K_S_c4']+c[3]) * \
        c[19]*c[3]/(c[3]+c[4]+10 ** (-6))*I8
    v[8] = ADM1_Instance.Model_Parameters['k_m_c4']*c[4] / \
        (ADM1_Instance.Model_Parameters['K_S_c4']+c[4]) * \
        c[19]*c[4]/(c[4]+c[3]+10 ** (-6))*I9
    v[9] = ADM1_Instance.Model_Parameters['k_m_pr']*c[5] / \
        (ADM1_Instance.Model_Parameters['K_S_pro']+c[5])*c[20]*I10
    v[10] = ADM1_Instance.Model_Parameters['k_m_ac']*c[6] / \
        (ADM1_Instance.Model_Parameters['K_S_ac']+c[6])*c[21]*I11
    v[11] = ADM1_Instance.Model_Parameters['k_m_h2']*c[7] / \
        (ADM1_Instance.Model_Parameters['K_S_h2']+c[7])*c[22]*I12
    v[12] = ADM1_Instance.Model_Parameters['k_dec_X_su']*c[16]
    v[13] = ADM1_Instance.Model_Parameters['k_dec_X_aa']*c[17]
    v[14] = ADM1_Instance.Model_Parameters['k_dec_X_fa']*c[18]
    v[15] = ADM1_Instance.Model_Parameters['k_dec_X_c4']*c[19]
    v[16] = ADM1_Instance.Model_Parameters['k_dec_X_pro']*c[20]
    v[17] = ADM1_Instance.Model_Parameters['k_dec_X_ac']*c[21]
    v[18] = ADM1_Instance.Model_Parameters['k_dec_X_h2']*c[22]
    v[19] = ADM1_Instance.Model_Parameters['k_A_B_va'] * \
        (c[27] * (ADM1_Instance.Model_Parameters['K_a_va'] + c[26]) -
         ADM1_Instance.Model_Parameters['K_a_va'] * c[3])
    v[20] = ADM1_Instance.Model_Parameters['k_A_B_bu'] * \
        (c[28] * (ADM1_Instance.Model_Parameters['K_a_bu'] + c[26]) -
         ADM1_Instance.Model_Parameters['K_a_bu'] * c[4])
    v[21] = ADM1_Instance.Model_Parameters['k_A_B_pro'] * \
        (c[29] * (ADM1_Instance.Model_Parameters['K_a_pro'] + c[26]) -
         ADM1_Instance.Model_Parameters['K_a_pro'] * c[5])
    v[22] = ADM1_Instance.Model_Parameters['k_A_B_ac'] * \
        (c[30] * (ADM1_Instance.Model_Parameters['K_a_ac'] + c[26]) -
         ADM1_Instance.Model_Parameters['K_a_ac'] * c[6])
    v[23] = ADM1_Instance.Model_Parameters['k_A_B_co2'] * \
        (c[31] * (ADM1_Instance.Model_Parameters['K_a_co2'] + c[26]) -
         ADM1_Instance.Model_Parameters['K_a_co2'] * c[9])
    v[24] = ADM1_Instance.Model_Parameters['k_A_B_IN'] * \
        (c[33] * (ADM1_Instance.Model_Parameters['K_a_IN'] + c[26]) -
         ADM1_Instance.Model_Parameters['K_a_IN'] * c[10])
    p_gas_h2 = c[35] * ADM1_Instance.Base_Parameters["R"] * \
        ADM1_Instance.Base_Parameters["T_op"] / 16
    p_gas_ch4 = c[36] * ADM1_Instance.Base_Parameters["R"] * \
        ADM1_Instance.Base_Parameters["T_op"] / 64
    p_gas_co2 = c[37] * ADM1_Instance.Base_Parameters["R"] * \
        ADM1_Instance.Base_Parameters["T_op"]
    p_gas_h2o = 0.0313 * \
        np.exp(5290 *
               (1 / ADM1_Instance.Base_Parameters["T_base"] - 1 / ADM1_Instance.Base_Parameters["T_op"]))
    P_gas = p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o
    q_gas = max(
        0, (ADM1_Instance.Model_Parameters['k_p'] * (P_gas - ADM1_Instance.Base_Parameters['P_atm'])))
    v[25] = ADM1_Instance.Model_Parameters['k_L_a'] * \
        (c[7] - 16 * ADM1_Instance.Model_Parameters['K_H_h2'] * p_gas_h2)
    v[26] = ADM1_Instance.Model_Parameters['k_L_a'] * \
        (c[8] - 64 * ADM1_Instance.Model_Parameters['K_H_ch4'] * p_gas_ch4)
    v[27] = ADM1_Instance.Model_Parameters['k_L_a'] * \
        (c[32] - ADM1_Instance.Model_Parameters['K_H_co2'] * p_gas_co2)
    dCdt = np.matmul(ADM1_Instance.S, v)
    phi = c[24]+c[34]-c[31] - (c[30] / 64) - (c[29] /
                                              112) - (c[28] / 160) - (c[27] / 208) - c[25]
    c[26] = (-1 * phi / 2) + \
        (0.5 * np.sqrt(phi**2 + 4 * ADM1_Instance.Model_Parameters['K_w']))
    dCdt[0: 35] = dCdt[0: 35]+ADM1_Instance.Base_Parameters['q_in'] / \
        ADM1_Instance.Base_Parameters["V_liq"] * \
        (ADM1_Instance.Inlet_Conditions[0: 35]-c[0:35].reshape(-1, 1))
    dCdt[35:] = dCdt[35:]+q_gas/ADM1_Instance.Base_Parameters["V_gas"] * \
        (ADM1_Instance.Inlet_Conditions[35:]-c[35:].reshape(-1, 1))
    dCdt[[26, 32, 34], 0] = 0
    if ADM1_Instance.Switch == "DAE":
        dCdt[7] = 0
        dCdt[27: 32] = 0
        dCdt[33] = 0
    return dCdt[:, 0]


def Build_Modified_ADM1_Stoiciometric_Matrix(Base_Parameters: dict, Model_Parameters: dict, Reactions: list, Species: list)->np.ndarray:
    """ 
    This function builds the stoichiometric matrix for the modified ADM Model.
        
        Model Parameters (dict): a dictionary which contains model parameters
        Base_Parameters (dict): a dictionary which contains base paramters
        Initial Conditions (dict): a dictionary containing inlet conditions for all species
        Inlet Conditions (dict): a dictionary containing inlet conditions for all species
        Reactions (list): a list containing all of the reaction names
        Species (list): a list containing all species
        
    Returns:
        np.ndarray: Returns an matrix of stochiometic values.
    """
    S = np.zeros((len(Species), len(Reactions)))
    S[list(map(Species.index, ["TSS", "X_ch", "X_pr", "X_li", "X_I"])),
      Reactions.index('TSS_Disintegration')] = [-1, Model_Parameters['f_ch_TSS'], Model_Parameters['f_pr_TSS'], Model_Parameters['f_li_TSS'], Model_Parameters['f_xI_TSS']]
    S[list(map(Species.index, ["TDS", "X_ch", "X_pr", "X_li", "S_I"])), Reactions.index('TDS_Disintegration')] = [-1,
                                                                                                                  Model_Parameters['f_ch_TDS'], Model_Parameters['f_pr_TDS'], Model_Parameters['f_li_TDS'], Model_Parameters['f_sI_TDS']]
    S[list(map(Species.index, ["X_ch", "S_su"])),
      Reactions.index('Hydrolysis carbohydrates')] = [-1, 1]
    S[list(map(Species.index, ["X_pr", "S_aa"])),
      Reactions.index('Hydrolysis proteins')] = [-1, 1]
    S[list(map(Species.index, ["X_li", "S_fa"])),
      Reactions.index('Hydrolysis lipids')] = [-1, 1]

    f_IC_su = -(-Model_Parameters['C_su'] +
                (1-Model_Parameters['Y_su'])*Model_Parameters['f_pro_su']*Model_Parameters['C_pro'] +
                (1-Model_Parameters['Y_su'])*Model_Parameters['f_et_su']*Model_Parameters['C_et'] +
                (1-Model_Parameters['Y_su'])*Model_Parameters['f_lac_su']*Model_Parameters['C_lac'] +
                (1-Model_Parameters['Y_su'])*Model_Parameters['f_ac_su']*Model_Parameters['C_ac'] +
                Model_Parameters['Y_su']*Model_Parameters['C_bac'])

    S[list(map(Species.index, ["S_su", "S_pro", "S_et", "S_lac", "S_ac", "S_IN", "S_IC", "X_su"])),
      Reactions.index('Uptake of sugars')] = [-1,
                                              (1-Model_Parameters['Y_su']) *
                                              Model_Parameters['f_pro_su'],
                                              (1-Model_Parameters['Y_su']) *
                                              Model_Parameters['f_et_su'],
                                              (1-Model_Parameters['Y_su']) *
                                              Model_Parameters['f_lac_su'],
                                              (1-Model_Parameters['Y_su']) *
                                              Model_Parameters['f_ac_su'],
                                              -Model_Parameters['N_bac']*Model_Parameters['Y_su'],
                                              f_IC_su,
                                              Model_Parameters['Y_su']]

    f_IC_aa = -(-Model_Parameters['C_aa'] +
                (1-Model_Parameters['Y_aa'])*Model_Parameters['f_pro_aa']*Model_Parameters['C_pro'] +
                (1-Model_Parameters['Y_aa'])*Model_Parameters['f_et_aa']*Model_Parameters['C_et'] +
                (1-Model_Parameters['Y_aa'])*Model_Parameters['f_lac_aa']*Model_Parameters['C_lac'] +
                (1-Model_Parameters['Y_aa'])*Model_Parameters['f_ac_aa']*Model_Parameters['C_ac'] +
                (1-Model_Parameters['Y_aa'])*Model_Parameters['C_bac'])

    S[list(map(Species.index, ["S_aa", "S_pro", "S_et", "S_lac", "S_ac", "S_IN", "S_IC", "X_aa"])),
      Reactions.index('Uptake of amino acids')] = [-1,
                                                   (1-Model_Parameters['Y_aa']) *
                                                   Model_Parameters['f_pro_aa'],
                                                   (1-Model_Parameters['Y_aa']) *
                                                   Model_Parameters['f_et_aa'],
                                                   (1-Model_Parameters['Y_aa']) *
                                                   Model_Parameters['f_lac_aa'],
                                                   (1-Model_Parameters['Y_aa']) *
                                                   Model_Parameters['f_ac_aa'],
                                                   Model_Parameters['N_aa']-Model_Parameters['Y_aa'] *
                                                   Model_Parameters['N_bac'],
                                                   f_IC_aa,
                                                   Model_Parameters['Y_aa']]
    f_IC_fa = -(-Model_Parameters['C_fa'] +
                (1-Model_Parameters['Y_fa'])*Model_Parameters['f_pro_fa']*Model_Parameters['C_pro'] +
                (1-Model_Parameters['Y_fa'])*Model_Parameters['f_et_fa']*Model_Parameters['C_et'] +
                (1-Model_Parameters['Y_fa'])*Model_Parameters['f_lac_fa']*Model_Parameters['C_lac'] +
                (1-Model_Parameters['Y_fa'])*Model_Parameters['f_ac_fa']*Model_Parameters['C_ac'] +
                (1-Model_Parameters['Y_fa'])*Model_Parameters['C_bac'])

    S[list(map(Species.index, ["S_fa", "S_pro", "S_et", "S_lac", "S_ac", "S_IN", "S_IC", "X_fa"])),
      Reactions.index('Uptake of LCFA')] = [-1,
                                            (1-Model_Parameters['Y_fa']) *
                                            Model_Parameters['f_pro_fa'],
                                            (1-Model_Parameters['Y_fa']) *
                                            Model_Parameters['f_et_fa'],
                                            (1-Model_Parameters['Y_fa']) *
                                            Model_Parameters['f_lac_fa'],
                                            (1-Model_Parameters['Y_fa']) *
                                            Model_Parameters['f_ac_fa'],
                                            -Model_Parameters['Y_fa'] *
                                            Model_Parameters['N_bac'],
                                            f_IC_fa,
                                            Model_Parameters['Y_fa']]

    f_IC_ac_et = -(-Model_Parameters['C_ac'] +
                   (1-Model_Parameters['Y_ac_et'])*Model_Parameters['f_et_ac']*Model_Parameters['C_et'] +
                   (1-Model_Parameters['Y_ac_et']) *
                   Model_Parameters['f_bu_ac']*Model_Parameters['C_bu'] +
                   (1-Model_Parameters['Y_ac_et'])*Model_Parameters['C_bac'])

    f_IC_ac_lac = -(-Model_Parameters['C_ac'] +
                    (1-Model_Parameters['Y_ac_lac'])*Model_Parameters['f_lac_ac']*Model_Parameters['C_lac'] +
                    (1-Model_Parameters['Y_ac_lac']) *
                    Model_Parameters['f_bu_ac']*Model_Parameters['C_bu'] +
                    (1-Model_Parameters['Y_ac_lac'])*Model_Parameters['C_bac'])

    S[list(map(Species.index, ["S_ac", "S_et", "S_bu", "S_IN", "S_IC", "S_h2", "X_ac_et"])),
      Reactions.index('Uptake of acetate_et')] = [-1,
                                                  (1-Model_Parameters['Y_ac_et']) *
                                                  Model_Parameters['f_et_ac'],
                                                  (1-Model_Parameters['Y_ac']) *
                                                  Model_Parameters['f_bu_ac'],
                                                  f_IC_ac_et,
                                                  -Model_Parameters['Y_ac_et'] *
                                                  Model_Parameters['N_bac'],
                                                  (1-Model_Parameters['Y_ac_et']) *
                                                  Model_Parameters['f_h2_ac'],
                                                  Model_Parameters['Y_ac_et']]

    S[list(map(Species.index, ["S_ac", "S_lac", "S_bu", "S_IN", "S_IC", "S_h2", "X_ac_lac"])),
        Reactions.index('Uptake of acetate_lac')] = [-1,
                                                     (1-Model_Parameters['Y_ac_lac']) *
                                                     Model_Parameters['f_lac_ac'],
                                                     (1-Model_Parameters['Y_ac_lac']) *
                                                     Model_Parameters['f_bu_ac'],
                                                     f_IC_ac_lac,
                                                     -Model_Parameters['Y_ac_lac'] *
                                                     Model_Parameters['N_bac'],
                                                     (1-Model_Parameters['Y_ac_lac']) *
                                                     Model_Parameters['f_h2_ac'],
                                                     Model_Parameters['Y_ac_lac']]

    f_IC_pro_et = -(-Model_Parameters['C_pro'] +
                    (1-Model_Parameters['Y_pro_et'])*Model_Parameters['f_et_pro']*Model_Parameters['C_et'] +
                    (1-Model_Parameters['Y_pro_et'])*Model_Parameters['f_va_pro']*Model_Parameters['C_va'] +
                    (1-Model_Parameters['Y_pro_et'])*Model_Parameters['C_bac'])

    f_IC_pro_lac = -(-Model_Parameters['C_pro'] +
                     (1-Model_Parameters['Y_pro_lac'])*Model_Parameters['f_lac_pro']*Model_Parameters['C_lac'] +
                     (1-Model_Parameters['Y_pro_lac'])*Model_Parameters['f_va_pro']*Model_Parameters['C_va'] +
                     (1-Model_Parameters['Y_pro_lac'])*Model_Parameters['C_bac'])

    S[list(map(Species.index, ["S_pro", "S_et", "S_va", "S_IN", "S_IC", "S_h2", "X_chain_et"])),
      Reactions.index('Uptake of propionate_et')] = [-1,
                                                     (1-Model_Parameters['Y_pro_et']) *
                                                     Model_Parameters['f_et_pro'],
                                                     (1-Model_Parameters['Y_pro_et']) *
                                                     Model_Parameters['f_va_pro'],
                                                     f_IC_pro_et,
                                                     -Model_Parameters['Y_pro_et'] *
                                                     Model_Parameters['N_bac'],
                                                     (1-Model_Parameters['Y_pro_et']) *
                                                     Model_Parameters['f_h2_pro'],
                                                     Model_Parameters['Y_chain_et_pro']]

    S[list(map(Species.index, ["S_pro", "S_lac", "S_va", "S_IN", "S_IC", "S_h2", "X_chain_lac"])),
        Reactions.index('Uptake of propionate_lac')] = [-1,
                                                        (1-Model_Parameters['Y_pro_lac']) *
                                                        Model_Parameters['f_lac_pro'],
                                                        (1-Model_Parameters['Y_pro_lac']) *
                                                        Model_Parameters['f_va_pro'],
                                                        f_IC_pro_lac,
                                                        -Model_Parameters['Y_pro_lac'] *
                                                        Model_Parameters['N_bac'],
                                                        (1-Model_Parameters['Y_pro_lac']) *
                                                        Model_Parameters['f_h2_pro'],
                                                        Model_Parameters['Y_chain_lac_pro']]

    f_IC_bu_et = -(-Model_Parameters['C_bu'] +
                   (1-Model_Parameters['Y_bu_et'])*Model_Parameters['f_et_bu']*Model_Parameters['C_et'] +
                   (1-Model_Parameters['Y_bu_et'])*Model_Parameters['f_cap_bu']*Model_Parameters['C_cap'] +
                   (1-Model_Parameters['Y_bu_et'])*Model_Parameters['C_bac'])

    f_IC_bu_lac = -(-Model_Parameters['C_bu'] +
                    (1-Model_Parameters['Y_bu_lac'])*Model_Parameters['f_lac_bu']*Model_Parameters['C_lac'] +
                    (1-Model_Parameters['Y_bu_lac'])*Model_Parameters['f_cap_bu']*Model_Parameters['C_cap'] +
                    (1-Model_Parameters['Y_bu_lac'])*Model_Parameters['C_bac'])

    S[list(map(Species.index, ["S_bu", "S_et", "S_cap", "S_IN", "S_IC", "S_h2", "X_chain_et"])),
        Reactions.index('Uptake of butyrate_et')] = [-1,
                                                     (1-Model_Parameters['Y_bu_et']
                                                      ) * Model_Parameters['f_et_bu'],
                                                     (1-Model_Parameters['Y_bu_et']) *
                                                     Model_Parameters['f_cap_bu'],
                                                     f_IC_bu_et,
                                                     -Model_Parameters['Y_bu_et'] * Model_Parameters['N_bac'],
                                                     (1-Model_Parameters['Y_bu_et']
                                                      )*Model_Parameters['f_h2_bu'],
                                                     Model_Parameters['Y_bu_et']]

    S[list(map(Species.index, ["S_bu", "S_lac", "S_cap", "S_IN", "S_IC", "S_h2", "X_chain_lac"])),
        Reactions.index('Uptake of butyrate_lac')] = [-1,
                                                      (1-Model_Parameters['Y_bu_lac']) *
                                                      Model_Parameters['f_lac_bu'],
                                                      (1-Model_Parameters['Y_bu_lac']) *
                                                      Model_Parameters['f_cap_bu'],
                                                      f_IC_bu_lac,
                                                      -Model_Parameters['Y_bu_lac'] *
                                                      Model_Parameters['N_bac'],
                                                      (1-Model_Parameters['Y_bu_lac']
                                                       )*Model_Parameters['f_h2_bu'],
                                                      Model_Parameters['Y_bu_lac']]

    S[list(map(Species.index, ["S_va", "S_ac", "X_VFA_deg"])),
        Reactions.index('Uptake of valerate')] = [-1,
                                                  (1-Model_Parameters['Y_va']),
                                                  Model_Parameters['Y_va']]

    S[list(map(Species.index, ["S_cap", "S_ac", "X_VFA_deg"])),
        Reactions.index('Uptake of caproate')] = [-1,
                                                  (1 -
                                                   Model_Parameters['Y_cap']),
                                                  Model_Parameters['Y_cap']]
    f_IC_Me_ach2 = 0
    S[list(map(Species.index, ["S_h2", "S_ac", "S_ch4", "X_Me_ac", 'S_IC'])),
        Reactions.index('Methanogenessis from acetate and h2')] = [-1,
                                                                   (1 - Model_Parameters['Y_h2_ac']
                                                                    )*Model_Parameters['f_ac_h2'],
                                                                   (1 -
                                                                       Model_Parameters['Y_Me_ac']),
                                                                   Model_Parameters['Y_Me_ac'],
                                                                   f_IC_Me_ach2]

    f_IC_Me_CO2h2 = -(Model_Parameters['Y_Me_CO2']*Model_Parameters['C_ch4'] +
                      Model_Parameters['Y_Me_h2']*Model_Parameters['C_bac'])
    
    S[list(map(Species.index, ["S_h2", "S_ch4", "X_Me_CO2", 'S_IC'])),
        Reactions.index('Methanogenessis from CO2 and h2')] = [-1,
                                                               (1 -
                                                                Model_Parameters['Y_h2_CO2']),
                                                               (Model_Parameters['Y_Me_CO2']),
                                                               f_IC_Me_CO2h2]
    
    f_IC_et_ox=-(-Model_Parameters['C_et'] +
                    (1-Model_Parameters['Y_ac_et_ox'])*Model_Parameters['C_bac']
                    +Model_Parameters['Y_ac_et_ox']*Model_Parameters['C_ac'])

    S[list(map(Species.index, ["S_et", "X_et","S_ac","S_IC"])),
        Reactions.index('Uptake of ethanol')] = [-1,1-Model_Parameters['Y_ac_et_ox'],Model_Parameters['Y_ac_et_ox'],f_IC_et_ox]

    
    
    f_IC_lac_ox=-(-Model_Parameters['C_lac'] +
                (1-Model_Parameters['Y_pro_lac_ox'])*Model_Parameters['C_bac']
                +Model_Parameters['Y_pro_lac_ox']*Model_Parameters['C_pro'])
    
    S[list(map(Species.index, ["S_lac", "X_lac","S_pro","S_IC"])),
        Reactions.index('Uptake of lactate')] = [-1, 1-Model_Parameters['Y_pro_lac_ox'],Model_Parameters['Y_pro_lac_ox'],f_IC_lac_ox]

    S[list(map(Species.index, ["X_su", "TSS"])),
        Reactions.index('Decay of Xsu')] = [-1, 1]

    S[list(map(Species.index, ["X_aa", "TSS"])),
        Reactions.index('Decay of Xaa')] = [-1, 1]

    S[list(map(Species.index, ["X_fa", "TSS"])),
        Reactions.index('Decay of Xfa')] = [-1, 1]

    S[list(map(Species.index, ["X_ac_et", "TSS"])),
        Reactions.index('Decay of X_ac_et')] = [-1, 1]

    S[list(map(Species.index, ["X_ac_lac", "TSS"])),
        Reactions.index('Decay of X_ac_lac')] = [-1, 1]

    S[list(map(Species.index, ["X_chain_et", "TSS"])),
        Reactions.index('Decay of X_chain_et')] = [-1, 1]

    S[list(map(Species.index, ["X_chain_lac", "TSS"])),
        Reactions.index('Decay of X_chain_lac')] = [-1, 1]

    S[list(map(Species.index, ["X_VFA_deg", "TSS"])),
        Reactions.index('Decay of X_VFA_deg')] = [-1, 1]

    S[list(map(Species.index, ["X_Me_ac", "TSS"])),
        Reactions.index('Decay of X_Me_ac')] = [-1, 1]

    S[list(map(Species.index, ["X_Me_CO2", "TSS"])),
        Reactions.index('Decay of X_Me_CO2')] = [-1, 1]

    S[list(map(Species.index, ["S_va_ion"])),
        Reactions.index('Acid Base Equilibrium (Va)')] = [-1]

    S[list(map(Species.index, ["S_bu_ion"])),
        Reactions.index('Acid Base Equilibrium (Bu)')] = [-1]

    S[list(map(Species.index, ["S_pro_ion"])),
        Reactions.index('Acid Base Equilibrium (Pro)')] = [-1]

    S[list(map(Species.index, ["S_cap_ion"])),
        Reactions.index('Acid Base Equilibrium (Cap)')] = [-1]

    S[list(map(Species.index, ["S_lac_ion"])),
        Reactions.index('Acid Base Equilibrium (Lac)')] = [-1]

    S[list(map(Species.index, ["S_ac_ion"])),
        Reactions.index('Acid Base Equilibrium (Ac)')] = [-1]

    # S[list(map(Species.index, ["S_CO2", "S_hco3_ion"])),  # I don't think this is right، should look at the reaction in ADM1
    #     Reactions.index('Acid Base Equilibrium (CO2)')] = [-1, 1]

    # S[list(map(Species.index, ["S_nh3", "S_nh4_ion"])),
    #     Reactions.index('Acid Base Equilibrium (IN)')] = [-1, 1]  # I don't think this is right، should look at the reaction in ADM1

    S[list(map(Species.index, ["S_h2", "S_gas_h2"])),
        Reactions.index('Gas Transfer H2')] = [-Base_Parameters['V_liq']/Base_Parameters['V_gas'], 1]
    S[list(map(Species.index, ["S_ch4", "S_gas_ch4"])),
        Reactions.index('Gas Transfer CH4')] = [-Base_Parameters['V_liq']/Base_Parameters['V_gas'], 1]
    S[list(map(Species.index, ["S_co2", "S_gas_co2"])),
        Reactions.index('Gas Transfer CO2')] = [-Base_Parameters['V_liq']/Base_Parameters['V_gas'], 1]
    return S


def Modified_ADM1_ODE_Sys(t: float, c: np.ndarray, Model: Model)-> np.ndarray:
    """
    This function is used to build the ODEs of the modified ADM1 model.
    
    Args:
        t (float):a matrix of zeros to be filled
        c (np.ndarray): an array of concentrations to be filled
        Model (Model): The model to calculate ODE with

    Returns:
        np.ndarray: The output is dCdt, the change of concentration with respect to time. 
    """
    c[c<0]=0
    c[Model.Species.index('S_nh4_ion')] = c[Model.Species.index(
        'S_IN')] - c[Model.Species.index('S_nh3')]
    c[Model.Species.index('S_co2')] = c[Model.Species.index(
        'S_IC')] - c[Model.Species.index('S_hco3_ion')]
    I_pH_aa = (Model.Model_Parameters["K_pH_aa"] ** Model.Model_Parameters['nn_aa'])/(np.power(
        c[Model.Species.index('S_H_ion')], Model.Model_Parameters['nn_aa']) + np.power(Model.Model_Parameters["K_pH_aa"], Model.Model_Parameters['nn_aa']))

    I_pH_ac = (Model.Model_Parameters['K_pH_ac'] ** Model.Model_Parameters["n_ac"])/(
        c[Model.Species.index('S_H_ion')] ** Model.Model_Parameters['n_ac'] + Model.Model_Parameters['K_pH_ac'] ** Model.Model_Parameters['n_ac'])

    I_pH_pro = (Model.Model_Parameters['K_pH_pro'] ** Model.Model_Parameters["n_pro"])/(
        c[Model.Species.index('S_H_ion')] ** Model.Model_Parameters['n_pro'] + Model.Model_Parameters['K_pH_pro'] ** Model.Model_Parameters['n_pro'])

    I_pH_bu = (Model.Model_Parameters['K_pH_bu'] ** Model.Model_Parameters["n_bu"])/(
        c[Model.Species.index('S_H_ion')] ** Model.Model_Parameters['n_bu'] + Model.Model_Parameters['K_pH_bu'] ** Model.Model_Parameters['n_bu'])

    I_pH_va = (Model.Model_Parameters['K_pH_va'] ** Model.Model_Parameters["n_va"])/(
        c[Model.Species.index('S_H_ion')] ** Model.Model_Parameters['n_va'] + Model.Model_Parameters['K_pH_va'] ** Model.Model_Parameters['n_va'])

    I_pH_cap = (Model.Model_Parameters['K_pH_cap'] ** Model.Model_Parameters["n_cap"])/(
        c[Model.Species.index('S_H_ion')] ** Model.Model_Parameters['n_cap'] + Model.Model_Parameters['K_pH_cap'] ** Model.Model_Parameters['n_cap'])

    I_pH_h2 = (Model.Model_Parameters['K_pH_h2']**Model.Model_Parameters['n_h2'])/(
        c[Model.Species.index('S_H_ion')] ** Model.Model_Parameters['n_h2'] + Model.Model_Parameters['K_pH_h2']**Model.Model_Parameters['n_h2'])

    I_IN_lim = 1 / \
        (1+(Model.Model_Parameters['K_S_IN'] / (c[Model.Species.index('S_IN')]+10**-9)))

    I_h2_fa = 1 / (1+(c[Model.Species.index('S_h2')] /
                   (Model.Model_Parameters['K_I_h2_fa']+10**-9)))

    I_h2_c4 = 1 / (1+(c[Model.Species.index('S_h2')] /
                   (Model.Model_Parameters['K_I_h2_c4']+10**-9)))

    I_h2_pro = (1/(1+(c[Model.Species.index('S_h2')] /
                (Model.Model_Parameters['K_I_h2_pro']+10**-9))))

    I_nh3 = 1/(1+(c[Model.Species.index('S_nh3')] /
               (Model.Model_Parameters['K_I_nh3']+10**-9)))

    I_h2_oxidation=(1/(1+(c[Model.Species.index('S_h2')] /
                (Model.Model_Parameters['K_I_h2_ox']+10**-9))))

    I5 = (I_pH_aa * I_IN_lim)
    I6 = I5.copy()
    I7 = (I_pH_aa * I_IN_lim * I_h2_fa)
    I8 = (I_pH_aa * I_IN_lim * I_h2_c4)
    I9 = I8.copy()
    I10 = (I_pH_pro * I_IN_lim * I_h2_pro)
    I11 = (I_pH_ac * I_IN_lim * I_nh3)
    I12 = (I_pH_h2 * I_IN_lim)
    I13 = (I_pH_cap * I_IN_lim * I_h2_c4)
    I14 = (I_pH_bu * I_IN_lim * I_h2_c4)
    I15 = (I_pH_va * I_IN_lim * I_h2_c4)
    I16 = I_IN_lim * I_nh3*I_pH_aa*I_h2_oxidation

    v = np.zeros((len(Model.Reactions), 1))

    v[Model.Reactions.index(
        'TSS_Disintegration')] = Model.Model_Parameters["k_dis_TSS"]*c[Model.Species.index('TSS')]

    v[Model.Reactions.index(
        'TDS_Disintegration')] = Model.Model_Parameters["k_dis_TDS"]*c[Model.Species.index('TDS')]

    v[Model.Reactions.index('Hydrolysis carbohydrates')
      ] = Model.Model_Parameters['k_hyd_ch']*c[Model.Species.index('X_ch')]

    v[Model.Reactions.index('Hydrolysis proteins')
      ] = Model.Model_Parameters['k_hyd_pr']*c[Model.Species.index('X_pr')]

    v[Model.Reactions.index('Hydrolysis lipids')
      ] = Model.Model_Parameters['k_hyd_li']*c[Model.Species.index('X_li')]

    v[Model.Reactions.index('Uptake of sugars')] = Model.Model_Parameters['k_m_su']*c[Model.Species.index('S_su')] / \
        (Model.Model_Parameters['K_S_su']+c[Model.Species.index('S_su')]
         )*c[Model.Species.index('X_su')]*I5

    v[Model.Reactions.index('Uptake of amino acids')] = Model.Model_Parameters['k_m_aa']*c[Model.Species.index('S_aa')] / \
        (Model.Model_Parameters['K_S_aa']+c[Model.Species.index('S_aa')]
         )*c[Model.Species.index('X_aa')]*I6

    v[Model.Reactions.index('Uptake of LCFA')] = Model.Model_Parameters['k_m_fa']*c[Model.Species.index('S_fa')] / \
        (Model.Model_Parameters['K_S_fa'] +
         c[Model.Species.index('S_fa')])*c[Model.Species.index('X_fa')]*I7

    v[Model.Reactions.index('Uptake of acetate_et')] = Model.Model_Parameters['k_m_ac']*c[Model.Species.index('S_ac')]*c[Model.Species.index('S_et')] / \
        (Model.Model_Parameters['K_S_ac']*c[Model.Species.index('S_ac')]+Model.Model_Parameters['K_S_ac_et']*c[Model.Species.index('S_et')]+c[Model.Species.index('S_ac')]*c[Model.Species.index('S_et')]+10**-9
         )*c[Model.Species.index('X_ac_et')]*I11

    v[Model.Reactions.index('Uptake of acetate_lac')] = Model.Model_Parameters['k_m_ac']*c[Model.Species.index('S_ac')]*c[Model.Species.index('S_lac')] / \
        (Model.Model_Parameters['K_S_ac']*c[Model.Species.index('S_ac')]+Model.Model_Parameters['K_S_ac_lac']*c[Model.Species.index('S_lac')]+c[Model.Species.index('S_ac')]*c[Model.Species.index('S_lac')]+10**-9
         )*c[Model.Species.index('X_ac_lac')]*I11

    v[Model.Reactions.index('Uptake of propionate_et')] = Model.Model_Parameters['k_m_pro']*c[Model.Species.index('S_pro')]*c[Model.Species.index('S_et')] / \
        (Model.Model_Parameters['K_S_pro']*c[Model.Species.index('S_pro')]+Model.Model_Parameters['K_S_pro_et']*c[Model.Species.index('S_et')]+c[Model.Species.index('S_pro')]*c[Model.Species.index('S_et')]+10**-9
         )*c[Model.Species.index('X_chain_et')]*I10

    v[Model.Reactions.index('Uptake of propionate_lac')] = Model.Model_Parameters['k_m_pro']*c[Model.Species.index('S_pro')]*c[Model.Species.index('S_lac')] / \
        (Model.Model_Parameters['K_S_pro']*c[Model.Species.index('S_pro')]+Model.Model_Parameters['K_S_pro_lac']*c[Model.Species.index('S_lac')]+c[Model.Species.index('S_pro')]*c[Model.Species.index('S_lac')]+10**-9
         )*c[Model.Species.index('X_chain_lac')]*I10

    v[Model.Reactions.index('Uptake of butyrate_et')] = Model.Model_Parameters['k_m_bu']*c[Model.Species.index('S_bu')]*c[Model.Species.index('S_et')] / \
        (Model.Model_Parameters['K_S_bu']*c[Model.Species.index('S_bu')]+Model.Model_Parameters['K_S_bu_et']*c[Model.Species.index('S_et')]+c[Model.Species.index('S_bu')]*c[Model.Species.index('S_et')]+10**-9
         )*c[Model.Species.index('X_chain_et')]*I14

    v[Model.Reactions.index('Uptake of butyrate_lac')] = Model.Model_Parameters['k_m_bu']*c[Model.Species.index('S_bu')]*c[Model.Species.index('S_lac')] / \
        (Model.Model_Parameters['K_S_bu']*c[Model.Species.index('S_bu')]+Model.Model_Parameters['K_S_bu_lac']*c[Model.Species.index('S_lac')]+c[Model.Species.index('S_bu')]*c[Model.Species.index('S_lac')]+10**-9
         )*c[Model.Species.index('X_chain_lac')]*I14

    v[Model.Reactions.index('Uptake of valerate')] = Model.Model_Parameters['k_m_va']*c[Model.Species.index('S_va')] / \
        (Model.Model_Parameters['K_S_va']+c[Model.Species.index('S_va')]
         )*c[Model.Species.index('X_VFA_deg')]*I15

    v[Model.Reactions.index('Uptake of caproate')] = Model.Model_Parameters['k_m_cap']*c[Model.Species.index('S_cap')] / \
        (Model.Model_Parameters['K_S_cap']+c[Model.Species.index('S_cap')]
         )*c[Model.Species.index('X_VFA_deg')]*I13

    v[Model.Reactions.index('Methanogenessis from acetate and h2')] = Model.Model_Parameters['k_m_h2_Me_ac']*c[Model.Species.index('S_h2')]*c[Model.Species.index('S_ac')] / \
        (Model.Model_Parameters['K_S_h2_Me_ac']*c[Model.Species.index('S_h2')]+Model.Model_Parameters['K_S_ac_Me']*c[Model.Species.index(
            'S_ac')]+c[Model.Species.index('S_ac')]*c[Model.Species.index('S_h2')]+10**-9)*c[Model.Species.index('X_Me_ac')]*I12

    v[Model.Reactions.index('Methanogenessis from CO2 and h2')] = Model.Model_Parameters['k_m_h2_Me_CO2']*c[Model.Species.index('S_h2')]*c[Model.Species.index('S_co2')] / \
        (Model.Model_Parameters['K_S_h2_Me_CO2']*c[Model.Species.index('S_h2')]+Model.Model_Parameters['K_S_CO2_Me']*c[Model.Species.index(
            'S_co2')]+c[Model.Species.index('S_co2')]*c[Model.Species.index('S_h2')]+10**-9)*c[Model.Species.index('X_Me_CO2')]*I12


    v[Model.Reactions.index('Uptake of ethanol')] = Model.Model_Parameters['k_m_et']*c[Model.Species.index('S_et')] / \
        (Model.Model_Parameters['K_S_et']+c[Model.Species.index('S_et')]
         )*c[Model.Species.index("X_et")]*I16

    v[Model.Reactions.index('Uptake of lactate')] = Model.Model_Parameters['k_m_lac']*c[Model.Species.index('S_lac')] / \
        (Model.Model_Parameters['K_S_lac']+c[Model.Species.index('S_lac')]
         )*c[Model.Species.index('X_lac')]*I16

    v[Model.Reactions.index(
        'Decay of Xsu')] = Model.Model_Parameters['k_dec_X_su']*c[Model.Species.index('X_su')]

    v[Model.Reactions.index(
        'Decay of Xaa')] = Model.Model_Parameters['k_dec_X_aa']*c[Model.Species.index('X_aa')]

    v[Model.Reactions.index(
        'Decay of Xfa')] = Model.Model_Parameters['k_dec_X_fa']*c[Model.Species.index('X_fa')]

    v[Model.Reactions.index(
        'Decay of X_ac_et')] = Model.Model_Parameters['k_dec_X_ac']*c[Model.Species.index('X_ac_et')]

    v[Model.Reactions.index(
        'Decay of X_ac_lac')] = Model.Model_Parameters['k_dec_X_ac']*c[Model.Species.index('X_ac_lac')]

    v[Model.Reactions.index(
        'Decay of X_chain_et')] = Model.Model_Parameters['k_dec_X_chain_et']*c[Model.Species.index('X_chain_et')]

    v[Model.Reactions.index('Decay of X_chain_lac')
      ] = Model.Model_Parameters['k_dec_X_chain_lac']*c[Model.Species.index('X_chain_lac')]

    v[Model.Reactions.index(
        'Decay of X_VFA_deg')] = Model.Model_Parameters['k_dec_X_VFA_deg']*c[Model.Species.index('X_VFA_deg')]

    v[Model.Reactions.index(
        'Decay of X_Me_ac')] = Model.Model_Parameters['k_dec_X_Me_ac']*c[Model.Species.index('X_Me_ac')]

    v[Model.Reactions.index(
        'Decay of X_Me_CO2')] = Model.Model_Parameters['k_dec_X_Me_CO2']*c[Model.Species.index('X_Me_CO2')]

    v[Model.Reactions.index(
        'Decay of Xet')] = Model.Model_Parameters['k_dec_X_et']*c[Model.Species.index('X_et')]

    v[Model.Reactions.index(
        'Decay of Xlac')] = Model.Model_Parameters['k_dec_X_lac']*c[Model.Species.index('X_lac')]

    v[Model.Reactions.index('Acid Base Equilibrium (Va)')] = Model.Model_Parameters['k_A_B_va'] * \
        (c[Model.Species.index('S_va_ion')] * (Model.Model_Parameters['K_a_va'] + c[Model.Species.index('S_H_ion')]) -
         Model.Model_Parameters['K_a_va'] * c[Model.Species.index('S_va')])

    v[Model.Reactions.index('Acid Base Equilibrium (Bu)')] = Model.Model_Parameters['k_A_B_bu'] * \
        (c[Model.Species.index('S_bu_ion')] * (Model.Model_Parameters['K_a_bu'] + c[Model.Species.index('S_H_ion')]) -
         Model.Model_Parameters['K_a_bu'] * c[Model.Species.index('S_bu')])

    v[Model.Reactions.index('Acid Base Equilibrium (Pro)')] = Model.Model_Parameters['k_A_B_pro'] * \
        (c[Model.Species.index('S_pro_ion')] * (Model.Model_Parameters['K_a_pro'] + c[Model.Species.index('S_H_ion')]) -
         Model.Model_Parameters['K_a_pro'] * c[Model.Species.index('S_pro')])

    v[Model.Reactions.index('Acid Base Equilibrium (Cap)')] = Model.Model_Parameters['k_A_B_cap'] * \
        (c[Model.Species.index('S_cap_ion')] * (Model.Model_Parameters['K_a_cap'] + c[Model.Species.index('S_H_ion')]) -
         Model.Model_Parameters['K_a_cap'] * c[Model.Species.index('S_cap')])

    v[Model.Reactions.index('Acid Base Equilibrium (Lac)')] = Model.Model_Parameters['k_A_B_lac'] * \
        (c[Model.Species.index('S_lac_ion')] * (Model.Model_Parameters['K_a_lac'] + c[Model.Species.index('S_H_ion')]) -
         Model.Model_Parameters['K_a_lac'] * c[Model.Species.index('S_lac')])

    v[Model.Reactions.index('Acid Base Equilibrium (Ac)')] = Model.Model_Parameters['k_A_B_ac'] * \
        (c[Model.Species.index('S_ac_ion')] * (Model.Model_Parameters['K_a_ac'] + c[Model.Species.index('S_H_ion')]) -
         Model.Model_Parameters['K_a_ac'] * c[Model.Species.index('S_ac')])

    v[Model.Reactions.index('Acid Base Equilibrium (CO2)')] = Model.Model_Parameters['k_A_B_co2'] * \
        (c[Model.Species.index('S_hco3_ion')] * (Model.Model_Parameters['K_a_co2'] + c[Model.Species.index('S_H_ion')]) -
         Model.Model_Parameters['K_a_co2'] * c[Model.Species.index('S_IC')])

    v[Model.Reactions.index('Acid Base Equilibrium (In)')] = Model.Model_Parameters['k_A_B_IN'] * \
        (c[Model.Species.index('S_nh3')] * (Model.Model_Parameters['K_a_IN'] + c[Model.Species.index('S_H_ion')]) -
         Model.Model_Parameters['K_a_IN'] * c[Model.Species.index('S_IC')])

    p_gas_h2 = c[Model.Species.index('S_gas_h2')] * Model.Base_Parameters["R"] * \
        Model.Base_Parameters["T_op"] / 16
    p_gas_ch4 = c[Model.Species.index('S_gas_ch4')] * Model.Base_Parameters["R"] * \
        Model.Base_Parameters["T_op"] / 64
    p_gas_co2 = c[Model.Species.index('S_gas_co2')] * Model.Base_Parameters["R"] * \
        Model.Base_Parameters["T_op"]
    p_gas_h2o = 0.0313 * \
        np.exp(5290 *
               (1 / Model.Base_Parameters["T_base"] - 1 / Model.Base_Parameters["T_op"]))
    P_gas = p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o
    q_gas = max(
        0, (Model.Model_Parameters['k_p'] * (P_gas - Model.Base_Parameters['P_atm'])))
    v[Model.Reactions.index('Gas Transfer H2')] = Model.Model_Parameters['k_L_a'] * \
        (c[Model.Species.index('S_h2')] - 16 *
         Model.Model_Parameters['K_H_h2'] * p_gas_h2)

    v[Model.Reactions.index('Gas Transfer CH4')] = max(0,Model.Model_Parameters['k_L_a'] * \
        (c[Model.Species.index('S_ch4')] - 64 *
         Model.Model_Parameters['K_H_ch4'] * p_gas_ch4))
    v[Model.Reactions.index('Gas Transfer CO2')] = max(0,Model.Model_Parameters['k_L_a'] * \
        (c[Model.Species.index('S_co2')] -
         Model.Model_Parameters['K_H_co2'] * p_gas_co2))

    dCdt = np.matmul(Model.S, v)

    phi = c[Model.Species.index('S_cation')]+c[Model.Species.index('S_nh4_ion')]-c[Model.Species.index('S_hco3_ion')]-(c[Model.Species.index('S_lac_ion')] / 88) - (c[Model.Species.index('S_ac_ion')] / 64) - (c[Model.Species.index('S_pro_ion')] /
                                                                                                                                                                     112) - (c[Model.Species.index('S_bu_ion')] / 160)-(c[Model.Species.index('S_cap_ion')] / 230) - (c[Model.Species.index('S_va_ion')] / 208) - c[Model.Species.index('S_anion')]
    if 'S_H_ion' in Model.Control_States.keys():
        c[Model.Species.index('S_H_ion')]=Model.Control_States['S_H_ion']
    else:
        c[Model.Species.index('S_H_ion')] = (-1 * phi / 2) + \
        (0.5 * np.sqrt(phi**2 + 4 * Model.Model_Parameters['K_w']))

    dCdt[0: Model.Species.__len__()-3] = dCdt[0: Model.Species.__len__()-3]+Model.Base_Parameters['q_in'] / \
        Model.Base_Parameters["V_liq"] * \
        (Model.Inlet_Conditions[0: Model.Species.__len__(
        )-3]-c[0: Model.Species.__len__()-3].reshape(-1, 1))

    dCdt[Model.Species.__len__()-3:] = dCdt[Model.Species.__len__()-3:]+q_gas/Model.Base_Parameters["V_gas"] * \
        (Model.Inlet_Conditions[Model.Species.__len__() -
         3:]-c[Model.Species.__len__()-3:].reshape(-1, 1))

    dCdt[[Model.Species.index('S_H_ion'), Model.Species.index(
        'S_co2'), Model.Species.index('S_nh4_ion')], 0] = 0

    if Model.Switch == "DAE":
        dCdt[Model.Species.index('S_h2')] = 0

        dCdt[Model.Species.index('S_va_ion'):Model.Species.index('S_co2')] = 0

        dCdt[Model.Species.index('S_nh3')] = 0
    
        c[Model.Species.index('S_va_ion')]=Model.Model_Parameters['K_a_va']/(Model.Model_Parameters['K_a_va']+c[Model.Species.index('S_H_ion')])*c[Model.Species.index('S_va')]
    
        c[Model.Species.index('S_bu_ion')]=Model.Model_Parameters['K_a_bu']/(Model.Model_Parameters['K_a_bu']+c[Model.Species.index('S_H_ion')])*c[Model.Species.index('S_bu')]
    
        c[Model.Species.index('S_pro_ion')]=Model.Model_Parameters['K_a_pro']/(Model.Model_Parameters['K_a_pro']+c[Model.Species.index('S_H_ion')])*c[Model.Species.index('S_pro')]
    
        c[Model.Species.index('S_cap_ion')]=Model.Model_Parameters['K_a_cap']/(Model.Model_Parameters['K_a_cap']+c[Model.Species.index('S_H_ion')])*c[Model.Species.index('S_cap')]
    
        c[Model.Species.index('S_ac_ion')]=Model.Model_Parameters['K_a_ac']/(Model.Model_Parameters['K_a_ac']+c[Model.Species.index('S_H_ion')])*c[Model.Species.index('S_ac')]
        
        c[Model.Species.index('S_lac_ion')]=Model.Model_Parameters['K_a_lac']/(Model.Model_Parameters['K_a_lac']+c[Model.Species.index('S_H_ion')])*c[Model.Species.index('S_lac')]    

    if Model.Control_States.keys():
        for state in Model.Control_States.keys():
            c[Model.Species.index(state)]=Model.Control_States[state]
            dCdt[Model.Species.index(state)]=0
    Model.Info["Fluxes"].append(v)
    return dCdt[:, 0]


if __name__ == "__main__":
   # Report = os.path.join(Main_Dir, "..", "Reports",
                         # "ADM_From_Alignment_JSON_Output.json")
   # with open(Report, 'r') as j:
        #Report = json.load(j)

   #adm1 = Model(Model_Parameters, Base_Parameters, Initial_Conditions,
   #             Inlet_Conditions, Reactions, Species, ADM1_ODE_Sys, Build_ADM1_Stoiciometric_Matrix, Name="ADM1", Switch="DAE")
   #Sol = adm1.Solve_Model(
   #    (0, 20), adm1.Initial_Conditions[:, 0], np.linspace(0, 20, 100))

    # adm1.Dash_App(Sol)
    with open("../Optimization/Tuned_Params.json", 'r') as f:
        mp=json.load(f)
    with open('/Users/parsaghadermarzi/Desktop/ADM_Mapping.json', 'r') as f:
        MR=json.load(f)
    with open('/Users/parsaghadermarzi/Desktop/ADToolbox/Database/Modified_ADM/Modified_ADM_Base_Parameters.json', 'r') as f:
        BP=json.load(f)
    with open('/Users/parsaghadermarzi/Desktop/ADToolbox/Database/Modified_ADM/Modified_ADM_Initial_Conditions.json', 'r') as f:
        IC=json.load(f)
    with open('/Users/parsaghadermarzi/Desktop/ADToolbox/Database/Modified_ADM/Modified_ADM_Inlet_Conditions.json', 'r') as f:
        InC=json.load(f)
    with open('/Users/parsaghadermarzi/Desktop/ADToolbox/Database/Modified_ADM/Modified_ADM_Reactions.json', 'r') as f:
        r=json.load(f)
    with open('/Users/parsaghadermarzi/Desktop/ADToolbox/Database/Modified_ADM/Modified_ADM_Species.json', 'r') as f:
        s=json.load(f)
    
    mod_adm1 = Model(mp, BP, IC, InC, r,
                     s, Modified_ADM1_ODE_Sys, Build_Modified_ADM1_Stoiciometric_Matrix,Control_States={"S_H_ion":0},Name="Modified_ADM1", Switch="DAE",Metagenome_Report=MR)
    Sol_mod_adm1 = mod_adm1.Solve_Model(mod_adm1.Initial_Conditions[:, 0], np.linspace(0,30, 10000))
    mod_adm1.Dash_App(Sol_mod_adm1)