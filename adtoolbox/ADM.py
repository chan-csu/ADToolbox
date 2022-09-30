from typing import Callable
import plotly
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
import Configs

### Note ###
# The following code is a modified version of the code from the PyADM1 package
# It is extensively based on PyADM1, and I would like to thank the author of PyADM1 for this work
# As far as implementing the original ADM1 in python goes, I still consider this as a modification from
# the PyADM1 code.

# ----------


RT = Reaction_Toolkit(reaction_db=Configs.Reaction_Toolkit().reaction_db)

class _Fake_Sol:
    def __init__(self, y,t):
        self.y = y
        self.t=t



class Model:

    """General Class For a Model
    Args:
        Model Parameters (dict): a dictionary which contains model parameters
        base_parameters (dict): a dictionary which contains base paramters
        Initial Conditions (dict): a dictionary containing inlet conditions for all species
        Inlet Conditions (dict): a dictionary containing inlet conditions for all species
        reactions (list): a list containing all types of reactions
        species (list): a list containing all species
        ODE_System (function): a function which solves inital value ODEs
        Build Stoiciometric Matrix (function): a function which outputs a matrix for of all reactions
        
    Returns:
        adtoolbox.ADM.Model: returns a model instance for downstream purposes.
    """
    


    def __init__(self, model_parameters: dict, base_parameters: dict, initial_conditions: dict, inlet_conditions:dict, reactions: list, species: list, ode_system:Callable, build_stoichiometric_matrix:Callable,control_state:dict={},metagenome_report=None, name="ADM1", switch="DAE",simulation_time=30):
        self.model_parameters = model_parameters
        self.base_parameters = base_parameters
        

        for items in control_state.keys():
            initial_conditions[items]
            initial_conditions[items]=control_state[items]
        self.control_state=control_state


        self.inlet_conditions = np.array(
            [inlet_conditions[i+"_in"] for i in species])[:, np.newaxis]
        self.reactions = reactions
        self.species = species
        self.initial_conditions = np.array(
            [initial_conditions[i] for i in species])[:, np.newaxis]
        self._ic=initial_conditions
        self._inc=inlet_conditions
        self.switch = switch
        self.name = name
        self.metagenome_report = metagenome_report
        self.S = build_stoichiometric_matrix(
            base_parameters, model_parameters, reactions, species)
        self.ode_system = ode_system
        self.sim_time=simulation_time

    # def Build_COBRA_Model(self, save_model=True):
    #     model = cobra.Model(self.Name)

    #     for i in range(len(self.reactions)):
    #         Temp_Rxn = cobra.Reaction(
    #             self.reactions[i].replace(" ", "_"), name=self.reactions[i].replace(" ", "_"))
    #         Temp_Mets = np.where(self.S[:, i] != 0)
    #         Met_Dict = {}
    #         for met in Temp_Mets[0]:
    #             Metabolite = cobra.Metabolite(self.species[met].replace(" ", "_"),
    #                                           name=self.species[met].replace(" ", "_"), compartment="Model")
    #             Met_Dict[Metabolite] = self.S[met, i]
    #         Temp_Rxn.add_metabolites(Met_Dict)
    #         # print(Temp_Rxn.reaction)
    #         model.add_reaction(Temp_Rxn)

        # if save_model:
        #     cobra.io.save_json_model(model, self.Name+'.json')

    def solve_model(self, y0: np.ndarray, t_eval: np.ndarray, method="LSODA", switch="DAF")->scipy.integrate._ivp.ivp.OdeResult:
        """
        Function to solve the model. 
        Examples:
            >>> from ADM import Model
            >>> import numpy as np
            >>> reactions=['rxn1','rxn2']
            >>> species=['a','b','c']
            >>> Initial_Conditions={'a':.001,'b':.002,'c':.003}
            >>> inlet_conditions={'a_in':.001,'b_in':.002,'c_in':.003}
            >>> model_parameters={'k1':0.001,'k2':0.002}
            >>> base_parameters={'T':0.1}
            >>> def Build_Stoiciometric_Matrix(base_parameters,model_parameters,reactions,species):
            ...    S = np.zeros((len(species), len(reactions)))
            ...    S[[0,1],0]=[-1,0.001]
            ...    S[[1,2],1]=[-5,1]
            ...    return S
            >>> def ODE_System(t,c,Model1):
            ...    v = np.zeros((len(Model1.reactions), 1))
            ...    v[0]=Model1.model_parameters['k1']*c[0]*Model1.base_parameters['T']/1000
            ...    v[1]=Model1.model_parameters['k2']*c[1]/1000
            ...    dCdt=np.matmul(Model1.S,v)
            ...    return dCdt[:, 0]
            >>> m= Model(model_parameters,base_parameters,Initial_Conditions,inlet_conditions,reactions,species,ODE_System,Build_Stoiciometric_Matrix)
            >>> m.Solve_Model((0,.1),m.Initial_Conditions[:, 0],np.linspace(0,0.1,10),method='RK45')['status']==0
            True
        
        Args:
            t_span (tuple): The range the time will include
            y0 (np.ndarray): The initial value at time = 0
            T_eval (np.ndarray): The time steps that will be evaluated
        
        Returns:
            scipy.integrate._ivp.ivp.OdeResult: Returns the results of the simulation being run and gives optimized paramters.
        """
        self.info={"Fluxes":[]}
        try:
            c = scipy.integrate.solve_ivp(
                self.ode_system, (0,self.sim_time), y0, t_eval=t_eval, method=method, args=[self])
        
        except Exception as e:
            print("Could not solve model, setting C to a very large value")
            c=_Fake_Sol(np.ones((y0.shape[0],t_eval.__len__()))*1e10,t_eval)
       
        return c

    
       #C = scipy.integrate.solve_ivp(
       #        self.ODE_System, t_span, y0, t_eval=T_eval, method=method, args=[self])
       #
       #return C

    def plot(self, Sol: scipy.integrate._ivp.ivp.OdeResult, type: str = "Line")-> None:
        """ A function which returns a plot of the solution from the ODE
        """
        solution = {
            't': Sol.t,
        }
        for i in range(len(self.species)):
            solution[self.species[i]] = Sol.y[i, :]
        sol_df = pd.DataFrame(solution)

        if type == "Line":
            fig = px.line(sol_df, x="t", y=sol_df.columns,
                          title="Concentration of species",text=sol_df.columns,textposition="bottom center")
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

        elif type == "Sankey":
            ### Maybe add a sankey plot here later
            pass

    def dash_app(self, sol: scipy.integrate._ivp.ivp.OdeResult, type: str = "Line")->None:
        """A method that creates the dash web app"""

        with open(os.path.join(PKG_DATA,"Modified_ADM_Map.json"),'rb') as f:
            escher_map=json.load(f)
        with open(os.path.join(PKG_DATA,"Modified_ADM_Model.json"),'rb') as f:
            cobra_model=json.load(f)


        app = Dash(__name__)
        colors = {
            'background': '#659dbd',
            'text': '#3e4444'
        }
        

        solution = {
            't': sol.t,
        }
        for i in range(len(self.species)):
            solution[self.species[i]] = sol.y[i, :]
        sol_df = pd.DataFrame(solution)

        if type == "Line":
            fig = px.line(sol_df, x="t", y=sol_df.columns,
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

        with_report=[html.Div(style={'backgroundColor': 'rgb(200, 200, 200)', 'height':'300px', "margin-top":'-50px'}, children=[html.H1(
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

            html.H2(f"{self.name} Concentration Plot", style={
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

            html.H3("base_parameters", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }),
            
            dash_table.DataTable(
            id='base_parameters',
            columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self.base_parameters.keys())],
            data=pd.DataFrame(self.base_parameters,index=[0]).to_dict('records'),
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
            html.H3("model_parameters", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }),



            dash_table.DataTable(
            id='model_parameters',
            columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self.model_parameters.keys())],
            data=pd.DataFrame(self.model_parameters,index=[0]).to_dict('records'),
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
            columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self._ic.keys())],
            data=pd.DataFrame(self._ic,index=[0]).to_dict('records'),
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
            html.H3("inlet_conditions", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }),
            dash_table.DataTable(
            id='inlet_conditions',
            columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self._inc.keys())],
            data=pd.DataFrame(self._inc,index=[0]).to_dict('records'),
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

                        html.H2("Escher Map", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }) ,
            
            dcc.Dropdown(["Show Map","Hide Map"],
                         self.reactions[0], style={"width": "300px","font-size":25}, id="Drop_Down_Escher"),
            html.Br(),
            html.Div(children=None,id="Escher_"),
            html.Div(children=None,id="Escher"),

            html.H2("Microbial Association", style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'font-size': '20px',
                    "BackgroundColor": "#fbeec1",
                    "font-family": "Trebuchet MS",
                    }) ,
            

            dcc.Dropdown(self.reactions,
                         self.reactions[0], style={"width": "300px","font-size":25}, id="Drop_Down") ,

            dcc.Graph(figure=None, id="Annotation_Graph", style={
            "height": "650px"}),


            
         ]
        without_report=with_report[:-3]
        
        
        
        app.layout = html.Div(with_report) if self.metagenome_report else html.Div(without_report)





        @app.callback(Output(component_id="Escher_", component_property='children'), Input(component_id="Drop_Down_Escher", component_property='value'))
        def escher_wrapper(drop_down_escher):
            if drop_down_escher=="Show Map":
                Labels={}
                for i in range(0,len(sol.t),int(len(sol.t)/20)):
                    Labels[i]={'label':str(int(sol.t[i])),'style':{'color': '#77b0b1'}}
                Labels[len(sol.t)-1]=self.sim_time
                return [html.Br(),html.H2("Time (Day)",style={'textAlign': 'center'}),dcc.Slider(0, len(sol.t),len(sol.t)/20,value=0,id="Escher_Slider",marks=Labels),html.Br()]

        @app.callback(Output(component_id="Escher", component_property='children'), Input(component_id="Drop_Down_Escher", component_property='value'),
        Input(component_id="Escher_Slider", component_property='value'))        
        def draw_escher(drop_down_escher,escher_slider):
            rxn_data={}
            for ind,i in enumerate(self.reactions):
                rxn_data[i.replace(" ","_")]=self.info["Fluxes"][int(escher_slider)][:][ind]
            
            if drop_down_escher=="Show Map":
                return [dash_escher.DashEscher(mapData=escher_map,modelData=cobra_model,
            options={
             'reaction_data':rxn_data
            }
            ,height='1000px',
        width='100%')
             ]

        @app.callback(Output(component_id="Annotation_Graph", component_property='figure'), Input(component_id='Drop_Down', component_property='value'))
        def update_graph_fig(input_value):
            reactions = self.reactions
            id_list = []
            label_list = []
            parents_list = []
            label_list.append(input_value)
            id_list.append("")
            parents_list.append(id_list[0])

            genome_list = self.metagenome_report[input_value]
            for item in genome_list.keys():
                label_list.append(item)
                id_list.append(item)
                parents_list.append(input_value)

            for item in genome_list.keys():
                C = 0
                if self.metagenome_report[input_value][item].__len__() > 0:
                    for ecs in self.metagenome_report[input_value][item].keys():
                        label_list.append(ecs)
                        id_list.append(item+" : "+ecs)
                        parents_list.append(item)
                        label_list.append("<br>".join(
                            self.metagenome_report[input_value][item][ecs]))
                        id_list.append(item+" : "+ecs+" : " + str(C))
                        parents_list.append(item+" : "+ecs)
                        C += 1

                else:
                    label_list.append(
                        f"No Relevant EC Numbers were found for {item.split(';')[-1]} !")
                    parents_list.append(item)

            fig = go.Figure(go.Treemap(
                ids=id_list,
                labels=label_list,
                parents=parents_list,
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
                    Input(component_id='base_parameters', component_property='data'),
                    Input(component_id='model_parameters', component_property='data'),
                    Input(component_id='Initial_Conditions', component_property='data'),
                    Input(component_id='inlet_conditions', component_property='data'))
        def update_graph_fig(base_parameters: dict, model_parameters:dict, initial_conditions: dict, inlet_conditions: dict,prevent_initial_call=True)->plotly.graph_objects.Figure:
            """A functions that imports base_parameters, model_parameters, Initial_Conditions, inlet_conditions into the plot.
            """
            if len(self.control_state.keys()):
                for i in self.control_state.keys():
                    self.control_state[i]=initial_conditions[0][i]

            self.base_parameters = base_parameters[0]
            self.model_parameters = model_parameters[0]
            self.initial_conditions = np.array(
            [initial_conditions[0][i] for i in self.species])[:, np.newaxis]
            self.inlet_conditions = np.array(
            [inlet_conditions[0][i+"_in"] for i in self.species])[:, np.newaxis]
            update_sol = self.solve_model(
                         self.initial_conditions[:, 0], np.linspace(0, 30, 10000))


            solution = {
                    't': update_sol.t,
                        }
            for i in range(len(self.species)):
                solution[self.species[i]] = update_sol.y[i, :]
            sol_df = pd.DataFrame(solution)

            fig = px.line(sol_df, x="t", y=sol_df.columns,
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

    def sol_to_df(self,Address: str)->None:
        """Converts the matrix to a pandas data frame then to a csv"""

        df = pd.DataFrame(self.S, columns=self.reactions, index=self.species)
        df.to_csv(Address, header=True,
                  index=True)

    def csv_report(self,sol: scipy.integrate._ivp.ivp.OdeResult ,address: str)->None:
        """Converts the results to a pandas data frame then to a csv"""
        df = pd.DataFrame(sol.y, columns=sol.t, index=self.species)
        df.to_csv(os.path.join(address,self.name+"_Report.csv"), header=True,
                  index=True)


def build_adm1_stoiciometric_Matrix(base_parameters: dict, model_parameters: dict, reactons: list, species:list)-> np.ndarray:
    """This function builds the stoichiometric matrix for the original ADM Model.
    No testing is done.
    """
    if type(base_parameters)!= dict and type(model_parameters)!= dict:
        raise TypeError("base_parameters and model_parameters need to be dictionary")
    if type(reactons)!= list and type(species)!= list:
        raise TypeError("reactions and species must be list")


    S = np.zeros((len(species), len(reactons)))
    S[0, [1, 3, 4]] = [1, (1-model_parameters["f_fa_li"]), - 1]
    S[1, [2, 5]] = [1, -1]
    S[2, [3, 6]] = [(model_parameters["f_fa_li"]), - 1]
    S[3, [5, 7]] = [(1-model_parameters['Y_aa']) *
                    model_parameters['f_va_aa'], - 1]
    S[4, [4, 5, 8]] = [(1-model_parameters['Y_su'])*model_parameters['f_bu_su'],
                       (1-model_parameters['Y_aa'])*model_parameters["f_bu_aa"], - 1]
    S[5, [4, 5, 7, 9]] = [(1-model_parameters["Y_su"])*model_parameters['f_pro_su'],
                          (1-model_parameters['Y_aa'])*model_parameters["f_pro_aa"], (1 - model_parameters['Y_c4'])*0.54, -1]
    S[6, [4, 5, 6, 7, 8, 9, 10]] = [(1-model_parameters['Y_su'])*model_parameters['f_ac_su'],
                                    (1-model_parameters['Y_aa']) *
                                    model_parameters['f_ac_aa'],
                                    (1-model_parameters['Y_fa'])*0.7,
                                    (1-model_parameters['Y_c4'])*0.31,
                                    (1-model_parameters['Y_c4'])*0.8,
                                    (1-model_parameters['Y_pro'])*0.57,
                                    -1]
    S[7, [4, 5, 6, 7, 8, 9, 11, 25]] = [(1-model_parameters['Y_su'])*model_parameters['f_h2_su'],
                                        (1-model_parameters['Y_aa']) *
                                        model_parameters['f_h2_aa'],
                                        (1-model_parameters['Y_fa'])*0.3,
                                        (1-model_parameters['Y_c4'])*0.15,
                                        (1-model_parameters['Y_c4'])*0.2,
                                        (1-model_parameters['Y_pro'])*0.43,
                                        -1,
                                        -1]
    S[8, [10, 11, 26]] = [(1-model_parameters['Y_ac']),
                          (1-model_parameters['Y_h2']),
                          -1]
    s_1 = (-1 * model_parameters['C_xc'] + model_parameters['f_sI_xc'] * model_parameters['C_sI'] + model_parameters['f_ch_xc'] * model_parameters['C_ch'] +
           model_parameters['f_pr_xc'] * model_parameters['C_pr'] + model_parameters['f_li_xc'] * model_parameters['C_li'] + model_parameters['f_xI_xc'] * model_parameters['C_xI'])
    s_2 = (-1 * model_parameters['C_ch'] + model_parameters['C_su'])
    s_3 = (-1 * model_parameters['C_pr'] + model_parameters['C_aa'])
    s_4 = (-1 * model_parameters['C_li'] + (1 - model_parameters['f_fa_li']) *
           model_parameters['C_su'] + model_parameters['f_fa_li'] * model_parameters['C_fa'])
    s_5 = (-1 * model_parameters['C_su'] + (1 - model_parameters['Y_su']) * (model_parameters['f_bu_su'] * model_parameters['C_bu'] + model_parameters['f_pro_su']
                                                                             * model_parameters['C_pro'] + model_parameters['f_ac_su'] * model_parameters['C_ac']) + model_parameters['Y_su'] * model_parameters['C_bac'])
    s_6 = (-1 * model_parameters['C_aa'] + (1 - model_parameters['Y_aa']) * (model_parameters['f_va_aa'] * model_parameters['C_va'] + model_parameters['f_bu_aa'] * model_parameters['C_bu'] +
                                                                             model_parameters['f_pro_aa'] * model_parameters['C_pro'] + model_parameters['f_ac_aa'] * model_parameters['C_ac']) + model_parameters['Y_aa'] * model_parameters['C_bac'])
    s_7 = (-1 * model_parameters['C_fa'] + (1 - model_parameters['Y_fa']) * 0.7 *
           model_parameters['C_ac'] + model_parameters['Y_fa'] * model_parameters['C_bac'])
    s_8 = (-1 * model_parameters['C_va'] + (1 - model_parameters['Y_c4']) * 0.54 * model_parameters['C_pro'] + (
        1 - model_parameters['Y_c4']) * 0.31 * model_parameters['C_ac'] + model_parameters['Y_c4'] * model_parameters['C_bac'])
    s_9 = (-1 * model_parameters['C_bu'] + (1 - model_parameters['Y_c4']) * 0.8 *
           model_parameters['C_ac'] + model_parameters['Y_c4'] * model_parameters['C_bac'])
    s_10 = (-1 * model_parameters['C_pro'] + (1 - model_parameters['Y_pro']) * 0.57 *
            model_parameters['C_ac'] + model_parameters['Y_pro'] * model_parameters['C_bac'])
    s_11 = (-1 * model_parameters['C_ac'] + (1 - model_parameters['Y_ac']) *
            model_parameters['C_ch4'] + model_parameters['Y_ac'] * model_parameters['C_bac'])
    s_12 = ((1 - model_parameters['Y_h2']) * model_parameters['C_ch4'] +
            model_parameters['Y_h2'] * model_parameters['C_bac'])
    s_13 = (-1 * model_parameters['C_bac'] + model_parameters['C_xc'])
    S[9, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 27]] = [-s_1, -s_2, -s_3, -s_4, -
                                                                                    s_5, -s_6, -s_7, -s_8, -s_9, -s_10, -s_11, -s_12, -s_13, -s_13, -s_13, -s_13, -s_13, -s_13, -s_13, -1]
    S[10, [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]] = [model_parameters['N_xc']-model_parameters['f_xI_xc']*model_parameters['N_I']-model_parameters['f_sI_xc']*model_parameters['N_I']-model_parameters['f_pr_xc']*model_parameters['N_aa'],
                                                                        -model_parameters['Y_su']*model_parameters['N_bac'],
                                                                        model_parameters['N_aa']-model_parameters['Y_aa'] *
                                                                        model_parameters['N_bac'],
                                                                        -model_parameters['Y_fa']*model_parameters['N_bac'],
                                                                        -model_parameters['Y_c4']*model_parameters['N_bac'],
                                                                        -model_parameters['Y_c4']*model_parameters['N_bac'],
                                                                        -model_parameters['Y_pro']*model_parameters['N_bac'],
                                                                        -model_parameters['Y_ac']*model_parameters['N_bac'],
                                                                        -model_parameters['Y_h2']*model_parameters['N_bac'],
                                                                        model_parameters['N_bac'] -
                                                                        model_parameters['N_xc'],
                                                                        model_parameters['N_bac'] -
                                                                        model_parameters['N_xc'],
                                                                        model_parameters['N_bac'] -
                                                                        model_parameters['N_xc'],
                                                                        model_parameters['N_bac'] -
                                                                        model_parameters['N_xc'],
                                                                        model_parameters['N_bac'] -
                                                                        model_parameters['N_xc'],
                                                                        model_parameters['N_bac'] -
                                                                        model_parameters['N_xc'],
                                                                        model_parameters['N_bac']-model_parameters['N_xc']]
    S[11, 0] = model_parameters['f_sI_xc']
    S[12, [0, 12, 13, 14, 15, 16, 17, 18]] = [-1, 1, 1, 1, 1, 1, 1, 1]
    S[13, [0, 1]] = [model_parameters['f_ch_xc'], -1]
    S[14, [0, 2]] = [model_parameters['f_pr_xc'], -1]
    S[15, [0, 3]] = [model_parameters['f_li_xc'], -1]
    S[16, [4, 12]] = [model_parameters['Y_su'], -1]
    S[17, [5, 13]] = [model_parameters['Y_aa'], -1]
    S[18, [6, 14]] = [model_parameters['Y_fa'], -1]
    S[19, [7, 8, 15]] = [model_parameters['Y_c4'], model_parameters['Y_c4'], -1]
    S[20, [9, 16]] = [model_parameters['Y_pro'], -1]
    S[21, [10, 17]] = [model_parameters['Y_ac'], -1]
    S[22, [11, 18]] = [model_parameters['Y_h2'], -1]
    S[23, 0] = model_parameters['f_xI_xc']
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
    S[35, 25] = base_parameters['V_liq']/base_parameters['V_gas']
    S[36, 26] = base_parameters['V_liq']/base_parameters['V_gas']
    S[37, 27] = base_parameters['V_liq']/base_parameters['V_gas']
    return S


def adm1_ode_sys(t: float, c: np.ndarray, adm1_instance:Model)-> np.ndarray:
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
    I_pH_aa = (adm1_instance.model_parameters["K_pH_aa"] ** adm1_instance.model_parameters['nn_aa'])/(np.power(
        c[26], adm1_instance.model_parameters['nn_aa']) + np.power(adm1_instance.model_parameters["K_pH_aa"], adm1_instance.model_parameters['nn_aa']))
    I_pH_ac = (adm1_instance.model_parameters['K_pH_ac'] ** adm1_instance.model_parameters["n_ac"])/(
        c[26] ** adm1_instance.model_parameters['n_ac'] + adm1_instance.model_parameters['K_pH_ac'] ** adm1_instance.model_parameters['n_ac'])
    I_pH_h2 = (adm1_instance.model_parameters['K_pH_h2']**adm1_instance.model_parameters['n_h2'])/(
        c[26] ** adm1_instance.model_parameters['n_h2'] + adm1_instance.model_parameters['K_pH_h2']**adm1_instance.model_parameters['n_h2'])
    I_IN_lim = 1 / (1+(adm1_instance.model_parameters['K_S_IN'] / c[10]))
    I_h2_fa = 1 / (1+(c[7] / adm1_instance.model_parameters['K_I_h2_fa']))
    I_h2_c4 = 1 / (1+(c[7]/adm1_instance.model_parameters['K_I_h2_c4']))
    I_h2_pro = (1/(1+(c[7]/adm1_instance.model_parameters['K_I_h2_pro'])))
    I_nh3 = 1/(1+(c[33]/adm1_instance.model_parameters['K_I_nh3']))
    I5 = (I_pH_aa * I_IN_lim)
    I6 = np.copy(I5)
    I7 = (I_pH_aa * I_IN_lim * I_h2_fa)
    I8 = (I_pH_aa * I_IN_lim * I_h2_c4)
    I9 = np.copy(I8)
    I10 = (I_pH_aa * I_IN_lim * I_h2_pro)
    I11 = (I_pH_ac * I_IN_lim * I_nh3)
    I12 = (I_pH_h2 * I_IN_lim)
    v = np.zeros((len(adm1_instance.reactions), 1))
    v[0] = adm1_instance.model_parameters["k_dis"]*c[12]
    v[1] = adm1_instance.model_parameters['k_hyd_ch']*c[13]
    v[2] = adm1_instance.model_parameters['k_hyd_pr']*c[14]
    v[3] = adm1_instance.model_parameters['k_hyd_li']*c[15]
    v[4] = adm1_instance.model_parameters['k_m_su']*c[0] / \
        (adm1_instance.model_parameters['K_S_su']+c[0])*c[16]*I5
    v[5] = adm1_instance.model_parameters['k_m_aa']*c[1] / \
        (adm1_instance.model_parameters['K_S_aa']+c[1])*c[17]*I6
    v[6] = adm1_instance.model_parameters['k_m_fa']*c[2] / \
        (adm1_instance.model_parameters['K_S_fa']+c[2])*c[18]*I7
    v[7] = adm1_instance.model_parameters['k_m_c4']*c[3] / \
        (adm1_instance.model_parameters['K_S_c4']+c[3]) * \
        c[19]*c[3]/(c[3]+c[4]+10 ** (-6))*I8
    v[8] = adm1_instance.model_parameters['k_m_c4']*c[4] / \
        (adm1_instance.model_parameters['K_S_c4']+c[4]) * \
        c[19]*c[4]/(c[4]+c[3]+10 ** (-6))*I9
    v[9] = adm1_instance.model_parameters['k_m_pr']*c[5] / \
        (adm1_instance.model_parameters['K_S_pro']+c[5])*c[20]*I10
    v[10] = adm1_instance.model_parameters['k_m_ac']*c[6] / \
        (adm1_instance.model_parameters['K_S_ac']+c[6])*c[21]*I11
    v[11] = adm1_instance.model_parameters['k_m_h2']*c[7] / \
        (adm1_instance.model_parameters['K_S_h2']+c[7])*c[22]*I12
    v[12] = adm1_instance.model_parameters['k_dec_X_su']*c[16]
    v[13] = adm1_instance.model_parameters['k_dec_X_aa']*c[17]
    v[14] = adm1_instance.model_parameters['k_dec_X_fa']*c[18]
    v[15] = adm1_instance.model_parameters['k_dec_X_c4']*c[19]
    v[16] = adm1_instance.model_parameters['k_dec_X_pro']*c[20]
    v[17] = adm1_instance.model_parameters['k_dec_X_ac']*c[21]
    v[18] = adm1_instance.model_parameters['k_dec_X_h2']*c[22]
    v[19] = adm1_instance.model_parameters['k_A_B_va'] * \
        (c[27] * (adm1_instance.model_parameters['K_a_va'] + c[26]) -
         adm1_instance.model_parameters['K_a_va'] * c[3])
    v[20] = adm1_instance.model_parameters['k_A_B_bu'] * \
        (c[28] * (adm1_instance.model_parameters['K_a_bu'] + c[26]) -
         adm1_instance.model_parameters['K_a_bu'] * c[4])
    v[21] = adm1_instance.model_parameters['k_A_B_pro'] * \
        (c[29] * (adm1_instance.model_parameters['K_a_pro'] + c[26]) -
         adm1_instance.model_parameters['K_a_pro'] * c[5])
    v[22] = adm1_instance.model_parameters['k_A_B_ac'] * \
        (c[30] * (adm1_instance.model_parameters['K_a_ac'] + c[26]) -
         adm1_instance.model_parameters['K_a_ac'] * c[6])
    v[23] = adm1_instance.model_parameters['k_A_B_co2'] * \
        (c[31] * (adm1_instance.model_parameters['K_a_co2'] + c[26]) -
         adm1_instance.model_parameters['K_a_co2'] * c[9])
    v[24] = adm1_instance.model_parameters['k_A_B_IN'] * \
        (c[33] * (adm1_instance.model_parameters['K_a_IN'] + c[26]) -
         adm1_instance.model_parameters['K_a_IN'] * c[10])
    p_gas_h2 = c[35] * adm1_instance.base_parameters["R"] * \
        adm1_instance.base_parameters["T_op"] / 16
    p_gas_ch4 = c[36] * adm1_instance.base_parameters["R"] * \
        adm1_instance.base_parameters["T_op"] / 64
    p_gas_co2 = c[37] * adm1_instance.base_parameters["R"] * \
        adm1_instance.base_parameters["T_op"]
    p_gas_h2o = 0.0313 * \
        np.exp(5290 *
               (1 / adm1_instance.base_parameters["T_base"] - 1 / adm1_instance.base_parameters["T_op"]))
    P_gas = p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o
    q_gas = max(
        0, (adm1_instance.model_parameters['k_p'] * (P_gas - adm1_instance.base_parameters['P_atm'])))
    v[25] = adm1_instance.model_parameters['k_L_a'] * \
        (c[7] - 16 * adm1_instance.model_parameters['K_H_h2'] * p_gas_h2)
    v[26] = adm1_instance.model_parameters['k_L_a'] * \
        (c[8] - 64 * adm1_instance.model_parameters['K_H_ch4'] * p_gas_ch4)
    v[27] = adm1_instance.model_parameters['k_L_a'] * \
        (c[32] - adm1_instance.model_parameters['K_H_co2'] * p_gas_co2)
    dCdt = np.matmul(adm1_instance.S, v)
    phi = c[24]+c[34]-c[31] - (c[30] / 64) - (c[29] /
                                              112) - (c[28] / 160) - (c[27] / 208) - c[25]
    c[26] = (-1 * phi / 2) + \
        (0.5 * np.sqrt(phi**2 + 4 * adm1_instance.model_parameters['K_w']))
    dCdt[0: 35] = dCdt[0: 35]+adm1_instance.base_parameters['q_in'] / \
        adm1_instance.base_parameters["V_liq"] * \
        (adm1_instance.inlet_conditions[0: 35]-c[0:35].reshape(-1, 1))
    dCdt[35:] = dCdt[35:]+q_gas/adm1_instance.base_parameters["V_gas"] * \
        (adm1_instance.inlet_conditions[35:]-c[35:].reshape(-1, 1))
    dCdt[[26, 32, 34], 0] = 0
    if adm1_instance.switch == "DAE":
        dCdt[7] = 0
        dCdt[27: 32] = 0
        dCdt[33] = 0
    return dCdt[:, 0]


def build_modified_adm_stoichiometric_matrix(base_parameters: dict, model_parameters: dict, reactions: list, species: list)->np.ndarray:
    """ 
    This function builds the stoichiometric matrix for the modified ADM Model.
        
        Model Parameters (dict): a dictionary which contains model parameters
        base_parameters (dict): a dictionary which contains base paramters
        Initial Conditions (dict): a dictionary containing inlet conditions for all species
        Inlet Conditions (dict): a dictionary containing inlet conditions for all species
        reactions (list): a list containing all of the reaction names
        species (list): a list containing all species
        
    Returns:
        np.ndarray: Returns an matrix of stochiometic values.
    """
    S = np.zeros((len(species), len(reactions)))
    S[list(map(species.index, ["TSS", "X_ch", "X_pr", "X_li", "X_I"])),
      reactions.index('TSS_Disintegration')] = [-1, model_parameters['f_ch_TSS'], model_parameters['f_pr_TSS'], model_parameters['f_li_TSS'], model_parameters['f_xI_TSS']]
    S[list(map(species.index, ["TDS", "X_ch", "X_pr", "X_li", "S_I"])), reactions.index('TDS_Disintegration')] = [-1,
                                                                                                                  model_parameters['f_ch_TDS'], model_parameters['f_pr_TDS'], model_parameters['f_li_TDS'], model_parameters['f_sI_TDS']]
    S[list(map(species.index, ["X_ch", "S_su"])),
      reactions.index('Hydrolysis carbohydrates')] = [-1, 1]
    S[list(map(species.index, ["X_pr", "S_aa"])),
      reactions.index('Hydrolysis proteins')] = [-1, 1]
    S[list(map(species.index, ["X_li", "S_fa"])),
      reactions.index('Hydrolysis lipids')] = [-1, 1]

    f_IC_su = -(-model_parameters['C_su'] +
                (1-model_parameters['Y_su'])*model_parameters['f_pro_su']*model_parameters['C_pro'] +
                (1-model_parameters['Y_su'])*model_parameters['f_et_su']*model_parameters['C_et'] +
                (1-model_parameters['Y_su'])*model_parameters['f_lac_su']*model_parameters['C_lac'] +
                (1-model_parameters['Y_su'])*model_parameters['f_ac_su']*model_parameters['C_ac'] +
                model_parameters['Y_su']*model_parameters['C_bac'])

    S[list(map(species.index, ["S_su", "S_pro", "S_et", "S_lac", "S_ac", "S_IN", "S_IC", "X_su"])),
      reactions.index('Uptake of sugars')] = [-1,
                                              (1-model_parameters['Y_su']) *
                                              model_parameters['f_pro_su'],
                                              (1-model_parameters['Y_su']) *
                                              model_parameters['f_et_su'],
                                              (1-model_parameters['Y_su']) *
                                              model_parameters['f_lac_su'],
                                              (1-model_parameters['Y_su']) *
                                              model_parameters['f_ac_su'],
                                              -model_parameters['N_bac']*model_parameters['Y_su'],
                                              f_IC_su,
                                              model_parameters['Y_su']]

    f_IC_aa = -(-model_parameters['C_aa'] +
                (1-model_parameters['Y_aa'])*model_parameters['f_pro_aa']*model_parameters['C_pro'] +
                (1-model_parameters['Y_aa'])*model_parameters['f_et_aa']*model_parameters['C_et'] +
                (1-model_parameters['Y_aa'])*model_parameters['f_lac_aa']*model_parameters['C_lac'] +
                (1-model_parameters['Y_aa'])*model_parameters['f_ac_aa']*model_parameters['C_ac'] +
                (1-model_parameters['Y_aa'])*model_parameters['C_bac'])

    S[list(map(species.index, ["S_aa", "S_pro", "S_et", "S_lac", "S_ac", "S_IN", "S_IC", "X_aa"])),
      reactions.index('Uptake of amino acids')] = [-1,
                                                   (1-model_parameters['Y_aa']) *
                                                   model_parameters['f_pro_aa'],
                                                   (1-model_parameters['Y_aa']) *
                                                   model_parameters['f_et_aa'],
                                                   (1-model_parameters['Y_aa']) *
                                                   model_parameters['f_lac_aa'],
                                                   (1-model_parameters['Y_aa']) *
                                                   model_parameters['f_ac_aa'],
                                                   model_parameters['N_aa']-model_parameters['Y_aa'] *
                                                   model_parameters['N_bac'],
                                                   f_IC_aa,
                                                   model_parameters['Y_aa']]
    f_IC_fa = -(-model_parameters['C_fa'] +
                (1-model_parameters['Y_fa'])*model_parameters['f_pro_fa']*model_parameters['C_pro'] +
                (1-model_parameters['Y_fa'])*model_parameters['f_et_fa']*model_parameters['C_et'] +
                (1-model_parameters['Y_fa'])*model_parameters['f_lac_fa']*model_parameters['C_lac'] +
                (1-model_parameters['Y_fa'])*model_parameters['f_ac_fa']*model_parameters['C_ac'] +
                (1-model_parameters['Y_fa'])*model_parameters['C_bac'])

    S[list(map(species.index, ["S_fa", "S_pro", "S_et", "S_lac", "S_ac", "S_IN", "S_IC", "X_fa"])),
      reactions.index('Uptake of LCFA')] = [-1,
                                            (1-model_parameters['Y_fa']) *
                                            model_parameters['f_pro_fa'],
                                            (1-model_parameters['Y_fa']) *
                                            model_parameters['f_et_fa'],
                                            (1-model_parameters['Y_fa']) *
                                            model_parameters['f_lac_fa'],
                                            (1-model_parameters['Y_fa']) *
                                            model_parameters['f_ac_fa'],
                                            -model_parameters['Y_fa'] *
                                            model_parameters['N_bac'],
                                            f_IC_fa,
                                            model_parameters['Y_fa']]

    f_IC_ac_et = -(-model_parameters['C_ac'] +
                   (1-model_parameters['Y_ac_et'])*model_parameters['f_et_ac']*model_parameters['C_et'] +
                   (1-model_parameters['Y_ac_et']) *
                   model_parameters['f_bu_ac']*model_parameters['C_bu'] +
                   (1-model_parameters['Y_ac_et'])*model_parameters['C_bac'])

    f_IC_ac_lac = -(-model_parameters['C_ac'] +
                    (1-model_parameters['Y_ac_lac'])*model_parameters['f_lac_ac']*model_parameters['C_lac'] +
                    (1-model_parameters['Y_ac_lac']) *
                    model_parameters['f_bu_ac']*model_parameters['C_bu'] +
                    (1-model_parameters['Y_ac_lac'])*model_parameters['C_bac'])

    S[list(map(species.index, ["S_ac", "S_et", "S_bu", "S_IN", "S_IC", "S_h2", "X_ac_et"])),
      reactions.index('Uptake of acetate_et')] = [-1,
                                                  (1-model_parameters['Y_ac_et']) *
                                                  model_parameters['f_et_ac'],
                                                  (1-model_parameters['Y_ac']) *
                                                  model_parameters['f_bu_ac'],
                                                  f_IC_ac_et,
                                                  -model_parameters['Y_ac_et'] *
                                                  model_parameters['N_bac'],
                                                  (1-model_parameters['Y_ac_et']) *
                                                  model_parameters['f_h2_ac'],
                                                  model_parameters['Y_ac_et']]

    S[list(map(species.index, ["S_ac", "S_lac", "S_bu", "S_IN", "S_IC", "S_h2", "X_ac_lac"])),
        reactions.index('Uptake of acetate_lac')] = [-1,
                                                     (1-model_parameters['Y_ac_lac']) *
                                                     model_parameters['f_lac_ac'],
                                                     (1-model_parameters['Y_ac_lac']) *
                                                     model_parameters['f_bu_ac'],
                                                     f_IC_ac_lac,
                                                     -model_parameters['Y_ac_lac'] *
                                                     model_parameters['N_bac'],
                                                     (1-model_parameters['Y_ac_lac']) *
                                                     model_parameters['f_h2_ac'],
                                                     model_parameters['Y_ac_lac']]

    f_IC_pro_et = -(-model_parameters['C_pro'] +
                    (1-model_parameters['Y_pro_et'])*model_parameters['f_et_pro']*model_parameters['C_et'] +
                    (1-model_parameters['Y_pro_et'])*model_parameters['f_va_pro']*model_parameters['C_va'] +
                    (1-model_parameters['Y_pro_et'])*model_parameters['C_bac'])

    f_IC_pro_lac = -(-model_parameters['C_pro'] +
                     (1-model_parameters['Y_pro_lac'])*model_parameters['f_lac_pro']*model_parameters['C_lac'] +
                     (1-model_parameters['Y_pro_lac'])*model_parameters['f_va_pro']*model_parameters['C_va'] +
                     (1-model_parameters['Y_pro_lac'])*model_parameters['C_bac'])

    S[list(map(species.index, ["S_pro", "S_et", "S_va", "S_IN", "S_IC", "S_h2", "X_chain_et"])),
      reactions.index('Uptake of propionate_et')] = [-1,
                                                     (1-model_parameters['Y_pro_et']) *
                                                     model_parameters['f_et_pro'],
                                                     (1-model_parameters['Y_pro_et']) *
                                                     model_parameters['f_va_pro'],
                                                     f_IC_pro_et,
                                                     -model_parameters['Y_pro_et'] *
                                                     model_parameters['N_bac'],
                                                     (1-model_parameters['Y_pro_et']) *
                                                     model_parameters['f_h2_pro'],
                                                     model_parameters['Y_chain_et_pro']]

    S[list(map(species.index, ["S_pro", "S_lac", "S_va", "S_IN", "S_IC", "S_h2", "X_chain_lac"])),
        reactions.index('Uptake of propionate_lac')] = [-1,
                                                        (1-model_parameters['Y_pro_lac']) *
                                                        model_parameters['f_lac_pro'],
                                                        (1-model_parameters['Y_pro_lac']) *
                                                        model_parameters['f_va_pro'],
                                                        f_IC_pro_lac,
                                                        -model_parameters['Y_pro_lac'] *
                                                        model_parameters['N_bac'],
                                                        (1-model_parameters['Y_pro_lac']) *
                                                        model_parameters['f_h2_pro'],
                                                        model_parameters['Y_chain_lac_pro']]

    f_IC_bu_et = -(-model_parameters['C_bu'] +
                   (1-model_parameters['Y_bu_et'])*model_parameters['f_et_bu']*model_parameters['C_et'] +
                   (1-model_parameters['Y_bu_et'])*model_parameters['f_cap_bu']*model_parameters['C_cap'] +
                   (1-model_parameters['Y_bu_et'])*model_parameters['C_bac'])

    f_IC_bu_lac = -(-model_parameters['C_bu'] +
                    (1-model_parameters['Y_bu_lac'])*model_parameters['f_lac_bu']*model_parameters['C_lac'] +
                    (1-model_parameters['Y_bu_lac'])*model_parameters['f_cap_bu']*model_parameters['C_cap'] +
                    (1-model_parameters['Y_bu_lac'])*model_parameters['C_bac'])

    S[list(map(species.index, ["S_bu", "S_et", "S_cap", "S_IN", "S_IC", "S_h2", "X_chain_et"])),
        reactions.index('Uptake of butyrate_et')] = [-1,
                                                     (1-model_parameters['Y_bu_et']
                                                      ) * model_parameters['f_et_bu'],
                                                     (1-model_parameters['Y_bu_et']) *
                                                     model_parameters['f_cap_bu'],
                                                     f_IC_bu_et,
                                                     -model_parameters['Y_bu_et'] * model_parameters['N_bac'],
                                                     (1-model_parameters['Y_bu_et']
                                                      )*model_parameters['f_h2_bu'],
                                                     model_parameters['Y_bu_et']]

    S[list(map(species.index, ["S_bu", "S_lac", "S_cap", "S_IN", "S_IC", "S_h2", "X_chain_lac"])),
        reactions.index('Uptake of butyrate_lac')] = [-1,
                                                      (1-model_parameters['Y_bu_lac']) *
                                                      model_parameters['f_lac_bu'],
                                                      (1-model_parameters['Y_bu_lac']) *
                                                      model_parameters['f_cap_bu'],
                                                      f_IC_bu_lac,
                                                      -model_parameters['Y_bu_lac'] *
                                                      model_parameters['N_bac'],
                                                      (1-model_parameters['Y_bu_lac']
                                                       )*model_parameters['f_h2_bu'],
                                                      model_parameters['Y_bu_lac']]

    S[list(map(species.index, ["S_va", "S_ac", "X_VFA_deg"])),
        reactions.index('Uptake of valerate')] = [-1,
                                                  (1-model_parameters['Y_va']),
                                                  model_parameters['Y_va']]

    S[list(map(species.index, ["S_cap", "S_ac", "X_VFA_deg"])),
        reactions.index('Uptake of caproate')] = [-1,
                                                  (1 -
                                                   model_parameters['Y_cap']),
                                                  model_parameters['Y_cap']]
    f_IC_Me_ach2 = 0
    S[list(map(species.index, ["S_h2", "S_ac", "S_ch4", "X_Me_ac", 'S_IC'])),
        reactions.index('Methanogenessis from acetate and h2')] = [-1,
                                                                   (1 - model_parameters['Y_h2_ac']
                                                                    )*model_parameters['f_ac_h2'],
                                                                   (1 -
                                                                       model_parameters['Y_Me_ac']),
                                                                   model_parameters['Y_Me_ac'],
                                                                   f_IC_Me_ach2]

    f_IC_Me_CO2h2 = -(model_parameters['Y_Me_CO2']*model_parameters['C_ch4'] +
                      model_parameters['Y_Me_h2']*model_parameters['C_bac'])
    
    S[list(map(species.index, ["S_h2", "S_ch4", "X_Me_CO2", 'S_IC'])),
        reactions.index('Methanogenessis from CO2 and h2')] = [-1,
                                                               (1 -
                                                                model_parameters['Y_h2_CO2']),
                                                               (model_parameters['Y_Me_CO2']),
                                                               f_IC_Me_CO2h2]
    
    f_IC_et_ox=-(-model_parameters['C_et'] +
                    (1-model_parameters['Y_ac_et_ox'])*model_parameters['C_bac']
                    +model_parameters['Y_ac_et_ox']*model_parameters['C_ac'])

    S[list(map(species.index, ["S_et", "X_et","S_ac","S_IC"])),
        reactions.index('Uptake of ethanol')] = [-1,1-model_parameters['Y_ac_et_ox'],model_parameters['Y_ac_et_ox'],f_IC_et_ox]

    
    
    f_IC_lac_ox=-(-model_parameters['C_lac'] +
                (1-model_parameters['Y_pro_lac_ox'])*model_parameters['C_bac']
                +model_parameters['Y_pro_lac_ox']*model_parameters['C_pro'])
    
    S[list(map(species.index, ["S_lac", "X_lac","S_pro","S_IC"])),
        reactions.index('Uptake of lactate')] = [-1, 1-model_parameters['Y_pro_lac_ox'],model_parameters['Y_pro_lac_ox'],f_IC_lac_ox]

    S[list(map(species.index, ["X_su", "TSS"])),
        reactions.index('Decay of Xsu')] = [-1, 1]

    S[list(map(species.index, ["X_aa", "TSS"])),
        reactions.index('Decay of Xaa')] = [-1, 1]

    S[list(map(species.index, ["X_fa", "TSS"])),
        reactions.index('Decay of Xfa')] = [-1, 1]

    S[list(map(species.index, ["X_ac_et", "TSS"])),
        reactions.index('Decay of X_ac_et')] = [-1, 1]

    S[list(map(species.index, ["X_ac_lac", "TSS"])),
        reactions.index('Decay of X_ac_lac')] = [-1, 1]

    S[list(map(species.index, ["X_chain_et", "TSS"])),
        reactions.index('Decay of X_chain_et')] = [-1, 1]

    S[list(map(species.index, ["X_chain_lac", "TSS"])),
        reactions.index('Decay of X_chain_lac')] = [-1, 1]

    S[list(map(species.index, ["X_VFA_deg", "TSS"])),
        reactions.index('Decay of X_VFA_deg')] = [-1, 1]

    S[list(map(species.index, ["X_Me_ac", "TSS"])),
        reactions.index('Decay of X_Me_ac')] = [-1, 1]

    S[list(map(species.index, ["X_Me_CO2", "TSS"])),
        reactions.index('Decay of X_Me_CO2')] = [-1, 1]

    S[list(map(species.index, ["S_va_ion"])),
        reactions.index('Acid Base Equilibrium (Va)')] = [-1]

    S[list(map(species.index, ["S_bu_ion"])),
        reactions.index('Acid Base Equilibrium (Bu)')] = [-1]

    S[list(map(species.index, ["S_pro_ion"])),
        reactions.index('Acid Base Equilibrium (Pro)')] = [-1]

    S[list(map(species.index, ["S_cap_ion"])),
        reactions.index('Acid Base Equilibrium (Cap)')] = [-1]

    S[list(map(species.index, ["S_lac_ion"])),
        reactions.index('Acid Base Equilibrium (Lac)')] = [-1]

    S[list(map(species.index, ["S_ac_ion"])),
        reactions.index('Acid Base Equilibrium (Ac)')] = [-1]

    # S[list(map(species.index, ["S_CO2", "S_hco3_ion"])),  # I don't think this is right، should look at the reaction in ADM1
    #     reactions.index('Acid Base Equilibrium (CO2)')] = [-1, 1]

    # S[list(map(species.index, ["S_nh3", "S_nh4_ion"])),
    #     reactions.index('Acid Base Equilibrium (IN)')] = [-1, 1]  # I don't think this is right، should look at the reaction in ADM1

    S[list(map(species.index, ["S_h2", "S_gas_h2"])),
        reactions.index('Gas Transfer H2')] = [-base_parameters['V_liq']/base_parameters['V_gas'], 1]
    S[list(map(species.index, ["S_ch4", "S_gas_ch4"])),
        reactions.index('Gas Transfer CH4')] = [-base_parameters['V_liq']/base_parameters['V_gas'], 1]
    S[list(map(species.index, ["S_co2", "S_gas_co2"])),
        reactions.index('Gas Transfer CO2')] = [-base_parameters['V_liq']/base_parameters['V_gas'], 1]
    return S


def modified_adm_ode_sys(t: float, c: np.ndarray, model: Model)-> np.ndarray:
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
    c[model.species.index('S_nh4_ion')] = c[model.species.index(
        'S_IN')] - c[model.species.index('S_nh3')]
    c[model.species.index('S_co2')] = c[model.species.index(
        'S_IC')] - c[model.species.index('S_hco3_ion')]
    I_pH_aa = (model.model_parameters["K_pH_aa"] ** model.model_parameters['nn_aa'])/(np.power(
        c[model.species.index('S_H_ion')], model.model_parameters['nn_aa']) + np.power(model.model_parameters["K_pH_aa"], model.model_parameters['nn_aa']))

    I_pH_ac = (model.model_parameters['K_pH_ac'] ** model.model_parameters["n_ac"])/(
        c[model.species.index('S_H_ion')] ** model.model_parameters['n_ac'] + model.model_parameters['K_pH_ac'] ** model.model_parameters['n_ac'])

    I_pH_pro = (model.model_parameters['K_pH_pro'] ** model.model_parameters["n_pro"])/(
        c[model.species.index('S_H_ion')] ** model.model_parameters['n_pro'] + model.model_parameters['K_pH_pro'] ** model.model_parameters['n_pro'])

    I_pH_bu = (model.model_parameters['K_pH_bu'] ** model.model_parameters["n_bu"])/(
        c[model.species.index('S_H_ion')] ** model.model_parameters['n_bu'] + model.model_parameters['K_pH_bu'] ** model.model_parameters['n_bu'])

    I_pH_va = (model.model_parameters['K_pH_va'] ** model.model_parameters["n_va"])/(
        c[model.species.index('S_H_ion')] ** model.model_parameters['n_va'] + model.model_parameters['K_pH_va'] ** model.model_parameters['n_va'])

    I_pH_cap = (model.model_parameters['K_pH_cap'] ** model.model_parameters["n_cap"])/(
        c[model.species.index('S_H_ion')] ** model.model_parameters['n_cap'] + model.model_parameters['K_pH_cap'] ** model.model_parameters['n_cap'])

    I_pH_h2 = (model.model_parameters['K_pH_h2']**model.model_parameters['n_h2'])/(
        c[model.species.index('S_H_ion')] ** model.model_parameters['n_h2'] + model.model_parameters['K_pH_h2']**model.model_parameters['n_h2'])

    I_IN_lim = 1 / \
        (1+(model.model_parameters['K_S_IN'] / (c[model.species.index('S_IN')]+10**-9)))

    I_h2_fa = 1 / (1+(c[model.species.index('S_h2')] /
                   (model.model_parameters['K_I_h2_fa']+10**-9)))

    I_h2_c4 = 1 / (1+(c[model.species.index('S_h2')] /
                   (model.model_parameters['K_I_h2_c4']+10**-9)))

    I_h2_pro = (1/(1+(c[model.species.index('S_h2')] /
                (model.model_parameters['K_I_h2_pro']+10**-9))))

    I_nh3 = 1/(1+(c[model.species.index('S_nh3')] /
               (model.model_parameters['K_I_nh3']+10**-9)))

    I_h2_oxidation=(1/(1+(c[model.species.index('S_h2')] /
                (model.model_parameters['K_I_h2_ox']+10**-9))))

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

    v = np.zeros((len(model.reactions), 1))

    v[model.reactions.index(
        'TSS_Disintegration')] = model.model_parameters["k_dis_TSS"]*c[model.species.index('TSS')]

    v[model.reactions.index(
        'TDS_Disintegration')] = model.model_parameters["k_dis_TDS"]*c[model.species.index('TDS')]

    v[model.reactions.index('Hydrolysis carbohydrates')
      ] = model.model_parameters['k_hyd_ch']*c[model.species.index('X_ch')]

    v[model.reactions.index('Hydrolysis proteins')
      ] = model.model_parameters['k_hyd_pr']*c[model.species.index('X_pr')]

    v[model.reactions.index('Hydrolysis lipids')
      ] = model.model_parameters['k_hyd_li']*c[model.species.index('X_li')]

    v[model.reactions.index('Uptake of sugars')] = model.model_parameters['k_m_su']*c[model.species.index('S_su')] / \
        (model.model_parameters['K_S_su']+c[model.species.index('S_su')]
         )*c[model.species.index('X_su')]*I5

    v[model.reactions.index('Uptake of amino acids')] = model.model_parameters['k_m_aa']*c[model.species.index('S_aa')] / \
        (model.model_parameters['K_S_aa']+c[model.species.index('S_aa')]
         )*c[model.species.index('X_aa')]*I6

    v[model.reactions.index('Uptake of LCFA')] = model.model_parameters['k_m_fa']*c[model.species.index('S_fa')] / \
        (model.model_parameters['K_S_fa'] +
         c[model.species.index('S_fa')])*c[model.species.index('X_fa')]*I7

    v[model.reactions.index('Uptake of acetate_et')] = model.model_parameters['k_m_ac']*c[model.species.index('S_ac')]*c[model.species.index('S_et')] / \
        (model.model_parameters['K_S_ac']*c[model.species.index('S_ac')]+model.model_parameters['K_S_ac_et']*c[model.species.index('S_et')]+c[model.species.index('S_ac')]*c[model.species.index('S_et')]+10**-9
         )*c[model.species.index('X_ac_et')]*I11

    v[model.reactions.index('Uptake of acetate_lac')] = model.model_parameters['k_m_ac']*c[model.species.index('S_ac')]*c[model.species.index('S_lac')] / \
        (model.model_parameters['K_S_ac']*c[model.species.index('S_ac')]+model.model_parameters['K_S_ac_lac']*c[model.species.index('S_lac')]+c[model.species.index('S_ac')]*c[model.species.index('S_lac')]+10**-9
         )*c[model.species.index('X_ac_lac')]*I11

    v[model.reactions.index('Uptake of propionate_et')] = model.model_parameters['k_m_pro']*c[model.species.index('S_pro')]*c[model.species.index('S_et')] / \
        (model.model_parameters['K_S_pro']*c[model.species.index('S_pro')]+model.model_parameters['K_S_pro_et']*c[model.species.index('S_et')]+c[model.species.index('S_pro')]*c[model.species.index('S_et')]+10**-9
         )*c[model.species.index('X_chain_et')]*I10

    v[model.reactions.index('Uptake of propionate_lac')] = model.model_parameters['k_m_pro']*c[model.species.index('S_pro')]*c[model.species.index('S_lac')] / \
        (model.model_parameters['K_S_pro']*c[model.species.index('S_pro')]+model.model_parameters['K_S_pro_lac']*c[model.species.index('S_lac')]+c[model.species.index('S_pro')]*c[model.species.index('S_lac')]+10**-9
         )*c[model.species.index('X_chain_lac')]*I10

    v[model.reactions.index('Uptake of butyrate_et')] = model.model_parameters['k_m_bu']*c[model.species.index('S_bu')]*c[model.species.index('S_et')] / \
        (model.model_parameters['K_S_bu']*c[model.species.index('S_bu')]+model.model_parameters['K_S_bu_et']*c[model.species.index('S_et')]+c[model.species.index('S_bu')]*c[model.species.index('S_et')]+10**-9
         )*c[model.species.index('X_chain_et')]*I14

    v[model.reactions.index('Uptake of butyrate_lac')] = model.model_parameters['k_m_bu']*c[model.species.index('S_bu')]*c[model.species.index('S_lac')] / \
        (model.model_parameters['K_S_bu']*c[model.species.index('S_bu')]+model.model_parameters['K_S_bu_lac']*c[model.species.index('S_lac')]+c[model.species.index('S_bu')]*c[model.species.index('S_lac')]+10**-9
         )*c[model.species.index('X_chain_lac')]*I14

    v[model.reactions.index('Uptake of valerate')] = model.model_parameters['k_m_va']*c[model.species.index('S_va')] / \
        (model.model_parameters['K_S_va']+c[model.species.index('S_va')]
         )*c[model.species.index('X_VFA_deg')]*I15

    v[model.reactions.index('Uptake of caproate')] = model.model_parameters['k_m_cap']*c[model.species.index('S_cap')] / \
        (model.model_parameters['K_S_cap']+c[model.species.index('S_cap')]
         )*c[model.species.index('X_VFA_deg')]*I13

    v[model.reactions.index('Methanogenessis from acetate and h2')] = model.model_parameters['k_m_h2_Me_ac']*c[model.species.index('S_h2')]*c[model.species.index('S_ac')] / \
        (model.model_parameters['K_S_h2_Me_ac']*c[model.species.index('S_h2')]+model.model_parameters['K_S_ac_Me']*c[model.species.index(
            'S_ac')]+c[model.species.index('S_ac')]*c[model.species.index('S_h2')]+10**-9)*c[model.species.index('X_Me_ac')]*I12

    v[model.reactions.index('Methanogenessis from CO2 and h2')] = model.model_parameters['k_m_h2_Me_CO2']*c[model.species.index('S_h2')]*c[model.species.index('S_co2')] / \
        (model.model_parameters['K_S_h2_Me_CO2']*c[model.species.index('S_h2')]+model.model_parameters['K_S_CO2_Me']*c[model.species.index(
            'S_co2')]+c[model.species.index('S_co2')]*c[model.species.index('S_h2')]+10**-9)*c[model.species.index('X_Me_CO2')]*I12


    v[model.reactions.index('Uptake of ethanol')] = model.model_parameters['k_m_et']*c[model.species.index('S_et')] / \
        (model.model_parameters['K_S_et']+c[model.species.index('S_et')]
         )*c[model.species.index("X_et")]*I16

    v[model.reactions.index('Uptake of lactate')] = model.model_parameters['k_m_lac']*c[model.species.index('S_lac')] / \
        (model.model_parameters['K_S_lac']+c[model.species.index('S_lac')]
         )*c[model.species.index('X_lac')]*I16

    v[model.reactions.index(
        'Decay of Xsu')] = model.model_parameters['k_dec_X_su']*c[model.species.index('X_su')]

    v[model.reactions.index(
        'Decay of Xaa')] = model.model_parameters['k_dec_X_aa']*c[model.species.index('X_aa')]

    v[model.reactions.index(
        'Decay of Xfa')] = model.model_parameters['k_dec_X_fa']*c[model.species.index('X_fa')]

    v[model.reactions.index(
        'Decay of X_ac_et')] = model.model_parameters['k_dec_X_ac']*c[model.species.index('X_ac_et')]

    v[model.reactions.index(
        'Decay of X_ac_lac')] = model.model_parameters['k_dec_X_ac']*c[model.species.index('X_ac_lac')]

    v[model.reactions.index(
        'Decay of X_chain_et')] = model.model_parameters['k_dec_X_chain_et']*c[model.species.index('X_chain_et')]

    v[model.reactions.index('Decay of X_chain_lac')
      ] = model.model_parameters['k_dec_X_chain_lac']*c[model.species.index('X_chain_lac')]

    v[model.reactions.index(
        'Decay of X_VFA_deg')] = model.model_parameters['k_dec_X_VFA_deg']*c[model.species.index('X_VFA_deg')]

    v[model.reactions.index(
        'Decay of X_Me_ac')] = model.model_parameters['k_dec_X_Me_ac']*c[model.species.index('X_Me_ac')]

    v[model.reactions.index(
        'Decay of X_Me_CO2')] = model.model_parameters['k_dec_X_Me_CO2']*c[model.species.index('X_Me_CO2')]

    v[model.reactions.index(
        'Decay of Xet')] = model.model_parameters['k_dec_X_et']*c[model.species.index('X_et')]

    v[model.reactions.index(
        'Decay of Xlac')] = model.model_parameters['k_dec_X_lac']*c[model.species.index('X_lac')]

    v[model.reactions.index('Acid Base Equilibrium (Va)')] = model.model_parameters['k_A_B_va'] * \
        (c[model.species.index('S_va_ion')] * (model.model_parameters['K_a_va'] + c[model.species.index('S_H_ion')]) -
         model.model_parameters['K_a_va'] * c[model.species.index('S_va')])

    v[model.reactions.index('Acid Base Equilibrium (Bu)')] = model.model_parameters['k_A_B_bu'] * \
        (c[model.species.index('S_bu_ion')] * (model.model_parameters['K_a_bu'] + c[model.species.index('S_H_ion')]) -
         model.model_parameters['K_a_bu'] * c[model.species.index('S_bu')])

    v[model.reactions.index('Acid Base Equilibrium (Pro)')] = model.model_parameters['k_A_B_pro'] * \
        (c[model.species.index('S_pro_ion')] * (model.model_parameters['K_a_pro'] + c[model.species.index('S_H_ion')]) -
         model.model_parameters['K_a_pro'] * c[model.species.index('S_pro')])

    v[model.reactions.index('Acid Base Equilibrium (Cap)')] = model.model_parameters['k_A_B_cap'] * \
        (c[model.species.index('S_cap_ion')] * (model.model_parameters['K_a_cap'] + c[model.species.index('S_H_ion')]) -
         model.model_parameters['K_a_cap'] * c[model.species.index('S_cap')])

    v[model.reactions.index('Acid Base Equilibrium (Lac)')] = model.model_parameters['k_A_B_lac'] * \
        (c[model.species.index('S_lac_ion')] * (model.model_parameters['K_a_lac'] + c[model.species.index('S_H_ion')]) -
         model.model_parameters['K_a_lac'] * c[model.species.index('S_lac')])

    v[model.reactions.index('Acid Base Equilibrium (Ac)')] = model.model_parameters['k_A_B_ac'] * \
        (c[model.species.index('S_ac_ion')] * (model.model_parameters['K_a_ac'] + c[model.species.index('S_H_ion')]) -
         model.model_parameters['K_a_ac'] * c[model.species.index('S_ac')])

    v[model.reactions.index('Acid Base Equilibrium (CO2)')] = model.model_parameters['k_A_B_co2'] * \
        (c[model.species.index('S_hco3_ion')] * (model.model_parameters['K_a_co2'] + c[model.species.index('S_H_ion')]) -
         model.model_parameters['K_a_co2'] * c[model.species.index('S_IC')])

    v[model.reactions.index('Acid Base Equilibrium (In)')] = model.model_parameters['k_A_B_IN'] * \
        (c[model.species.index('S_nh3')] * (model.model_parameters['K_a_IN'] + c[model.species.index('S_H_ion')]) -
         model.model_parameters['K_a_IN'] * c[model.species.index('S_IC')])

    p_gas_h2 = c[model.species.index('S_gas_h2')] * model.base_parameters["R"] * \
        model.base_parameters["T_op"] / 16
    p_gas_ch4 = c[model.species.index('S_gas_ch4')] * model.base_parameters["R"] * \
        model.base_parameters["T_op"] / 64
    p_gas_co2 = c[model.species.index('S_gas_co2')] * model.base_parameters["R"] * \
        model.base_parameters["T_op"]
    p_gas_h2o = 0.0313 * \
        np.exp(5290 *
               (1 / model.base_parameters["T_base"] - 1 / model.base_parameters["T_op"]))
    P_gas = p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o
    q_gas = max(
        0, (model.model_parameters['k_p'] * (P_gas - model.base_parameters['P_atm'])))
    v[model.reactions.index('Gas Transfer H2')] = model.model_parameters['k_L_a'] * \
        (c[model.species.index('S_h2')] - 16 *
         model.model_parameters['K_H_h2'] * p_gas_h2)

    v[model.reactions.index('Gas Transfer CH4')] = max(0,model.model_parameters['k_L_a'] * \
        (c[model.species.index('S_ch4')] - 64 *
         model.model_parameters['K_H_ch4'] * p_gas_ch4))
    v[model.reactions.index('Gas Transfer CO2')] = max(0,model.model_parameters['k_L_a'] * \
        (c[model.species.index('S_co2')] -
         model.model_parameters['K_H_co2'] * p_gas_co2))

    dCdt = np.matmul(model.S, v)

    phi = c[model.species.index('S_cation')]+c[model.species.index('S_nh4_ion')]-c[model.species.index('S_hco3_ion')]-(c[model.species.index('S_lac_ion')] / 88) - (c[model.species.index('S_ac_ion')] / 64) - (c[model.species.index('S_pro_ion')] /
                                                                                                                                                                     112) - (c[model.species.index('S_bu_ion')] / 160)-(c[model.species.index('S_cap_ion')] / 230) - (c[model.species.index('S_va_ion')] / 208) - c[model.species.index('S_anion')]
    if 'S_H_ion' in model.control_state.keys():
        c[model.species.index('S_H_ion')]=model.control_state['S_H_ion']
    else:
        c[model.species.index('S_H_ion')] = (-1 * phi / 2) + \
        (0.5 * np.sqrt(phi**2 + 4 * model.model_parameters['K_w']))

    dCdt[0: model.species.__len__()-3] = dCdt[0: model.species.__len__()-3]+model.base_parameters['q_in'] / \
        model.base_parameters["V_liq"] * \
        (model.inlet_conditions[0: model.species.__len__(
        )-3]-c[0: model.species.__len__()-3].reshape(-1, 1))

    dCdt[model.species.__len__()-3:] = dCdt[model.species.__len__()-3:]+q_gas/model.base_parameters["V_gas"] * \
        (model.inlet_conditions[model.species.__len__() -
         3:]-c[model.species.__len__()-3:].reshape(-1, 1))

    dCdt[[model.species.index('S_H_ion'), model.species.index(
        'S_co2'), model.species.index('S_nh4_ion')], 0] = 0

    if model.switch == "DAE":
        dCdt[model.species.index('S_h2')] = 0

        dCdt[model.species.index('S_va_ion'):model.species.index('S_co2')] = 0

        dCdt[model.species.index('S_nh3')] = 0
    
        c[model.species.index('S_va_ion')]=model.model_parameters['K_a_va']/(model.model_parameters['K_a_va']+c[model.species.index('S_H_ion')])*c[model.species.index('S_va')]
    
        c[model.species.index('S_bu_ion')]=model.model_parameters['K_a_bu']/(model.model_parameters['K_a_bu']+c[model.species.index('S_H_ion')])*c[model.species.index('S_bu')]
    
        c[model.species.index('S_pro_ion')]=model.model_parameters['K_a_pro']/(model.model_parameters['K_a_pro']+c[model.species.index('S_H_ion')])*c[model.species.index('S_pro')]
    
        c[model.species.index('S_cap_ion')]=model.model_parameters['K_a_cap']/(model.model_parameters['K_a_cap']+c[model.species.index('S_H_ion')])*c[model.species.index('S_cap')]
    
        c[model.species.index('S_ac_ion')]=model.model_parameters['K_a_ac']/(model.model_parameters['K_a_ac']+c[model.species.index('S_H_ion')])*c[model.species.index('S_ac')]
        
        c[model.species.index('S_lac_ion')]=model.model_parameters['K_a_lac']/(model.model_parameters['K_a_lac']+c[model.species.index('S_H_ion')])*c[model.species.index('S_lac')]    

    if model.control_state.keys():
        for state in model.control_state.keys():
            c[model.species.index(state)]=model.control_state[state]
            dCdt[model.species.index(state)]=0
    model.info["Fluxes"].append(v)
    return dCdt[:, 0]


if __name__ == "__main__":
   # Report = os.path.join(Main_Dir, "..", "Reports",
                         # "ADM_From_Alignment_JSON_Output.json")
   # with open(Report, 'r') as j:
        #Report = json.load(j)

   #adm1 = Model(model_parameters, base_parameters, Initial_Conditions,
   #             inlet_conditions, reactions, species, ADM1_ODE_Sys, Build_ADM1_Stoiciometric_Matrix, Name="ADM1", Switch="DAE")
   #Sol = adm1.Solve_Model(
   #    (0, 20), adm1.Initial_Conditions[:, 0], np.linspace(0, 20, 100))

    # adm1.Dash_App(Sol)
    with open("../Optimization/Tuned_Params.json", 'r') as f:
        mp=json.load(f)
    with open('/Users/parsaghadermarzi/Desktop/ADM_Mapping.json', 'r') as f:
        MR=json.load(f)
    with open('/Users/parsaghadermarzi/Desktop/ADToolbox/Database/Modified_ADM/Modified_ADM_base_parameters.json', 'r') as f:
        BP=json.load(f)
    with open('/Users/parsaghadermarzi/Desktop/ADToolbox/Database/Modified_ADM/Modified_ADM_Initial_Conditions.json', 'r') as f:
        IC=json.load(f)
    with open('/Users/parsaghadermarzi/Desktop/ADToolbox/Database/Modified_ADM/Modified_ADM_inlet_conditions.json', 'r') as f:
        InC=json.load(f)
    with open('/Users/parsaghadermarzi/Desktop/ADToolbox/Database/Modified_ADM/Modified_ADM_reactions.json', 'r') as f:
        r=json.load(f)
    with open('/Users/parsaghadermarzi/Desktop/ADToolbox/Database/Modified_ADM/Modified_ADM_species.json', 'r') as f:
        s=json.load(f)
    
    mod_adm1 = Model(mp, BP, IC, InC, r,
                     s, modified_adm_ode_sys, build_modified_adm_stoichiometric_matrix,control_state={"S_H_ion":0},name="Modified_ADM1", switch="DAE",metagenome_report=MR)
    Sol_mod_adm1 = mod_adm1.solve_model(mod_adm1.initial_conditions[:, 0], np.linspace(0,30, 10000))
    mod_adm1.dash_app(Sol_mod_adm1)