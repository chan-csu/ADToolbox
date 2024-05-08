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
from core import Database as Database
from core import SeedDB as SeedDB
from core import Feed
import pandas as pd
import dash_bootstrap_components as dbc
import utils
from adtoolbox import Main_Dir,PKG_DATA
import dash_escher
import configs
import time

### Note ###
# The following code is a modified version of the code from the PyADM1 package
# It is extensively based on PyADM1, and I would like to thank the author of PyADM1 for this work
# As far as implementing the original ADM1 in python goes, I still consider this as a modification from
# the PyADM1 code.
# ----------

DEFAULT_FEED=Feed(name="Default Feed",
                    carbohydrates=10,
                    proteins=20, 
                    lipids=20, 
                    si=30,
                    xi=50,
                    tss=80)
RT = SeedDB(config=configs.Database())

class _Fake_Sol:
    def __init__(self, y,t):
        self.y = y
        self.t=t

class Model:

    """Any kinetic model could be an instance of this class.
    Args:
        model_parameters (dict): a dictionary which contains model parameters
        base_parameters (dict): a dictionary which contains base paramters
        initial_conditions (dict): a dictionary containing inlet conditions for all species
        inlet_conditions (dict): a dictionary containing inlet conditions for all species
        feed (Feed): a Feed instance which contains the feed information
        reactions (list): a list containing all types of reactions
        species (list): a list containing all species
        ode_system (Callable): a callable which outputs the ODE system compatible with Scipy.integrate.solve_ivp
        build_stoichiometric_matrix(Callable): a callable which builds the stoichiometric matrix
        control_state (dict, optional): a dictionary containing the states that are desired to be constant. Defaults to {}.

        
        
    Returns:
        Model: returns a model instance for downstream purposes.
    """
    def __init__(self, 
                 model_parameters: dict,
                 base_parameters: dict,
                 initial_conditions: dict,
                 inlet_conditions:dict,
                 feed:Feed,
                 reactions: list, 
                 species: list, 
                 ode_system:Callable, 
                 build_stoichiometric_matrix:Callable,
                 control_state:dict={},
                 name:str="ADM", 
                 switch:str="DAE",
                 simulation_time:float=30,
                 time_limit:float=-1):
        
        self.model_parameters = model_parameters
        self.base_parameters = base_parameters
        self.feed=feed
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
        self.build_stoichiometric_matrix = build_stoichiometric_matrix
        self.ode_system = ode_system
        self.sim_time=simulation_time
        self.time_limit=time_limit
        self.nitrogen_limited=False

    @property
    def s(self):
        """Returns the stoichiometric matrix of a model"""
        return self.build_stoichiometric_matrix(
            base_parameters=self.base_parameters,model_parameters= self.model_parameters,reactions= self.reactions,species= self.species,feed=self.feed, nitrogen_limited=self.nitrogen_limited)

    def update_parameters(self, 
                        model_parameters: dict|None=None,
                        base_parameters:  dict|None=None,
                        initial_conditions: dict|None=None,
                        inlet_conditions: dict|None=None)->None:
        """
        This method updates the parameters of the model. Each argument can be a dictionary containing the parameters to be updated.
        NOTE: It is important to note that you have to separate different kind parameters.
        Args:
            model_parameters (dict): a dictionary which contains the model parameters to be updated as keys and their values as values.
            base_parameters (dict): a dictionary which contains the base parameters to be updated as keys and their values as values.
            initial_conditions (dict): a dictionary containing the initial conditions to be updated as keys and their values as values.
            inlet_conditions (dict): a dictionary containing the inlet conditions to be updated as keys and their values as values.

        Returns:
            None: This method does not return anything.
        """
        if model_parameters is not None:
            self.model_parameters.update(model_parameters)
        if base_parameters is not None:
            self.base_parameters.update(base_parameters)
        if initial_conditions is not None:
            for k,v in initial_conditions.items():
                self.initial_conditions[self.species.index(k)]=v
        if inlet_conditions is not None:
            for k,v in inlet_conditions.items():
                self.inlet_conditions[self.species.index(k)]=v
            

    

    def solve_model(self, t_eval: np.ndarray, method="BDF")->scipy.integrate._ivp.ivp.OdeResult:
        """
        Function to solve the model. 
        Examples:
            >>> import numpy as np
            >>> reactions=['rxn1','rxn2']
            >>> species=['a','b','c']
            >>> initial_conditions={'a':.001,'b':.002,'c':.003}
            >>> inlet_conditions={'a_in':.001,'b_in':.002,'c_in':.003}
            >>> model_parameters={'k1':0.001,'k2':0.002}
            >>> base_parameters={'T':0.1}
            >>> feed=Feed(10,20,20,20)
            >>> def build_stoiciometric_matrix(base_parameters,model_parameters,reactions,species):
            ...    s = np.zeros((len(species), len(reactions)))
            ...    s[[0,1],0]=[-1,0.001]
            ...    s[[1,2],1]=[-5,1]
            ...    return s
            >>> def ode_system(t,c,Model1):
            ...    v = np.zeros((len(Model1.reactions), 1))
            ...    v[0]=Model1.model_parameters['k1']*c[0]*Model1.base_parameters['T']/1000
            ...    v[1]=Model1.model_parameters['k2']*c[1]/1000
            ...    dCdt=np.matmul(Model1.S,v)
            ...    return dCdt[:, 0]
            >>> m= Model(model_parameters,base_parameters,initial_conditions,inlet_conditions,reactions,species,ODE_System,Build_Stoiciometric_Matrix)
            >>> m.solve_model((0,.1),np.linspace(0,0.1,10),method='RK45')['status']==0
            True
        
        Args:
            t_eval (np.ndarray): Time points at which the solution is reported
            method (str, optional): The method used to solve the ODE. Defaults to "BDF".
        
        Returns:
            scipy.integrate._ivp.ivp.OdeResult: Returns the results of the simulation being run and gives optimized paramters.
        """
        self.info={"Fluxes":[]}
        y0=self.initial_conditions[:, 0]
        try:
            self._be_time=time.time()
            c = scipy.integrate.solve_ivp(self.ode_system, (0,self.sim_time), y0, t_eval=t_eval, method=method, args=[self],rtol=1e-6)
            if not c.success:
                raise Exception
        except Exception as e:
            print("Could not solve model, setting C to a very large value")
            c=_Fake_Sol(np.ones((y0.shape[0],len(t_eval)))*1e10,t_eval)
       
        return c

    
       #C = scipy.integrate.solve_ivp(
       #        self.ODE_System, t_span, y0, t_eval=T_eval, method=method, args=[self])
       #
       #return C


        
    def plot(self, Sol: scipy.integrate._ivp.ivp.OdeResult, type: str = "Line")-> go.Figure:
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
                          title="Concentration of species")
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

        elif type == "Sankey":
            ### Maybe add a sankey plot here later
            pass
        
        return fig
    
        

    def dash_app(self, sol: scipy.integrate._ivp.ivp.OdeResult,
                 escher_map:str|None=os.path.join(PKG_DATA,"Modified_ADM_Map.json"),
                 cobra_model:str|None=os.path.join(PKG_DATA,"Modified_ADM_Model.json"),
                 **kwargs)->None:
        """A method that creates the dash web app for a model based on an ODE solution.
        
        Examples:
            >>> import numpy as np
            >>> reactions=['rxn1','rxn2']
            >>> species=['a','b','c']
            >>> initial_conditions={'a':.001,'b':.002,'c':.003}
            >>> inlet_conditions={'a_in':.001,'b_in':.002,'c_in':.003}
            >>> model_parameters={'k1':0.001,'k2':0.002}
            >>> base_parameters={'T':0.1}
            >>> feed=Feed(10,20,20,20)
            >>> def build_stoiciometric_matrix(base_parameters,model_parameters,reactions,species):
            ...    s = np.zeros((len(species), len(reactions)))
            ...    s[[0,1],0]=[-1,0.001]
            ...    s[[1,2],1]=[-5,1]
            ...    return s
            >>> def ode_system(t,c,Model1):
            ...    v = np.zeros((len(Model1.reactions), 1))
            ...    v[0]=Model1.model_parameters['k1']*c[0]*Model1.base_parameters['T']/1000
            ...    v[1]=Model1.model_parameters['k2']*c[1]/1000
            ...    dCdt=np.matmul(Model1.S,v)
            ...    return dCdt[:, 0]
            >>> m= Model(model_parameters,base_parameters,initial_conditions,inlet_conditions,reactions,species,ODE_System,Build_Stoiciometric_Matrix)
            >>> m.solve_model((0,.1),np.linspace(0,0.1,10),method='RK45')['status']==0
            True
            >>> m.dash_app(m.solve_model(np.linspace(0,30,1000)))
        
        Args:
            sol (scipy.integrate._ivp.ivp.OdeResult): The solution of the ODE system. This should be the output of the solve_model method.

        Returns:
            None: This method does not return anything.
        
        
        """
        if escher_map is not None:
            with open(escher_map,'rb') as f:
                escher_map=json.load(f)
        if cobra_model is not None:
            with open(cobra_model,'rb') as f:
                cobra_model=json.load(f)

        app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY])
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

        
        fig = px.line(sol_df, x="t", y=sol_df.columns,
                      title="Concentration of species")
        fig.update_layout(
        title={
        'y': 0.95,
        'x': 0.5,
        "font_size": 30,
        'xanchor': 'center',
        'yanchor': 'top'},
        legend=dict(font=dict(size= 20),),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
            )
        fig.update_xaxes(
        title={
        "text": "Time (Days)",
        "font_size": 25,
            },
             tickfont_size=20,
        linecolor='grey',
        gridcolor='grey',
            )
        fig.update_yaxes(
        title={
        "text": "Concentrations (kg COD/m^3)",
        "font_size": 25,
         },
        tickfont_size=20,
        linecolor='grey',
        gridcolor='grey',
        
            )
        fig.update_traces(line=dict(width=3))

        styles={
            'table_width': '95%',
            'padding-left': '20px',
            'container_width': '85%'
        }
        page=[dbc.Container(
                        html.H1("ADToolbox Web Interface",style={"font-size":"70px", "padding-top":"50px"}),className="text-white bg-primary",style={"height":"300px","text-align": "center"}, fluid=True),
                        dbc.Container([dbc.Row(
                                    [dbc.Card([
                                        html.H2(f"{self.name} Concentration Plot", style={
                                            'textAlign': 'left',
                                            'color': colors['text'],
                                            'font-size': '15',
                                            'padding-top': '50px',
                                            'padding-bottom': '20px',
                                            'padding-left': styles['padding-left'] },
                                             className="card-title"),
                                        dcc.Graph(figure=fig, id='Concentrations_Line_Plot',
                                                style={
                                                "height":"600px",
                                                "padding-left": styles['padding-left'],
                                                'background-color': 'rgba(0,0,0,0)'}
                                                ),],className='bg-light'),

                                    dbc.Card([html.H3("Base Parameters", style={
                                        'textAlign': 'left',
                                        'color': colors['text'],
                                        'font-size': '15',
                                        'padding-top': '50px',
                                        'padding-bottom': '20px',
                                        'padding-left': styles['padding-left']
                                        }),
                                        dash_table.DataTable(
                                        id='base_parameters',
                                        columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self.base_parameters.keys())],
                                        data=pd.DataFrame(self.base_parameters,index=[0]).to_dict('records'),
                                        editable=True,
                                        style_table={'overflowX': 'scroll', 'padding-left': '20px','padding-bottom':'30px', 'width': styles['table_width']},
                                        style_header={
                                        'color': 'black',
                                        'font-size': '30px',
                                            },
                                        style_data={
                                        'backgroundColor': 'rgb(250, 250, 250)',
                                        'color': 'black',
                                        'font-size': '25px'}),],className="bg-light"),

                                    dbc.Card([html.H3("Model Parameters", style={
                                        'textAlign': 'left',
                                        'color': colors['text'],
                                        'font-size': '15',
                                        'padding-top': '50px',
                                        'padding-bottom': '20px',
                                        'padding-left': styles['padding-left']
                                        }),
                                        dash_table.DataTable(
                                        id='model_parameters',
                                        columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self.model_parameters.keys())],
                                        data=pd.DataFrame(self.model_parameters,index=[0]).to_dict('records'),
                                        editable=True,
                                        style_table={'overflowX': 'scroll', 'padding-left': '20px','padding-bottom':'30px', 'width': styles['table_width']},
                                        style_header={
                                        'color': 'black',
                                        'font-size': '30px',
                                            },
                                        style_data={
                                        'backgroundColor': 'rgb(250, 250, 250)',
                                        'color': 'black',
                                        'font-size': '25px'}),],className="bg-light"),
                                    
                                    dbc.Card([html.H3("Initial Conditions", style={
                                        'textAlign': 'left',
                                        'color': colors['text'],
                                        'font-size': '15',
                                        'padding-top': '50px',
                                        'padding-bottom': '20px',
                                        'padding-left': styles['padding-left']
                                        }),
                                        dash_table.DataTable(
                                        id='initial_conditions',
                                        columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self._ic.keys())],
                                        data=pd.DataFrame(self._ic,index=[0]).to_dict('records'),
                                        editable=True,
                                        style_table={'overflowX': 'scroll', 'padding-left': '20px','padding-bottom':'30px', 'width': styles['table_width']},
                                        style_header={
                                        'color': 'black',
                                        'font-size': '30px',
                                            },
                                        style_data={
                                        'backgroundColor': 'rgb(250, 250, 250)',
                                        'color': 'black',
                                        'font-size': '25px'}),],className="bg-light"),

                                    dbc.Card([html.H3("Inlet Conditions", style={
                                        'textAlign': 'left',
                                        'color': colors['text'],
                                        'font-size': '15',
                                        'padding-top': '50px',
                                        'padding-bottom': '20px',
                                        'padding-left': styles['padding-left']
                                        }),
                                        dash_table.DataTable(
                                        id='inlet_conditions',
                                        columns=[{"name": i, "id": i,"type":"numeric"} for i in list(self._inc.keys())],
                                        data=pd.DataFrame(self._inc,index=[0]).to_dict('records'),
                                        editable=True,
                                        style_table={'overflowX': 'scroll', 'padding-left': '20px','padding-bottom':'30px', 'width': styles['table_width']},
                                        style_header={
                                        'color': 'black',
                                        'font-size': '30px',
                                            },
                                        style_data={
                                        'backgroundColor': 'rgb(250, 250, 250)',
                                        'color': 'black',
                                        'font-size': '25px'}),],className="bg-light"),
                                        ],className="bg-light")],fluid=True,className="bg-light",style={"width": styles['container_width']}),
                                    dbc.Container([dbc.Row(
                                    [
                                    html.H2("Escher Map", style={
                                    'textAlign': 'left',
                                    'color': colors['text'],
                                    'font-size': '15',
                                    'padding-top': '20px',
                                    'padding-bottom': '20px',
                                    'padding-left': styles['padding-left']
                                    }) ,
            
                                    dcc.Dropdown(["Show Map","Hide Map"],
                                     self.reactions[0], style={"width": "300px","font-size":25,'padding-left':'2-px'}, id="Drop_Down_Escher"),
                                    html.Div(children=None,id="Escher_",style={"height": "100px",'padding-buttom':'20px'}),
                                    ])], fluid=True,className="bg-light pb-3",style={"width": styles['container_width']}),
            dbc.Container(html.Div(children=None,id="Escher",style={'align':'center'}),fluid=True,className="bg-light pb-3",style={"width": styles['container_width']}),
        ]
        if escher_map is None:
            page.pop(-1)
            page.pop(-1)
            page.pop(-1)

        
        app.layout = html.Div(page)

        @app.callback(Output(component_id="Escher_", component_property='children'), Input(component_id="Drop_Down_Escher", component_property='value'))
        def escher_wrapper(drop_down_escher):
            print("drop_down_escher")
            if drop_down_escher=="Show Map":
                Labels={}
                for i in range(0,self.sim_time,int(self.sim_time/20)):
                    Labels[i]={'label':str(i),'style':{'color': '#77b0b1'}}
                Labels[self.sim_time]=self.sim_time
                return [html.H2("Time (Day)",style={'textAlign': 'center'}),dcc.Slider(0,self.sim_time,int(self.sim_time/20),value=0,id="Escher_Slider",marks=None,tooltip={"placement": "bottom", "always_visible": True})]

        @app.callback(Output(component_id="Escher", component_property='children'), Input(component_id="Drop_Down_Escher", component_property='value'),
        Input(component_id="Escher_Slider", component_property='value'),prevent_initial_call=True)        
        def draw_escher(drop_down_escher,escher_slider):
            rxn_data={}
            self.ode_system(0,sol.y[:,int(sol.y.shape[1]/self.sim_time*escher_slider)],self)
            fluxes=self.info["Fluxes"]
            for ind,i in enumerate(self.reactions):
                rxn_data[i.replace(" ","_")]= fluxes[ind]
            if kwargs.get('min_flux',None):
                min_scale={ 'type': 'value','value':kwargs.get('min_flux') , 'color': 'red','size':10 }
            else:
                min_scale={ 'type': 'min' , 'color': 'red','size':10 }
            if kwargs.get('max_flux',None):
                max_scale={ 'type': 'value','value':kwargs.get('max_flux') , 'color': 'green','size':10 }
            else:
                max_scale={ 'type': 'max', 'color': 'green','size':10 }
                
            if drop_down_escher=="Show Map":
                return [dash_escher.DashEscher(mapData=escher_map,modelData=cobra_model,
            options={
             'reaction_data':rxn_data,
             'enable_keys':False,
             'reaction_scale':[min_scale,max_scale],
            }
            ,height='1000px',
        width='100%')
             ]
        @app.callback(Output(component_id='Concentrations_Line_Plot', component_property='figure'),
                    Input(component_id='base_parameters', component_property='data'),
                    Input(component_id='model_parameters', component_property='data'),
                    Input(component_id='initial_conditions', component_property='data'),
                    Input(component_id='inlet_conditions', component_property='data'),
                    prevent_initial_call=True
                    )
        def update_graph_fig(base_parameters: dict, model_parameters:dict, initial_conditions: dict, inlet_conditions: dict)->plotly.graph_objects.Figure:
            
            if len(self.control_state.keys()):
                for i in self.control_state.keys():
                    self.control_state[i]=initial_conditions[0][i]
            if len(base_parameters):
                self.base_parameters = base_parameters[0]
            if len(model_parameters):
                self.model_parameters = model_parameters[0]
            self.initial_conditions = np.array(
            [initial_conditions[0][i] for i in self.species])[:, np.newaxis]
            self.inlet_conditions = np.array(
            [inlet_conditions[0][i+"_in"] for i in self.species])[:, np.newaxis]
            update_sol = self.solve_model(np.linspace(0, self.sim_time, 10000))

            sol=update_sol
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
            legend=dict(font=dict(size= 20),),
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",

                )
            fig.update_xaxes(
            title={
            "text": "Time (Days)",
            "font_size": 25,
                },
                 tickfont_size=20,
            linecolor='grey',
            gridcolor='grey',
                )
            fig.update_yaxes(
            title={
            "text": "Concentrations (kg COD/m^3)",
            "font_size": 25,
             },
            tickfont_size=20,
            linecolor='grey',
            gridcolor='grey',

            
                )
            fig.update_traces(line=dict(width=3))
            return fig
            


        app.run_server(**kwargs)

    def csv_report(self,sol: scipy.integrate._ivp.ivp.OdeResult ,address: str)->None:
        """Converts the results to a pandas data frame then to a csv"""
        df = pd.DataFrame(sol.y, columns=sol.t, index=self.species)
        df.to_csv(os.path.join(address,self.name+"_Report.csv"), header=True,
                  index=True)
        
    def copy(self):
        """Returns a copy of the model"""
        return type(self)(model_parameters=self.model_parameters.copy(),
                          base_parameters=self.base_parameters.copy(),
                          initial_conditions=self._ic.copy(),
                          inlet_conditions=self._inc.copy(),
                          feed=self.feed,
                          reactions=self.reactions.copy(),
                          species=self.species.copy(),
                          ode_system=self.ode_system,
                          build_stoichiometric_matrix=self.build_stoichiometric_matrix,
                          control_state=self.control_state.copy(),
                          name=self.name,
                          switch=self.switch,
                          time_limit=self.time_limit,
                          simulation_time=self.sim_time)
    
    def build_cobra_model(self,address:str=None):
        """This method builds a cobra model from an instance of Model. One particular use
        of such models is to build an escher map from the model.
        Args:
            address (str, optional): The address to save the model. Defaults to None.
        """
        try:
            import cobra
        except ImportError:
            raise ImportError("CobraPy is not installed, please install it to use this function")
        model = cobra.Model(self.name)
        for reaction in self.reactions:
            temp_reaction = cobra.Reaction(reaction.replace(" ", "_"), name=reaction.replace(" ", "_"))
            temp_mets = np.where(self.s[:, self.reactions.index(reaction)] != 0)
            met_dict = {}
            for met in temp_mets[0]:
                metabolite = cobra.Metabolite(self.species[met].replace(" ", "_"),
                                              name=self.species[met].replace(" ", "_"), compartment="Model")
                met_dict[metabolite] = self.s[met, self.reactions.index(reaction)]
            temp_reaction.add_metabolites(met_dict)
            model.add_reactions([temp_reaction])
        if address:
            cobra.io.save_json_model(model, address)
        return model



def build_adm1_stoiciometric_matrix(base_parameters: dict, model_parameters: dict, reactons: list, species:list,feed:Feed,nitrogen_limited:bool=False)-> np.ndarray:
    """This function builds the stoichiometric matrix for the ADM1 Model.
    Args:
        base_parameters (dict): a dictionary containing the base parameters
        model_parameters (dict): a dictionary containing the model parameters
        reactons (list): a list containing all reactions
        species (list): a list containing all species
        feed (Feed): a Feed instance which contains the feed information
        nitrogen_limited (bool, optional): A boolean which indicates whether the model is nitrogen limited. Defaults to False.
    
    Returns:
        np.ndarray: Returns the stoichiometric matrix of the ADM1 model.
    """

    S = np.zeros((len(species), len(reactons)))
    S[0, [1, 3, 4]] = [1, (1-model_parameters["f_fa_li"]), - 1]
    S[1, [2, 5]] = [1, -1]
    S[2, [3, 6]] = [(model_parameters["f_fa_li"]), - 1]
    Y_aa=0 if nitrogen_limited else model_parameters['Y_aa']
    S[3, [5, 7]] = [(1-Y_aa) *
                    model_parameters['f_va_aa'], - 1]
    Y_su=0 if nitrogen_limited else model_parameters['Y_su']
    S[4, [4, 5, 8]] = [(1-Y_su)*model_parameters['f_bu_su'],
                       (1-Y_aa)*model_parameters["f_bu_aa"], - 1]
    S[5, [4, 5, 7, 9]] = [(1-model_parameters["Y_su"])*model_parameters['f_pro_su'],
                          (1-Y_aa)*model_parameters["f_pro_aa"], (1 - model_parameters['Y_c4'])*0.54, -1]
    
    Y_fa=0 if nitrogen_limited else model_parameters['Y_fa'] 
    S[6, [4, 5, 6, 7, 8, 9, 10]] = [(1-Y_su)*model_parameters['f_ac_su'],
                                    (1-Y_aa) *
                                    model_parameters['f_ac_aa'],
                                    (1-Y_fa)*0.7,
                                    (1-model_parameters['Y_c4'])*0.31,
                                    (1-model_parameters['Y_c4'])*0.8,
                                    (1-model_parameters['Y_pro'])*0.57,
                                    -1]
    S[7, [4, 5, 6, 7, 8, 9, 11, 25]] = [(1-Y_su)*model_parameters['f_h2_su'],
                                        (1-Y_aa) *
                                        model_parameters['f_h2_aa'],
                                        (1-Y_fa)*0.3,
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
    s_5 = (-1 * model_parameters['C_su'] + (1 - Y_su) * (model_parameters['f_bu_su'] * model_parameters['C_bu'] + model_parameters['f_pro_su']
                                                                             * model_parameters['C_pro'] + model_parameters['f_ac_su'] * model_parameters['C_ac']) + Y_su * model_parameters['C_bac'])
    s_6 = (-1 * model_parameters['C_aa'] + (1 - Y_aa) * (model_parameters['f_va_aa'] * model_parameters['C_va'] + model_parameters['f_bu_aa'] * model_parameters['C_bu'] +
                                                                             model_parameters['f_pro_aa'] * model_parameters['C_pro'] + model_parameters['f_ac_aa'] * model_parameters['C_ac']) + Y_aa * model_parameters['C_bac'])
    s_7 = (-1 * model_parameters['C_fa'] + (1 - Y_fa) * 0.7 *
           model_parameters['C_ac'] + Y_fa * model_parameters['C_bac'])
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
                                                                        -Y_su*model_parameters['N_bac'],
                                                                        model_parameters['N_aa']-Y_aa *
                                                                        model_parameters['N_bac'],
                                                                        -Y_fa*model_parameters['N_bac'],
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
    S[16, [4, 12]] = [Y_su, -1]
    S[17, [5, 13]] = [Y_aa, -1]
    S[18, [6, 14]] = [Y_fa, -1]
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


def adm1_ode_sys(t: float, c: np.ndarray, model:Model)-> np.ndarray:
    """ The ODE system for the original ADM.
        No testing is done.

        Args:
            t (float):a matrix of zeros to be filled
            c (np.ndarray): an array of concentrations to be filled
            Model (Model): The an instance of Model to calculate ODE with

        Returns:
            np.ndarray: The output is dCdt, the change of concentration with respect to time.
    """
    c[34] = c[10] - c[33]
    c[32] = c[9] - c[31]
    I_pH_aa = (model.model_parameters["K_pH_aa"] ** model.model_parameters['nn_aa'])/(np.power(
        c[26], model.model_parameters['nn_aa']) + np.power(model.model_parameters["K_pH_aa"], model.model_parameters['nn_aa']))
    I_pH_ac = (model.model_parameters['K_pH_ac'] ** model.model_parameters["n_ac"])/(
        c[26] ** model.model_parameters['n_ac'] + model.model_parameters['K_pH_ac'] ** model.model_parameters['n_ac'])
    I_pH_h2 = (model.model_parameters['K_pH_h2']**model.model_parameters['n_h2'])/(
        c[26] ** model.model_parameters['n_h2'] + model.model_parameters['K_pH_h2']**model.model_parameters['n_h2'])
    I_IN_lim = 1 / (1+(model.model_parameters['K_S_IN'] / c[10]))
    I_h2_fa = 1 / (1+(c[7] / model.model_parameters['K_I_h2_fa']))
    I_h2_c4 = 1 / (1+(c[7]/model.model_parameters['K_I_h2_c4']))
    I_h2_pro = (1/(1+(c[7]/model.model_parameters['K_I_h2_pro'])))
    I_nh3 = 1/(1+(c[33]/model.model_parameters['K_I_nh3']))
    I5 = (I_pH_aa * I_IN_lim)
    I6 = np.copy(I5)
    I7 = (I_pH_aa * I_IN_lim * I_h2_fa)
    I8 = (I_pH_aa * I_IN_lim * I_h2_c4)
    I9 = np.copy(I8)
    I10 = (I_pH_aa * I_IN_lim * I_h2_pro)
    I11 = (I_pH_ac * I_IN_lim * I_nh3)
    I12 = (I_pH_h2 * I_IN_lim)
    v = np.zeros((len(model.reactions), 1))
    v[0] = model.model_parameters["k_dis"]*c[12]
    
    v[1] = model.model_parameters['k_hyd_ch']*c[13]
    v[2] = model.model_parameters['k_hyd_pr']*c[14]
    v[3] = model.model_parameters['k_hyd_li']*c[15]
    
    v[4] = model.model_parameters['k_m_su']*c[0] / \
(model.model_parameters['K_S_su']+c[0])*c[16]*I5
    v[5] = model.model_parameters['k_m_aa']*c[1] / \
        (model.model_parameters['K_S_aa']+c[1])*c[17]*I6
    v[6] = model.model_parameters['k_m_fa']*c[2] / \
        (model.model_parameters['K_S_fa']+c[2])*c[18]*I7
    v[7] = model.model_parameters['k_m_c4']*c[3] / \
        (model.model_parameters['K_S_c4']+c[3]) * \
        c[19]*c[3]/(c[3]+c[4]+10 ** (-6))*I8
    v[8] = model.model_parameters['k_m_c4']*c[4] / \
        (model.model_parameters['K_S_c4']+c[4]) * \
        c[19]*c[4]/(c[4]+c[3]+10 ** (-6))*I9
    v[9] = model.model_parameters['k_m_pr']*c[5] / \
        (model.model_parameters['K_S_pro']+c[5])*c[20]*I10
    v[10] = model.model_parameters['k_m_ac']*c[6] / \
        (model.model_parameters['K_S_ac']+c[6])*c[21]*I11
    v[11] = model.model_parameters['k_m_h2']*c[7] / \
        (model.model_parameters['K_S_h2']+c[7])*c[22]*I12
    v[12] = model.model_parameters['k_dec_X_su']*c[16]
    v[13] = model.model_parameters['k_dec_X_aa']*c[17]
    v[14] = model.model_parameters['k_dec_X_fa']*c[18]
    v[15] = model.model_parameters['k_dec_X_c4']*c[19]
    v[16] = model.model_parameters['k_dec_X_pro']*c[20]
    v[17] = model.model_parameters['k_dec_X_ac']*c[21]
    v[18] = model.model_parameters['k_dec_X_h2']*c[22]
    v[19] = model.model_parameters['k_A_B_va'] * \
        (c[27] * (model.model_parameters['K_a_va'] + c[26]) -
         model.model_parameters['K_a_va'] * c[3])
    v[20] = model.model_parameters['k_A_B_bu'] * \
        (c[28] * (model.model_parameters['K_a_bu'] + c[26]) -
         model.model_parameters['K_a_bu'] * c[4])
    v[21] = model.model_parameters['k_A_B_pro'] * \
        (c[29] * (model.model_parameters['K_a_pro'] + c[26]) -
         model.model_parameters['K_a_pro'] * c[5])
    v[22] = model.model_parameters['k_A_B_ac'] * \
        (c[30] * (model.model_parameters['K_a_ac'] + c[26]) -
         model.model_parameters['K_a_ac'] * c[6])
    v[23] = model.model_parameters['k_A_B_co2'] * \
        (c[31] * (model.model_parameters['K_a_co2'] + c[26]) -
         model.model_parameters['K_a_co2'] * c[9])
    v[24] = model.model_parameters['k_A_B_IN'] * \
        (c[33] * (model.model_parameters['K_a_IN'] + c[26]) -
         model.model_parameters['K_a_IN'] * c[10])
    p_gas_h2 = c[35] * model.base_parameters["R"] * \
        model.base_parameters["T_op"] / 16
    p_gas_ch4 = c[36] * model.base_parameters["R"] * \
        model.base_parameters["T_op"] / 64
    p_gas_co2 = c[37] * model.base_parameters["R"] * \
        model.base_parameters["T_op"]
    p_gas_h2o = 0.0313 * \
        np.exp(5290 *
               (1 / model.base_parameters["T_base"] - 1 / model.base_parameters["T_op"]))
    P_gas = p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o
    q_gas = max(
        0, (model.model_parameters['k_p'] * (P_gas - model.base_parameters['P_atm'])))
    v[25] = model.model_parameters['k_L_a'] * \
        (c[7] - 16 * model.model_parameters['K_H_h2'] * p_gas_h2)
    v[26] = model.model_parameters['k_L_a'] * \
        (c[8] - 64 * model.model_parameters['K_H_ch4'] * p_gas_ch4)
    v[27] = model.model_parameters['k_L_a'] * \
        (c[32] - model.model_parameters['K_H_co2'] * p_gas_co2)
    dCdt = np.matmul(model.s, v)
    
    if c[model.species.index('S_IN')]<0.01:
        model.nitrogen_limited=True
    else:
        model.nitrogen_limited=False
    
    phi = c[24]+c[34]-c[31] - (c[30] / 64) - (c[29] / 112) - (c[28] / 160) - (c[27] / 208) - c[25]
    c[26] = (-1 * phi / 2) + (0.5 * np.sqrt(phi**2 + 4 * model.model_parameters['K_w']))
    
    dCdt[0: 35] = dCdt[0: 35]+model.base_parameters['q_in'] / model.base_parameters["V_liq"] * \
        (model.inlet_conditions[0: 35]-c[0:35].reshape(-1, 1))
    
        
    dCdt[35:] = dCdt[35:]+q_gas/model.base_parameters["V_gas"] * (model.inlet_conditions[35:]-c[35:].reshape(-1, 1))
    dCdt[[26, 32, 34], 0] = 0
    if model.switch == "DAE":
        dCdt[7] = 0
        dCdt[27: 32] = 0
        dCdt[33] = 0
    
    if model.control_state.keys():
        for state in model.control_state.keys():
            c[model.species.index(state)]=model.control_state[state]
            dCdt[model.species.index(state)]=0
    
    
    return dCdt[:, 0]


def build_e_adm_2_stoichiometric_matrix(base_parameters: dict,
                                             model_parameters: dict,
                                             reactions: list,
                                             species: list,
                                             feed:Feed,
                                             nitrogen_limited:bool=False)->np.ndarray:
    """ 
    This function builds the stoichiometric matrix for e-ADM2 Model.
        
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
      reactions.index('TSS_Disintegration')] = [-1,feed.ch_tss, feed.prot_tss, feed.lip_tss, feed.xi_tss]
    S[list(map(species.index, ["TDS", "X_ch", "X_pr", "X_li", "S_I"])), reactions.index('TDS_Disintegration')] = [-1,
                                                                                                                  feed.ch_tds, feed.prot_tds, feed.lip_tds, feed.si_tds]
    S[list(map(species.index, ["X_ch", "S_su"])),
      reactions.index('Hydrolysis carbohydrates')] = [-1, 1]
    S[list(map(species.index, ["X_pr", "S_aa"])),
      reactions.index('Hydrolysis proteins')] = [-1, 1]
    S[list(map(species.index, ["X_li", "S_fa"])),
      reactions.index('Hydrolysis lipids')] = [-1, 1]
    
    Y_su=0 if nitrogen_limited else model_parameters['Y_su']
    f_ac_su=1-model_parameters['f_pro_su']-model_parameters['f_et_su']-model_parameters['f_lac_su']
    f_IC_su = -(-model_parameters['C_su'] +
                (1-Y_su)*model_parameters['f_pro_su']*model_parameters['C_pro'] +
                (1-Y_su)*model_parameters['f_et_su']*model_parameters['C_et'] +
                (1-Y_su)*model_parameters['f_lac_su']*model_parameters['C_lac'] +
                (1-Y_su)*f_ac_su*model_parameters['C_ac'] +
                Y_su*model_parameters['C_bac'])


    S[list(map(species.index, ["S_su", "S_pro", "S_et", "S_lac", "S_ac", "S_IN", "S_IC", "X_su"])),
      reactions.index('Uptake of sugars')] = [-1,
                                              (1-Y_su) * model_parameters['f_pro_su'],
                                              (1-Y_su) * model_parameters['f_et_su'],
                                              (1-Y_su) * model_parameters['f_lac_su'],
                                              (1-Y_su) * f_ac_su,
                                              -model_parameters['N_bac']*Y_su,
                                              f_IC_su,
                                              Y_su]
      
    Y_aa=0 if nitrogen_limited else model_parameters['Y_aa']
    f_ac_aa=1-model_parameters['f_pro_aa']-model_parameters['f_et_aa']-model_parameters['f_lac_aa']
    f_IC_aa = -(-model_parameters['C_aa'] +
                (1-Y_aa)*model_parameters['f_pro_aa']*model_parameters['C_pro'] +
                (1-Y_aa)*model_parameters['f_et_aa']*model_parameters['C_et'] +
                (1-Y_aa)*model_parameters['f_lac_aa']*model_parameters['C_lac'] +
                (1-Y_aa)*f_ac_aa*model_parameters['C_ac'] +
                Y_aa*model_parameters['C_bac'])



    S[list(map(species.index, ["S_aa", "S_pro", "S_et", "S_lac", "S_ac", "S_IN", "S_IC", "X_aa"])),
      reactions.index('Uptake of amino acids')] = [-1,
                                                   (1-Y_aa) * model_parameters['f_pro_aa'],
                                                   (1-Y_aa) * model_parameters['f_et_aa'],
                                                   (1-Y_aa) * model_parameters['f_lac_aa'],
                                                   (1-Y_aa) * f_ac_aa,
                                                   model_parameters['N_aa']-Y_aa * model_parameters['N_bac'],
                                                   f_IC_aa,
                                                   Y_aa]
      
    Y_fa=0 if nitrogen_limited else model_parameters['Y_fa']
    f_ac_fa=1-model_parameters['f_pro_fa']-model_parameters['f_et_fa']-model_parameters['f_lac_fa']
    f_IC_fa = -(-model_parameters['C_fa']+
                (1-Y_fa)*model_parameters['f_pro_fa']*model_parameters['C_pro'] +
                (1-Y_fa)*model_parameters['f_et_fa']*model_parameters['C_et'] +
                (1-Y_fa)*model_parameters['f_lac_fa']*model_parameters['C_lac'] +
                (1-Y_fa)*f_ac_fa*model_parameters['C_ac'] +
                Y_fa*model_parameters['C_bac'])
    # if f_IC_fa<0:
    #     raise ValueError("f_IC_fa is negative") 

    S[list(map(species.index, ["S_fa", "S_pro", "S_et", "S_lac", "S_ac", "S_IN", "S_IC", "X_fa"])),
      reactions.index('Uptake of LCFA')] = [-1,
                                            (1-Y_fa) * model_parameters['f_pro_fa'],
                                            (1-Y_fa) * model_parameters['f_et_fa'],
                                            (1-Y_fa) * model_parameters['f_lac_fa'],
                                            (1-Y_fa) * f_ac_fa,
                                            -Y_fa * model_parameters['N_bac'],
                                            f_IC_fa,
                                            Y_fa]
    if any([f_ac_fa<0,f_ac_aa<0,f_ac_su<0]):
        raise ValueError("f_ac is negative")
    Y_ac_et=0 if nitrogen_limited else model_parameters['Y_ac_et']
    Y_ac_lac=0 if nitrogen_limited else model_parameters['Y_ac_lac']
    f_IC_ac_et = -(-model_parameters['C_ac'] +
                    model_parameters['f_et_ac']*model_parameters['C_et'] +
                   (1-model_parameters['f_et_ac']-Y_ac_et) * model_parameters['f_bu_ac']*model_parameters['C_bu'] +
                   Y_ac_et*model_parameters['C_bac'])

    f_IC_ac_lac = -(-model_parameters['C_ac'] +
                    model_parameters['f_lac_ac']*model_parameters['C_lac'] +
                    (1-model_parameters['f_lac_ac']-Y_ac_lac) * model_parameters['f_bu_ac']*model_parameters['C_bu'] +
                    Y_ac_lac*model_parameters['C_bac'])


    S[list(map(species.index, ["S_ac", "S_et", "S_bu", "S_IN", "S_IC", "S_h2", "X_ac_et"])),
      reactions.index('Uptake of acetate_et')] = [-1,
                                                  model_parameters['f_et_ac'],
                                                  (1- model_parameters['f_et_ac']-model_parameters['Y_ac']) * model_parameters['f_bu_ac'],
                                                  -Y_ac_et * model_parameters['N_bac'],
                                                  f_IC_ac_et,
                                                  (1- model_parameters['f_et_ac']-Y_ac_et) * (1-model_parameters['f_bu_ac']),
                                                  Y_ac_et]

    S[list(map(species.index, ["S_ac", "S_lac", "S_bu", "S_IN", "S_IC", "S_h2", "X_ac_lac"])),
        reactions.index('Uptake of acetate_lac')] = [-1,
                                                    model_parameters['f_lac_ac'],
                                                     (1-model_parameters['f_lac_ac']-Y_ac_lac) * model_parameters['f_bu_ac'],
                                                     -Y_ac_lac * model_parameters['N_bac'],
                                                     f_IC_ac_lac,
                                                     (1-model_parameters['f_lac_ac']-Y_ac_lac) * (1-model_parameters['f_bu_ac']),
                                                     Y_ac_lac]
    
    Y_pro_et=0 if nitrogen_limited else model_parameters['Y_pro_et']
    Y_pro_lac=0 if nitrogen_limited else model_parameters['Y_pro_et']
    
    f_IC_pro_et = -(-model_parameters['C_pro'] +
                    model_parameters['f_et_pro']*model_parameters['C_et'] +
                    (1-model_parameters['f_et_pro']-Y_pro_et)*model_parameters['f_va_pro']*model_parameters['C_va'] +
                    (Y_pro_et)*model_parameters['C_bac'])

    f_IC_pro_lac = -(-model_parameters['C_pro'] +
                     model_parameters['f_lac_pro']*model_parameters['C_lac'] +
                     (1-model_parameters['f_lac_pro']-Y_pro_lac)*model_parameters['f_va_pro']*model_parameters['C_va'] +
                     (Y_pro_lac)*model_parameters['C_bac'])
    


    
    S[list(map(species.index, ["S_pro", "S_et", "S_va", "S_IN", "S_IC", "S_h2", "X_chain_et"])),
      reactions.index('Uptake of propionate_et')] = [-1,
                                                    model_parameters['f_et_pro'],
                                                     (1-model_parameters['f_et_pro']-Y_pro_et) * model_parameters['f_va_pro'],
                                                     -Y_pro_et *  model_parameters['N_bac'],
                                                     f_IC_pro_et,
                                                     (1-model_parameters['f_et_pro']-Y_pro_et) * (1-model_parameters['f_va_pro']),
                                                     model_parameters['Y_chain_et_pro']]

    S[list(map(species.index, ["S_pro", "S_lac", "S_va", "S_IN", "S_IC", "S_h2", "X_chain_lac"])),
        reactions.index('Uptake of propionate_lac')] = [-1,
                                                        model_parameters['f_lac_pro'],
                                                        (1-model_parameters['f_lac_pro']-Y_pro_lac) * model_parameters['f_va_pro'],
                                                        -Y_pro_lac * model_parameters['N_bac'],
                                                        f_IC_pro_lac,
                                                        (1-model_parameters['f_lac_pro']-Y_pro_lac) * (1-model_parameters['f_va_pro']),
                                                        model_parameters['Y_chain_lac_pro']]

    Y_bu_et=0 if nitrogen_limited else model_parameters['Y_bu_et']
    Y_bu_lac=0 if nitrogen_limited else model_parameters['Y_bu_lac']
    f_IC_bu_et = -(-model_parameters['C_bu'] +
                    model_parameters['f_et_bu']*model_parameters['C_et'] +
                   (1-model_parameters['f_et_bu']-Y_bu_et)*model_parameters['f_cap_bu']*model_parameters['C_cap'] +
                   (Y_bu_et)*model_parameters['C_bac'])

    f_IC_bu_lac = -(-model_parameters['C_bu'] +
                    model_parameters['f_lac_bu']*model_parameters['C_lac'] +
                    (1-model_parameters['f_lac_bu']-Y_bu_lac)*model_parameters['f_cap_bu']*model_parameters['C_cap'] +
                    (Y_bu_lac)*model_parameters['C_bac'])
    

    S[list(map(species.index, ["S_bu", "S_et", "S_cap", "S_IN", "S_IC", "S_h2", "X_chain_et"])),
        reactions.index('Uptake of butyrate_et')] = [-1,
                                                     model_parameters['f_et_bu'],
                                                     (1-model_parameters['f_et_bu']-Y_bu_et) * model_parameters['f_cap_bu'],
                                                     -Y_bu_et * model_parameters['N_bac'],
                                                     f_IC_bu_et,
                                                     (1-model_parameters['f_et_bu']-Y_bu_et)*(1-model_parameters['f_cap_bu']),
                                                     Y_bu_et]

    S[list(map(species.index, ["S_bu", "S_lac", "S_cap", "S_IN", "S_IC", "S_h2", "X_chain_lac"])),
        reactions.index('Uptake of butyrate_lac')] = [-1,
                                                      model_parameters['f_lac_bu'],
                                                      (1- model_parameters['f_lac_bu']-Y_bu_lac) * model_parameters['f_cap_bu'],
                                                      -Y_bu_lac *model_parameters['N_bac'],
                                                      f_IC_bu_lac,
                                                      (1- model_parameters['f_lac_bu']-Y_bu_lac)*(1-model_parameters['f_cap_bu']),
                                                      Y_bu_lac]


    Y_va=0 if nitrogen_limited else model_parameters['Y_va']
                
    S[list(map(species.index, ["S_va", "S_pro", "X_VFA_deg"])),
        reactions.index('Uptake of valerate')] = [-1,
                                                  (1-Y_va),
                                                  Y_va,
                                                  ]

    Y_cap=0 if nitrogen_limited else model_parameters['Y_cap']
    S[list(map(species.index, ["S_cap", "S_ac", "X_VFA_deg"])),
        reactions.index('Uptake of caproate')] = [-1,
                                                  (1 - Y_cap),
                                                  Y_cap]
        
    Y_cap=0 if nitrogen_limited else model_parameters['Y_bu']
    S[list(map(species.index, ["S_bu", "S_ac", "X_VFA_deg"])),
        reactions.index('Uptake of butyrate')] = [-1,
                                                  (1 - Y_cap),
                                                  Y_cap]
        
        
    
    Y_Me_ac=0 if nitrogen_limited else model_parameters["Y_Me_ac"]
    f_IC_Me_ach2 =0
    S[list(map(species.index, ["S_gas_h2", "S_ac", "S_ch4", "X_Me_ac", 'S_IC', 'S_IN'])),
        reactions.index('Methanogenessis from acetate and h2')] = [-1,
                                                                   model_parameters['f_ac_h2'],
                                                                   (1 +model_parameters['f_ac_h2']- Y_Me_ac),
                                                                   Y_Me_ac,
                                                                   f_IC_Me_ach2,
                                                                    -Y_Me_ac *model_parameters['N_bac']
                                                                   ]
    
    Y_Me_CO2=0 if nitrogen_limited else model_parameters["Y_Me_CO2"]

    
    S[list(map(species.index, ["S_gas_h2", "S_gas_ch4", "X_Me_CO2", 'S_gas_co2',"S_IN"])),
        reactions.index('Methanogenessis from CO2 and h2')] = [-1,
                                                               (1 -model_parameters['f_co2_ch4']- Y_Me_CO2),
                                                               (Y_Me_CO2),
                                                               model_parameters['f_co2_ch4'],
                                                                -Y_Me_CO2 *model_parameters['N_bac']
                                                                ]
    
    
    
    Y_ac_et_ox=0 if nitrogen_limited else model_parameters["Y_ac_et_ox"]
    f_IC_et_ox=-(-model_parameters['C_et'] +
                    (1-Y_ac_et_ox)*model_parameters['C_bac']
                    +Y_ac_et_ox*model_parameters['C_ac'])

    S[list(map(species.index, ["S_et", "X_et","S_ac","S_IC"])),
        reactions.index('Uptake of ethanol')] = [-1,1-Y_ac_et_ox,Y_ac_et_ox,f_IC_et_ox]

    
    Y_pro_lac_ox=0 if nitrogen_limited else model_parameters['Y_pro_lac_ox']
    f_IC_lac_ox=-(-model_parameters['C_lac'] +
                (1-Y_pro_lac_ox)*model_parameters['C_bac']
                +Y_pro_lac_ox*model_parameters['C_pro'])
    
    S[list(map(species.index, ["S_lac" ,"S_pro","X_lac","S_IC"])),
        reactions.index('Uptake of lactate')] = [-1, 1-Y_pro_lac_ox,Y_pro_lac_ox,f_IC_lac_ox]

    S[list(map(species.index, ["X_su", "TSS","S_IN","S_IC"])),
        reactions.index('Decay of Xsu')] = [-1, 1,model_parameters['N_bac'],model_parameters['C_bac']]

    S[list(map(species.index, ["X_aa", "TSS","S_IN","S_IC"])),
        reactions.index('Decay of Xaa')] = [-1, 1,model_parameters['N_bac'],model_parameters['C_bac']]

    S[list(map(species.index, ["X_fa", "TSS","S_IN","S_IC"])),
        reactions.index('Decay of Xfa')] = [-1, 1,model_parameters['N_bac'],model_parameters['C_bac']]

    S[list(map(species.index, ["X_ac_et", "TSS","S_IN","S_IC"])),
        reactions.index('Decay of X_ac_et')] = [-1, 1,model_parameters['N_bac'],model_parameters['C_bac']]

    S[list(map(species.index, ["X_ac_lac", "TSS","S_IN","S_IC"])),
        reactions.index('Decay of X_ac_lac')] = [-1, 1,model_parameters['N_bac'],model_parameters['C_bac']]

    S[list(map(species.index, ["X_chain_et", "TSS", "S_IN","S_IC"])),
        reactions.index('Decay of X_chain_et')] = [-1, 1,model_parameters['N_bac'],model_parameters['C_bac']]

    S[list(map(species.index, ["X_chain_lac", "TSS", "S_IN","S_IC"])),
        reactions.index('Decay of X_chain_lac')] = [-1, 1,model_parameters['N_bac'],model_parameters['C_bac']]

    S[list(map(species.index, ["X_VFA_deg", "TSS", "S_IN","S_IC"])),
        reactions.index('Decay of X_VFA_deg')] = [-1, 1,model_parameters['N_bac'],model_parameters['C_bac']]

    S[list(map(species.index, ["X_Me_ac", "TSS", "S_IN","S_IC"])),
        reactions.index('Decay of X_Me_ac')] = [-1, 1,model_parameters['N_bac'],model_parameters['C_bac']]

    S[list(map(species.index, ["X_Me_CO2", "TSS", "S_IN","S_IC"])),
        reactions.index('Decay of X_Me_CO2')] = [-1, 1,model_parameters['N_bac'],model_parameters['C_bac']]

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

    S[list(map(species.index, ["S_hco3_ion"])),  # I don't think this is right، should look at the reaction in ADM1
        reactions.index('Acid Base Equilibrium (CO2)')] = [-1]

    S[list(map(species.index, ["S_nh3", "S_nh4_ion"])),
        reactions.index('Acid Base Equilibrium (In)')] = [-1, 1]  # I don't think this is right، should look at the reaction in ADM1

    S[list(map(species.index, ["S_h2", "S_gas_h2"])),
        reactions.index('Gas Transfer H2')] = [-base_parameters['V_liq']/base_parameters['V_gas'], 1]
    S[list(map(species.index, ["S_ch4", "S_gas_ch4"])),
        reactions.index('Gas Transfer CH4')] = [-base_parameters['V_liq']/base_parameters['V_gas'], 1]
    S[list(map(species.index, ["S_co2", "S_gas_co2"])),
        reactions.index('Gas Transfer CO2')] = [-base_parameters['V_liq']/base_parameters['V_gas'], 1]
    
    return S


def build_e_adm_stoiciometric_matrix(base_parameters: dict,
                                     model_parameters: dict,
                                     reactions: list,
                                     species: list,
                                     feed:Feed,
                                     nitrogen_limited:bool=False)->np.ndarray:
    """ 
    This function builds the stoichiometric matrix for the e_ADM Model.
        
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
      reactions.index('TSS_Disintegration')] = [-1,feed.ch_tss, feed.prot_tss, feed.lip_tss,feed.xi_tss]
    S[list(map(species.index, ["TDS", "X_ch", "X_pr", "X_li", "S_I"])), reactions.index('TDS_Disintegration')] = [-1,
                                                                                                            feed.ch_tds, feed.prot_tds, feed.lip_tds, feed.si_tds]
    S[list(map(species.index, ["X_ch", "S_su"])),reactions.index('Hydrolysis carbohydrates')] = [-1, 1]
    S[list(map(species.index, ["X_pr", "S_aa"])),reactions.index('Hydrolysis proteins')] = [-1, 1]
    S[list(map(species.index, ["X_li", "S_fa"])),reactions.index('Hydrolysis lipids')] = [-1, 1]

    f_IC_su_et=-(-model_parameters['C_su']+
               (1-model_parameters['Y_su_et']) * model_parameters['C_et']+
               (1-model_parameters['Y_su_et']) * model_parameters['C_bac']
              )
    
    f_IC_su_lac=-(-model_parameters['C_su']+
              (1-model_parameters['Y_su_lac']) * model_parameters['C_lac']+
              (1-model_parameters['Y_su_lac']) * model_parameters['C_bac']
              )
    
    f_IC_su_ac=-(-model_parameters['C_su']+
               (1-model_parameters['Y_su_ac']) * model_parameters['C_ac']+
               (1-model_parameters['Y_su_ac']) * model_parameters['C_bac']
              )
    f_IC_su_pro=-(-model_parameters['C_su']+
                (1-model_parameters['Y_su_pro']) * model_parameters['C_pro']+
                (1-model_parameters['Y_su_pro']) * model_parameters['C_bac']
              )
    
    S[list(map(species.index, ["S_su","S_et","S_IN","S_IC","X_su"])),
     reactions.index('Su_to_et')] = [-1,
                                         (1-model_parameters['Y_su_et']),
                                           -model_parameters['N_bac']* model_parameters['Y_su_et'],
                                            f_IC_su_et,
                                            model_parameters['Y_su_et']]
    
    S[list(map(species.index, ["S_su","S_lac","S_IN","S_IC","X_su"])),
     reactions.index('Su_to_lac')] = [-1,
                                                (1-model_parameters['Y_su_lac']),
                                                -model_parameters['N_bac']* model_parameters['Y_su_lac'],
                                                f_IC_su_lac,
                                                model_parameters['Y_su_lac']]
    
    S[list(map(species.index, ["S_su","S_ac","S_IN","S_IC","X_su"])),
      reactions.index('Su_to_ac')] = [-1,
                                    (1-model_parameters['Y_su_ac']),
                                    -model_parameters['N_bac']* model_parameters['Y_su_ac'],
                                    f_IC_su_ac,
                                    model_parameters['Y_su_ac']]
      
    S[list(map(species.index, ["S_su","S_pro","S_IN","S_IC","X_su"])),
            reactions.index('Su_to_pro')] = [-1,
                                            (1-model_parameters['Y_su_pro']),
                                              -model_parameters['N_bac']* model_parameters['Y_su_pro'],
                                              f_IC_su_pro,
                                              model_parameters['Y_su_pro']]

    f_IC_aa_lac=-(-model_parameters['C_aa']+
              (1-model_parameters['Y_aa_lac']) * model_parameters['C_lac']+
              (1-model_parameters['Y_aa_lac']) * model_parameters['C_bac']
              )
    
    f_IC_aa_ac=-(-model_parameters['C_aa']+
              (1-model_parameters['Y_aa_ac']) * model_parameters['C_ac']+
              (1-model_parameters['Y_aa_ac']) * model_parameters['C_bac']
              )
    
    f_IC_aa_pro=-(-model_parameters['C_aa']+
              (1-model_parameters['Y_aa_pro']) * model_parameters['C_pro']+
              (1-model_parameters['Y_aa_pro']) * model_parameters['C_bac']
              )


    S[list(map(species.index, ["S_aa","S_lac","S_IN", "S_IC", "X_aa"])),
      reactions.index('aas_to_lac')] = [-1,
                                             (1-model_parameters['Y_aa_lac']),
                                             model_parameters['N_aa']- model_parameters['Y_aa_lac'] * model_parameters['N_bac'],
                                             f_IC_aa_lac,
                                             model_parameters['Y_aa_lac']]
    
    S[list(map(species.index, ["S_aa","S_pro","S_IN", "S_IC", "X_aa"])),
      reactions.index('aas_to_pro')] = [-1,
                                             (1-model_parameters['Y_aa_pro']),
                                             model_parameters['N_aa']- model_parameters['Y_aa_pro'] * model_parameters['N_bac'],
                                             f_IC_aa_pro,
                                             model_parameters['Y_aa_pro']]
      
    
    S[list(map(species.index, ["S_aa","S_ac","S_IN", "S_IC", "X_aa"])),
      reactions.index('aas_to_ac')] = [-1,
                                             (1-model_parameters['Y_aa_ac']),
                                             model_parameters['N_aa']- model_parameters['Y_aa_ac'] * model_parameters['N_bac'],
                                             f_IC_aa_ac,
                                             model_parameters['Y_aa_ac']]
      
    Y_fa=0 if nitrogen_limited else model_parameters['Y_fa']
    f_IC_fa = -(-model_parameters['C_fa'] +
                (1-Y_fa)*model_parameters['f_pro_fa']*model_parameters['C_pro'] +
                (1-Y_fa)*model_parameters['f_ac_fa']*model_parameters['C_ac'] +
                (1-Y_fa)*model_parameters['C_bac'])

    S[list(map(species.index, ["S_fa", "S_pro", "S_ac", "S_IN", "S_IC", "X_fa"])),
      reactions.index('Uptake of LCFA')] = [-1,
                                            (1-Y_fa) * model_parameters['f_pro_fa'],
                                            (1-Y_fa) * model_parameters['f_ac_fa'],
                                              -Y_fa * model_parameters['N_bac'],
                                              f_IC_fa,
                                              Y_fa]
#HERE
    Y_ac_et=0 if nitrogen_limited else model_parameters['Y_ac_et']
    Y_ac_lac=0 if nitrogen_limited else model_parameters['Y_ac_lac']
    f_IC_ac_et = -((-1-(1-Y_ac_et) * model_parameters['f_et_ac'])*model_parameters['C_ac'] +
                   (1-Y_ac_et)* model_parameters['f_et_ac']*model_parameters['C_et'] +
                   (1-Y_ac_et) * model_parameters['f_bu_ac']*model_parameters['C_bu'] +
                   (1-Y_ac_et)* model_parameters['C_bac'])  

    f_IC_ac_lac = -((-1-(1-Y_ac_lac) * model_parameters['f_lac_ac'])*model_parameters['C_ac'] +
                    (1-Y_ac_lac)* model_parameters['f_lac_ac']* model_parameters['C_lac'] +
                    (1-Y_ac_lac)* model_parameters['f_bu_ac']* model_parameters['C_bu'] +
                    (1-Y_ac_lac)* model_parameters['C_bac'])

    S[list(map(species.index, ["S_ac", "S_et", "S_bu", "S_IN", "S_IC", "S_h2", "X_ac_et"])),
      reactions.index('Uptake of acetate_et')] = [-1-(1-Y_ac_et) * model_parameters['f_et_ac'],
                                                  (1-Y_ac_et) * model_parameters['f_et_ac'],
                                                  (1-model_parameters['Y_ac']) * model_parameters['f_bu_ac'],
                                                  -Y_ac_et * model_parameters['N_bac'],
                                                  f_IC_ac_et,
                                                  (1-Y_ac_et) * (1-model_parameters['f_bu_ac']),
                                                  Y_ac_et]

    S[list(map(species.index, ["S_ac", "S_lac", "S_bu", "S_IN", "S_IC", "S_h2", "X_ac_lac"])),
        reactions.index('Uptake of acetate_lac')] = [-1-(1-Y_ac_lac) * model_parameters['f_lac_ac'],
                                                     (1-Y_ac_lac) * model_parameters['f_lac_ac'],
                                                     (1-Y_ac_lac) * model_parameters['f_bu_ac'],
                                                     -Y_ac_lac * model_parameters['N_bac'], 
                                                     f_IC_ac_lac,
                                                     (1-Y_ac_lac) * (1-model_parameters['f_bu_ac']),
                                                     Y_ac_lac]

    Y_pro_et=0 if nitrogen_limited else model_parameters['Y_pro_et']
    Y_pro_lac=0 if nitrogen_limited else model_parameters['Y_pro_lac']
    
    f_IC_pro_et = -((-1-(1-Y_pro_et) * model_parameters['f_et_pro'])*model_parameters['C_pro'] +
                    (1-Y_pro_et)*model_parameters['f_et_pro']*model_parameters['C_et'] +
                    (1-Y_pro_et)*model_parameters['f_va_pro']*model_parameters['C_va'] +
                    (1-Y_pro_et)*model_parameters['C_bac'])

    f_IC_pro_lac = -((-1-(1-Y_pro_lac) * model_parameters['f_lac_pro'])*model_parameters['C_pro'] +
                     (1-Y_pro_lac)*model_parameters['f_lac_pro']*model_parameters['C_lac'] +
                     (1-Y_pro_lac)*model_parameters['f_va_pro']*model_parameters['C_va'] +
                     (1-Y_pro_lac)*model_parameters['C_bac'])

    S[list(map(species.index, ["S_pro", "S_et", "S_va","S_IC","S_IN","S_h2", "X_chain_et"])),
      reactions.index('Uptake of propionate_et')] = [-1-(1-model_parameters['Y_chain_et_pro']) * model_parameters['f_et_pro'],
                                                     (1-model_parameters['Y_chain_et_pro']) * model_parameters['f_et_pro'],
                                                     (1-model_parameters['Y_chain_et_pro']) * model_parameters['f_va_pro'],
                                                     f_IC_pro_et,
                                                     -model_parameters['Y_chain_et_pro'] * model_parameters['N_bac'],
                                                     (1-model_parameters['Y_chain_et_pro']) * (1-model_parameters['f_va_pro']),
                                                     model_parameters['Y_chain_et_pro']]

    S[list(map(species.index, ["S_pro", "S_lac", "S_va", "S_IC", "S_IN", "S_h2", "X_chain_lac"])),
        reactions.index('Uptake of propionate_lac')] = [-1-(1-model_parameters['Y_chain_lac_pro']) * model_parameters['f_lac_pro'],
                                                        (1-model_parameters['Y_chain_lac_pro']) * model_parameters['f_lac_pro'],
                                                        (1-model_parameters['Y_chain_lac_pro']) * model_parameters['f_va_pro'],
                                                        f_IC_pro_lac,
                                                        -model_parameters['Y_chain_lac_pro'] * model_parameters['N_bac'],
                                                        (1-model_parameters['Y_chain_lac_pro']) * (1-model_parameters['f_va_pro']),
                                                        model_parameters['Y_chain_lac_pro']]

    Y_bu_et=0 if nitrogen_limited else model_parameters['Y_bu_et']
    Y_pro_lac=0 if nitrogen_limited else model_parameters['Y_pro_lac']
    
    f_IC_bu_et = -((-1-(1-Y_bu_et) * model_parameters['f_et_bu'])*model_parameters['C_bu'] +
                   (1-Y_bu_et)*model_parameters['f_et_bu']*model_parameters['C_et'] +
                   (1-Y_bu_et)*model_parameters['f_cap_bu']*model_parameters['C_cap'] +
                   (1-Y_bu_et)*model_parameters['C_bac'])

    f_IC_bu_lac = -((-1-(1-Y_pro_lac) * model_parameters['f_lac_bu'])*model_parameters['C_bu'] +
                    (1-Y_pro_lac)*model_parameters['f_lac_bu']*model_parameters['C_lac'] +
                    (1-Y_pro_lac)*model_parameters['f_cap_bu']*model_parameters['C_cap'] +
                    (1-Y_pro_lac)*model_parameters['C_bac'])

    S[list(map(species.index, ["S_bu", "S_et", "S_cap", "S_IC", "S_IN", "S_h2", "X_chain_et"])),
        reactions.index('Uptake of butyrate_et')] = [-1-(1-Y_bu_et) * model_parameters['f_et_bu'],
                                                     (1-Y_bu_et) * model_parameters['f_et_bu'],
                                                     (1-Y_bu_et) * model_parameters['f_cap_bu'],
                                                     f_IC_bu_et,
                                                     -Y_bu_et * model_parameters['N_bac'],
                                                     (1-Y_bu_et)*(1-model_parameters['f_cap_bu']),
                                                     Y_bu_et]

    S[list(map(species.index, ["S_bu", "S_lac", "S_cap", "S_IC", "S_IN", "S_h2", "X_chain_lac"])),
        reactions.index('Uptake of butyrate_lac')] = [-1-(1-Y_pro_lac) * model_parameters['f_lac_bu'],
                                                      (1-Y_pro_lac) * model_parameters['f_lac_bu'],
                                                      (1-Y_pro_lac) * model_parameters['f_cap_bu'],
                                                      f_IC_bu_lac,
                                                      -Y_pro_lac * model_parameters['N_bac'],
                                                      (1-Y_pro_lac)*(1-model_parameters['f_cap_bu']),
                                                      Y_pro_lac]

    Y_va=0 if nitrogen_limited else model_parameters['Y_va']
    Y_cap=0 if nitrogen_limited else model_parameters['Y_cap']
    S[list(map(species.index, ["S_va", "S_pro", "X_VFA_deg"])),
        reactions.index('Uptake of valerate')] = [-1,
                                                  (1-Y_va),
                                                  Y_va]

    S[list(map(species.index, ["S_cap", "S_ac", "X_VFA_deg"])),
        reactions.index('Uptake of caproate')] = [-1,
                                                  (1 - Y_cap),
                                                  Y_cap]
    
    S[list(map(species.index, ["S_bu", "S_ac", "X_VFA_deg"])),
        reactions.index('Uptake of butyrate')] = [-1,
                                                  (1 - model_parameters['Y_bu']),
                                                  model_parameters['Y_bu']]

    Y_Me_ac=0 if nitrogen_limited else model_parameters["Y_Me_ac"]
    f_IC_Me_ach2 =0
    f_IC_Me_ach2 = -((1 - model_parameters['Y_h2_ac'])*model_parameters['f_ac_h2']*model_parameters['C_ac']+
                     (1 -Y_Me_ac)*model_parameters['C_ch4']+
                     Y_Me_ac*model_parameters['C_bac'])
                     
        
    
    S[list(map(species.index, ["S_h2", "S_ac", "S_ch4", "X_Me_ac", 'S_IC'])),
        reactions.index('Methanogenessis from acetate and h2')] = [-1-(1 - model_parameters['Y_h2_ac'])*model_parameters['f_ac_h2'],
                    (1 - model_parameters['Y_h2_ac'])*model_parameters['f_ac_h2'],
                    (1 -model_parameters['Y_h2_ac']),
                    model_parameters['Y_h2_ac'],
                    f_IC_Me_ach2]

    f_IC_Me_CO2h2 = -(model_parameters['Y_h2_CO2']*model_parameters['C_ch4'] +
                      model_parameters['Y_h2_CO2']*model_parameters['C_bac'])
    
    S[list(map(species.index, ["S_h2", "S_ch4", "X_Me_CO2", 'S_IC'])),
        reactions.index('Methanogenessis from CO2 and h2')] = [-1,
                                                               (1 - model_parameters['Y_h2_CO2']),
                                                               (model_parameters['Y_h2_CO2']),
                                                               f_IC_Me_CO2h2]

    Y_ac_et_ox=0 if nitrogen_limited else model_parameters["Y_ac_et_ox"]
    
    f_IC_et_ox=-(-model_parameters['C_et'] +
                    (1-Y_ac_et_ox)*model_parameters['C_bac']
                    +Y_ac_et_ox*model_parameters['C_ac'])

    S[list(map(species.index, ["S_et", "X_et","S_ac","S_IC"])),
        reactions.index('Uptake of ethanol')] = [-1,Y_ac_et_ox,(1-Y_ac_et_ox),f_IC_et_ox]

    
    Y_pro_lac_ox=0 if nitrogen_limited else model_parameters['Y_pro_lac_ox']
    f_IC_lac_ox=-(-model_parameters['C_lac'] +
                (1-Y_pro_lac_ox)*model_parameters['C_bac']
                +Y_pro_lac_ox*model_parameters['C_pro'])
    
    S[list(map(species.index, ["S_lac", "X_lac","S_pro","S_IC"])),
        reactions.index('Uptake of lactate')] = [-1, Y_pro_lac_ox,(1-Y_pro_lac_ox),f_IC_lac_ox]

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

    S[list(map(species.index, ["S_va_ion","S_va"])),
        reactions.index('Acid Base Equilibrium (Va)')] = [-1,1]

    S[list(map(species.index, ["S_bu_ion","S_bu"])),
        reactions.index('Acid Base Equilibrium (Bu)')] = [-1,1]

    S[list(map(species.index, ["S_pro_ion","S_pro"])),
        reactions.index('Acid Base Equilibrium (Pro)')] = [-1,1]

    S[list(map(species.index, ["S_cap_ion","S_cap"])),
        reactions.index('Acid Base Equilibrium (Cap)')] = [-1,1]

    S[list(map(species.index, ["S_lac_ion","S_lac"])),
        reactions.index('Acid Base Equilibrium (Lac)')] = [-1,1]

    S[list(map(species.index, ["S_ac_ion","S_ac"])),
        reactions.index('Acid Base Equilibrium (Ac)')] = [-1,1]

    S[list(map(species.index, ["S_co2", "S_hco3_ion"])),  # I don't think this is right، should look at the reaction in ADM1
        reactions.index('Acid Base Equilibrium (CO2)')] = [-1, 1]

    S[list(map(species.index, ["S_nh3", "S_nh4_ion"])),
        reactions.index('Acid Base Equilibrium (In)')] = [-1, 1]  # I don't think this is right، should look at the reaction in ADM1

    S[list(map(species.index, ["S_h2", "S_gas_h2"])),
        reactions.index('Gas Transfer H2')] = [-base_parameters['V_liq']/base_parameters['V_gas'], 1]
    S[list(map(species.index, ["S_ch4", "S_gas_ch4"])),
        reactions.index('Gas Transfer CH4')] = [-base_parameters['V_liq']/base_parameters['V_gas'], 1]
    S[list(map(species.index, ["S_co2", "S_gas_co2"])),
        reactions.index('Gas Transfer CO2')] = [-base_parameters['V_liq']/base_parameters['V_gas'], 1]
    return S

def e_adm_2_ode_sys(t: float, c: np.ndarray, model: Model)-> np.ndarray:
    """
    This function is used to build the ODEs of the e-adm2 model.
    
    Args:
        t (float):a matrix of zeros to be filled
        c (np.ndarray): an array of concentrations to be filled
        Model (Model): The model to calculate ODE with

    Returns:
        np.ndarray: The output is dCdt, the change of concentration with respect to time. 
    """
    ### Initialize the ion concentrations
    # if t==0:
    if t==0:
        c[model.species.index('S_va_ion')]=model.model_parameters['K_a_va']/(model.model_parameters['K_a_va']+c[model.species.index('S_H_ion')])*c[model.species.index('S_va')]
        c[model.species.index('S_bu_ion')]=model.model_parameters['K_a_bu']/(model.model_parameters['K_a_bu']+c[model.species.index('S_H_ion')])*c[model.species.index('S_bu')]
        c[model.species.index('S_pro_ion')]=model.model_parameters['K_a_pro']/(model.model_parameters['K_a_pro']+c[model.species.index('S_H_ion')])*c[model.species.index('S_pro')]
        c[model.species.index('S_cap_ion')]=model.model_parameters['K_a_cap']/(model.model_parameters['K_a_cap']+c[model.species.index('S_H_ion')])*c[model.species.index('S_cap')]
        c[model.species.index('S_ac_ion')]=model.model_parameters['K_a_ac']/(model.model_parameters['K_a_ac']+c[model.species.index('S_H_ion')])*c[model.species.index('S_ac')]
        c[model.species.index('S_lac_ion')]=model.model_parameters['K_a_lac']/(model.model_parameters['K_a_lac']+c[model.species.index('S_H_ion')])*c[model.species.index('S_lac')]    
        c[model.species.index('S_hco3_ion')] = c[model.species.index('S_IC')] - c[model.species.index('S_co2')]
        phi=(model.model_parameters['K_w']/c[model.species.index('S_H_ion')]-c[model.species.index('S_H_ion')])
        c[model.species.index('S_anion')] = c[model.species.index('S_cation')]+c[model.species.index('S_nh4_ion')]-c[model.species.index('S_hco3_ion')]-(c[model.species.index('S_lac_ion')] / 88) - (c[model.species.index('S_ac_ion')] / 64) - (c[model.species.index('S_pro_ion')] /
                                                                                                                                                                     112) - (c[model.species.index('S_bu_ion')] / 160)-(c[model.species.index('S_cap_ion')] / 230) - (c[model.species.index('S_va_ion')] / 208)-phi

    c[model.species.index('S_hco3_ion')] = model.model_parameters['K_a_co2'] * c[model.species.index('S_IC')]/(model.model_parameters['K_a_co2'] + c[model.species.index('S_H_ion')])
    c[model.species.index('S_nh4_ion')]=  model.model_parameters['K_b_nh3'] * c[model.species.index('S_IN')]/(model.model_parameters['K_b_nh3'] + model.base_parameters['K_W'] / c[model.species.index('S_H_ion')])
    
    c[model.species.index('S_co2')]= c[model.species.index('S_IC')] -  c[model.species.index('S_hco3_ion')]
    c[model.species.index('S_nh3')]= c[model.species.index('S_IN')] - c[model.species.index('S_nh4_ion')]
        
    if (time.time()-model._be_time )>model.time_limit and model.time_limit!=-1:
        raise Exception("Time limit exceeded")

        
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
    
    I_IN_lim = 1 / (1+(c[model.species.index('S_IN')] / (model.model_parameters['K_S_IN']+10**-9)))
    
    I_h2_fa = 1 /  (1+(c[model.species.index('S_h2')] /(model.model_parameters['K_I_h2_fa']+10**-9)))

    I_h2_c4 = 1 /  (1+(c[model.species.index('S_h2')] / (model.model_parameters['K_I_h2_c4']+10**-9)))

    I_h2_pro = 1/  (1+(c[model.species.index('S_h2')] / (model.model_parameters['K_I_h2_pro']+10**-9)))

    I_nh3 =    1/  (1+(c[model.species.index('S_nh3')] / (model.model_parameters['K_I_nh3']+10**-9)))

    I_h2_oxidation=1/(1+(c[model.species.index('S_h2')] / (model.model_parameters['K_I_h2_ox']+10**-9)))

    I5 =    max(0,(I_pH_aa * I_IN_lim))
    I6 =    max(0,I5)
    I7 =    max(0,(I_pH_aa * I_IN_lim * I_h2_fa))
    I8 =    max(0,(I_pH_aa * I_IN_lim * I_h2_c4))
    I9 =    max(0,I8)
    I10 =   max(0,(I_pH_pro * I_IN_lim * I_h2_pro))
    I11 =   max(0,(I_pH_ac * I_IN_lim * I_nh3))
    I12 =   max(0,(I_pH_h2 * I_IN_lim))
    I13 =   max(0,(I_pH_cap * I_IN_lim * I_h2_c4))
    I14 =   max(0,(I_pH_bu * I_IN_lim * I_h2_c4))
    I15 =   max(0,(I_pH_va * I_IN_lim * I_h2_c4))
    I16 =   max(0,I_IN_lim * I_nh3*I_pH_aa*I_h2_oxidation)

    v = np.zeros((len(model.reactions), 1))

    v[model.reactions.index('TSS_Disintegration')] = model.model_parameters["k_dis_TSS"]*c[model.species.index('TSS')]

    v[model.reactions.index('TDS_Disintegration')] = model.model_parameters["k_dis_TDS"]*c[model.species.index('TDS')]

    v[model.reactions.index('Hydrolysis carbohydrates')] = model.model_parameters['k_hyd_ch']*c[model.species.index('X_ch')]

    v[model.reactions.index('Hydrolysis proteins')] = model.model_parameters['k_hyd_pr']*c[model.species.index('X_pr')]

    v[model.reactions.index('Hydrolysis lipids')] = model.model_parameters['k_hyd_li']*c[model.species.index('X_li')]

    v[model.reactions.index('Uptake of sugars')] = model.model_parameters['k_m_su']*c[model.species.index('S_su')] / \
        (model.model_parameters['K_S_su']+c[model.species.index('S_su')])*c[model.species.index('X_su')]*I5

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
    
    v[model.reactions.index('Uptake of butyrate')] = model.model_parameters['k_m_bu_deg']*c[model.species.index('S_bu')] / \
        (model.model_parameters['K_S_bu']+c[model.species.index('S_bu')]
         )*c[model.species.index('X_VFA_deg')]*I13

    v[model.reactions.index('Methanogenessis from acetate and h2')] = model.model_parameters['k_m_h2_Me_ac']*c[model.species.index('S_gas_h2')]*c[model.species.index('S_ac')] / \
        (model.model_parameters['K_S_h2_Me_ac']*c[model.species.index('S_gas_h2')]+model.model_parameters['K_S_ac_Me']*c[model.species.index(
            'S_ac')]+c[model.species.index('S_ac')]*c[model.species.index('S_gas_h2')]+10**-9)*c[model.species.index('X_Me_ac')]*I12

    v[model.reactions.index('Methanogenessis from CO2 and h2')] = model.model_parameters['k_m_h2_Me_CO2']*c[model.species.index('S_gas_h2')]*c[model.species.index('S_gas_co2')] / \
        (model.model_parameters['K_S_h2_Me_CO2']*c[model.species.index('S_gas_h2')]+model.model_parameters['K_S_CO2_Me']*c[model.species.index(
            'S_gas_co2')]+c[model.species.index('S_gas_co2')]*c[model.species.index('S_gas_h2')]+10**-9)*c[model.species.index('X_Me_CO2')]*I12


    v[model.reactions.index('Uptake of ethanol')] = model.model_parameters['k_m_et']*c[model.species.index('S_et')] / \
        (model.model_parameters['K_S_et']+c[model.species.index('S_et')]
         )*c[model.species.index("X_et")]*I16

    v[model.reactions.index('Uptake of lactate')] = model.model_parameters['k_m_lac']*c[model.species.index('S_lac')] / \
        (model.model_parameters['K_S_lac']+c[model.species.index('S_lac')]
         )*c[model.species.index('X_lac')]*I16

    v[model.reactions.index('Decay of Xsu')] = model.model_parameters['k_dec_X_su']*c[model.species.index('X_su')]
    v[model.reactions.index('Decay of Xaa')] = model.model_parameters['k_dec_X_aa']*c[model.species.index('X_aa')]
    v[model.reactions.index('Decay of Xfa')] = model.model_parameters['k_dec_X_fa']*c[model.species.index('X_fa')]
    v[model.reactions.index('Decay of X_ac_et')] = model.model_parameters['k_dec_X_ac']*c[model.species.index('X_ac_et')]
    v[model.reactions.index('Decay of X_ac_lac')] = model.model_parameters['k_dec_X_ac']*c[model.species.index('X_ac_lac')]
    v[model.reactions.index('Decay of X_chain_et')] = model.model_parameters['k_dec_X_chain_et']*c[model.species.index('X_chain_et')]
    v[model.reactions.index('Decay of X_chain_lac')] = model.model_parameters['k_dec_X_chain_lac']*c[model.species.index('X_chain_lac')]
    v[model.reactions.index('Decay of X_VFA_deg')] = model.model_parameters['k_dec_X_VFA_deg']*c[model.species.index('X_VFA_deg')]
    v[model.reactions.index('Decay of X_Me_ac')] = model.model_parameters['k_dec_X_Me_ac']*c[model.species.index('X_Me_ac')]
    v[model.reactions.index('Decay of X_Me_CO2')] = model.model_parameters['k_dec_X_Me_CO2']*c[model.species.index('X_Me_CO2')]
    v[model.reactions.index('Decay of Xet')] = model.model_parameters['k_dec_X_et']*c[model.species.index('X_et')]
    v[model.reactions.index('Decay of Xlac')] = model.model_parameters['k_dec_X_lac']*c[model.species.index('X_lac')]
    
    
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

    
    p_gas_h2 = c[model.species.index('S_gas_h2')] * model.base_parameters["R"] * model.base_parameters["T_op"] / 16
    p_gas_ch4 = c[model.species.index('S_gas_ch4')] * model.base_parameters["R"] * model.base_parameters["T_op"] / 64
    p_gas_co2 = c[model.species.index('S_gas_co2')] * model.base_parameters["R"] * model.base_parameters["T_op"]
    p_gas_h2o = 0.0313 * np.exp(5290 *(1 / model.base_parameters["T_base"] - 1 / model.base_parameters["T_op"]))
    
    P_gas = p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o
    
    q_gas = max(0, (model.model_parameters['k_p'] * (P_gas - model.base_parameters['P_atm'])))
    
    v[model.reactions.index('Gas Transfer H2')] = max(0,model.model_parameters['k_L_a'] * (c[model.species.index('S_h2')] - 16 *model.model_parameters['K_H_h2'] * p_gas_h2))
    v[model.reactions.index('Gas Transfer CH4')] = max(0,model.model_parameters['k_L_a'] * (c[model.species.index('S_ch4')] - 64 * model.model_parameters['K_H_ch4'] * p_gas_ch4))
    v[model.reactions.index('Gas Transfer CO2')] = max(0,model.model_parameters['k_L_a'] * (c[model.species.index('S_co2')] - model.model_parameters['K_H_co2'] * p_gas_co2))

    if c[model.species.index('S_IN')]<0.01:
        model.nitrogen_limited=True
    else:
        model.nitrogen_limited=False
        
    dCdt = np.matmul(model.s, v)
    phi = c[model.species.index('S_cation')]+c[model.species.index('S_nh4_ion')]-c[model.species.index('S_hco3_ion')]-(c[model.species.index('S_lac_ion')] / 88) - \
    (c[model.species.index('S_ac_ion')] / 64) - (c[model.species.index('S_pro_ion')] / 112) - (c[model.species.index('S_bu_ion')] / 160)-(c[model.species.index('S_cap_ion')] / 230) - (c[model.species.index('S_va_ion')] / 208) - c[model.species.index('S_anion')]
    
    if 'S_H_ion' in model.control_state.keys():
        c[model.species.index('S_H_ion')]=model.control_state['S_H_ion']
    else:
        c[model.species.index('S_H_ion')] = (-1 * phi / 2) + (0.5 * np.sqrt(phi**2 + 4 * model.model_parameters['K_w']))

    dCdt[0: len(model.species)-3] = dCdt[0: len(model.species)-3]+model.base_parameters['q_in'] / model.base_parameters["V_liq"] * (model.inlet_conditions[0: len(model.species)-3]-c[0: len(model.species)-3].reshape(-1, 1))

    dCdt[len(model.species)-3:] = dCdt[len(model.species)-3:]+q_gas/model.base_parameters["V_gas"] * (model.inlet_conditions[len(model.species)-3:]-c[len(model.species)-3:].reshape(-1, 1))

    dCdt[[model.species.index('S_H_ion'), model.species.index('S_co2'), model.species.index('S_nh4_ion')], 0] = 0
    
    if model.switch == "DAE":
        dCdt[model.species.index('S_va_ion'):model.species.index('S_co2')] = 0
        dCdt[model.species.index('S_nh3')] = 0
        c[model.species.index('S_va_ion')]=model.model_parameters['K_a_va']/(model.model_parameters['K_a_va']+c[model.species.index('S_H_ion')])*c[model.species.index('S_va')]
        c[model.species.index('S_bu_ion')]=model.model_parameters['K_a_bu']/(model.model_parameters['K_a_bu']+c[model.species.index('S_H_ion')])*c[model.species.index('S_bu')]
        c[model.species.index('S_pro_ion')]=model.model_parameters['K_a_pro']/(model.model_parameters['K_a_pro']+c[model.species.index('S_H_ion')])*c[model.species.index('S_pro')]
        c[model.species.index('S_cap_ion')]=model.model_parameters['K_a_cap']/(model.model_parameters['K_a_cap']+c[model.species.index('S_H_ion')])*c[model.species.index('S_cap')]
        c[model.species.index('S_ac_ion')]=model.model_parameters['K_a_ac']/(model.model_parameters['K_a_ac']+c[model.species.index('S_H_ion')])*c[model.species.index('S_ac')]
        c[model.species.index('S_lac_ion')]=model.model_parameters['K_a_lac']/(model.model_parameters['K_a_lac']+c[model.species.index('S_H_ion')])*c[model.species.index('S_lac')]    
        c[model.species.index('S_hco3_ion')] = c[model.species.index('S_IC')] - c[model.species.index('S_co2')]


  

    if model.control_state.keys():
        for state in model.control_state.keys():
            c[model.species.index(state)]=model.control_state[state]
            dCdt[model.species.index(state)]=0
    
    model.info["Fluxes"]=v
    return dCdt[:, 0]

def e_adm_ode_sys(t: float, c: np.ndarray, model: Model)-> np.ndarray:
    """
    This function is used to build the ODEs of the e_adm model.
    
    Args:
        t (float):a matrix of zeros to be filled
        c (np.ndarray): an array of concentrations to be filled
        Model (Model): The model to calculate ODE with

    Returns:
        np.ndarray: The output is dCdt, the change of concentration with respect to time. 
    """
    c[c<0]=0
    c[model.species.index('S_H_ion')]=0.000001
    if model.switch == "DAE":
        
        c[model.species.index('S_va_ion')]=model.model_parameters['K_a_va']/(model.model_parameters['K_a_va']+c[model.species.index('S_H_ion')])*c[model.species.index('S_va')]
    
        c[model.species.index('S_bu_ion')]=model.model_parameters['K_a_bu']/(model.model_parameters['K_a_bu']+c[model.species.index('S_H_ion')])*c[model.species.index('S_bu')]
    
        c[model.species.index('S_pro_ion')]=model.model_parameters['K_a_pro']/(model.model_parameters['K_a_pro']+c[model.species.index('S_H_ion')])*c[model.species.index('S_pro')]
    
        c[model.species.index('S_cap_ion')]=model.model_parameters['K_a_cap']/(model.model_parameters['K_a_cap']+c[model.species.index('S_H_ion')])*c[model.species.index('S_cap')]
    
        c[model.species.index('S_ac_ion')]=model.model_parameters['K_a_ac']/(model.model_parameters['K_a_ac']+c[model.species.index('S_H_ion')])*c[model.species.index('S_ac')]
        
        c[model.species.index('S_lac_ion')]=model.model_parameters['K_a_lac']/(model.model_parameters['K_a_lac']+c[model.species.index('S_H_ion')])*c[model.species.index('S_lac')]
    else: 
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

    # I5 = (I_pH_aa * I_IN_lim)
    # I6 = I5.copy()
    # I7 = (I_pH_aa * I_IN_lim * I_h2_fa)
    # I8 = (I_pH_aa * I_IN_lim * I_h2_c4)
    # I9 = I8.copy()
    # I10 = (I_pH_pro * I_IN_lim * I_h2_pro)
    # I11 = (I_pH_ac * I_IN_lim * I_nh3)
    # I12 = (I_pH_h2 * I_IN_lim)
    # I13 = (I_pH_cap * I_IN_lim * I_h2_c4)
    # I14 = (I_pH_bu * I_IN_lim * I_h2_c4)
    # I15 = (I_pH_va * I_IN_lim * I_h2_c4)
    # I16 = I_IN_lim * I_nh3*I_pH_aa*I_h2_oxidation
    I5  = 1
    I6  = 1
    I7  = 1
    I8  = 1
    I9  = 1 #one
    I10 = 1
    I11 = 1
    I12 = 1
    I13 = 1
    I14 = 1
    I15 = 1
    I16 = 1

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

    v[model.reactions.index('Su_to_et')] = model.model_parameters['k_m_su_et']*c[model.species.index('S_su')] / \
        (model.model_parameters['K_S_su_et']+c[model.species.index('S_su')])*c[model.species.index('X_su')]*I5
    
    v[model.reactions.index('Su_to_lac')] = model.model_parameters['k_m_su_lac']*c[model.species.index('S_su')] / \
        (model.model_parameters['K_S_su_lac']+c[model.species.index('S_su')]
         )*c[model.species.index('X_su')]/(c[model.species.index('X_su')]+model.model_parameters['K_S_X_su_lac'])*I5    

    v[model.reactions.index('Su_to_ac')] = model.model_parameters['k_m_su_ac']*c[model.species.index('S_su')] / \
        (model.model_parameters['K_S_su_ac']+c[model.species.index('S_su')]
         )*c[model.species.index('X_su')]*I5
    
    v[model.reactions.index('Su_to_pro')] = model.model_parameters['k_m_su_pro']*c[model.species.index('S_su')] / \
        (model.model_parameters['K_S_su_pro']+c[model.species.index('S_su')]
         )*c[model.species.index('X_su')]/(c[model.species.index('X_su')]+model.model_parameters['K_S_X_su_pro'])*I5        
        
    
    v[model.reactions.index('aas_to_lac')] = model.model_parameters['k_m_aa_lac']*c[model.species.index('S_aa')] / \
        (model.model_parameters['K_S_aa_lac']+c[model.species.index('S_aa')]
         )*c[model.species.index('X_aa')]/(c[model.species.index('X_aa')]+model.model_parameters['K_S_X_aa_lac'])*I6

    v[model.reactions.index('aas_to_pro')] = model.model_parameters['k_m_aa_pro']*c[model.species.index('S_aa')] / \
        (model.model_parameters['K_S_aa_pro']+c[model.species.index('S_aa')]
         )*c[model.species.index('X_aa')]/(c[model.species.index('X_aa')]+model.model_parameters['K_S_X_aa_pro'])*I6
    
    v[model.reactions.index('aas_to_ac')] = model.model_parameters['k_m_aa_ac']*c[model.species.index('S_aa')] / \
        (model.model_parameters['K_S_aa_ac']+c[model.species.index('S_aa')]
         )*c[model.species.index('X_aa')]/(c[model.species.index('X_aa')]+model.model_parameters['K_S_X_aa_ac'])*I6

    v[model.reactions.index('Uptake of LCFA')] = model.model_parameters['k_m_fa']*c[model.species.index('S_fa')] / \
        (model.model_parameters['K_S_fa'] +
         c[model.species.index('S_fa')])*c[model.species.index('X_fa')]/(c[model.species.index('X_fa')]+model.model_parameters['K_S_X_fa'])*I7

    v[model.reactions.index('Uptake of acetate_et')] = model.model_parameters['k_m_ac_et']*c[model.species.index('S_ac')]*c[model.species.index('S_et')] / \
        (model.model_parameters['K_S_ac']*c[model.species.index('S_ac')]+model.model_parameters['K_S_ac_et']*c[model.species.index('S_et')]+c[model.species.index('S_ac')]*c[model.species.index('S_et')]+10**-9
         )*c[model.species.index('X_ac_et')]/(c[model.species.index('X_ac_et')]+model.model_parameters['K_S_X_ac_et'])*I11

    v[model.reactions.index('Uptake of acetate_lac')] = model.model_parameters['k_m_ac_lac']*c[model.species.index('S_ac')]*c[model.species.index('S_lac')] / \
        (model.model_parameters['K_S_ac']*c[model.species.index('S_ac')]+model.model_parameters['K_S_ac_lac']*c[model.species.index('S_lac')]+c[model.species.index('S_ac')]*c[model.species.index('S_lac')]+10**-9
         )*c[model.species.index('X_ac_lac')]/(c[model.species.index('X_ac_lac')]+model.model_parameters['K_S_X_ac_lac'])*I11

    v[model.reactions.index('Uptake of propionate_et')] = model.model_parameters['k_m_pro_et']*c[model.species.index('S_pro')]*c[model.species.index('S_et')] / \
        (model.model_parameters['K_S_pro']*c[model.species.index('S_pro')]+model.model_parameters['K_S_pro_et']*c[model.species.index('S_et')]+c[model.species.index('S_pro')]*c[model.species.index('S_et')]+10**-9
         )*c[model.species.index('X_chain_et')]/(c[model.species.index('X_chain_et')]+model.model_parameters['K_S_X_chain_et'])*I10

    v[model.reactions.index('Uptake of propionate_lac')] = model.model_parameters['k_m_pro_lac']*c[model.species.index('S_pro')]*c[model.species.index('S_lac')] / \
        (model.model_parameters['K_S_pro']*c[model.species.index('S_pro')]+model.model_parameters['K_S_pro_lac']*c[model.species.index('S_lac')]+c[model.species.index('S_pro')]*c[model.species.index('S_lac')]+10**-9
         )*c[model.species.index('X_chain_lac')]/(c[model.species.index('X_chain_lac')]+model.model_parameters['K_S_X_chain_lac'])*I10

    v[model.reactions.index('Uptake of butyrate_et')] = model.model_parameters['k_m_bu_et']*c[model.species.index('S_bu')]*c[model.species.index('S_et')] / \
        (model.model_parameters['K_S_bu']*c[model.species.index('S_bu')]+model.model_parameters['K_S_bu_et']*c[model.species.index('S_et')]+c[model.species.index('S_bu')]*c[model.species.index('S_et')]+10**-9
         )*c[model.species.index('X_chain_et')]*I14

    v[model.reactions.index('Uptake of butyrate_lac')] = model.model_parameters['k_m_bu_lac']*c[model.species.index('S_bu')]*c[model.species.index('S_lac')] / \
        (model.model_parameters['K_S_bu']*c[model.species.index('S_bu')]+model.model_parameters['K_S_bu_lac']*c[model.species.index('S_lac')]+c[model.species.index('S_bu')]*c[model.species.index('S_lac')]+10**-9
         )*c[model.species.index('X_chain_lac')]/(c[model.species.index('X_chain_lac')]+model.model_parameters['K_S_X_chain_lac'])*I14
    
    v[model.reactions.index('Uptake of butyrate')] = model.model_parameters['k_m_bu']*c[model.species.index('S_bu')]/ \
        (model.model_parameters['K_S_bu']+c[model.species.index('S_bu')])*c[model.species.index('X_VFA_deg')]/(c[model.species.index('X_VFA_deg')]+model.model_parameters['K_S_X_VFA_deg'])*I14
    
    v[model.reactions.index('Uptake of valerate')] = model.model_parameters['k_m_va']*c[model.species.index('S_va')] / \
        (model.model_parameters['K_S_va']+c[model.species.index('S_va')])*c[model.species.index('X_VFA_deg')]/(c[model.species.index('X_VFA_deg')]+model.model_parameters['K_S_X_VFA_deg'])*I15

    v[model.reactions.index('Uptake of caproate')] = model.model_parameters['k_m_cap']*c[model.species.index('S_cap')] / \
        (model.model_parameters['K_S_cap']+c[model.species.index('S_cap')])*c[model.species.index('X_VFA_deg')]/(c[model.species.index('X_VFA_deg')]+model.model_parameters['K_S_X_VFA_deg'])*I13

    v[model.reactions.index('Methanogenessis from acetate and h2')] = model.model_parameters['k_m_h2_Me_ac']*c[model.species.index('S_h2')]*c[model.species.index('S_ac')] / \
        (model.model_parameters['K_S_h2_Me_ac']*c[model.species.index('S_h2')]+model.model_parameters['K_S_ac_Me']*c[model.species.index(
            'S_ac')]+c[model.species.index('S_ac')]*c[model.species.index('S_h2')]+10**-9)*c[model.species.index('X_Me_ac')]*I12

    v[model.reactions.index('Methanogenessis from CO2 and h2')] = model.model_parameters['k_m_h2_Me_CO2']*c[model.species.index('S_h2')]*c[model.species.index('S_co2')] / \
        (model.model_parameters['K_S_h2_Me_CO2']*c[model.species.index('S_h2')]+model.model_parameters['K_S_CO2_Me']*c[model.species.index(
            'S_co2')]+c[model.species.index('S_co2')]*c[model.species.index('S_h2')]+10**-9)*c[model.species.index('X_Me_CO2')]*I12


    v[model.reactions.index('Uptake of ethanol')] = model.model_parameters['k_m_et']*c[model.species.index('S_et')] / \
        (model.model_parameters['K_S_et']+c[model.species.index('S_et')]
         )*c[model.species.index("X_et")]/(c[model.species.index("X_et")]+model.model_parameters['K_S_X_et'])*I16

    v[model.reactions.index('Uptake of lactate')] = model.model_parameters['k_m_lac']*c[model.species.index('S_lac')] / \
        (model.model_parameters['K_S_lac']+c[model.species.index('S_lac')]
         )*c[model.species.index('X_lac')]/(c[model.species.index('X_lac')]+model.model_parameters['K_S_X_lac'])*I16

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

    dCdt = np.matmul(model.s, v)

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
    
    if c[model.species.index('S_IN')]<0.01:
        model.nitrogen_limited=True
    else:
        model.nitrogen_limited=False
    
    if model.switch == "DAE":
        # dCdt[model.species.index('S_h2')] = 0

        dCdt[model.species.index('S_va_ion'):model.species.index('S_co2')] = 0

        dCdt[model.species.index('S_nh3')] = 0
    
    if model.control_state.keys():
        for state in model.control_state.keys():
            c[model.species.index(state)]=model.control_state[state]
            dCdt[model.species.index(state)]=0
    model.info["Fluxes"]=v
    return dCdt[:, 0]



if __name__ == "__main__":
    
    local_params=dict(
        model_parameters="/Users/parsaghadermarzi/Desktop/Academics/Projects/Anaerobic_Digestion_Modeling/16s_study/Database/ADM_Parameters/Modified_ADM_Model_Parameters.json",
        base_parameters="/Users/parsaghadermarzi/Desktop/Academics/Projects/Anaerobic_Digestion_Modeling/16s_study/Database/ADM_Parameters/Modified_ADM_Base_Parameters.json",
        initial_conditions="/Users/parsaghadermarzi/Desktop/Academics/Projects/Anaerobic_Digestion_Modeling/16s_study/Database/ADM_Parameters/Modified_ADM_Initial_Conditions.json",
        inlet_conditions="/Users/parsaghadermarzi/Desktop/Academics/Projects/Anaerobic_Digestion_Modeling/16s_study/Database/ADM_Parameters/Modified_ADM_Inlet_Conditions.json",
        species="/Users/parsaghadermarzi/Desktop/Academics/Projects/Anaerobic_Digestion_Modeling/16s_study/Database/ADM_Parameters/Modified_ADM_Species.json",
        reactions="/Users/parsaghadermarzi/Desktop/Academics/Projects/Anaerobic_Digestion_Modeling/16s_study/Database/ADM_Parameters/Modified_ADM_Reactions.json",
    )
    params=utils.load_multiple_json_files(local_params)
    # params_adm1=utils.load_multiple_json_files(configs.ADM1_LOCAL)
    # pd.DataFrame(params.reactions)
    # model=Model(model_parameters=params.model_parameters,
    #             base_parameters=params.base_parameters,
    #             initial_conditions=params.initial_conditions,
    #             build_stoichiometric_matrix=build_adm1_stoiciometric_matrix,
    #             ode_system=adm1_ode_sys,
    #             inlet_conditions=params.inlet_conditions,
    #             species=params.species,
    #             reactions=params.reactions,
    #             feed=Feed)
    # model.solve_model(t_eval=np.linspace(0, 30, 1000))
    params.initial_conditions['S_cation']=0.1
    mod_adm1 = Model(model_parameters=params.model_parameters,
                    base_parameters=params.base_parameters, 
                             initial_conditions=params.initial_conditions,
                             inlet_conditions=params.inlet_conditions,
                             feed=DEFAULT_FEED,
                             reactions=params.reactions,
                            species=params.species,
                            ode_system=e_adm_2_ode_sys,
                            build_stoichiometric_matrix=build_e_adm_2_stoichiometric_matrix,
                            control_state={
                             
                                },
                            name="Modified_ADM1",
                            switch="DAE",
                            metagenome_report=None)
    mod_adm1.control_state['S_H_ion']=0.000001
    mod_adm1.time_limit=100000
    
    mod_adm1.dash_app(mod_adm1.solve_model(t_eval=np.linspace(0, 30, 1000)))
    
    
    
