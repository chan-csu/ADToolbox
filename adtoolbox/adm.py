from typing import Callable
import plotly
import numpy as np
import scipy.optimize
import scipy.integrate
import json
import os
import plotly.express as px
import plotly.graph_objects as go
from core import Database as Database
from core import SeedDB as SeedDB
from core import Feed
import utils
from adtoolbox import PKG_DATA
import configs
import time
import warnings
import polars as pl

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
class _Fake_Sol:
    def __init__(self, y,t):
        self.y = y
        self.t=t


def _require_model_parameters(model_parameters: dict, required: list[str], model_name: str) -> None:
    missing = sorted(set(required) - set(model_parameters))
    if missing:
        raise ValueError(
            f"{model_name} model_parameters is missing required keys: {', '.join(missing)}"
        )


def _monod_limitation(concentration: float, half_saturation: float, eps: float = 10**-9) -> float:
    return concentration / (half_saturation + concentration + eps)


class Model:

    """Any kinetic model could be an instance of this class.
    Args:
        model_parameters (dict): a dictionary which contains model parameters
        base_parameters (dict): a dictionary which contains base paramters
        initial_conditions (dict): a dictionary containing initial conditions for all species
        inlet_conditions (dict): a dictionary containing inlet conditions for all species
        feed (Feed): a Feed instance which contains the feed information
        reactions (list): a list containing all types of reactions
        species (list): a list containing all species
        ode_system (Callable): a callable which outputs the ODE system compatible with Scipy.integrate.solve_ivp
        build_stoichiometric_matrix(Callable): a callable which builds the stoichiometric matrix
        control_state (dict, optional): a dictionary containing the states that are desired to be constant. Defaults to {}.
        name (str, optional): a name for the model. Defaults to "ADM".
        switch (str, optional): whether acid/base states are solved algebraically ("DAE") or kinetically. Defaults to "DAE".
        simulation_time (float, optional): default simulation duration in days. Defaults to 30.
        time_limit (float, optional): wall-clock limit for a solve, or -1 for no limit. Defaults to -1.
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
        self._s_cache={}   # stoichiometric matrix cache, keyed by nitrogen_limited

    @property
    def s(self):
        """Stoichiometric matrix, cached per solve.

        The matrix depends only on the parameters, feed, and the binary
        ``nitrogen_limited`` flag — never on the state or time. Rebuilding it on
        every ODE right-hand-side call (hundreds per solve) dominated runtime, so
        we memoise by ``nitrogen_limited`` (at most two variants). The cache is
        reset whenever parameters change (``update_parameters``) or a new solve
        starts (``solve_model``), so it can never go stale within a solve."""
        cache = self.__dict__.setdefault("_s_cache", {})
        key = self.nitrogen_limited
        matrix = cache.get(key)
        if matrix is None:
            matrix = self.build_stoichiometric_matrix(
                base_parameters=self.base_parameters, model_parameters=self.model_parameters,
                reactions=self.reactions, species=self.species, feed=self.feed,
                nitrogen_limited=key)
            cache[key] = matrix
        return matrix

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
        # model_parameters / base_parameters / feed feed the stoichiometric
        # matrix; invalidate its cache so the next access rebuilds it.
        self._s_cache={}

    def fit(self, train_data, search_space, *, optimizer="scipy",
            parameter_target="model_parameters", random_state=None, **optimize_kwargs):
        """Calibrate this model's parameters to experimental data.

        This gives every model type one homogeneous fitting call. Mechanistic
        models fit their ``search_space`` with a black-box optimizer here;
        :class:`NeuralADM` overrides this with gradient descent. Callers write
        ``model.fit(train_data, ...)`` regardless of the model type.

        Args:
            train_data: Iterable of :class:`core.Experiment`.
            search_space: ``{parameter: (low, high)}`` bounds to search.
            optimizer: ``"scipy"`` (differential evolution), ``"genetic"``, or ``"surrogate"``.
            parameter_target: which parameter group to fit. Defaults to ``"model_parameters"``.
            random_state: optional seed for reproducibility.
            **optimize_kwargs: forwarded to the optimizer's ``optimize`` (e.g. ``maxiter``,
                ``popsize``, ``workers``).

        Returns:
            Model: ``self``, with the fitted parameters applied.
        """
        import optimize as _optimize
        optimizers = {"scipy": _optimize.ScipyOptimizer,
                      "genetic": _optimize.GeneticOptimizer,
                      "surrogate": _optimize.SurrogateOptimizer}
        if optimizer not in optimizers:
            raise ValueError(f"optimizer must be one of {sorted(optimizers)}")
        opt = optimizers[optimizer](base_model=self, train_data=list(train_data),
                                    search_space=search_space,
                                    parameter_target=parameter_target,
                                    random_state=random_state)
        opt.optimize(**optimize_kwargs)
        self.update_parameters(**{parameter_target: opt.optimized_parameters})
        return self


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
        self._s_cache={}   # fresh matrix cache for this solve (params/feed fixed within it)
        y0=self.initial_conditions[:, 0]
        self._be_time=time.time()
        # Only numerical instability is caught here: the integrator returns
        # success=False when it cannot converge (e.g. a stiff blow-up), in which
        # case we hand back a large sentinel solution so callers such as the
        # parameter optimizer can penalise that parameter set. Any other
        # exception (a genuine bug, a bad configuration, a KeyError in the ODE,
        # ...) is left to propagate and stop the flow rather than being masked.
        c = scipy.integrate.solve_ivp(self.ode_system, (0,self.sim_time), y0, t_eval=t_eval, method=method, args=[self],rtol=1e-6)
        if not c.success:
            warnings.warn(
                f"solve_model: integration did not converge ({c.message}); "
                "returning a sentinel solution (numerical instability).",
                RuntimeWarning,
            )
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
        sol_df = pl.DataFrame(solution)

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

        try:
            from dash import Dash, dcc, html, Input, Output, dash_table
            import dash_bootstrap_components as dbc
            import dash_escher
        except ImportError as exc:
            raise ImportError(
                "Dash reports require optional dashboard dependencies. "
                "Install them with `pip install adtoolbox[dashboard]`."
            ) from exc

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
        sol_df = pl.DataFrame(solution)

        
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
                                        data=[self.base_parameters.copy()],
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
                                        data=[self.model_parameters.copy()],
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
                                        data=[self._ic.copy()],
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
                                        data=[self._inc.copy()],
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
            sol_df = pl.DataFrame(solution)

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
        """Write the model solution to CSV."""
        rows = []
        for species, values in zip(self.species, sol.y):
            row = {"species": species}
            row.update({str(time_point): value for time_point, value in zip(sol.t, values)})
            rows.append(row)
        pl.DataFrame(rows).write_csv(os.path.join(address,self.name+"_Report.csv"))
        
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

    def save(self, path: str | os.PathLike) -> str:
        """Serialize everything needed to reproduce this model to one JSON file.

        Captures the current parameters, initial/inlet conditions, feed, control
        states, and solver settings. The two callables (``ode_system`` and
        ``build_stoichiometric_matrix``) are stored by name and re-resolved from
        the :mod:`adtoolbox.adm` module on load, so only the packaged model
        variants (ADM1 / e-ADM) round-trip -- custom callables are not supported.

        Args:
            path: Destination ``.json`` file.

        Returns:
            str: The path written.
        """
        species = list(self.species)
        payload = {
            "format": "adtoolbox-model/1",
            "name": self.name,
            "switch": self.switch,
            "simulation_time": self.sim_time,
            "time_limit": self.time_limit,
            "ode_system": self.ode_system.__name__,
            "build_stoichiometric_matrix": self.build_stoichiometric_matrix.__name__,
            "control_state": {k: float(v) for k, v in self.control_state.items()},
            "model_parameters": self.model_parameters,
            "base_parameters": self.base_parameters,
            "reactions": list(self.reactions),
            "species": species,
            "initial_conditions": {s: float(self.initial_conditions[i, 0]) for i, s in enumerate(species)},
            "inlet_conditions": {s + "_in": float(self.inlet_conditions[i, 0]) for i, s in enumerate(species)},
            "feed": {
                "name": self.feed.name, "carbohydrates": self.feed.carbohydrates,
                "lipids": self.feed.lipids, "proteins": self.feed.proteins,
                "tss": self.feed.tss, "si": self.feed.si, "xi": self.feed.xi,
                "reference": getattr(self.feed, "reference", ""),
            },
        }
        with open(path, "w") as f:
            json.dump(payload, f, indent=1, default=float)
        return str(path)

    @classmethod
    def load(cls, path: str | os.PathLike) -> "Model":
        """Reconstruct a model previously written by :meth:`save`.

        Args:
            path: A JSON file produced by :meth:`save`.

        Returns:
            Model: A ready-to-solve model with the saved state restored.

        Raises:
            ValueError: If a stored callable cannot be resolved from
                :mod:`adtoolbox.adm` (e.g. a custom, non-packaged model).
        """
        with open(path) as f:
            payload = json.load(f)

        def _resolve(fn_name):
            fn = globals().get(fn_name)
            if not callable(fn):
                raise ValueError(
                    f"Cannot resolve callable {fn_name!r} from adtoolbox.adm; "
                    "Model.load only supports the packaged model variants."
                )
            return fn

        feed = payload["feed"]
        return cls(
            model_parameters=payload["model_parameters"],
            base_parameters=payload["base_parameters"],
            initial_conditions=payload["initial_conditions"],
            inlet_conditions=payload["inlet_conditions"],
            feed=Feed(name=feed["name"], carbohydrates=feed["carbohydrates"],
                      lipids=feed["lipids"], proteins=feed["proteins"],
                      tss=feed["tss"], si=feed["si"], xi=feed["xi"],
                      reference=feed.get("reference", "")),
            reactions=payload["reactions"],
            species=payload["species"],
            ode_system=_resolve(payload["ode_system"]),
            build_stoichiometric_matrix=_resolve(payload["build_stoichiometric_matrix"]),
            control_state=payload.get("control_state", {}),
            name=payload.get("name", "ADM"),
            switch=payload.get("switch", "DAE"),
            simulation_time=payload.get("simulation_time", 30),
            time_limit=payload.get("time_limit", -1),
        )

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



def build_adm1_stoichiometric_matrix(base_parameters: dict,
                                     model_parameters: dict,
                                     reactions: list,
                                     species: list,
                                     feed: Feed,
                                     nitrogen_limited: bool = False)-> np.ndarray:
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

    S = np.zeros((len(species), len(reactions)))
    S[species.index('S_su'), list(map(reactions.index, ['Hydrolysis carbohydrates', 'Hydrolysis of lipids', 'Uptake of sugars']))] = [1, (1-model_parameters["f_fa_li"]), - 1]
    S[species.index('S_aa'), list(map(reactions.index, ['Hydrolysis of proteins', 'Uptake of amino acids']))] = [1, -1]
    S[species.index('S_fa'), list(map(reactions.index, ['Hydrolysis of lipids', 'Uptake of LCFA']))] = [(model_parameters["f_fa_li"]), - 1]
    Y_aa=0 if nitrogen_limited else model_parameters['Y_aa']
    S[species.index('S_va'), list(map(reactions.index, ['Uptake of amino acids', 'Uptake of valerate']))] = [(1-Y_aa) *
                    model_parameters['f_va_aa'], - 1]
    Y_su=0 if nitrogen_limited else model_parameters['Y_su']
    S[species.index('S_bu'), list(map(reactions.index, ['Uptake of sugars', 'Uptake of amino acids', 'Uptake of butyrate']))] = [(1-Y_su)*model_parameters['f_bu_su'],
                       (1-Y_aa)*model_parameters["f_bu_aa"], - 1]
    S[species.index('S_pro'), list(map(reactions.index, ['Uptake of sugars', 'Uptake of amino acids', 'Uptake of valerate', 'Uptake of propionate']))] = [(1-model_parameters["Y_su"])*model_parameters['f_pro_su'],
                          (1-Y_aa)*model_parameters["f_pro_aa"], (1 - model_parameters['Y_c4'])*0.54, -1]
    
    Y_fa=0 if nitrogen_limited else model_parameters['Y_fa'] 
    S[species.index('S_ac'), list(map(reactions.index, ['Uptake of sugars', 'Uptake of amino acids', 'Uptake of LCFA', 'Uptake of valerate', 'Uptake of butyrate', 'Uptake of propionate', 'Uptake of acetate']))] = [(1-Y_su)*model_parameters['f_ac_su'],
                                    (1-Y_aa) *
                                    model_parameters['f_ac_aa'],
                                    (1-Y_fa)*0.7,
                                    (1-model_parameters['Y_c4'])*0.31,
                                    (1-model_parameters['Y_c4'])*0.8,
                                    (1-model_parameters['Y_pro'])*0.57,
                                    -1]
    S[species.index('S_h2'), list(map(reactions.index, ['Uptake of sugars', 'Uptake of amino acids', 'Uptake of LCFA', 'Uptake of valerate', 'Uptake of butyrate', 'Uptake of propionate', 'Uptake of Hydrogen', 'Gas Transfer H2']))] = [(1-Y_su)*model_parameters['f_h2_su'],
                                        (1-Y_aa) *
                                        model_parameters['f_h2_aa'],
                                        (1-Y_fa)*0.3,
                                        (1-model_parameters['Y_c4'])*0.15,
                                        (1-model_parameters['Y_c4'])*0.2,
                                        (1-model_parameters['Y_pro'])*0.43,
                                        -1,
                                        -1]
    S[species.index('S_ch4'), list(map(reactions.index, ['Uptake of acetate', 'Uptake of Hydrogen', 'Gas Transfer CH4']))] = [(1-model_parameters['Y_ac']),
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
    S[species.index('S_IC'), list(map(reactions.index, ['Disintegration', 'Hydrolysis carbohydrates', 'Hydrolysis of proteins', 'Hydrolysis of lipids', 'Uptake of sugars', 'Uptake of amino acids', 'Uptake of LCFA', 'Uptake of valerate', 'Uptake of butyrate', 'Uptake of propionate', 'Uptake of acetate', 'Uptake of Hydrogen', 'Decay of Xsu', 'Decay of Xaa', 'Decay of Xfa', 'Decay of Xc4', 'Decay of Xpro', 'Decay of Xac', 'Decay of Xh2', 'Gas Transfer CO2']))] = [-s_1, -s_2, -s_3, -s_4, -
                                                                                    s_5, -s_6, -s_7, -s_8, -s_9, -s_10, -s_11, -s_12, -s_13, -s_13, -s_13, -s_13, -s_13, -s_13, -s_13, -1]
    S[species.index('S_IN'), list(map(reactions.index, ['Disintegration', 'Uptake of sugars', 'Uptake of amino acids', 'Uptake of LCFA', 'Uptake of valerate', 'Uptake of butyrate', 'Uptake of propionate', 'Uptake of acetate', 'Uptake of Hydrogen', 'Decay of Xsu', 'Decay of Xaa', 'Decay of Xfa', 'Decay of Xc4', 'Decay of Xpro', 'Decay of Xac', 'Decay of Xh2']))] = [model_parameters['N_xc']-model_parameters['f_xI_xc']*model_parameters['N_I']-model_parameters['f_sI_xc']*model_parameters['N_I']-model_parameters['f_pr_xc']*model_parameters['N_aa'],
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
    S[species.index('S_I'), reactions.index('Disintegration')] = model_parameters['f_sI_xc']
    S[species.index('X_xc'), list(map(reactions.index, ['Disintegration', 'Decay of Xsu', 'Decay of Xaa', 'Decay of Xfa', 'Decay of Xc4', 'Decay of Xpro', 'Decay of Xac', 'Decay of Xh2']))] = [-1, 1, 1, 1, 1, 1, 1, 1]
    S[species.index('X_ch'), list(map(reactions.index, ['Disintegration', 'Hydrolysis carbohydrates']))] = [model_parameters['f_ch_xc'], -1]
    S[species.index('X_pr'), list(map(reactions.index, ['Disintegration', 'Hydrolysis of proteins']))] = [model_parameters['f_pr_xc'], -1]
    S[species.index('X_li'), list(map(reactions.index, ['Disintegration', 'Hydrolysis of lipids']))] = [model_parameters['f_li_xc'], -1]
    S[species.index('X_su'), list(map(reactions.index, ['Uptake of sugars', 'Decay of Xsu']))] = [Y_su, -1]
    S[species.index('X_aa'), list(map(reactions.index, ['Uptake of amino acids', 'Decay of Xaa']))] = [Y_aa, -1]
    S[species.index('X_fa'), list(map(reactions.index, ['Uptake of LCFA', 'Decay of Xfa']))] = [Y_fa, -1]
    S[species.index('X_c4'), list(map(reactions.index, ['Uptake of valerate', 'Uptake of butyrate', 'Decay of Xc4']))] = [model_parameters['Y_c4'], model_parameters['Y_c4'], -1]
    S[species.index('X_pro'), list(map(reactions.index, ['Uptake of propionate', 'Decay of Xpro']))] = [model_parameters['Y_pro'], -1]
    S[species.index('X_ac'), list(map(reactions.index, ['Uptake of acetate', 'Decay of Xac']))] = [model_parameters['Y_ac'], -1]
    S[species.index('X_h2'), list(map(reactions.index, ['Uptake of Hydrogen', 'Decay of Xh2']))] = [model_parameters['Y_h2'], -1]
    S[species.index('X_I'), reactions.index('Disintegration')] = model_parameters['f_xI_xc']
    S[species.index('S_cation'), :] = 0
    S[species.index('S_anion'), :] = 0
    S[species.index('S_H_ion'), :] = 0
    S[species.index('S_va_ion'), reactions.index('Acid Base Equilibrium (Va)')] = -1
    S[species.index('S_bu_ion'), reactions.index('Acid Base Equilibrium (Bu)')] = -1
    S[species.index('S_pro_ion'), reactions.index('Acid Base Equilibrium (Pro)')] = -1
    S[species.index('S_ac_ion'), reactions.index('Acid Base Equilibrium (Ac)')] = -1
    S[species.index('S_hco3_ion'), reactions.index('Acid Base Equilibrium (CO2)')] = -1
    S[species.index('S_co2'), :] = 0
    S[species.index('S_nh3'), reactions.index('Acid Base Equilibrium (In)')] = -1
    S[species.index('S_nh4_ion'), :] = 0
    S[species.index('S_gas_h2'), reactions.index('Gas Transfer H2')] = base_parameters['V_liq']/base_parameters['V_gas']
    S[species.index('S_gas_ch4'), reactions.index('Gas Transfer CH4')] = base_parameters['V_liq']/base_parameters['V_gas']
    S[species.index('S_gas_co2'), reactions.index('Gas Transfer CO2')] = base_parameters['V_liq']/base_parameters['V_gas']
    return S


build_adm1_stoiciometric_matrix = build_adm1_stoichiometric_matrix


def adm1_ode_sys(t: float, c: np.ndarray, model:Model)-> np.ndarray:
    """ The ODE system for the original ADM.
        No testing is done.

        Args:
            t (float):a matrix of zeros to be filled
            c (np.ndarray): an array of concentrations to be filled
            model (Model): An instance of Model to calculate the ODE with

        Returns:
            np.ndarray: The output is dCdt, the change of concentration with respect to time.
    """
    species_index = model.species.index
    reaction_index = model.reactions.index

    c[species_index('S_nh4_ion')] = c[species_index('S_IN')] - c[species_index('S_nh3')]
    c[species_index('S_co2')] = c[species_index('S_IC')] - c[species_index('S_hco3_ion')]
    I_pH_aa = (model.model_parameters["K_pH_aa"] ** model.model_parameters['nn_aa'])/(np.power(
        c[species_index('S_H_ion')], model.model_parameters['nn_aa']) + np.power(model.model_parameters["K_pH_aa"], model.model_parameters['nn_aa']))
    I_pH_ac = (model.model_parameters['K_pH_ac'] ** model.model_parameters["n_ac"])/(
        c[species_index('S_H_ion')] ** model.model_parameters['n_ac'] + model.model_parameters['K_pH_ac'] ** model.model_parameters['n_ac'])
    I_pH_h2 = (model.model_parameters['K_pH_h2']**model.model_parameters['n_h2'])/(
        c[species_index('S_H_ion')] ** model.model_parameters['n_h2'] + model.model_parameters['K_pH_h2']**model.model_parameters['n_h2'])
    I_IN_lim = 1 / (1+(model.model_parameters['K_S_IN'] / c[species_index('S_IN')]))
    I_h2_fa = 1 / (1+(c[species_index('S_h2')] / model.model_parameters['K_I_h2_fa']))
    I_h2_c4 = 1 / (1+(c[species_index('S_h2')]/model.model_parameters['K_I_h2_c4']))
    I_h2_pro = (1/(1+(c[species_index('S_h2')]/model.model_parameters['K_I_h2_pro'])))
    I_nh3 = 1/(1+(c[species_index('S_nh3')]/model.model_parameters['K_I_nh3']))
    I5 = (I_pH_aa * I_IN_lim)
    I6 = np.copy(I5)
    I7 = (I_pH_aa * I_IN_lim * I_h2_fa)
    I8 = (I_pH_aa * I_IN_lim * I_h2_c4)
    I9 = np.copy(I8)
    I10 = (I_pH_aa * I_IN_lim * I_h2_pro)
    I11 = (I_pH_ac * I_IN_lim * I_nh3)
    I12 = (I_pH_h2 * I_IN_lim)
    v = np.zeros((len(model.reactions), 1))
    v[reaction_index('Disintegration')] = model.model_parameters["k_dis"]*c[species_index('X_xc')]
    
    v[reaction_index('Hydrolysis carbohydrates')] = model.model_parameters['k_hyd_ch']*c[species_index('X_ch')]
    v[reaction_index('Hydrolysis of proteins')] = model.model_parameters['k_hyd_pr']*c[species_index('X_pr')]
    v[reaction_index('Hydrolysis of lipids')] = model.model_parameters['k_hyd_li']*c[species_index('X_li')]
    
    v[reaction_index('Uptake of sugars')] = model.model_parameters['k_m_su']*c[species_index('S_su')] / \
(model.model_parameters['K_S_su']+c[species_index('S_su')])*c[species_index('X_su')]*I5
    v[reaction_index('Uptake of amino acids')] = model.model_parameters['k_m_aa']*c[species_index('S_aa')] / \
        (model.model_parameters['K_S_aa']+c[species_index('S_aa')])*c[species_index('X_aa')]*I6
    v[reaction_index('Uptake of LCFA')] = model.model_parameters['k_m_fa']*c[species_index('S_fa')] / \
        (model.model_parameters['K_S_fa']+c[species_index('S_fa')])*c[species_index('X_fa')]*I7
    v[reaction_index('Uptake of valerate')] = model.model_parameters['k_m_c4']*c[species_index('S_va')] / \
        (model.model_parameters['K_S_c4']+c[species_index('S_va')]) * \
        c[species_index('X_c4')]*c[species_index('S_va')]/(c[species_index('S_va')]+c[species_index('S_bu')]+10 ** (-6))*I8
    v[reaction_index('Uptake of butyrate')] = model.model_parameters['k_m_c4']*c[species_index('S_bu')] / \
        (model.model_parameters['K_S_c4']+c[species_index('S_bu')]) * \
        c[species_index('X_c4')]*c[species_index('S_bu')]/(c[species_index('S_bu')]+c[species_index('S_va')]+10 ** (-6))*I9
    v[reaction_index('Uptake of propionate')] = model.model_parameters['k_m_pr']*c[species_index('S_pro')] / \
        (model.model_parameters['K_S_pro']+c[species_index('S_pro')])*c[species_index('X_pro')]*I10
    v[reaction_index('Uptake of acetate')] = model.model_parameters['k_m_ac']*c[species_index('S_ac')] / \
        (model.model_parameters['K_S_ac']+c[species_index('S_ac')])*c[species_index('X_ac')]*I11
    v[reaction_index('Uptake of Hydrogen')] = model.model_parameters['k_m_h2']*c[species_index('S_h2')] / \
        (model.model_parameters['K_S_h2']+c[species_index('S_h2')])*c[species_index('X_h2')]*I12
    v[reaction_index('Decay of Xsu')] = model.model_parameters['k_dec_X_su']*c[species_index('X_su')]
    v[reaction_index('Decay of Xaa')] = model.model_parameters['k_dec_X_aa']*c[species_index('X_aa')]
    v[reaction_index('Decay of Xfa')] = model.model_parameters['k_dec_X_fa']*c[species_index('X_fa')]
    v[reaction_index('Decay of Xc4')] = model.model_parameters['k_dec_X_c4']*c[species_index('X_c4')]
    v[reaction_index('Decay of Xpro')] = model.model_parameters['k_dec_X_pro']*c[species_index('X_pro')]
    v[reaction_index('Decay of Xac')] = model.model_parameters['k_dec_X_ac']*c[species_index('X_ac')]
    v[reaction_index('Decay of Xh2')] = model.model_parameters['k_dec_X_h2']*c[species_index('X_h2')]
    v[reaction_index('Acid Base Equilibrium (Va)')] = model.model_parameters['k_A_B_va'] * \
        (c[species_index('S_va_ion')] * (model.model_parameters['K_a_va'] + c[species_index('S_H_ion')]) -
         model.model_parameters['K_a_va'] * c[species_index('S_va')])
    v[reaction_index('Acid Base Equilibrium (Bu)')] = model.model_parameters['k_A_B_bu'] * \
        (c[species_index('S_bu_ion')] * (model.model_parameters['K_a_bu'] + c[species_index('S_H_ion')]) -
         model.model_parameters['K_a_bu'] * c[species_index('S_bu')])
    v[reaction_index('Acid Base Equilibrium (Pro)')] = model.model_parameters['k_A_B_pro'] * \
        (c[species_index('S_pro_ion')] * (model.model_parameters['K_a_pro'] + c[species_index('S_H_ion')]) -
         model.model_parameters['K_a_pro'] * c[species_index('S_pro')])
    v[reaction_index('Acid Base Equilibrium (Ac)')] = model.model_parameters['k_A_B_ac'] * \
        (c[species_index('S_ac_ion')] * (model.model_parameters['K_a_ac'] + c[species_index('S_H_ion')]) -
         model.model_parameters['K_a_ac'] * c[species_index('S_ac')])
    v[reaction_index('Acid Base Equilibrium (CO2)')] = model.model_parameters['k_A_B_co2'] * \
        (c[species_index('S_hco3_ion')] * (model.model_parameters['K_a_co2'] + c[species_index('S_H_ion')]) -
         model.model_parameters['K_a_co2'] * c[species_index('S_IC')])
    v[reaction_index('Acid Base Equilibrium (In)')] = model.model_parameters['k_A_B_IN'] * \
        (c[species_index('S_nh3')] * (model.model_parameters['K_a_IN'] + c[species_index('S_H_ion')]) -
         model.model_parameters['K_a_IN'] * c[species_index('S_IN')])
    p_gas_h2 = c[species_index('S_gas_h2')] * model.base_parameters["R"] * \
        model.base_parameters["T_op"] / 16
    p_gas_ch4 = c[species_index('S_gas_ch4')] * model.base_parameters["R"] * \
        model.base_parameters["T_op"] / 64
    p_gas_co2 = c[species_index('S_gas_co2')] * model.base_parameters["R"] * \
        model.base_parameters["T_op"]
    p_gas_h2o = 0.0313 * \
        np.exp(5290 *
               (1 / model.base_parameters["T_base"] - 1 / model.base_parameters["T_op"]))
    P_gas = p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o
    q_gas = max(
        0, (model.model_parameters['k_p'] * (P_gas - model.base_parameters['P_atm'])))
    v[reaction_index('Gas Transfer H2')] = model.model_parameters['k_L_a'] * \
        (c[species_index('S_h2')] - 16 * model.model_parameters['K_H_h2'] * p_gas_h2)
    v[reaction_index('Gas Transfer CH4')] = model.model_parameters['k_L_a'] * \
        (c[species_index('S_ch4')] - 64 * model.model_parameters['K_H_ch4'] * p_gas_ch4)
    v[reaction_index('Gas Transfer CO2')] = model.model_parameters['k_L_a'] * \
        (c[species_index('S_co2')] - model.model_parameters['K_H_co2'] * p_gas_co2)
    dCdt = np.matmul(model.s, v)
    
    if c[species_index('S_IN')]<0.01:
        model.nitrogen_limited=True
    else:
        model.nitrogen_limited=False
    
    phi = c[species_index('S_cation')]+c[species_index('S_nh4_ion')]-c[species_index('S_hco3_ion')] - (c[species_index('S_ac_ion')] / 64) - (c[species_index('S_pro_ion')] / 112) - (c[species_index('S_bu_ion')] / 160) - (c[species_index('S_va_ion')] / 208) - c[species_index('S_anion')]
    c[species_index('S_H_ion')] = (-1 * phi / 2) + (0.5 * np.sqrt(phi**2 + 4 * model.model_parameters['K_w']))
    
    gas_start = species_index('S_gas_h2')
    dCdt[0: gas_start] = dCdt[0: gas_start]+model.base_parameters['q_in'] / model.base_parameters["V_liq"] * \
        (model.inlet_conditions[0: gas_start]-c[0:gas_start].reshape(-1, 1))
    
        
    dCdt[gas_start:] = dCdt[gas_start:]+q_gas/model.base_parameters["V_gas"] * (model.inlet_conditions[gas_start:]-c[gas_start:].reshape(-1, 1))
    dCdt[[species_index('S_H_ion'), species_index('S_co2'), species_index('S_nh4_ion')], 0] = 0
    if model.switch == "DAE":
        # NOTE: dissolved S_h2 is left dynamic (previously frozen here, which
        # discarded fermentation H2 and starved hydrogenotrophic methanogenesis).
        dCdt[species_index('S_va_ion'): species_index('S_co2')] = 0
        dCdt[species_index('S_nh3')] = 0
    
    if model.control_state.keys():
        for state in model.control_state.keys():
            c[model.species.index(state)]=model.control_state[state]
            dCdt[model.species.index(state)]=0
    
    
    return dCdt[:, 0]


def build_e_adm_stoichiometric_matrix(base_parameters: dict,
                                      model_parameters: dict,
                                      reactions: list,
                                      species: list,
                                      feed:Feed,
                                      nitrogen_limited:bool=False)->np.ndarray:
    """ 
    This function builds the stoichiometric matrix for the e-ADM model.
        
        Model Parameters (dict): a dictionary which contains model parameters
        base_parameters (dict): a dictionary which contains base paramters
        Initial Conditions (dict): a dictionary containing inlet conditions for all species
        Inlet Conditions (dict): a dictionary containing inlet conditions for all species
        reactions (list): a list containing all of the reaction names
        species (list): a list containing all species
        
    Returns:
        np.ndarray: Returns an matrix of stochiometic values.
    """
    _require_model_parameters(
        model_parameters,
        [
            'Y_su', 'Y_aa', 'Y_fa', 'Y_ac_et', 'Y_ac_lac', 'Y_pro_et',
            'Y_pro_lac', 'Y_bu_et', 'Y_bu_lac', 'Y_va', 'Y_cap', 'Y_bu',
            'Y_Me_ac', 'Y_Me_CO2', 'Y_ac_et_ox', 'Y_pro_lac_ox',
        ],
        "e-ADM",
    )
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
                                                  (1- model_parameters['f_et_ac']-Y_ac_et) * model_parameters['f_bu_ac'],
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
    Y_pro_lac=0 if nitrogen_limited else model_parameters['Y_pro_lac']
    
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
                                                     Y_pro_et]

    S[list(map(species.index, ["S_pro", "S_lac", "S_va", "S_IN", "S_IC", "S_h2", "X_chain_lac"])),
        reactions.index('Uptake of propionate_lac')] = [-1,
                                                        model_parameters['f_lac_pro'],
                                                        (1-model_parameters['f_lac_pro']-Y_pro_lac) * model_parameters['f_va_pro'],
                                                        -Y_pro_lac * model_parameters['N_bac'],
                                                        f_IC_pro_lac,
                                                        (1-model_parameters['f_lac_pro']-Y_pro_lac) * (1-model_parameters['f_va_pro']),
                                                        Y_pro_lac]

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
        
    Y_bu=0 if nitrogen_limited else model_parameters['Y_bu']
    S[list(map(species.index, ["S_bu", "S_ac", "X_VFA_deg"])),
        reactions.index('Uptake of butyrate')] = [-1,
                                                  (1 - Y_bu),
                                                  Y_bu]
        
        
    
    Y_Me_ac=0 if nitrogen_limited else model_parameters["Y_Me_ac"]
    # Acetoclastic methanogenesis is a disproportionation of acetate alone
    # (CH3COO- + H2O -> CH4 + HCO3-); it does not consume H2. COD closes:
    # -1 acetate -> (1-Y) CH4 + Y biomass, CO2 released via the S_IC balance.
    f_IC_Me_ach2 = -(-1*model_parameters['C_ac'] +
                     (1 - Y_Me_ac)*model_parameters['C_ch4'] +
                     Y_Me_ac*model_parameters['C_bac'])
    S[list(map(species.index, ["S_ac", "S_ch4", "X_Me_ac", 'S_IC', 'S_IN'])),
        reactions.index('Methanogenessis from acetate and h2')] = [-1,
                                                                   (1 - Y_Me_ac),
                                                                   Y_Me_ac,
                                                                   f_IC_Me_ach2,
                                                                    -Y_Me_ac *model_parameters['N_bac']
                                                                   ]
    
    Y_Me_CO2=0 if nitrogen_limited else model_parameters["Y_Me_CO2"]
    # Hydrogenotrophic methanogenesis acts on DISSOLVED H2 and produces
    # DISSOLVED CH4 (which then reaches the headspace via Gas Transfer CH4).
    # COD closes: -1 H2 -> (1-Y) CH4 + Y biomass; CO2 consumed via S_IC.
    f_IC_Me_co2 = -((1 - Y_Me_CO2)*model_parameters['C_ch4'] +
                    Y_Me_CO2*model_parameters['C_bac'])
    S[list(map(species.index, ["S_h2", "S_ch4", "X_Me_CO2", 'S_IC',"S_IN"])),
        reactions.index('Methanogenessis from CO2 and h2')] = [-1,
                                                               (1 - Y_Me_CO2),
                                                               (Y_Me_CO2),
                                                               f_IC_Me_co2,
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
        reactions.index('Gas Transfer H2')] = [-1, base_parameters['V_liq']/base_parameters['V_gas']]
    S[list(map(species.index, ["S_ch4", "S_gas_ch4"])),
        reactions.index('Gas Transfer CH4')] = [-1, base_parameters['V_liq']/base_parameters['V_gas']]
    S[list(map(species.index, ["S_co2", "S_gas_co2"])),
        reactions.index('Gas Transfer CO2')] = [-1, base_parameters['V_liq']/base_parameters['V_gas']]
    
    return S


def e_adm_ode_sys(t: float, c: np.ndarray, model: Model)-> np.ndarray:
    """
    This function is used to build the ODEs of the e-ADM model.
    
    Args:
        t (float):a matrix of zeros to be filled
        c (np.ndarray): an array of concentrations to be filled
        model (Model): An instance of Model to calculate the ODE with

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
    
    I_IN_lim = _monod_limitation(
        c[model.species.index('S_IN')],
        model.model_parameters['K_S_IN'],
    )
    
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

    # Acetoclastic: Monod on dissolved acetate only (no H2 dependence).
    v[model.reactions.index('Methanogenessis from acetate and h2')] = model.model_parameters['k_m_h2_Me_ac']*c[model.species.index('S_ac')] / \
        (model.model_parameters['K_S_ac_Me']+c[model.species.index('S_ac')]+10**-9)*c[model.species.index('X_Me_ac')]*I12

    # Hydrogenotrophic: Monod on dissolved H2 (CO2 assumed non-limiting).
    v[model.reactions.index('Methanogenessis from CO2 and h2')] = model.model_parameters['k_m_h2_Me_CO2']*c[model.species.index('S_h2')] / \
        (model.model_parameters['K_S_h2_Me_CO2']+c[model.species.index('S_h2')]+10**-9)*c[model.species.index('X_Me_CO2')]*I12


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
         model.model_parameters['K_a_IN'] * c[model.species.index('S_IN')])

    
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


# ============================================================================
# Data-driven models that share the Model interface (option 1: COD-balanced)
# ----------------------------------------------------------------------------
# A neural network outputs non-negative reaction *rates*; concentrations follow
# dc/dt = S @ rate with e-ADM's COD-balanced stoichiometric matrix S. Because
# every biochemical column of S sums to zero in COD, total COD is conserved for
# ANY network output -- the net only learns the (uncertain) rate laws. These
# models expose the same interface as Model (species/solve_model/copy/fit) and
# are trained by gradient descent, so `model.fit(train_data, ...)` works for
# both mechanistic and neural models (the fit internals differ, the call does not).
# ============================================================================

DEFAULT_OBSERVED = ["S_ac", "S_pro", "S_bu", "S_va", "S_cap"]
_NEURAL_BIOMASS = ["X_su", "X_aa", "X_fa", "X_ac_et", "X_lac", "X_et",
                   "X_pr", "X_ch", "X_li", "X_Me_ac", "X_Me_CO2"]


def _require_torch():
    try:
        import torch
        return torch
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "Neural ADM models require PyTorch. Install it with `pip install torch`."
        ) from exc


def _neural_features(ic: dict, t: float, observed) -> list:
    """Condition features fed to the rate network at time t (see _NEURAL_BIOMASS)."""
    import math
    pH = -math.log10(ic.get("S_H_ion", 10 ** -6.5))
    cod = float(ic.get("TSS", 0.0)) + float(ic.get("TDS", 0.0))
    bio = [float(ic.get(b, 0.0)) for b in _NEURAL_BIOMASS]
    v0 = [float(ic.get(v, 0.0)) for v in observed]
    return [float(t), pH, cod, sum(bio), *bio, *v0]


def _cumulative_trapz_np(rates, t):
    dt = np.diff(t)
    incr = 0.5 * (rates[1:] + rates[:-1]) * dt[:, None]
    cum = np.zeros_like(rates)
    cum[1:] = np.cumsum(incr, axis=0)
    return cum


def _build_rate_net(backbone, in_dim, out_dim, hidden=32, seed=0):
    """Return a torch module mapping (batch, time, in_dim) -> non-negative rates."""
    torch = _require_torch()
    import torch.nn as nn
    import torch.nn.functional as F
    torch.manual_seed(seed)

    class MLP(nn.Module):
        def __init__(self):
            super().__init__()
            self.net = nn.Sequential(nn.Linear(in_dim, hidden), nn.SiLU(),
                                     nn.Linear(hidden, hidden), nn.SiLU(),
                                     nn.Linear(hidden, out_dim))
        def forward(self, x):
            return F.softplus(self.net(x))

    class LSTMNet(nn.Module):
        def __init__(self):
            super().__init__()
            self.rnn = nn.LSTM(in_dim, hidden, batch_first=True)
            self.head = nn.Linear(hidden, out_dim)
        def forward(self, x):
            return F.softplus(self.head(self.rnn(x)[0]))

    class TransformerNet(nn.Module):
        def __init__(self):
            super().__init__()
            self.inp = nn.Linear(in_dim, hidden)
            layer = nn.TransformerEncoderLayer(hidden, 2, 2 * hidden, batch_first=True, dropout=0.0)
            self.enc = nn.TransformerEncoder(layer, 1)
            self.head = nn.Linear(hidden, out_dim)
        def forward(self, x):
            return F.softplus(self.head(self.enc(self.inp(x))))

    b = backbone.lower()
    if b == "mlp":
        net = MLP()
    elif b == "lstm":
        net = LSTMNet()
    elif b in ("transformer", "tfm"):
        net = TransformerNet()
    else:
        raise ValueError(f"backbone must be 'mlp', 'lstm', or 'transformer'; got {backbone!r}")
    # Start from near-zero rates: with softplus(bias<<0) the model begins pinned at
    # the initial condition and grows rates as needed -- essential for a stable,
    # non-diverging integration through S.
    for module in reversed(list(net.modules())):
        if isinstance(module, nn.Linear):
            nn.init.zeros_(module.weight)
            nn.init.constant_(module.bias, -4.0)
            break
    return net


def _neural_ode_sys(t, c, model):
    """Open-loop neural RHS (dc/dt = S @ rate(t)); available for scipy integration.

    ``solve_model`` uses the cumulative-rate rollout directly so that training and
    inference share identical math, but this keeps a Model-compatible ODE callable.
    """
    return np.asarray(model.s, float) @ model._rate_at(float(t))


class NeuralADM(Model):
    """A COD-balanced, data-driven model that shares the :class:`Model` interface.

    The network (``backbone`` = ``"lstm"``, ``"transformer"``, or ``"mlp"``) emits
    non-negative reaction rates; concentrations follow ``dc/dt = S @ rate`` using
    e-ADM's COD-balanced stoichiometric matrix, so total COD is conserved by
    construction. Fit with :meth:`fit` (gradient descent); solve with the usual
    :meth:`solve_model`. Build one from an existing e-ADM model with
    :meth:`from_model`.
    """

    def __init__(self, model_parameters, base_parameters, initial_conditions, inlet_conditions,
                 feed, reactions, species, ode_system=None, build_stoichiometric_matrix=None,
                 control_state={}, name="NeuralADM", switch="DAE", simulation_time=30, time_limit=-1,
                 *, backbone="lstm", observed=None, hidden=32, rate_net=None, seed=0):
        super().__init__(model_parameters=model_parameters, base_parameters=base_parameters,
                         initial_conditions=initial_conditions, inlet_conditions=inlet_conditions,
                         feed=feed, reactions=reactions, species=species,
                         ode_system=(ode_system or _neural_ode_sys),
                         build_stoichiometric_matrix=(build_stoichiometric_matrix or build_e_adm_stoichiometric_matrix),
                         control_state=control_state, name=name, switch=switch,
                         simulation_time=simulation_time, time_limit=time_limit)
        self.backbone = backbone
        self.hidden = hidden
        self.observed = list(observed) if observed else list(DEFAULT_OBSERVED)
        self._feat_dim = 4 + len(_NEURAL_BIOMASS) + len(self.observed)
        self.rate_net = rate_net if rate_net is not None else _build_rate_net(
            backbone, self._feat_dim, len(reactions), hidden, seed)
        self._feat_mean = None      # feature standardisation, learned in fit()
        self._feat_std = None

    @classmethod
    def from_model(cls, base_model, *, backbone="lstm", observed=None, hidden=32, seed=0):
        """Build a NeuralADM that reuses a fitted e-ADM model's structure (its S)."""
        return cls(model_parameters=base_model.model_parameters.copy(),
                   base_parameters=base_model.base_parameters.copy(),
                   initial_conditions=base_model._ic.copy(), inlet_conditions=base_model._inc.copy(),
                   feed=base_model.feed, reactions=base_model.reactions.copy(),
                   species=base_model.species.copy(), control_state=base_model.control_state.copy(),
                   simulation_time=base_model.sim_time, backbone=backbone, observed=observed,
                   hidden=hidden, seed=seed)

    def _features(self, t_eval):
        # build conditions from the LIVE initial_conditions array (update_parameters
        # writes there, not to self._ic), so features reflect the current experiment.
        ic = {s: float(self.initial_conditions[i, 0]) for i, s in enumerate(self.species)}
        return np.array([_neural_features(ic, float(t), self.observed) for t in t_eval], float)

    def _norm(self, feats):
        if self._feat_mean is None:
            return feats
        return (feats - self._feat_mean) / self._feat_std

    def _rate_at(self, t):
        torch = _require_torch()
        feats = self._norm(self._features([t]))
        with torch.no_grad():
            return self.rate_net(torch.tensor(feats, dtype=torch.float32).unsqueeze(0)).squeeze(0).numpy()[0]

    def solve_model(self, t_eval, method="BDF"):
        """Predict the COD-balanced trajectory at ``t_eval`` (same return as Model)."""
        torch = _require_torch()
        t_eval = np.asarray(t_eval, float)
        feats = self._norm(self._features(t_eval))
        with torch.no_grad():
            rates = self.rate_net(torch.tensor(feats, dtype=torch.float32).unsqueeze(0)).squeeze(0).numpy()
        cum = _cumulative_trapz_np(rates, t_eval)          # (n_t, n_rxn)
        S = np.asarray(self.s, float)                      # (n_species, n_rxn), COD-balanced
        y = self.initial_conditions[:, 0][:, None] + S @ cum.T
        self.info = {"Fluxes": rates}
        return _Fake_Sol(y=y, t=t_eval)

    def fit(self, train_data, *, epochs=400, lr=5e-3, weight_decay=1e-4, seed=0, verbose=False):
        """Train the rate network by gradient descent on observed VFA trajectories.

        Same call as :meth:`Model.fit` (``model.fit(train_data, ...)``); the internals
        are gradient descent instead of a black-box search. COD balance is structural
        (via S), so it is never part of the loss.
        """
        torch = _require_torch()
        import torch.nn as nn
        torch.manual_seed(seed)
        S = torch.tensor(np.asarray(self.s, float), dtype=torch.float32)
        vfa_rows = [self.species.index(v) for v in self.observed]
        S_vfa = S[vfa_rows]                                # (n_vfa, n_rxn)

        raw = []
        for e in train_data:
            ic = dict(self._ic); ic.update(e.initial_concentrations)
            obs = np.asarray(e.data, float)                # (n_t, n_vars) -- Experiment stores (time, var)
            for j, v in enumerate(e.variables):
                ic[v] = float(obs[0, j])
            t = np.asarray(e.time, float)
            feats = np.array([_neural_features(ic, float(ti), self.observed) for ti in t], float)
            c0 = np.array([ic[v] for v in self.observed], float)
            var_idx = [self.observed.index(v) for v in e.variables]
            raw.append((feats, t, obs, c0, var_idx))

        # standardise features (per-column) for stable optimisation
        allf = np.concatenate([r[0] for r in raw], axis=0)
        self._feat_mean = allf.mean(axis=0)
        self._feat_std = allf.std(axis=0) + 1e-8

        data = []
        for feats, t, obs, c0, var_idx in raw:
            feats = (feats - self._feat_mean) / self._feat_std
            data.append((torch.tensor(feats, dtype=torch.float32),
                         torch.tensor(t, dtype=torch.float32),
                         torch.tensor(obs, dtype=torch.float32),
                         torch.tensor(c0, dtype=torch.float32), var_idx))

        opt = torch.optim.Adam(self.rate_net.parameters(), lr=lr, weight_decay=weight_decay)
        lossf = nn.MSELoss()
        for epoch in range(epochs):
            opt.zero_grad(); total = 0.0
            for feats, t, obs, c0, var_idx in data:
                rates = self.rate_net(feats.unsqueeze(0)).squeeze(0)      # (n_t, n_rxn)
                dt = (t[1:] - t[:-1]).unsqueeze(1)
                incr = 0.5 * (rates[1:] + rates[:-1]) * dt
                cum = torch.zeros_like(rates)
                cum[1:] = torch.cumsum(incr, dim=0)
                pred = (c0 + cum @ S_vfa.T)[:, var_idx]                   # anchored at IC, COD-balanced
                total = total + lossf(pred, obs)
            loss = total / len(data)
            loss.backward(); opt.step()
            if verbose and epoch % 100 == 0:
                print(f"  epoch {epoch:4d}  loss={float(loss):.5f}")
        return self

    def copy(self):
        new = type(self)(model_parameters=self.model_parameters.copy(),
                         base_parameters=self.base_parameters.copy(),
                         initial_conditions=self._ic.copy(), inlet_conditions=self._inc.copy(),
                         feed=self.feed, reactions=self.reactions.copy(), species=self.species.copy(),
                         control_state=self.control_state.copy(), name=self.name, switch=self.switch,
                         time_limit=self.time_limit, simulation_time=self.sim_time,
                         backbone=self.backbone, observed=list(self.observed), hidden=self.hidden)
        new.rate_net.load_state_dict(self.rate_net.state_dict())
        new._feat_mean = None if self._feat_mean is None else self._feat_mean.copy()
        new._feat_std = None if self._feat_std is None else self._feat_std.copy()
        return new
