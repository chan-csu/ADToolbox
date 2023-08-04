import plotly.express as px
from dash import Dash, html, dcc
import plotly.graph_objects as go
from dash import Dash, dcc,ctx, html, Input, Output,dash_table,callback,no_update
from dash.dash_table.Format import Format, Scheme, Sign, Symbol
import pandas as pd
import dash_bootstrap_components as dbc
import json
from collections import OrderedDict
from sklearn.manifold import MDS,TSNE
import base64
from adtoolbox import core,adm,configs
import numpy as np

"""
samples file should be formatted like this:

{
    "sample1":
    {
        "src":"Path/to/figure",
        "data":{
            "X_ch":0.1,
            ....
                },
        "group":"group1"
    },
    "sample2":
    {
        "src":"Path/to/figure",
        "data":{
            "X_ch":0.1,
            ....
                },
        "group":"group2"
    },
    ...
}

"""


# Small molecule drugbank dataset
# Source: https://raw.githubusercontent.com/plotly/dash-sample-apps/main/apps/dash-drug-discovery/data/small_molecule_drugbank.csv'
# data_path = '/Users/parsaghadermarzi/Desktop/Academics/Projects/Anaerobic_Digestion_Modeling/experiment/datasets/small_molecule_drugbank.csv'

# df = pd.read_csv(data_path, header=0, index_col=0)

# fig = go.Figure(data=[
#     go.Scatter(
#         x=df["LOGP"],
#         y=df["PKA"],
#         mode="markers",
#         marker=dict(
#             colorscale='viridis',
#             color=df["MW"],
#             size=df["MW"],
#             colorbar={"title": "Molecular<br>Weight"},
#             line={"color": "#444"},
#             reversescale=True,
#             sizeref=45,
#             sizemode="diameter",
#             opacity=0.8,
#         )
#     )
# ])

# turn off native plotly.js hover effects - make sure to use
# hoverinfo="none" rather than "skip" which also halts events.
# fig.update_traces(hoverinfo="none", hovertemplate=None)

# fig.update_layout(
#     xaxis=dict(title='Log P'),
#     yaxis=dict(title='pkA'),
#     plot_bgcolor='rgba(255,255,255,0.1)'
# )
_compounds=["S_bu","S_ac","S_lac","S_et"]




TOTALL_MICROBIAL_COD=10
DURATION=5

def process_input(json_data,method="MDS",normalize=True,n_components=2):
    samples=list(json_data.keys())
    df=pd.DataFrame([json_data[i]["data"] for i in samples],index=samples)
    df=df.div(df.sum(axis=1), axis=0)
    df_copy=df.copy()
    X=df.to_numpy()
    if method=="MDS":
        mds=MDS(n_components=n_components,max_iter=10000)
        X=mds.fit_transform(X)
    elif method=="TSNE":
        tsne=TSNE(n_components=n_components,n_iter=10000)
        X=tsne.fit_transform(X)
    df=pd.DataFrame(X,index=df_copy.index,columns=["x","y"] if n_components==2 else ["x","y","z"])
    df["group"]=[json_data[i]["group"] for i in samples]
    df["src"]=[json_data[i]["src"] for i in samples]
    df_copy["group"]=df["group"]
    return df,df_copy


def main(json_file:str,method:str="MDS",normalize:bool=True,n_components:int=2):
    with open(json_file,"r") as f:
        json_data=json.load(f)

    df,df_copy=process_input(json_data,method=method,normalize=normalize,n_components=n_components)
    
    if n_components==2:
        fig=px.scatter(df,x="x",y="y",color="group",hover_name=process_input(json_data,method="MDS",normalize=True)[0].index,hover_data=["src"])
    else:
        fig=px.scatter_3d(df,x="x",y="y",z="z",color="group",hover_name=process_input(json_data,method="MDS",normalize=True)[0].index,hover_data=["src"])
        
    fig.update_traces(hoverinfo="none", hovertemplate=None)
    fig.update_traces(marker=dict(size=15,))
    #increase the font size 
    fig.update_layout(
        font=dict(
            size=19,
        )
    )
    # convert df_copy to tall format

    df_copy_=df_copy.reset_index().melt(id_vars=["index","group"],var_name="compound",value_name="concentration")

    fig_cod_portions=px.box(df_copy_,x="compound",y="concentration",color="group",)
    fig_cod_portions.update_layout(font=dict(size=19,))
    


    #### Modeling Panel ####

    base_feed=adm.DEFAULT_FEED

    with open(configs.Database().model_parameters,"r") as f:
        mp=json.load(f)
    with open(configs.Database().base_parameters,"r") as f:
        bp=json.load(f)
    with open(configs.Database().species,"r") as f:
        sp=json.load(f)
    with open(configs.Database().reactions,"r") as f:
        rxns=json.load(f)
    with open(configs.Database().initial_conditions,"r") as f:
        ic=json.load(f)
    with open(configs.Database().inlet_conditions,"r") as f:
        inc=json.load(f)


    base_model=adm.Model(model_parameters=mp,
                         base_parameters=bp,
                         initial_conditions=ic,
                         inlet_conditions=inc,
                         feed=base_feed,
                         species=sp,
                         reactions=rxns,
                         ode_system=adm.modified_adm_ode_sys,
                         build_stoichiometric_matrix=adm.build_modified_adm_stoichiometric_matrix,
                         )
    sol_df={}

    for i in json_data.keys():
        prep={"sample":i,
              "host":json_data[i]["group"],}
        ic_=ic.copy()
        for j in json_data[i]["data"].keys():
            ic_[j]=df_copy.loc[i,j]*TOTALL_MICROBIAL_COD
        model=adm.Model(model_parameters=mp,
                         base_parameters=bp,
                         initial_conditions=ic_,
                         inlet_conditions=inc,
                         feed=base_feed,
                         species=sp,
                         reactions=rxns,
                         ode_system=adm.modified_adm_ode_sys,
                         build_stoichiometric_matrix=adm.build_modified_adm_stoichiometric_matrix,
                         )
        sol=model.solve_model(np.linspace(0, DURATION, 1000))
        prep.update({"sol":sol.y})
        sol_df.update({i:prep})
    comp_inds=[sp.index(i) for i in _compounds]
    solution=[]
    for sample in sol_df.keys():
        for ind in comp_inds:
            solution.append({"sample":sample,
                             "host":sol_df[sample]["host"],
                             "compound":sp[ind],
                             "concentration":sol_df[sample]["sol"][ind,-1]})


    solution_df=pd.DataFrame(solution)
    conc_plot=px.box(solution_df,x="compound",y="concentration",color="host",points="all")
    ## add all points and make the font size bigger
    conc_plot.update_traces(marker=dict(size=6,),
                            pointpos=0,)
    conc_plot.update_layout(font=dict(size=19,))



    app = Dash(__name__,external_stylesheets=[dbc.themes.FLATLY])
    styles={
                'table_width': '95%',
                'padding-left': '20px',
                'container_width': '95%'
            }
    app.layout =html.Div([dbc.Container(
                            html.H1("ADToolbox Metagenomics Report",
                                    style={"font-size":"70px", "padding-top":"50px"}),
                                    className="text-white bg-primary",
                                    style={"height":"300px","text-align": "center"},
                                    fluid=True),
                            dbc.Container([
                                            dbc.Row(
                                                [dbc.Card([
                                                    html.H2("High-Level Relationship Between Samples",style={"font-size":"40px","padding-top":"5px"}),
                                                    dcc.Graph(figure=fig, 
                                                    id='MDS_Plot',
                                                    clear_on_unhover=True,
                                                    style={
                                                    "height":"800px",
                                                    "padding-top":"20px",
                                                    "padding-bottom":"20px",
                                                    "width":styles['container_width'],
                                                    "align-items": "center",
                                                    "display": "inline-block",}
                                                            ),
                                                    dcc.Tooltip(id="cod-tooltip"),],
                                                style={"align-items": "center"},
                                                class_name='bg-light'

                                                        ),
                                                 dbc.Card([
                                                    html.H2("Portion of COD Concentrations of ADM Microbial Species",style={"font-size":"40px","padding-top":"5px"}),
                                                    dcc.Graph(figure=fig_cod_portions, 
                                                    id='box_plot',
                                                    clear_on_unhover=True,
                                                    style={
                                                    "height":"800px",
                                                    "padding-top":"20px",
                                                    "padding-bottom":"20px",
                                                    "width":styles['container_width'],
                                                    "align-items": "center",
                                                    "display": "inline-block",}
                                                            )],
                                                style={"align-items": "center","padding-top":"10px" },
                                                class_name='bg-light'),

                                                dbc.Card([
                                                    dbc.Row([
                                                    dbc.Col([
                                                    html.H2("Select Species to Show",style={"width":"300px","font-size":"25px","padding-bottom":"5px","align-items":"center",},),
                                                    dcc.Checklist(sp, id="species",
                                                                  value=_compounds,
                                                                  style={"font-size":"22px",
                                                                         "padding-left":"5px",
                                                                         "padding-buttom":"100px",
                                                                         "height":"150px",
                                                                         "overflow-x":"scroll",
                                                                         "width":"250px",
                                                                          "background-color":"#f8f9fa",
                                                                          "color":"#808080"
                                                                         },

                                                                  inline=False,
                                                                  labelStyle={'display': 'block'}),],style={"padding-top":"5px",
                                                                                                            "height":"280px",
                                                                                                            "padding-bottom":"50px"}),
                                                    dbc.Col(
                                                        [
                                                        html.H2("Feed Composition",style={"width":"300px","font-size":"25px","padding-bottom":"5px","align-items":"center",},),
                                                        dash_table.DataTable(
                                                                id='feed_comp',
                                                                columns=[{"name": i, "id": i} for i in ["Lipid","Carbohydrates","Proteins"]],
                                                                data=[
                                                                    {"Lipid":10,"Carbohydrates":20,"Proteins":10}
                                                                ],
                                                                editable=True,
                                                                style_header={
                                                                    'backgroundColor': 'gray',
                                                                    'fontWeight': 'bold',
                                                                    "fontFamily": "system-ui",
                                                                    'color':'white',
                                                                    'font-size':'20px',
                                                                    'height': '60px',
                                                                    'width':'200px',
                                                                    'text-align':'center',

                                                                },
                                                                style_data={
                                                                    'backgroundColor': 'white',
                                                                    'color': 'black',
                                                                    'font-size':'20px',
                                                                    'height': '60px',
                                                                    'text-align':'center',
                                                                    'width':'150px',}
                                                                  ),
                                                        dbc.Spinner(html.Div(id="loading-output"),
                                                                    color="primary",type="grow",
                                                                    fullscreen=True,
                                                                    fullscreen_style={"background-color":"rgba(255,255,255,0.8)"}
                                                                    ),
                                                        ]
                                                        , style={"width":"150px"}),
                                                     dbc.Col([
                                                         html.H3("Feed Characteristics",style={"font-size":"25px","padding-bottom":"5px","align-items":"center",},),
                                                         dash_table.DataTable(
                                                                id='feed_char',
                                                                columns=[{"name": i, "id": i} for i in ["TSS","Total COD","Particulate Inerts","Soluble Inerts"]],

                                                                data=[
                                                                    {"TSS":10,"Total COD":20,"Particulate Inerts":10,"Soluble Inerts":10}
                                                                ],
                                                                editable=True,
                                                                style_header={
                                                                    'backgroundColor': 'gray',
                                                                    'fontWeight': 'bold',
                                                                    "fontFamily": "system-ui",
                                                                    'color':'white',
                                                                    'font-size':'20px',
                                                                    'height': '60px',
                                                                    'width':'200px',
                                                                    'text-align':'center',

                                                                },
                                                                style_data={
                                                                    'backgroundColor': 'white',
                                                                    'color': 'black',
                                                                    'font-size':'20px',
                                                                    'height': '60px',
                                                                    'text-align':'center',
                                                                    'width':'150px',}
                                                                  ),


                                                    ]),


                                                    ],class_name='text-white bg-primary',style={"width":"95%","padding-top":"25px","height":"250px" }),
                                                    dcc.Graph(figure=conc_plot,id="conc_plot",
                                                              style={"height":"800px","padding-top":"20px",
                                                                     "padding-bottom":"20px","width":styles['container_width'],
                                                                     "align-items": "center","display": "inline-block",})





                                                 ],style={"align-items": "center","padding-top":"5px" },class_name='bg-light'),

                                            ]),]),])





    @callback(
        Output("cod-tooltip", "show"),
        Output("cod-tooltip", "bbox"),
        Output("cod-tooltip", "children"),
        Input("MDS_Plot", "hoverData"),
    )
    def display_hover(hoverData):
        if hoverData is None:
            return False, no_update, no_update
        pt = hoverData["points"][0]
        bbox = pt["bbox"]
        img_src=base64.b64encode(open(pt["customdata"][0], 'rb').read()).decode('ascii')
        children = [
            html.Div([
                html.Img(src='data:image/png;base64,{}'.format(img_src), style={"width": "100%"}),
                # html.P(f"{form}"),
                # html.P(f"{desc}"),
            ], style={'width': '500px', 'white-space': 'normal'})
        ]

        return True, bbox, children


    @callback(
        Output("conc_plot", "figure"),
        Output("loading-output", "children"),
        Input("species", "value"),
        Input("feed_comp", "data"),
        Input("feed_char", "data")
    )
    def checkbox_update(value,data1,data2):
        if ctx.triggered[0]["prop_id"]=="feed_comp.data" or ctx.triggered[0]["prop_id"]=="feed_char.data":

            ic_["TSS"]=float(data2[0]["TSS"])
            ic_["TDS"]=float(data2[0]["Total COD"])-float(data2[0]["TSS"])
            feed=core.Feed(
                carbohydrates=float(data1[0]["Carbohydrates"]),
                proteins=float(data1[0]["Proteins"]),
                lipids=float(data1[0]["Lipid"]),
                tss=float(data2[0]["TSS"]),
                si=float(data2[0]["Soluble Inerts"]),
                xi=float(data2[0]["Particulate Inerts"]),)

            sol_df_={}
            for i in json_data.keys():
                prep={"sample":i,
                      "host":json_data[i]["group"],}
                for j in json_data[i]["data"].keys():
                    ic_[j]=df_copy.loc[i,j]*TOTALL_MICROBIAL_COD
                model=adm.Model(model_parameters=mp,
                                 base_parameters=bp,
                                 initial_conditions=ic_,
                                 inlet_conditions=inc,
                                 feed=feed,
                                 species=sp,
                                 reactions=rxns,
                                 ode_system=adm.modified_adm_ode_sys,
                                 build_stoichiometric_matrix=adm.build_modified_adm_stoichiometric_matrix,
                                 )
                sol=model.solve_model(np.linspace(0, DURATION, 1000))
                prep.update({"sol":sol.y})
                sol_df_.update({i:prep})
            comp_inds=[sp.index(i) for i in value]
            solution=[]
            for sample in sol_df_.keys():
                for ind in comp_inds:
                    solution.append({"sample":sample,
                                     "host":sol_df_[sample]["host"],
                                     "compound":sp[ind],
                                     "concentration":sol_df_[sample]["sol"][ind,-1]})


            solution_df=pd.DataFrame(solution)
            conc_plot=px.box(solution_df,x="compound",y="concentration",color="host",points="all")
            conc_plot.update_traces(marker=dict(size=6,),
                                    pointpos=0,)
            conc_plot.update_layout(font=dict(size=19,))
        else:
            inds=[sp.index(i) for i in value]
            solution=[]
            for sample in sol_df.keys():
                for ind in inds:
                    solution.append({"sample":sample,
                                     "host":sol_df[sample]["host"],
                                     "compound":sp[ind],
                                     "concentration":sol_df[sample]["sol"][ind,-1]})
            solution_df=pd.DataFrame(solution)
            conc_plot=px.box(solution_df,x="compound",y="concentration",color="host",points="all")
            ## add all points and make the font size bigger
            conc_plot.update_traces(marker=dict(size=6,),
                                    pointpos=0,)
            conc_plot.update_layout(font=dict(size=19,))
        return conc_plot, ""



    @callback(Output("loading", "children"), Input("feed_comp", "data"))
    def input_triggers_spinner(value):
        return value
    
    app.run_server()



if __name__ == "__main__":
    main("/Users/parsaghadermarzi/Desktop/Academics/Projects/Anaerobic_Digestion_Modeling/experiment/full_data.json",method="MDS",normalize=True,n_components=3)