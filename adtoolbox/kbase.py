"""This module inludes the functionalities to populate the ADToolbox kbase. It is based on sqlite3"""
import os
from __init__ import Main_Dir
import rich
import Configs
import dash
import dash.html as html
import pandas as pd
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc



def get_studies(Config=Configs.Kbase()):
    metagenomics_studies=pd.read_csv(Config.metagenomics_studies)
    experimental_data_references=pd.read_csv(Config.experimental_data_references)
    studies={"Name":[],"Type":[],"Reference":[],"Comments":[]}
    studies['Name'].extend(metagenomics_studies['Name'].tolist())
    studies['Type'].extend(["Metagenomics"]*metagenomics_studies.shape[0])
    studies['Reference'].extend(metagenomics_studies['Reference'].tolist())
    studies['Comments'].extend(metagenomics_studies['Comments'].tolist())
    studies['Name'].extend(experimental_data_references['Name'].tolist())
    studies['Type'].extend(["experimental"]*experimental_data_references.shape[0])
    studies['Reference'].extend(experimental_data_references['Reference'].tolist())
    studies['Comments'].extend(experimental_data_references['Comments'].tolist())
    return pd.DataFrame(studies)



def add_metagenomics_study(metagenomics_studies_table, name, type, microbiome,SRA_accession, reference,comments):
    metagenomics_studies_table=metagenomics_studies_table.append({"Name":name,"Type":type,"Microbiome":microbiome,"SRA_accession":SRA_accession,"Reference":reference,"Comments":comments},ignore_index=True)
    rich.print(f"[bold green]The metagenomics study '{name}' was added successfully!")
    return metagenomics_studies_table

def delete_metagenomics_study(metagenomics_studies_table, metagenomics_study_id):
    metagenomics_studies_table.drop(labels=metagenomics_study_id,inplace=True,axis=0)
    rich.print(f"[bold green]The metagenomics study with id '{metagenomics_study_id}' was deleted successfully!")
    return metagenomics_studies_table

def get_metagenomics_studies(Config=Configs.Kbase()):
    metagenomics_studies=pd.read_csv(Config.metagenomics_studies)
    return metagenomics_studies

def write_metagenomics_studies(metagenomics_studies_table,Config=Configs.Kbase()):
    metagenomics_studies_table.to_csv(Config.metagenomics_studies,index=False)

def create_metagenomics_study_table(Config=Configs.Kbase()):
    metagenomics_studies_table={"Name":[],"Type":[],"Microbiome":[],"SRA_accession":[],"Reference":[],"Comments":[]}
    metagenomics_studies_table=pd.DataFrame(metagenomics_studies_table)
    metagenomics_studies_table.to_csv(Config.metagenomics_studies,index=False)
    rich.print(f"[bold green]The metagenomics studies table was created successfully!")

def create_experimental_data_references_table(Config=Configs.Kbase()):
    experimental_data_references_table={"Name":[],"Type":[],"Reference":[],"Comments":[]}
    experimental_data_references_table=pd.DataFrame(experimental_data_references_table)
    experimental_data_references_table.to_csv(Config.experimental_data_references,index=False)
    rich.print(f"[bold green]The experimental data references table was created successfully!")

def _initialize_databases(Config=Configs.Kbase()):
    if not os.path.exists(Config.base_dir):
        os.makedirs(Config.base_dir)
    if not os.path.exists(Config.metagenomics_studies):
        create_metagenomics_study_table(Config)
    if not os.path.exists(Config.experimental_data_references):
        create_experimental_data_references_table(Config)
    
def dash_app(configs:Configs.Kbase,table:str='all') -> None:
    """Main function of the app."""
    # Create the app.
    app = dash.Dash(__name__,external_stylesheets=[dbc.themes.BOOTSTRAP])
    _initialize_databases(Configs.Kbase())

    if table=='all':
        studies=get_studies(configs)
        metagenomics_studies=get_metagenomics_studies(configs)
        if studies.shape[0]:
            studies=studies.to_json(orient='records')
        else:
            studies=[{"Index":'',"Name":'',"Type":'',"Reference":''}]
            
        if metagenomics_studies.shape[0]:
            metagenomics_studies=metagenomics_studies.to_json(orient='records')
        else:
            metagenomics_studies=[{"Index":'',"Name":'',"Study":'',"Type":'',"Microbiome":'',"SRA_accession":'',"Reference":'',"Comments":''}]
            

    app.layout = html.Div(children=[html.P("All Studies"),
                                    dash.dash_table.DataTable(studies,
                                    columns=[{"name": i, "id": i} for i in studies[0].keys()],
                                     id='studies_table',
                                     style_cell={'textAlign': 'left',
                                     'fontSize':15, 'font-family':'sans-serif',
                                     'maxWidth': '250px'},
                                     style_data={'whiteSpace': 'normal',
                                    'height': 'auto',},
                                    style_table={'overflowX': 'auto'},
                                    export_format='csv'),
                                    html.Br(),
                                    html.P("Metagenomics Studies"),
                                    dash.dash_table.DataTable(metagenomics_studies,
                                     id='metagenomics_studies_table',
                                        columns=[{"name": i, "id": i} for i in metagenomics_studies[0].keys()],
                                     style_cell={'textAlign': 'left',
                                     'fontSize':15, 'font-family':'sans-serif',
                                     'maxWidth': '250px'},
                                     style_data={'whiteSpace': 'normal',
                                                'height': 'auto'},
                                     style_table={'overflowX': 'auto'},
                                     export_format='csv',
                                    editable=True,
                                    row_deletable=True,
                                     ),
                                     html.Button('Add metagenomics study', id='add-metagenomics-study', n_clicks=0),
                                     html.Button('Submit metagenomics study', id='submit-metagenomics-study', n_clicks=0),
                                     html.Div(id='submission-status'),
    ])
    @app.callback(
    Output('metagenomics_studies_table', 'data'),
    Input('add-metagenomics-study', 'n_clicks'),
    State('metagenomics_studies_table', 'data'),
    State('metagenomics_studies_table', 'columns'),
    prevent_initial_call=True)
    def add_metag_row(n_clicks, rows, columns):
        if n_clicks > 0:
            rows.append({c['id']: '' for c in columns})
        return rows
    
    @app.callback(
    Output('submission-status', 'children'),
    Output('studies_table',"data"),
    Input('submit-metagenomics-study', 'n_clicks'),
    State('metagenomics_studies_table', 'data'),
    State('metagenomics_studies_table', 'columns'),
    prevent_initial_call=True
    )
    def submit_metag_row(n_clicks, rows, columns):
        if n_clicks > 0:

            write_metagenomics_studies(pd.DataFrame(rows),configs)
            studies=get_studies(configs)
            if studies.shape[0]:
                studies=studies.to_json(orient='records')
            else:
                studies=[{"Index":'',"Name":'',"Type":'',"Reference":''}]

            

        return html.Div([dbc.Alert(f"All of the changes were made successfully!", color="success")]),studies





    # @app.callback(
    # Output('adding-rows-graph', 'figure'),
    # Input('adding-rows-table', 'data'),
    # Input('adding-rows-table', 'columns'))
    # def display_output(rows, columns):
    #     return {
    #     'data': [{
    #         'type': 'heatmap',
    #         'z': [[row.get(c['id'], None) for c in columns] for row in rows],
    #         'x': [c['name'] for c in columns]
    #     }]
    # }
    app.run_server(debug=True)

dash_app(configs=Configs.Kbase(),table='all')