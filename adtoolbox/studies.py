"""This module inludes the functionalities to populate the ADToolbox kbase
TODO:
- kbase class design
"""
import os
from __init__ import Main_Dir
import rich
import configs
import dash
import dash.html as html
import pandas as pd
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc


def get_studies(config=configs.Kbase()):
    """Collects all of the studies in the databases
    
    args:
        config (configs.Kbase): The configuration file for studies database
    Returns:
        pd.DataFrame: The studies table
        
    """
    metagenomics_studies=pd.read_csv(config.metagenomics_studies,delimiter="\t")
    experimental_data_references=pd.read_csv(config.experimental_data_references,delimiter="\t")
    studies=pd.concat([metagenomics_studies,experimental_data_references],axis=0,join="inner")
    return studies



def add_metagenomics_study(metagenomics_studies_table:str, name:str, type:str, microbiome,SRA_accession:str, reference:str,comments:str):
    """Adds a metagenomics study to the database
    Args:
        metagenomics_studies_table (pd.DataFrame): The metagenomics studies table
        name (str): The name of the metagenomics study
        type (str): The type of the metagenomics study 16s or shotgun
        microbiome (str): The microbiome of the metagenomics study
        SRA_accession (str): The SRA accession of the metagenomics study
        reference (str): The reference of the metagenomics study
        comments (str): The comments of the metagenomics study
    Returns:
        pd.DataFrame: The updated metagenomics studies table
    """ 
    metagenomics_studies_table=metagenomics_studies_table.append({"Name":name,"Type":type,"Microbiome":microbiome,"SRA_accession":SRA_accession,"Reference":reference,"Comments":comments},ignore_index=True)
    rich.print(f"[bold green]The metagenomics study '{name}' was added successfully!")
    return metagenomics_studies_table

def delete_metagenomics_study(metagenomics_studies_table:pd.DataFrame, metagenomics_study_id:int):
    """Deletes a metagenomics study from the database
    Args:
        metagenomics_studies_table (pd.DataFrame): The metagenomics studies table
        metagenomics_study_id (str): The id of the metagenomics study to be deleted
    Returns:
        pd.DataFrame: The updated metagenomics studies table
        """
    metagenomics_studies_table.drop(labels=metagenomics_study_id,inplace=True,axis=0)
    rich.print(f"[bold green]The metagenomics study with id '{metagenomics_study_id}' was deleted successfully!")
    return metagenomics_studies_table

def get_metagenomics_studies(Config=configs.Kbase()):
    """Collects all of the metagenomics studies in the databases"""
    metagenomics_studies=pd.read_csv(Config.metagenomics_studies,delimiter="\t")
    return metagenomics_studies

def write_metagenomics_studies(metagenomics_studies_table:pd.DataFrame,Config=configs.Kbase()):
    metagenomics_studies_table.to_csv(Config.metagenomics_studies,index=False,sep="\t")

def create_metagenomics_study_table(Config=configs.Kbase()):
    metagenomics_studies_table={"Name":[],"Type":[],"Microbiome":[],"SRA_accession":[],"Reference":[],"Comments":[]}
    metagenomics_studies_table=pd.DataFrame(metagenomics_studies_table)
    metagenomics_studies_table.to_csv(Config.metagenomics_studies,index=False,sep="\t")
    rich.print(f"[bold green]The metagenomics studies table was created successfully!")

def create_experimental_data_references_table(Config=configs.Kbase()):
    experimental_data_references_table={"Name":[],"Type":[],"Reference":[],"Comments":[]}
    experimental_data_references_table=pd.DataFrame(experimental_data_references_table)
    experimental_data_references_table.to_csv(Config.experimental_data_references,index=False,sep="\t")
    rich.print(f"[bold green]The experimental data references table was created successfully!")

def _initialize_databases(Config=configs.Kbase()):
    if not os.path.exists(Config.base_dir):
        os.makedirs(Config.base_dir)
    if not os.path.exists(Config.metagenomics_studies):
        create_metagenomics_study_table(Config)
    if not os.path.exists(Config.experimental_data_references):
        create_experimental_data_references_table(Config)
    
def dash_app(configs:configs.Kbase,table:str='all') -> None:
    """Main function of the app."""
    # Create the app.
    app = dash.Dash(__name__,external_stylesheets=[dbc.themes.BOOTSTRAP])
    _initialize_databases(configs)

    if table=='all':
        studies=get_studies(configs)
        metagenomics_studies=get_metagenomics_studies(configs)
        if studies.shape[0]:
            studies=studies.to_dict(orient='records')
        else:
            studies=[{"Index":"","Name":"","Type":"","Reference":""}]
            
        if metagenomics_studies.shape[0]:
            metagenomics_studies=metagenomics_studies.to_dict(orient='records')
        else:
            metagenomics_studies=[{"Index":"","Name":"","Type":"","Microbiome":"","SRA_accession":"","Reference":"","Comments":""}]
            

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
                                     export_format='tsv',
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

            if rows:
                write_metagenomics_studies(pd.DataFrame(rows),configs)
            else:
                rows=[{"Index":"","Name":"","Type":"","Microbiome":"","SRA_accession":"","Reference":"","Comments":""}]
                write_metagenomics_studies(pd.DataFrame(columns=rows[0].keys()),configs)
            studies=get_studies(configs)
            if studies.shape[0]:
                studies=studies.to_dict(orient='records')
            else:
                studies=[{"Index":"","Name":"","Type":"","Reference":""}]

            
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


if __name__ == '__main__':
    get_studies(configs.Kbase())