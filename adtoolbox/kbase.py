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
    studies={"name":[],"type":[],"reference":[],"comments":[]}
    studies['name'].extend(metagenomics_studies['name'].tolist())
    studies['type'].extend(["metagenomics"]*metagenomics_studies.shape[0])
    studies['reference'].extend(metagenomics_studies['reference'].tolist())
    studies['comments'].extend(metagenomics_studies['comments'].tolist())
    studies['name'].extend(experimental_data_references['name'].tolist())
    studies['type'].extend(["experimental"]*experimental_data_references.shape[0])
    studies['reference'].extend(experimental_data_references['reference'].tolist())
    studies['comments'].extend(experimental_data_references['comments'].tolist())
    studies=pd.DataFrame(studies)




def _create_metagenoics_studies_table(address):
    create_metagenomics_studies_table_query = """ CREATE TABLE IF NOT EXISTS metagenomics_studies (
                                        metagenomics_study_id integer PRIMARY KEY,
                                        name text NOT NULL,
                                        type text NOT NULL,
                                        microbiome text NOT NULL,
                                        SRA_accession text NOT NULL,
                                        reference text NOT NULL,
                                        comments text NOT NULL,);"""
    # create a database connection
    if not os.path.exists(address.split("studies.sqlite")[0]):
        os.mkdir(address.split("studies.sqlite")[0])
    connection = connect_to_db(address)
    if connection is not None:
        # create metagenomics studies table
        execute_query(connection, create_metagenomics_studies_table_query,address)
    else:
        print("Error! cannot create the database connection.")

def add_metagenomics_study(address, name, type, microbiome,SRA_accession, reference,comments):
    add_metagenomics_study_query = f""" INSERT INTO metagenomics_studies(name,type, microbiome,SRA_accession, reference,comments)
                            VALUES("{name}", "{type}", "{microbiome}","{SRA_accession}", "{reference}","{comments}"); """
    # create a database connection
    connection = connect_to_db(address)
    if connection is not None:
        # add metagenomics study
        execute_query(connection, add_metagenomics_study_query,address)
        rich.print(f"[bold green]The metagenomics study '{name}' was added successfully!")
    else:
        print("Error! cannot create the database connection.")

def delete_metagenomics_study(address, metagenomics_study_id):
    study = f""" SELECT FROM metagenomics_studies WHERE metagenomics_study_id = {metagenomics_study_id}; """
    
    delete_metagenomics_study_query = f""" DELETE FROM metagenomics_studies WHERE metagenomics_study_id = {metagenomics_study_id}; """
    # create a database connection
    connection = connect_to_db(address)
    if connection is not None:
        # delete metagenomics study

        rich.print(f"[bold green]The metagenomics study with id '{metagenomics_study_id}' was deleted successfully!")
    else:
        print("Error! cannot create the database connection.")

def get_metagenomics_studies(address):
    get_metagenomics_studies_query = "SELECT * FROM metagenomics_studies"
    # create a database connection
    connection = connect_to_db(address)
    if connection is not None:
        # get metagenomics studies
        metagenomics_studies = read_query(connection, get_metagenomics_studies_query)
        return metagenomics_studies


def create_experimental_data_references_table(address):
    create_experimental_data_references_table_query = """ CREATE TABLE IF NOT EXISTS experimental_data_references (
                                        experimental_data_reference_id integer PRIMARY KEY,
                                        name text NOT NULL,
                                        type text NOT NULL,
                                        reference text NOT NULL,
                                        comments text NOT NULL
                                    ); """
    # create a database connection
    connection = connect_to_db(address)
    if connection is not None:
        # create experimental data references table
        execute_query(connection, create_experimental_data_references_table_query,address)
    else:
        print("Error! cannot create the database connection.")

def add_experimental_data_reference(address, name, type, reference, comments):
    add_experimental_data_reference_query = f""" INSERT INTO experimental_data_references(name, type, reference, comments)
                            VALUES('{name}', '{type}', '{reference}', '{comments}'); """
    # create a database connection
    connection = connect_to_db(address)
    if connection is not None:
        # add experimental data reference
        execute_query(connection, add_experimental_data_reference_query,address)
    else:
        print("Error! cannot create the database connection.")

def _get_last_study_id(address):
    get_last_study_id_query = "SELECT id FROM study ORDER BY id DESC LIMIT 1"
    # create a database connection

    connection = connect_to_db(address)

    if connection is not None:
        # get last study id
        last_study_id = read_query(connection, get_last_study_id_query)
        return last_study_id

def _delete_study(address, study_id):
    delete_study_query = f""" DELETE FROM study WHERE id = {study_id}; """
    # create a database connection
    connection = connect_to_db(address)
    if connection is not None:
        # delete study
        execute_query(connection)

def dash_app(configs:Configs.Kbase,table:str='all') -> None:
    """Main function of the app."""
    # Create the app.
    app = dash.Dash(__name__,external_stylesheets=[dbc.themes.BOOTSTRAP])
    studies_df=[]

    metagenomics_studies_table=[]
    experimental_data_references_table={"Index":[],
                                        "Name":[],
                                        "Study":[],
                                        "Type":[],
                                        "Reference":[],
                                        "Comments":[]}
    if table=='all':
        studies=get_studies(configs.studies)
        metagenomics_studies=get_metagenomics_studies(configs.studies)
        if studies:
            for study in studies:
                studies_df.append({"Index":study[0],"Name":study[1],"Type":study[2],"Reference":study[3]})
        else:
            studies_df.append({"Index":'',"Name":'',"Type":'',"Reference":''})
            
        if metagenomics_studies:
            for metagenomics_study in metagenomics_studies:
                metagenomics_studies_table.append({"Index":metagenomics_study[0],"Name":metagenomics_study[1],"Study":metagenomics_study[2],"Type":metagenomics_study[3],"Microbiome":metagenomics_study[4],"SRA_accession":metagenomics_study[5],"Reference":metagenomics_study[6],"Comments":metagenomics_study[7]})
        else:
            metagenomics_studies_table.append({"Index":'',"Name":'',"Study":'',"Type":'',"Microbiome":'',"SRA_accession":'',"Reference":'',"Comments":''})
             






    app.layout = html.Div(children=[html.P("All Studies"),
                                    dash.dash_table.DataTable(studies_df,
                                    columns=[{"name": i, "id": i} for i in studies_df[0].keys()],
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
                                    dash.dash_table.DataTable(metagenomics_studies_table,
                                     id='metagenomics_studies_table',
                                        columns=[{"name": i, "id": i} for i in metagenomics_studies_table[0].keys()],
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
    Output('metagenomics_studies_table', 'data'),
    Output('studies_table',"data"),
    Input('submit-metagenomics-study', 'n_clicks'),
    State('metagenomics_studies_table', 'data'),
    State('metagenomics_studies_table', 'columns'),
    metagenomics_studies,
    prevent_initial_call=True
    )
    def submit_metag_row(n_clicks, rows, columns,metagenomics_studies=metagenomics_studies):
        if n_clicks > 0:
            prev=set([i[0] for i in metagenomics_studies]) 
            cur=set([i["Index"] for i in rows])
            to_be_removed=prev-cur
            to_be_added=cur-prev
            if to_be_added:
                for i in to_be_removed:
                    delete_metagenomics_study(configs.studies,i)
            
            for index,item in enumerate(to_be_added):
                add_metagenomics_study(configs.studies,rows[index+num_metagenomics_studies]['Name'],rows[index+num_metagenomics_studies]['Type'],rows[index+num_metagenomics_studies]['Microbiome'],rows[index+num_metagenomics_studies]['SRA_accession'],rows[index+num_metagenomics_studies]['Reference'],rows[index+num_metagenomics_studies]['Comments'])

            
            
            studies=get_studies(configs.studies)
            metagenomics_studies=get_metagenomics_studies(configs.studies)
            # _studies_df=[]
            # _metagenomics_studies_table=[]
            
            # if not studies:
            #     studies=[]
            #     _create_study_table(configs.studies)
            for study in studies:
                _studies_df.append({"Index":study[0],"Name":study[1],"Type":study[2],"Reference":study[3]})
            
            if not metagenomics_studies:
                metagenomics_studies=[]
                _create_metagenoics_studies_table(configs.studies)

            for metagenomics_study in metagenomics_studies:
                _metagenomics_studies_table.append({"Index":metagenomics_study[0],"Name":metagenomics_study[1],"Study":metagenomics_study[2],"Type":metagenomics_study[3],"Microbiome":metagenomics_study[4],"SRA_accession":metagenomics_study[5],"Reference":metagenomics_study[6],"Comments":metagenomics_study[7]})
            
            num_metagenomics_studies=len(_metagenomics_studies_table)
            num_studies=len(studies_df)
            

        return studies_rows,metagenomics_rows,html.Div([dbc.Alert(f"{len(to_be_removed)} Studies were removed successfully!", color="success"),
                        dbc.Alert(f"{len(to_be_added)} Studies were added successfully!", color="success")])





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