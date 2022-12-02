"""This module inludes the functionalities to populate the ADToolbox kbase. It is based on sqlite3"""
import os
import sqlite3
from sqlite3 import Error
from __init__ import Main_Dir
import rich



def connect_to_db(db_file):
    """ create a database connection to the ADToolbox kbase """
    
    try:
        connect = sqlite3.connect(db_file)
    except:
        connect = None
    return connect
def execute_query(connection, query,address):
    cursor = connection.cursor()
    try:
        cursor.execute(query)
        connection.commit()
    
    except Error as e:
        print(f"The error '{e}' occurred")

def read_query(connection, query):
    
    try:
        cursor = connection.cursor()
        cursor.execute(query)
        result = cursor.fetchall()
        return result
    except Error as e:
        print(f"The error '{e}' occurred")

def _create_study_table(address):
    create_study_table_query = """ CREATE TABLE IF NOT EXISTS study (
                                        id integer PRIMARY KEY,
                                        name text NOT NULL,
                                        type text NOT NULL,
                                        reference text NOT NULL); """
    # create a database connection
    if not os.path.exists(address.split("studies.sqlite")[0]):
        os.mkdir(address.split("studies.sqlite")[0])
    
    connection = connect_to_db(address)
    if connection is not None:
        # create study table
        execute_query(connection, create_study_table_query,address)
        

def _add_study(address, name, type, reference):
    add_study_query = f""" INSERT INTO study(name, type, reference)
                            VALUES('{name}', '{type}', '{reference}'); """
    # create a database connection
    connection = connect_to_db(address)
    if connection is not None:
        # add study
        execute_query(connection, add_study_query,address)

def get_studies(address):
    get_studies_query = "SELECT * FROM study"
    # create a database connection
    connection = connect_to_db(address)
    if connection is not None:
        # get studies
        studies = read_query(connection, get_studies_query)
        return studies
    


def _create_metagenoics_studies_table(address):
    create_metagenomics_studies_table_query = """ CREATE TABLE IF NOT EXISTS metagenomics_studies (
                                        metagenomics_study_id integer PRIMARY KEY,
                                        name text NOT NULL,
                                        study_id integer NOT NULL,
                                        type text NOT NULL,
                                        microbiome text NOT NULL,
                                        SRA_accession text NOT NULL,
                                        reference text NOT NULL,
                                        comments text NOT NULL,
                                        FOREIGN KEY(study_id) REFERENCES study(id));"""
    # create a database connection
    if not os.path.exists(address.split("studies.sqlite")[0]):
        os.mkdir(address.split("studies.sqlite")[0])
    connection = connect_to_db(address)
    if connection is not None:
        # create metagenomics studies table
        execute_query(connection, create_metagenomics_studies_table_query,address)
    else:
        print("Error! cannot create the database connection.")

def add_metagenomics_study(address, name, study_id, type, microbiome,SRA_accession, reference,comments):
    add_metagenomics_study_query = f""" INSERT INTO metagenomics_studies(name, study_id,type, microbiome,SRA_accession, reference,comments)
                            VALUES("{name}", "{study_id}", "{type}", "{microbiome}","{SRA_accession}", "{reference}","{comments}"); """
    # create a database connection
    connection = connect_to_db(address)
    if connection is not None:
        # add metagenomics study
        execute_query(connection, add_metagenomics_study_query,address)
        rich.print(f"[bold green]The metagenomics study '{name}' was added successfully!")
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
                                        study_id integer NOT NULL,
                                        FOREIGN KEY (study_id) REFERENCES study (id),
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

def add_experimental_data_reference(address, name, study_id, type, reference, comments):
    add_experimental_data_reference_query = f""" INSERT INTO experimental_data_references(name, study_id, type, reference, comments)
                            VALUES('{name}', '{study_id}', '{type}', '{reference}', '{comments}'); """
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

