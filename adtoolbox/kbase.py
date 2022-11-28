"""This module inludes the functionalities to populate the ADToolbox kbase. It is based on sqlite3"""
import os
import sqlite3
from sqlite3 import Error
from __init__ import Main_Dir
import rich



def connect_to_db(db_file):
    """ create a database connection to the ADToolbox kbase """
    connect = None
    try:
        connect = sqlite3.connect(db_file)
    except Error as e:
        print(e)
    return connect

def execute_query(connection, query):
    cursor = connection.cursor()
    try:
        cursor.execute(query)
        connection.commit()
    except Error as e:
        print(f"The error '{e}' occurred")

def read_query(connection, query):
    cursor = connection.cursor()
    result = None
    try:
        cursor.execute(query)
        result = cursor.fetchall()
        return result
    except Error as e:
        print(f"The error '{e}' occurred")

def create_study_table(address):
    create_study_table_query = """ CREATE TABLE IF NOT EXISTS study (
                                        id integer PRIMARY KEY,
                                        name text NOT NULL,
                                        type text NOT NULL,
                                        reference text NOT NULL); """
    # create a database connection
    if not os.path.exists(address.split("Studies.sqlite")[0]):
        os.mkdir(address.split("Studies.sqlite")[0])
    
    connection = connect_to_db(address)
    if connection is not None:
        # create study table
        execute_query(connection, create_study_table_query,)
        rich.print("[bold green]The study database was created successfully!")
        
    else:
        print("Error! cannot create the database connection.")

def add_study(address, name, type, reference):
    add_study_query = f""" INSERT INTO study(name, type, reference)
                            VALUES('{name}', '{type}', '{reference}'); """
    # create a database connection
    connection = connect_to_db(address)
    if connection is not None:
        # add study
        execute_query(connection, add_study_query)
    else:
        print("Error! cannot create the database connection.")

def get_studies(address):
    get_studies_query = "SELECT * FROM study"
    # create a database connection
    connection = connect_to_db(address)
    if connection is not None:
        # get studies
        studies = read_query(connection, get_studies_query)
        return studies
    else:
        print("Error! cannot create the database connection.")

def create_metagenoics_studies_table(address):
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
        execute_query(connection, create_metagenomics_studies_table_query)
        rich.print("[bold green]The metagenomics studies database was created successfully!")
    else:
        print("Error! cannot create the database connection.")

def add_metagenomics_study(address, name, study_id, type, microbiome,SRA_accession, reference,comments):
    add_metagenomics_study_query = f""" INSERT INTO metagenomics_studies(name, study_id, type, microbiome,SRA_accession, reference,comments)
                            VALUES('{name}', '{study_id}', '{type}', '{microbiome}','{SRA_accession}', '{reference}','{comments}'); """
    # create a database connection
    connection = connect_to_db(address)
    if connection is not None:
        # add metagenomics study
        execute_query(connection, add_metagenomics_study_query)
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
    else:
        print("Error! cannot create the database connection.")

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
        execute_query(connection, create_experimental_data_references_table_query)
    else:
        print("Error! cannot create the database connection.")

def add_experimental_data_reference(address, name, study_id, type, reference, comments):
    add_experimental_data_reference_query = f""" INSERT INTO experimental_data_references(name, study_id, type, reference, comments)
                            VALUES('{name}', '{study_id}', '{type}', '{reference}', '{comments}'); """
    # create a database connection
    connection = connect_to_db(address)
    if connection is not None:
        # add experimental data reference
        execute_query(connection, add_experimental_data_reference_query)
    else:
        print("Error! cannot create the database connection.")

