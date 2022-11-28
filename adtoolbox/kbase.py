import sqlite3
from sqlite3 import Error
from __init__ import Main_Dir



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
        print("Query executed successfully")
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

