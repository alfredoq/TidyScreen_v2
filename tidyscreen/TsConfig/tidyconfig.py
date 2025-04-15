# Here I will include all actions aimed to modify the main TidyScreen database that is installed in the environment
import sys

class MainDbConfigs:

    def __init__(self):
        self.main_projects_db = sys.modules['tidyscreen_dbs'].__file__.replace('__init__.py','projects_database.db')

    def list_projects(self):
        pass