from sqlalchemy import create_engine
import yaml

class SQL_Connection:
    def __init__(self,config_name, config_file):
        with open(config_file,'r') as file:
            configs = yaml.safe_load(file)
        cfg = configs[config_name]
        self.user = cfg['user']
        self.password = cfg['password']
        self.host = cfg['host']
        self.port = cfg['port']
        self.database = cfg['database']

    def get_engine(self):
        engine = create_engine(f'mysql+mysqlconnector://{self.user}:{self.password}@{self.host}:{self.port}/{self.database}', echo=False)
        return engine