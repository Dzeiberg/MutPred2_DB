from sqlalchemy import create_engine
import yaml

def sql_connect(config_name, config_file):
    with open(config_file,'r') as file:
        configs = yaml.safe_load(file)
    connection_config = configs[config_name]
    engine = create_engine(f'mysql+mysqlconnector://{connection_config["user"]}:{connection_config["password"]}@{connection_config["host"]}:{connection_config["port"]}/{connection_config["database"]}', echo=False)
    return engine

