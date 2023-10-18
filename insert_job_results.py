import fire
from pathlib import Path
from models.sql_connection import SQL_Connection
from models.job_processor import Processor

def insert_job_results(job_path, db_config_file="sql_configs.yaml", db_config_name="Local"):
    job_path = Path(job_path)
    con = SQL_Connection(db_config_name, db_config_file)
    processor = Processor(con)
    processor.process_job(job_path)

if __name__ == "__main__":
    insert_job_results()