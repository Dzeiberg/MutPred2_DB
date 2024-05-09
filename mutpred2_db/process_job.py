from models.job_processor import Processor
from models.sql_connection import SQL_Connection, initialize_mechanisms
from fire import Fire
from sqlalchemy.exc import DatabaseError
from tqdm import tqdm
# tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)

class FileCursor:
    def __init__(self, filepath):
        self.file = open(filepath, 'w')

    def __exit__(self):
        super().__exit__()

    def execute(self, query, data=None):
        if data is not None:
            query = query.replace("%s","{}")
            w = query.format(*data)
        else:
            w = query
        self.file.write(w)
        self.file.write("\n")

class FileConnection:
    def __init__(self, filepath):
        self._cursor = FileCursor(filepath)

    @property
    def cursor(self):
        return self._cursor

    def commit(self):
        self._cursor.file.flush()


def process_job_list(sql_config_name : str, sql_config_file : str,
                    job_list_file : str|None=None,
                    job_path : str|None=None,**run_option_kwargs):
    if job_list_file is not None:
        with open(job_list_file,'r') as file:
            job_list = list(map(str.strip, file.readlines()))
    elif job_path is not None:
        job_list = [job_path, ]
    else:
        raise ValueError("Either job_list_file or job_path must be provided")
    sql_connection = SQL_Connection(sql_config_name, sql_config_file)
    cursor, cnx = sql_connection.open()
    initialize_mechanisms(cursor,cnx)
    def process_job(cursor, cnx, job_dir : str, sql_config_name : str, sql_config_file : str,**run_option_kwargs):
        job_processor = Processor(cursor, cnx)
        job_processor.process(job_dir, write_to_db=True,**run_option_kwargs)
    failed_jobs = []
    for job in tqdm(job_list):
        try:
            process_job(cursor, cnx, job, sql_config_name, sql_config_file, **run_option_kwargs)
        except DatabaseError as e:
            print(f"JOB_ERROR: Error processing {job}")
            failed_jobs.append(job)
            print(e)
            continue
    sql_connection.close(cursor,cnx)
    if len(failed_jobs) > 0:
        print("Failed jobs:")
        for job in failed_jobs:
            print(job)
    else:
        print("All jobs processed successfully")
if __name__ == "__main__":
    Fire(process_job_list)
