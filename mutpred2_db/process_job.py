from models.job_processor import Processor
from models.sql_connection import SQL_Connection, initialize_mechanisms
from fire import Fire
# from joblib import Parallel, delayed

from tqdm import tqdm
from functools import partialmethod
# tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)

class FileCursor:
    def __init__(self, filepath):
        self.file = open(filepath, 'w')

    def __exit__(self):
        # self.file.close()
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


def process_job_list(job_list_file : str, sql_config_name : str, sql_config_file : str,**run_option_kwargs):
    with open(job_list_file,'r') as file:
        job_list = list(map(str.strip, file.readlines()))
    sql_connection = SQL_Connection(sql_config_name, sql_config_file)
    cursor, cnx = sql_connection.open()
    # cnx = FileConnection("initialize.sql")
    # cursor = cnx._cursor
    initialize_mechanisms(cursor,cnx)
    # cursor.file.close()
    def process_job(job_dir : str, sql_config_name : str, sql_config_file : str,**run_option_kwargs):
        # print(f"Processing {job_dir}")
        cursor, cnx = sql_connection.open()
        # cnx = FileConnection(f"{job_dir}/insert.sql")
        # cursor = cnx._cursor
        job_processor = Processor(cursor, cnx)
        job_processor.process(job_dir, write_to_db=True,**run_option_kwargs)
        sql_connection.close(cursor,cnx)
        # cursor.file.close()
    [process_job(job,
                sql_config_name,
                sql_config_file,
                **run_option_kwargs) for job in tqdm(job_list)]
if __name__ == "__main__":
    Fire(process_job_list)


    def close(self):
        pass