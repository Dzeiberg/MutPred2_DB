import os
import fire
from pathlib import Path
from models.sql_connection import SQL_Connection
from models.job_processor import Processor
from tqdm import tqdm

def insert_job_results(jobs_filepath=None,jobs_dirs_base=None,job_path=None, db_config_file="sql_configs.yaml", db_config_name="Local"):
    """
    Insert MutPred2 results into database

    Requires either jobs_dirs_base or job_path

    Optional Arguments
    - jobs_filepath : path to file containing absolute paths to list of jobs to insert
    - jobs_dirs_base : path to directory containing folders for each job to be processed
    - job_path : path to job to be processed
    """
    if (jobs_dirs_base is not None and job_path is not None and jobs_filepath is not None) or \
        (jobs_dirs_base is None and job_path is None and jobs_filepath is None):
        raise ValueError("specify either jobs_dirs_base or job_path or jobs_filepath")
    elif jobs_filepath is not None:
        jobs = [l.strip() for l in open(jobs_filepath).readlines()]
        jobs = [j for j in jobs if os.path.isdir(j)]
    elif jobs_dirs_base is not None:
        print(f"processing jobs in directory {jobs_dirs_base}")
        jobs_dirs_base = Path(jobs_dirs_base)
        jobs = [jobs_dirs_base/j for j in os.listdir(jobs_dirs_base)]
    else:
        print(f"processing job {job_path}")
        jobs = [Path(job_path)]
    
    con = SQL_Connection(db_config_name, db_config_file)
    processor = Processor(con)
    for job in tqdm(jobs):
        processor.process_job(job)

if __name__ == "__main__":
    fire.Fire(insert_job_results)