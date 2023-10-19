# MutPred2_DB
> Tools for managing a MySQL database of pre-computed MutPred2 predictions

## Install Dependencies
```console
foo@bar:~$ python -m pip install -r requirements.txt
```

## Create MySQL Database schema

1. Create MySQL Database `MutPred2_DB`
2. Create `sql_configs.yaml` with database configuration information in base directory, e.g.
    ```yml
    Local:
        user: your_user
        host: localhost
        password: your_password
        database: MutPred2_DB
        port: 3306
    ```
3. Create Schema
    ```console
    foo@bar:~$ cd sql;
    foo@bar:~$ mysql --host=localhost --user=your_username --password=your_password  -e "building_tables.sql"
    ```
4. Seed Sequences and Sequence Mappings Tables
    ```console
    foo@bar:~$ python initialize_sequence_and_mapping_tables.py
    ```
## Updating Database with MutPred2 Results
### Single Job
```console
foo@bar:~$ python insert_job_results.py --job_path=~/path/to/job/directory --db_config_file=/path/to/sql_configs.yaml --db_config_name=Local
```
### Directory Containing Multiple Jobs
```console
foo@bar:~$ python insert_job_results.py --job_dirs_base=~/path/to/directory/of/multiple/jobs/ --db_config_file=/path/to/sql_configs.yaml --db_config_name=Local
```