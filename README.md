# MutPred2_DB
> Tools for managing a MySQL database of pre-computed MutPred2 predictions

## Install Dependencies
```console
foo@bar:~$ python -m pip install -r requirements.txt
```

## Create MySQL Database schema

1. Create MySQL Database `MutPred2_DB`
2. Create yaml file with database configuration information in base directory, e.g.
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

## Updating Database with MutPred2 Results
```console
foo@bar:~$ python insert_job_results --job_path=~/path/to/job/directory
```