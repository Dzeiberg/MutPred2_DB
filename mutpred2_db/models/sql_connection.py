import mysql.connector
import yaml
from .sequence import Sequence
from .variant import Variant
from .mechanism import Mechanism
from .feature_set import Features_Set
from typing import List
from tqdm import tqdm
from joblib import Parallel, delayed
import pandas as pd

class SQL_Connection(object):
    def __init__(self,config_name, config_file):
        with open(config_file,'r') as file:
            configs = yaml.safe_load(file)
        cfg = configs[config_name]
        self.user = cfg['user']
        self.password = cfg['password']
        self.host = cfg['host']
        self.port = cfg['port']
        self.database = cfg['database']
        self.pool = mysql.connector.pooling.MySQLConnectionPool(pool_name="my_pool", pool_size=5,
                                                            user=self.user,
                                                            database=self.database,
                                                            password=self.password,
                                                            host=self.host,
                                                            port=self.port)

    def open(self,):
        conn = self.pool.get_connection()
        cursor = conn.cursor()
        return cursor, conn

    def close(self,conn,cursor):
        cursor.close()
        conn.close()

    # def get_engine(self):
    #     self.cnx = mysql.connector.connect(user=self.user, database=self.database, password=self.password, host=self.host, port=self.port)
    #     self.cursor = self.cnx.cursor()

    def __exit__(self):
        # self.cursor.close()
        # self.cnx.close()
        self.pool.close()
        super().__exit__()

def write_sequence(cursor, cnx,sequence : Sequence) -> int:
    add_sequence = ("INSERT INTO Protein "
            "(seq_hash, sequence) "
            "VALUES (%s, %s)"
            "ON DUPLICATE KEY UPDATE seq_hash=seq_hash, sequence=sequence;")
    data_sequence = (sequence.seq_hash, sequence.seq)
    cursor.execute(add_sequence, data_sequence)
    sequence_number = 0#cursor.lastrowid
    cnx.commit()
    return sequence_number

def write_variant( cursor, cnx,variant : Variant, **kwargs) -> int:
    do_commit = kwargs.get("do_commit",True)
    add_variant = ("INSERT INTO Variant "
            "(variant_id, seq_hash, reference_aa, position, alternate_aa, score, option_id) "
            "VALUES (%s, %s, %s, %s, %s, %s, %s)"
            "ON DUPLICATE KEY UPDATE variant_id=variant_id, seq_hash=seq_hash, reference_aa=reference_aa, position=position, alternate_aa=alternate_aa, score=score, option_id=option_id;")
    data_variant = (variant.variant_id, variant.seq_hash, variant.reference_aa, variant.position, variant.alternate_aa, variant.mutpred_score, variant.option_id)
    try:
        _ = cursor.execute(add_variant, data_variant)
    except Exception as e:
        print(data_variant)
        raise e
    if do_commit:
        cnx.commit()
    # cursor.execute('SELECT LAST_INSERT_ID()')
    return 0#last_insert_id[0]

def write_variants(cursor, cnx,variants : List[Variant]) -> List[int]:
    variant_numbers = []
    for variant in tqdm(variants,desc="Writing variants",leave=False):
        variant_numbers.append(write_variant(cursor, cnx,variant, do_commit=True))
    return variant_numbers

def query_runoption(cursor, cnx,compute_homology_profile : bool,
                    use_predicted_conservation_scores : bool,
                    skip_psi_blast : bool,
                    p_value_threshold : float) -> int:
    # query = ("SELECT option_id FROM RunOption "
    #             "WHERE compute_homology_profile = %s AND use_predicted_conservation_scores = %s AND skip_psi_blast = %s AND p_value_threshold = %s")
    # data = (compute_homology_profile, use_predicted_conservation_scores, skip_psi_blast, p_value_threshold)
    # cursor.execute(query, data)
    # option_id = cursor.fetchone()
    if compute_homology_profile and use_predicted_conservation_scores and ~skip_psi_blast and p_value_threshold == 1:
        option_id = [1,]
    else:
        option_id = None
    if option_id is None:
        raise NotImplementedError("not supporting this option as I'm using the write to file, rather than write to db")
        return insert_runoption(cursor, cnx,compute_homology_profile, use_predicted_conservation_scores, skip_psi_blast, p_value_threshold)
    if len(option_id) > 1:
        raise ValueError("More than one option_id returned")
    return option_id[0]

def insert_runoption(cursor, cnx, compute_homology_profile : bool,
                    use_predicted_conservation_scores : bool,
                    skip_psi_blast : bool,
                    p_value_threshold : float) -> int:
    add_option = ("INSERT INTO RunOption "
                    "(compute_homology_profile, use_predicted_conservation_scores, skip_psi_blast, p_value_threshold) "
                    "VALUES (%s, %s, %s, %s);")
    data_option = (compute_homology_profile, use_predicted_conservation_scores, skip_psi_blast, p_value_threshold)
    cursor.execute(add_option, data_option)
    option_id = 0#cursor.lastrowid
    cnx.commit()
    return option_id

def write_mechanism(cursor, cnx, mechanism : Mechanism|pd.Series, **kwargs) -> int:
    do_commit = kwargs.get("do_commit",True)
    add_mechanism = ("INSERT INTO VariantMechanism "
                        "(variant_id, mechanism_id, mechanism_type, altered_position, score, pvalue, description) "
                        "VALUES (%s, %s, %s, %s, %s, %s, %s)"
                        "ON DUPLICATE KEY UPDATE variant_id=variant_id, mechanism_id=mechanism_id, mechanism_type=mechanism_type, altered_position=altered_position, score=score, pvalue=pvalue, description=description;")
    data_mechanism = (mechanism['variant_id'],mechanism['mechanism_id'], mechanism['mechanism_type'], mechanism['altered_position'], mechanism['score'], mechanism['pvalue'], mechanism['description'])
    try:
        cursor.execute(add_mechanism, data_mechanism)
    except Exception as e:
        print(add_mechanism)
        print(data_mechanism)
        raise e
    variant_mechanism_id = 0#cursor.lastrowid
    if do_commit:
        cnx.commit()
    return variant_mechanism_id

def write_mechanisms(cursor, cnx,mechanisms : List[Mechanism]) -> List[int]:
    variant_mechanism_ids = [write_mechanism(cursor, cnx,mechanism, do_commit=False) for mechanism in tqdm(map(lambda m: m.to_dict(), mechanisms),
                                                                                                        total=len(mechanisms),
                                                                                                        desc="Writing mechanisms",
                                                                                                        leave=False)]
    cnx.commit()
    return variant_mechanism_ids

def initialize_mechanisms(cursor, cnx):
    for idx, mechanism in enumerate(Mechanism.mechanism_order):
        add_mechanism = ("INSERT IGNORE INTO Mechanism "
                            "(mechanism_id, name) "
                            "VALUES (%s, %s);")
        data_mechanism = (idx, mechanism)
        cursor.execute(add_mechanism, data_mechanism)
    cnx.commit()

def write_feature_set(cursor, cnx,feature_set : Features_Set, feature_set_table, **kwargs) -> int:
    do_commit = kwargs.get("do_commit",True)
    column_names, column_values = zip(*list(feature_set.to_series().items()))
    column_str = ", ".join(column_names)
    s = ["%s",]
    add_feature_set = (f"INSERT INTO {feature_set_table} "
                        f"({column_str}) "
                        f"VALUES (%s, {', '.join(s * (len(column_values) - 1))})"
                        f"ON DUPLICATE KEY UPDATE {', '.join([f'{column_name}={column_name}' for column_name in column_names])};")
    try:
        cursor.execute(add_feature_set, column_values)
    except Exception as e:
        print(add_feature_set)
        print(column_values)
        raise e
    feature_set_id = 0#cursor.lastrowid
    if do_commit:
        cnx.commit()
    return feature_set_id

def write_feature_sets(cursor, cnx,feature_sets : List[Features_Set], feature_set_table) -> List[int]:
    feature_set_ids = []
    for feature_set in tqdm(feature_sets, desc=f"Writing {feature_set_table} sets",leave=False):
        feature_set_ids.append(write_feature_set(cursor, cnx, feature_set, feature_set_table, do_commit=False))
    cnx.commit()
    return feature_set_ids