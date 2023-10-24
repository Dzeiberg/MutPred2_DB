import fire
import pandas as pd
import gzip
import urllib.request
from Bio.SeqIO import parse
import os
from models.sequence import Sequence
from models.sequence_mapping import SequenceMapping
from models.sql_connection import SQL_Connection

def init_sequence_and_mapping_tables(db_config_file="sql_configs.yaml", db_config_name="Local"):
    con = SQL_Connection(db_config_name, db_config_file)
    mane_summary = pd.read_csv("https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.2/MANE.GRCh38.v1.2.summary.txt.gz",compression='gzip',delimiter='\t').set_index("RefSeq_prot")
    mane_sequences = {}
    response = urllib.request.urlretrieve("https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.2/MANE.GRCh38.v1.2.refseq_protein.faa.gz",'MANE.GRCh38.v1.2.refseq_protein.faa.gz')
    with gzip.open('MANE.GRCh38.v1.2.refseq_protein.faa.gz','rt') as handle:
        for record in parse(handle, 'fasta'):
            mane_sequences[record.id] = str(record.seq)
    os.remove("MANE.GRCh38.v1.2.refseq_protein.faa.gz")
    mane_summary = mane_summary.assign(sequence = mane_sequences)
    mappings = []
    for _,r in mane_summary.iterrows():
        seq = Sequence(r.sequence)
        mapping = SequenceMapping(seq, r.Ensembl_prot, r.Ensembl_Gene, r.Ensembl_nuc, r.symbol, True, True)
        mappings.append(mapping)
    mapping_df = pd.DataFrame.from_records([m.to_series() for m in mappings])
    sequence_df = pd.DataFrame.from_records([m.seq.to_series() for m in mappings])
    existing_mappings = pd.read_sql('SELECT * FROM sequence_mapping',con=con.get_engine())
    existing_mappings.canon = existing_mappings.canon.astype(bool)
    existing_mappings.MANE_select = existing_mappings.MANE_select.astype(bool)
    existing_sequences = pd.read_sql('SELECT * FROM sequences',con=con.get_engine())
    existing_mappings = set(map(tuple,existing_mappings.values))
    existing_sequences = set(map(tuple,existing_sequences.values))
    new_mappings = set(map(tuple,mapping_df.values))
    new_sequences = set(map(tuple,sequence_df.values))
    mappings_to_add = pd.DataFrame(list(new_mappings - existing_mappings),columns=mapping_df.columns)
    sequences_to_add = pd.DataFrame(list(new_sequences - existing_sequences),columns=sequence_df.columns)
    print(f"preparing to add {mappings_to_add.shape[0]} mappings")
    print(f"preparing to add {sequences_to_add.shape[0]} sequences")
    mappings_to_add.to_sql("sequence_mapping", con=con.get_engine(),if_exists='append',index=False)
    sequences_to_add.to_sql("sequences",con=con.get_engine(),if_exists='append',index=False)
    


if __name__ == "__main__":
    fire.Fire(init_sequence_and_mapping_tables)