import pandas as pd
import urllib
from Bio import SeqIO
import gzip
import hashlib
from fire import Fire

def fetch_mane_summary():
    # Fetch the MANE summary from the MANE website
    # and return it as a list of dictionaries
    summary = pd.read_csv('https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.3.summary.txt.gz', compression='gzip',sep='\t')
    column_mapping = {
        "Ensembl_Gene" : 'ensembl_gene_id',
        'Ensembl_nuc' : 'ensembl_nuc_id',
        'Ensembl_prot' : 'ensembl_prot_id',
        'symbol' : 'gene_symbol',
        'MANE_status' : 'MANE_select'
    }
    df = summary.loc[:, list(column_mapping.keys())]
    df = df.rename(columns=column_mapping)
    df = df.assign(MANE_select = df.MANE_select == 'MANE Select').set_index('ensembl_prot_id')
    return df

def fetch_fasta():
    # Fetch the MANE FASTA file from the MANE website
    # and return it as a list of dictionaries
    urllib.request.urlretrieve('https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.3.ensembl_protein.faa.gz', 'MANE.GRCh38.v1.3.ensembl_protein.faa.gz')
    protein_records = []
    ensembl_prot_id_to_hash = {}
    with gzip.open('MANE.GRCh38.v1.3.ensembl_protein.faa.gz', "rt") as handle:
        for r in SeqIO.parse(handle, "fasta"):
            seq = str(r.seq)
            seq_hash = hashlib.md5(seq.encode()).hexdigest()
            protein_records.append((seq, seq_hash))
            ensembl_prot_id_to_hash[r.id] = seq_hash
    return protein_records, ensembl_prot_id_to_hash

def create_insert_files(protein_output_filepath, sequence_mapping_output_filepath):
    # Create sql insert file to populate the proteins table
    protein_records,ensembl_prot_id_to_hash = fetch_fasta()
    with open(protein_output_filepath, 'w') as f:
        for record in protein_records:
            f.write(f"INSERT IGNORE INTO proteins (sequence, seq_hash) VALUES ('{record[0]}', '{record[1]}');\n")
    summary = fetch_mane_summary()
    summary = summary.assign(seq_hash = summary.index.map(ensembl_prot_id_to_hash))
    with open(sequence_mapping_output_filepath, 'w') as f:
        for index, row in summary.iterrows():
            f.write(f"INSERT INTO sequence_mapping (seq_hash, ensembl_gene_id, ensembl_nuc_id, ensembl_prot_id, gene_symbol, MANE_select) VALUES ('{row.seq_hash}', '{row['ensembl_gene_id']}', '{row['ensembl_nuc_id']}', '{index}', '{row['gene_symbol']}', {row['MANE_select']});\n")

if __name__ == '__main__':
    Fire(create_insert_files)