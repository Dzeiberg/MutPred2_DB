#!/usr/bin/env python
# coding: utf-8

import mysql.connector
import Bio.SeqIO
from tqdm.autonotebook import tqdm
import os
from itertools import chain
from tqdm.autonotebook import tqdm,trange
import pandas as pd

######## GLOBALS ########
phase2_root = "/Users/danielzeiberg/Documents/research/IGVF/phase2/"
#########################

""" mysql connection information """

mydb = mysql.connector.connect(
  host="localhost",
  user="root",
  password="",
  database="igvf"
)

mycursor = mydb.cursor()

""" create database schema """
mycursor.execute("CREATE database if not exists igvf")
mycursor.execute("CREATE TABLE if not exists proteins (ensembl_protein_id CHAR(32) PRIMARY KEY NOT NULL, sequence TEXT NOT NULL)")
mycursor.execute("CREATE TABLE if not exists run_options (option_id int PRIMARY KEY AUTO_INCREMENT, compute_homology_profile BIT(1), use_predicted_conservation_scores BIT(1), skip_psi_blast BIT(1), p_value_threshold DECIMAL(2,1), UNIQUE (compute_homology_profile, use_predicted_conservation_scores, skip_psi_blast, p_value_threshold))")
mycursor.execute("CREATE TABLE if not exists variant (ensembl_protein_id CHAR(32), reference_aa CHAR(1) NOT NULL, position int NOT NULL, alternate_aa CHAR(1) NOT NULL, score DECIMAL(4,3) NOT NULL, option_id int, PRIMARY KEY(ensembl_protein_id, position, reference_aa, alternate_aa), FOREIGN KEY (ensembl_protein_id) REFERENCES proteins(ensembl_protein_id), FOREIGN KEY (option_id) REFERENCES run_options(option_id))")

# insert the run options used to generate phase 2 results
phase2_options = (1,1,0,1)
mycursor.execute("INSERT IGNORE INTO run_options (compute_homology_profile, use_predicted_conservation_scores, skip_psi_blast, p_value_threshold) VALUES (%s, %s, %s, %s)", phase2_options)

"""  Functions to parse a result line from a mutpred2 output file """
def parse_line(line):
    """
    Parse a line from a mutpred2 output file

    Args:
    line (str): a line from a mutpred2 output file : assumed (minimum) format: <ensembl_protein_id>,<reference_aa><position><alternate_aa>,<score>
    """
    parts = line.split(",")
    ensemble_protein_id = ".".join(parts[0].split("_")[:2])
    reference_aa = parts[1][0]
    position = int(parts[1][1:-1])
    alternate_aa = parts[1][-1]
    score = float(parts[2])
    return (ensemble_protein_id, reference_aa, position, alternate_aa, score)


def parse_input(fasta_file):
    "Return the ensembl protein id and sequence from a fasta file assuemd to have a single record"
    records = list(Bio.SeqIO.parse(fasta_file, "fasta").records)
    assert len(records) == 1
    return (".".join(records[0].id.split("_")[:2]), str(records[0].seq))

"""Read into memory the input and output files from the phase 2 mutpred2 runs"""
inputs = []
errors = []
for job in tqdm(os.listdir(phase2_root)):
    try:
        protein_input = parse_input(os.path.join(phase2_root, job, "input.faa"))
    except (FileNotFoundError, NotADirectoryError):
        errors.append(f"cannot find {phase2_root}/{job}/input.faa")
        continue
    try:
        with open(os.path.join(phase2_root, job, "output.txt")) as f:
            lines = f.readlines()[1:]
    except (FileNotFoundError, NotADirectoryError):
        errors.append(f"cannot find {phase2_root}/{job}/output.txt")
        continue
    variant_results = list(map(parse_line,lines))
    inputs.append((protein_input, variant_results))

# set the max_allowed_packet to 500MB
mycursor.execute('SET GLOBAL max_allowed_packet=500*1024*1024')
# insert the protein sequences into the database
mycursor.executemany("INSERT IGNORE INTO proteins (ensembl_protein_id, sequence) VALUES (%s, %s)",list(set([x[0] for x in inputs])))
mydb.commit()

# add the run_option id to the variant records
processed_inputs = list(map(lambda tup: (*tup, 1), chain(*[x[1] for x in inputs])))
# insert the variant records into the database
batch_size = 25000
for start in trange(0,len(processed_inputs), batch_size):
    mycursor.executemany("INSERT IGNORE INTO variant (ensembl_protein_id, reference_aa, position, alternate_aa, score, option_id) VALUES (%s, %s, %s, %s, %s, %s)", processed_inputs[start:start+batch_size])
    mydb.commit()
