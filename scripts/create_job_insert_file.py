import Bio.SeqIO
import os
import hashlib

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

def create_insert_file(job_path, insert_filepath):
    """
    process the job at job_path and add the insert statements to insert_filepath
    """
    try:
        (ensembl_prot_id, seq) = parse_input(os.path.join(job_path, "input.faa"))
    except (FileNotFoundError, NotADirectoryError):
        print(f"cannot find {job_path}/input.faa")
        return
    seq_hash = hashlib.md5(seq.encode()).hexdigest()
    try:
        with open(os.path.join(job_path, "output.txt")) as f:
            lines = f.readlines()[1:]
    except (FileNotFoundError, NotADirectoryError):
        print(f"cannot find {job_path}/output.txt")
        return
    variant_results = list(map(parse_line,lines))
    with open(insert_filepath, "w") as f:
        f.write(f"INSERT IGNORE INTO proteins (seq_hash, sequence) VALUES ('{seq_hash}', '{seq}');\n")
        for variant in variant_results:
            f.write(f"INSERT IGNORE INTO variant (seq_hash, reference_aa, position, alternate_aa, score, option_id) VALUES ('{seq_hash}', '{variant[1]}', {variant[2]}, '{variant[3]}', {variant[4]}, 1);\n")

def process_many_jobs(jobs_dir,job_subset=None):
    "process all jobs in a directory"
    jobs = os.listdir(jobs_dir)
    for job in jobs:
        if job_subset is not None and job not in job_subset: continue
        job_pth = os.path.join(jobs_dir, job)
        create_insert_file(job_pth, os.path.join(job_pth,"insert_job.sql"))

if __name__ == "__main__":
    import sys
    subset = None
    if len(sys.argv) < 2:
        subset = set(sys.argv[2:])
    process_many_jobs(sys.argv[1], subset)