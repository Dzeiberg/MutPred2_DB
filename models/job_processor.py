from sequence import Sequence
from output_record import MutPred2Output
from pathlib import Path
from typing import List,Dict
from Bio.SeqIO import parse

class Processor:
    def __init__(self,sequence_objects : Dict[str,Sequence]):
        self.sequence_objects = sequence_objects

    def process_job(self,job_dir : Path) -> List[MutPred2Output.T]:
        sequence = str(next(iter(parse(job_dir/'input.faa', 'fasta').records)).seq)
        seq_hash = Sequence.get_sequence_hash(sequence)
        seq = self.sequence_objects[seq_hash]
