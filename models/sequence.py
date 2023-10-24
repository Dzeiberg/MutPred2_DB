import hashlib
import pandas as pd
from pathlib import Path
from Bio.SeqIO import parse
from sqlalchemy.engine import Engine

class Sequence:
    def __init__(self,seq):
        self.seq = seq
        h = hashlib.new('md5')
        h.update(self.seq.encode('utf-8'))
        self.seq_hash = Sequence.get_sequence_hash(seq)

    def to_series(self):
        return pd.Series({'seq_hash': self.seq_hash,
                          'sequence':self.seq})

    @staticmethod
    def get_sequence_hash(sequence : str):
        return hashlib.md5(sequence.encode(encoding='utf-8')).hexdigest()
