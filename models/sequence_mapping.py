from sequence import Sequence
import pandas as pd

class SequenceMapping:
    def __init__(self, seq: Sequence, ensp_id: str, ensg_id : str, enst_id : str, gene_symbol : str, canon : bool, MANE_select : bool):
        self.seq = seq
        self.ensp_id = ensp_id
        self.ensg_id = ensg_id
        self.enst_id = enst_id
        self.gene_symbol = gene_symbol
        self.canon = canon
        self.MANE_select = MANE_select
    def to_series(self):
        return pd.Series({'seq_hash' : self.seq.seq_hash,
                   'ensp_id' : self.ensp_id,
                   'ensg_id' : self.ensg_id,
                   'enst_id' : self.enst_id,
                   'gene_symbol' : self.gene_symbol,
                   'canon' : self.canon,
                    'MANE_select' : self.MANE_select})