# from .sequence_mapping import SequenceMapping
import pandas as pd
from Bio.PDB.Polypeptide import protein_letters_3to1
import hashlib
RESIDUES = set(protein_letters_3to1 .values())

class Variant:
    def __init__(self, seq_hash : str, substitution: str, mutpred_score: float, option_id: int):
        self.seq_hash = seq_hash
        self.option_id = option_id
        self.parse_substitution(substitution)
        hashstr = f"{self.seq_hash}{self.reference_aa}{self.position}{self.alternate_aa}{self.option_id}"
        self.variant_id = hashlib.md5(hashstr.encode(encoding='utf-8')).hexdigest()
        self.mutpred_score = mutpred_score

    def parse_substitution(self, substitution: str):
        self.reference_aa = substitution[0]
        assert self.reference_aa in RESIDUES, f"Invalid residue {self.reference_aa} for substitution {substitution}"
        self.position = int(substitution[1:-1])
        self.alternate_aa = substitution[-1]
        assert self.alternate_aa in RESIDUES, f"Invalid residue {self.alternate_aa} for substitution {substitution}"

    def to_series(self):
        return pd.Series({'seq_hash': self.seq_hash,
                          'reference_aa': self.reference_aa,
                          'position': self.position,
                          'alternate_aa': self.alternate_aa,
                          'score': self.mutpred_score,
                          'option_id': self.option_id})
