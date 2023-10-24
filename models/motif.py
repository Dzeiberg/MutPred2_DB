from .sequence_mapping import SequenceMapping
from .mutation import Mutation
from typing import List,TypeVar
import pandas as pd

class Motif:
    M = TypeVar('M', bound='Motif')
    def __init__(self,seq_mapping : SequenceMapping,mutation : Mutation,description : str, posterior : float, pvalue : float):
        self.seq_mapping = seq_mapping
        self.mutation = mutation
        self.description = description
        self.posterior = posterior
        self.pvalue = pvalue

    @classmethod
    def motifs_from_string(Motif, seq_mapping : SequenceMapping, mutation : Mutation, motif_str : str, posterior : float, pvalue : float) -> List[M]:
        motif_strings = motif_str.split(";")
        motifs = []
        for ms in motif_strings:
            if ms == "None": continue
            motifs.append(Motif(seq_mapping, mutation, ms, posterior, pvalue))
        return motifs
    
    def to_series(self):
        return pd.Series({"seq_hash" : self.seq_mapping.seq.seq_hash,"mutation" : self.mutation.mutation,
          "motif" : self.description, "posterior" : self.posterior, "pvalue" : self.pvalue})
