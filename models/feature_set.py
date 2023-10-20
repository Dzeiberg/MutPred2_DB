import pandas as pd
import abc
from .sequence_mapping import SequenceMapping
from .mutation import Mutation

class Features_Set(metaclass=abc.ABCMeta):
    __feature_order__ = []
    
    def __init__(self,mapping, mutation,feature_vec):
        self.mapping = mapping
        self.mutation = mutation
        self.feature_vec = feature_vec

    def to_series(self):
        d = dict(zip(self.__feature_order__, self.feature_vec))
        d["seq_hash"] = self.mapping.seq.seq_hash
        d['mutation'] = self.mutation.mutation
        return pd.Series(d)