import pandas as pd
import abc

class Features_Set(metaclass=abc.ABCMeta):
    __feature_order__ = []
    
    def __init__(self,variant_id,feature_vec):
        self.variant_id = variant_id
        self.feature_vec = feature_vec

    def to_series(self):
        d = dict(zip(self.__feature_order__, self.feature_vec))
        d["variant_id"] = self.variant_id
        return pd.Series(d)