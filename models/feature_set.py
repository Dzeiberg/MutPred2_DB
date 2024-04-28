import pandas as pd
import abc

class Features_Set(metaclass=abc.ABCMeta):
    __features_order__ = []
    
    def __init__(self,variant_id,runoption_id,feature_vec):
        self.variant_id = variant_id
        self.runoption_id = runoption_id
        self.feature_vec = feature_vec

    def to_series(self):
        d = dict(zip(self.__features_order__, self.feature_vec))
        d["variant_id"] = self.variant_id
        d["runoption_id"] = self.runoption_id
        return pd.Series(d)