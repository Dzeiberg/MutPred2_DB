from sequence_mapping import SequenceMapping
import pandas as pd

class Mutation:
    def __init__(self, mapping: SequenceMapping, mutation: str, mutpred_score: float):
        self.mapping = mapping
        self.mutation = mutation
        self.mutpred_score = mutpred_score
    def to_series(self):
        return pd.Series({'seq_hash': self.mapping.seq.seq_hash,
                          'mutation': self.mutation,
                          'mutpred_score': self.mutpred_score})
