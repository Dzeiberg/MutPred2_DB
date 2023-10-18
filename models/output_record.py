from .mechanism import Mechanism
from .mutation import Mutation
from typing import List, TypeVar
import pandas as pd

from numpy import ndarray

class MutPred2Output:
    T = TypeVar('T', bound='MutPred2Output')
    M = TypeVar('M', bound='Mechanism')

    __altered_set__ = set(('N_terminal_signal', 'Signal_helix',
                                        'C_terminal_signal', 'Signal_cleavage',
                                        'Intracellular_loop', 'Transmembrane_region',
                                        'Extracellular_loop', 'Non_transmembrane',
                                        'Coiled_coil', 'Calmodulin_binding',
                                        'DNA_binding', 'RNA_binding',
                                        'PPI_residue', 'PPI_hotspot',
                                        'MoRF', 'Stability'))
    
    __no_region_set__ = __altered_set__.union(set(('Helix', 'Strand',
                                                   'Loop', 'VSL2B_disorder',
                                                   'B_factor', 'Surface_accessibility')))
    
    __mechanism_order__ = ["VSL2B_disorder","B_factor","Surface_accessibility","Helix",
                          "Strand","Loop","N_terminal_signal","Signal_helix","C_terminal_signal",
                          "Signal_cleavage","Intracellular_loop","Transmembrane_region","Extracellular_loop","Non_transmembrane","Coiled_coil",
                          "Catalytic_site","Calmodulin_binding","DNA_binding","RNA_binding","PPI_residue","PPI_hotspot","MoRF",
                          "Allosteric_site","Cadmium_binding","Calcium_binding","Cobalt_binding","Copper_binding","Iron_binding",
                          "Magnesium_binding","Manganese_binding","Nickel_binding","Potassium_binding","Sodium_binding","Zinc_binding",
                          "Acetylation","ADP_ribosylation","Amidation","C_linked_glycosylation","Carboxylation","Disulfide_linkage",
                          "Farnesylation","Geranylgeranylation","GPI_anchor_amidation","Hydroxylation","Methylation","Myristoylation",
                          "N_linked_glycosylation","N_terminal_acetylation","O_linked_glycosylation","Palmitoylation","Phosphorylation",
                          "Proteolytic_cleavage","Pyrrolidone_carboxylic_acid","Sulfation","SUMOylation","Ubiquitylation","Motifs","Stability"]
    
    def __init__(self,mutation : Mutation, mechanisms : List[M]) -> None:
        self.mutation = mutation
        self.mechanisms = mechanisms

    def to_series(self,flatfile_info=True) -> pd.Series:
        if flatfile_info:
            data = {"RefSeq_nuc" : self.mutation.mapping.seq.refseq_nuc,
                     "RefSeq_prot" : self.mutation.mapping.seq.refseq_prot,
                     "Ensembl_gene" : self.mutation.mapping.seq.ensg_id,
                     "Ensembl_nuc" : self.mutation.mapping.seq.enst_id,
                    "Ensembl_prot" : self.mutation.mapping.seq.ensp_id,
                    "Substitution" : self.mutation.mutation,
                    "MutPred2_Score" : self.mutation.mutpred_score}
        else:
            data = {'seq_hash' : self.mutation.mapping.seq.seq_hash,
                    'mutation' : self.mutation.mutation}
        for mechanism in self.mechanisms:
            if mechanism.name not in MutPred2Output.__no_region_set__:
                data[mechanism.name + "_position"] = mechanism.position
            data[mechanism.name + "_posterior"] = mechanism.posterior
            data[mechanism.name + "_pvalue"] = mechanism.p_value
            data[mechanism.name + "_effect"] = mechanism.effect_type
        return pd.Series(data)

    @staticmethod
    def read_mechanisms(positions_pu : ndarray, pvals_pu : ndarray, scores_pu : ndarray, types_pu : ndarray) -> List[M]:
        mechanisms = []
        for mechanism, position_encoding, p_val, score, effect_num in zip(MutPred2Output.__mechanism_order__,
                                                                          positions_pu,
                                                                          pvals_pu,
                                                                          scores_pu,
                                                                          types_pu):
            if mechanism in MutPred2Output.__altered_set__:
                effect_type = 'altered'
            elif effect_num == 1:
                effect_type = "loss"
            else:
                effect_type = "gain"
            if mechanism in MutPred2Output.__no_region_set__:
                position = None
            else:
                position = position_encoding
            try:
                mechanisms.append(Mechanism(mechanism, effect_type, p_val, score, position))
            except TypeError as e:
                print(f"Failed to initialize mechanism {mechanism}")
                raise e
        return mechanisms

