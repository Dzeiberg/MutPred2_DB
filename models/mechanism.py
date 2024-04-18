import typing
import pandas as pd
class Mechanism:
    altered_set = set(('N_terminal_signal', 'Signal_helix',
                                        'C_terminal_signal', 'Signal_cleavage',
                                        'Intracellular_loop', 'Transmembrane_region',
                                        'Extracellular_loop', 'Non_transmembrane',
                                        'Coiled_coil', 'Calmodulin_binding',
                                        'DNA_binding', 'RNA_binding',
                                        'PPI_residue', 'PPI_hotspot',
                                        'MoRF', 'Stability'))
    
    no_region_set = altered_set.union(set(('Helix', 'Strand',
                                                'Loop', 'VSL2B_disorder',
                                                'B_factor', 'Surface_accessibility')))
    
    mechanism_order = ["VSL2B_disorder","B_factor","Surface_accessibility","Helix",
                        "Strand","Loop","N_terminal_signal","Signal_helix","C_terminal_signal",
                        "Signal_cleavage","Intracellular_loop","Transmembrane_region","Extracellular_loop","Non_transmembrane","Coiled_coil",
                        "Catalytic_site","Calmodulin_binding","DNA_binding","RNA_binding","PPI_residue","PPI_hotspot","MoRF",
                        "Allosteric_site","Cadmium_binding","Calcium_binding","Cobalt_binding","Copper_binding","Iron_binding",
                        "Magnesium_binding","Manganese_binding","Nickel_binding","Potassium_binding","Sodium_binding","Zinc_binding",
                        "Acetylation","ADP_ribosylation","Amidation","C_linked_glycosylation","Carboxylation","Disulfide_linkage",
                        "Farnesylation","Geranylgeranylation","GPI_anchor_amidation","Hydroxylation","Methylation","Myristoylation",
                        "N_linked_glycosylation","N_terminal_acetylation","O_linked_glycosylation","Palmitoylation","Phosphorylation",
                        "Proteolytic_cleavage","Pyrrolidone_carboxylic_acid","Sulfation","SUMOylation","Ubiquitylation","Motifs","Stability"]

    motif_index = [i for i,mechanism in enumerate(mechanism_order) if mechanism == "Motifs"][0]
    def __init__(self, mechanism_id : int,
                        variant_id : int,
                        altered_position : int,
                        pvalue : float,
                        score : float,
                        mechanism_type : int,
                        description : typing.Optional[str]=None,
                        ) -> None:
        self.mechanism_id = mechanism_id

        self.variant_id = variant_id
        self.process_mechanism_type(mechanism_type)
        self.score = score
        self.pvalue = pvalue
        self.process_alterred_position(altered_position)
        self.description = description

    def process_mechanism_type(self, mechanism_type : int) -> None:
        if Mechanism.mechanism_order[self.mechanism_id] in Mechanism.altered_set:
            self.mechanism_type = 'altered'
        elif mechanism_type == 1:
            self.mechanism_type = "loss"
        else:
            self.mechanism_type = "gain"

    def process_alterred_position(self, altered_position : int) -> None:
        if Mechanism.mechanism_order[self.mechanism_id] in Mechanism.no_region_set:
            self.position = None
        else:
            self.position = int(altered_position)

    def to_series(self):
        d = dict(variant_id = self.variant_id,
                mechanism_id=self.mechanism_id,
                mechanism_type=self.mechanism_type,
                altered_position=self.position,
                description=self.description,
                score=self.score,
                pvalue=self.pvalue)
        return pd.Series(d)
