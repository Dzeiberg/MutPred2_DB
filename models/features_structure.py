from .feature_set import Features_Set

class Features_Structure(Features_Set):
    __features_order__ = ["VSL2B_disorder_wild","VSL2B_disorder_mutant","VSL2B_disorder_loss","VSL2B_disorder_gain","B_factor_wild","B_factor_mutant","B_factor_loss",
    "B_factor_gain","Surface_accessibility_wild","Surface_accessibility_mutant","Surface_accessibility_loss","Surface_accessibility_gain","Helix_wild",
    "Helix_mutant","Helix_loss","Helix_gain","Strand_wild","Strand_mutant","Strand_loss","Strand_gain","Loop_wild","Loop_mutant","Loop_loss","Loop_gain",
    "N-terminal_signal_wild","N-terminal_signal_mutant","N-terminal_signal_loss","N-terminal_signal_gain","Signal_helix_wild","Signal_helix_mutant",
    "Signal_helix_loss","Signal_helix_gain","C-terminal_signal_wild","C-terminal_signal_mutant","C-terminal_signal_loss","C-terminal_signal_gain",
    "Signal_cleavage_wild","Signal_cleavage_mutant","Signal_cleavage_loss","Signal_cleavage_gain","Cytoplasmic_loop_wild","Cytoplasmic_loop_mutant",
    "Cytoplasmic_loop_loss","Cytoplasmic_loop_gain","Transmembrane_region_wild","Transmembrane_region_mutant","Transmembrane_region_loss",
    "Transmembrane_region_gain","Non-cytoplasmic_loop_wild","Non-cytoplasmic_loop_mutant","Non-cytoplasmic_loop_loss","Non-cytoplasmic_loop_gain",
    "Non_transmembrane_wild","Non_transmembrane_mutant","Non_transmembrane_loss","Non_transmembrane_gain","Coiled_coil_wild","Coiled_coil_mutant",
    "Coiled_coil_loss","Coiled_coil_gain","Stability_reliability","Stability_ddG","VSL2B_disorder_exist","B_factor_exist","Surface_accessibility_exist",
    "Helix_exist","Strand_exist","Loop_exist","N-terminal_signal_exist","Signal_helix_exist","C-terminal_signal_exist","Signal_cleavage_exist","Cytoplasmic_loop_exist",
    "Transmembrane_region_exist","Non-cytoplasmic_loop_exist","Non_transmembrane_exist","Coiled_coil_exist","Stability_reliability_exist","Stability_ddG_exist"]