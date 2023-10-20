import pandas as pd
from .feature_set import Features_Set

class Features_Sequence(Features_Set):
    __feature_order__ = ["AA_frequency_A","AA_frequency_C","AA_frequency_D","AA_frequency_E","AA_frequency_F","AA_frequency_G","AA_frequency_H","AA_frequency_I","AA_frequency_K","AA_frequency_L",
                    "AA_frequency_M","AA_frequency_N","AA_frequency_P","AA_frequency_Q","AA_frequency_R","AA_frequency_S","AA_frequency_T","AA_frequency_V","AA_frequency_W","AA_frequency_Y",
                    "Entropy","Beta_entropy_1.25","Beta_entropy_1.50","Beta_entropy_1.75","Net_charge","Total_charge","Aromatics","Charge_hydrophobicity","AA_frequency_A","AA_frequency_C",
                    "AA_frequency_D","AA_frequency_E","AA_frequency_F","AA_frequency_G","AA_frequency_H","AA_frequency_I","AA_frequency_K","AA_frequency_L","AA_frequency_M","AA_frequency_N",
                    "AA_frequency_P","AA_frequency_Q","AA_frequency_R","AA_frequency_S","AA_frequency_T","AA_frequency_V","AA_frequency_W","AA_frequency_Y","Entropy","Beta_entropy_1.25",
                    "Beta_entropy_1.50","Beta_entropy_1.75","Net_charge","Total_charge","Aromatics","Charge_hydrophobicity","AA_frequency_A","AA_frequency_C","AA_frequency_D","AA_frequency_E",
                    "AA_frequency_F","AA_frequency_G","AA_frequency_H","AA_frequency_I","AA_frequency_K","AA_frequency_L","AA_frequency_M","AA_frequency_N","AA_frequency_P","AA_frequency_Q",
                    "AA_frequency_R","AA_frequency_S","AA_frequency_T","AA_frequency_V","AA_frequency_W","AA_frequency_Y","Entropy","Beta_entropy_1.25","Beta_entropy_1.50","Beta_entropy_1.75",
                    "Net_charge","Total_charge","Aromatics","Charge_hydrophobicity","AA_frequency_A","AA_frequency_C","AA_frequency_D","AA_frequency_E","AA_frequency_F","AA_frequency_G",
                    "AA_frequency_H","AA_frequency_I","AA_frequency_K","AA_frequency_L","AA_frequency_M","AA_frequency_N","AA_frequency_P","AA_frequency_Q","AA_frequency_R","AA_frequency_S",
                    "AA_frequency_T","AA_frequency_V","AA_frequency_W","AA_frequency_Y","Entropy","Beta_entropy_1.25","Beta_entropy_1.50","Beta_entropy_1.75","Net_charge","Total_charge",
                    "Aromatics","Charge_hydrophobicity","Hydrophobicity_mean","Hydrophobicity_std","Hydrophobicity_max","Hydrophobic_moments_100_mean","Hydrophobic_moments_100_std",
                    "Hydrophobic_moments_100_max","Hydrophobic_moments_160_mean","Hydrophobic_moments_160_std","Hydrophobic_moments_160_max","Hydrophobic_moments_120_mean",
                    "Hydrophobic_moments_120_std","Hydrophobic_moments_120_max","Vihinen_flexibility_mean","Vihinen_flexibility_std","Vihinen_flexibility_max","Amino_acid_volume_mean",
                    "Amino_acid_volume_std","Amino_acid_volume_max","Hydrophobicity_mean","Hydrophobicity_std","Hydrophobicity_max","Hydrophobic_moments_100_mean","Hydrophobic_moments_100_std",
                    "Hydrophobic_moments_100_max","Hydrophobic_moments_160_mean","Hydrophobic_moments_160_std","Hydrophobic_moments_160_max","Hydrophobic_moments_120_mean",
                    "Hydrophobic_moments_120_std","Hydrophobic_moments_120_max","Vihinen_flexibility_mean","Vihinen_flexibility_std","Vihinen_flexibility_max","Amino_acid_volume_mean",
                    "Amino_acid_volume_std","Amino_acid_volume_max","Hydrophobicity_mean","Hydrophobicity_std","Hydrophobicity_max","Hydrophobic_moments_100_mean","Hydrophobic_moments_100_std",
                    "Hydrophobic_moments_100_max","Hydrophobic_moments_160_mean","Hydrophobic_moments_160_std","Hydrophobic_moments_160_max","Hydrophobic_moments_120_mean","Hydrophobic_moments_120_std",
                    "Hydrophobic_moments_120_max","Vihinen_flexibility_mean","Vihinen_flexibility_std","Vihinen_flexibility_max","Amino_acid_volume_mean","Amino_acid_volume_std",
                    "Amino_acid_volume_max","Hydrophobicity_mean","Hydrophobicity_std","Hydrophobicity_max","Hydrophobic_moments_100_mean","Hydrophobic_moments_100_std",
                    "Hydrophobic_moments_100_max","Hydrophobic_moments_160_mean","Hydrophobic_moments_160_std","Hydrophobic_moments_160_max","Hydrophobic_moments_120_mean",
                    "Hydrophobic_moments_120_std","Hydrophobic_moments_120_max","Vihinen_flexibility_mean","Vihinen_flexibility_std","Vihinen_flexibility_max","Amino_acid_volume_mean","Amino_acid_volume_std","Amino_acid_volume_max"]
