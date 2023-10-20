from .feature_set import Features_Set

class Features_Homology(Features_Set):
    __features_order__ = ["Homology_count_human_50","Homology_count_human_55","Homology_count_human_60","Homology_count_human_65",
                          "Homology_count_human_70","Homology_count_human_75","Homology_count_human_80","Homology_count_human_85",
                          "Homology_count_human_90","Homology_count_human_95","Homology_count_mouse_50","Homology_count_mouse_55",
                          "Homology_count_mouse_60","Homology_count_mouse_65","Homology_count_mouse_70","Homology_count_mouse_75",
                          "Homology_count_mouse_80","Homology_count_mouse_85","Homology_count_mouse_90","Homology_count_mouse_95"]