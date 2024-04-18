import pandas as pd

class RunOption:
    def __init__(self, compute_homology_profile : bool,
                        use_predicted_conservation_scores : bool,
                        skip_psi_blast : bool,
                        p_value_threshold : float) : 
        self.compute_homology_profile = compute_homology_profile
        self.use_predicted_conservation_scores = use_predicted_conservation_scores
        self.skip_psi_blast = skip_psi_blast
        self.p_value_threshold = p_value_threshold

    def to_series(self):
        return pd.Series({'compute_homology_profile': self.compute_homology_profile,
                          'use_predicted_conservation_scores': self.use_predicted_conservation_scores,
                          'skip_psi_blast': self.skip_psi_blast,
                          'p_value_threshold': self.p_value_threshold})