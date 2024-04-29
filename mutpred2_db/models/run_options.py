import pandas as pd

class RunOption:
    def __init__(self, compute_homology_profile : bool,
                        use_predicted_conservation_scores : bool,
                        skip_psi_blast : bool,
                        p_value_threshold : float) : 
        """
        Arguments used in computing MutPred2 results

        Parameters
        ----------
        compute_homology_profile : bool
            Whether to compute the homology profile of the sequence

        use_predicted_conservation_scores : bool
            Whether to use predicted conservation scores

        skip_psi_blast : bool
            Whether to skip PSI-BLAST

        p_value_threshold : float
            The threshold for the p-value

        """
        self.compute_homology_profile = compute_homology_profile
        self.use_predicted_conservation_scores = use_predicted_conservation_scores
        self.skip_psi_blast = skip_psi_blast
        self.p_value_threshold = p_value_threshold

    def to_series(self):
        return pd.Series({'compute_homology_profile': self.compute_homology_profile,
                          'use_predicted_conservation_scores': self.use_predicted_conservation_scores,
                          'skip_psi_blast': self.skip_psi_blast,
                          'p_value_threshold': self.p_value_threshold})