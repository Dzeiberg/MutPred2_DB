from .sequence import Sequence
from .variant import Variant
from .features_conservation import Features_Conservation
from .features_function import Features_Function
from .features_homology import Features_Homology
from .features_pssm import Features_PSSM
from .features_sequence import Features_Sequence
from .features_structure import Features_Structure
from .features_substitution import Features_Substitution
from .mechanism import Mechanism

from scipy.io import loadmat
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List,Dict
import re
import os
from typing import Iterable

class Processor:

    def process(self, job_dir : Path):
        sequence = self.make_sequence(job_dir)
        variants = self.make_variants(job_dir, sequence)
        mechanisms = self.make_mechanisms(job_dir, variants)
        (features_sequence, features_substitution,
                features_pssm, features_conservation,
                features_homology, features_structure,
                features_function) = self.make_features(job_dir, variants)
        return dict(sequence=sequence,
                    variants=variants,
                    features_sequence=features_sequence,
                    features_substitution=features_substitution,
                    features_pssm=features_pssm,
                    features_conservation=features_conservation,
                    features_homology=features_homology,
                    features_structure=features_structure,
                    features_function=features_function,
                    mechanisms=mechanisms)

    def make_sequence(self, job_dir : Path) -> str:
        # seq = self.read_mat_files(job_dir, pattern='.*.txt.sequences.mat',key_value='sequences').item().item()
        seq = loadmat(f"{job_dir}/output.txt.sequences.mat")['sequences'].item().item()
        assert isinstance(seq,str), "Sequence must be a string, not {}".format(type(seq))
        return Sequence(seq)

    def read_substitutions(self, job_dir : Path) -> List[str]:
        # subs = self.read_mat_files(job_dir, pattern='.*.txt.substitutions.mat',key_value='substitutions')
        subs = loadmat(f"{job_dir}/output.txt.substitutions.mat")['substitutions']
        substitutions = list(map(lambda i : i.item(),subs.item().ravel()))
        return substitutions

    def read_mutpred2_scores(self, job_dir : Path) -> Iterable[float]:
        try:
            scores = self.read_mat_files(job_dir, pattern='.*.txt.MutPred2Score_\d+.mat',key_value='S')
        except AttributeError:
            df = pd.read_csv(f"{job_dir}/output.txt")
            scores = df.loc[:,'MutPred2 score'].values.astype(float)
        return scores

    def make_variants(self, job_dir : Path, sequence : Sequence) -> Iterable[Variant]:
        substitutions = self.read_substitutions(job_dir)
        scores = self.read_mutpred2_scores(job_dir)
        variants = [Variant(sequence.seq_hash,sub,score) for sub,score in zip(substitutions,scores)]
        self.validate_variants(variants,sequence)
        return variants

    def validate_variants(self, variants : List[Variant], sequence : Sequence) -> None:
        for variant in variants:
            if sequence.seq[variant.position-1] != variant.reference_aa:
                raise ValueError(f"Reference amino acid {variant.reference_aa} at position {variant.position} does not match sequence {sequence.seq}")

    def read_mat_files(self, job_dir : Path, pattern : str, key_value : str):
        directory = os.listdir(job_dir)
        re_pattern = re.compile(pattern)
        files = sorted([f for f in directory if re_pattern.match(f)])
        data = pd.DataFrame({
            "Files": files
        })
        data["file_num"] = data["Files"].str.split(".").str[-2].str.split("_").str[-1].astype(int)
        files = data.sort_values(by="file_num")["Files"]
        if key_value == "motif_info" or key_value=="models":
            data = np.concatenate([loadmat(f"{job_dir}/{f}")[key_value].ravel() for f in files])

        else:
            data = np.concatenate([loadmat(f"{job_dir}/{f}")[key_value] for f in files])

        return data

    def read_mechanism_info(self, job_dir : Path) -> Dict[str,np.ndarray]:
        positions_pu = self.read_mat_files(job_dir=job_dir, pattern='.*.txt.positions_pu_\d+.mat', key_value='positions_pu')
        positions_pu = np.concatenate((positions_pu,
                                    np.ones((positions_pu.shape[0], 1)) * -1),
                                    axis=1)
        pvals_pu        = self.read_mat_files(job_dir=job_dir, pattern='.*.txt.prop_pvals_pu_\d+.mat', key_value='prop_pvals_pu')
        scores_pu       = self.read_mat_files(job_dir=job_dir, pattern='.*.txt.prop_scores_pu_\d+.mat', key_value='prop_scores_pu')
        prop_types_pu   = self.read_mat_files(job_dir=job_dir, pattern='.*.txt.prop_types_pu_\d+.mat', key_value='prop_types_pu')
        motif_info      = self.read_mat_files(job_dir=job_dir, pattern='.*.txt.motif_info_\d+.mat', key_value='motif_info')
        return dict(positions_pu=positions_pu,
                    pvals_pu=pvals_pu,
                    scores_pu=scores_pu,
                    prop_types_pu=prop_types_pu,
                    motif_info=motif_info)

    def make_features(self, job_dir : Path, variants : List[Variant]):
        features = self.read_mat_files(job_dir=job_dir, pattern='.*.txt.feats_\d+.mat', key_value='feats')
        features_conservation = []
        features_function = []
        features_homology = []
        features_pssm = []
        features_sequence = []
        features_structure = []
        features_substitution = []
        for i,(variant,feats) in enumerate(zip(variants, features)):
            features_sequence.append(Features_Sequence(variant.seq_hash, feats[:184]))
            features_substitution.append(Features_Substitution(variant.seq_hash, feats[184:630]))
            features_pssm.append(Features_PSSM(variant.seq_hash, feats[630:799]))
            features_conservation.append(Features_Conservation(variant.seq_hash, feats[799:1036]))
            features_homology.append(Features_Homology(variant.seq_hash, feats[1036:1056]))
            features_structure.append(Features_Structure(variant.seq_hash, feats[1056:1135]))
            features_function.append(Features_Function(variant.seq_hash, feats[1135:]))
        return (features_sequence, features_substitution,
                features_pssm, features_conservation,
                features_homology, features_structure,
                features_function)

    def make_mechanisms(self,job_dir : Path, variants : List[Variant]) -> List[Mechanism]:

        mechanism_info = self.read_mechanism_info(job_dir)
        mechanisms = []
        for i,variant in enumerate(variants):
            mechanism_info_i = {k:v[i] for k,v in mechanism_info.items()}
            motif = mechanism_info_i.get('motif_info')
            for mechanism_idx, (position, pvalue, score, prop_type) in enumerate(zip(*[mechanism_info_i.get(k) for k in ['positions_pu',
                                                                                                                'pvals_pu',
                                                                                                                'scores_pu',
                                                                                                                'prop_types_pu']])):
                mechanisms.append(Mechanism(mechanism_idx,
                                            variant.seq_hash,
                                            position,
                                            pvalue,
                                            score,
                                            prop_type,
                                            motif if mechanism_idx == Mechanism.motif_index else None))
        return mechanisms

    