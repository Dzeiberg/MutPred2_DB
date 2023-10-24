from .sequence_mapping import SequenceMapping
from .sequence import Sequence
from .mutation import Mutation
from .output_record import MutPred2Output
from .sql_connection import SQL_Connection
from .motif import Motif
from .features_conservation import Features_Conservation
from .features_function import Features_Function
from .features_homology import Features_Homology
from .features_pssm import Features_PSSM
from .features_sequence import Features_Sequence
from .features_structure import Features_Structure
from .features_substitution import Features_Substitution

from scipy.io import loadmat
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List,Dict
from Bio.SeqIO import parse


class Processor:
    def __init__(self, con : SQL_Connection):
        self.con = con
        self.init_mapping_objects()
        
    def init_mapping_objects(self):
        self.mapping_objects = {}
        seq_df = pd.read_sql("SELECT * FROM sequence_mapping INNER JOIN sequences ON sequence_mapping.seq_hash = sequences.seq_hash",con=self.con.get_engine())
        for _,r in seq_df.iterrows():
            sequence = Sequence(r.sequence)
            mapping = SequenceMapping(sequence, r.ensp_id, r.ensg_id, r.enst_id, r.gene_symbol,r.canon,r.MANE_select)
            self.mapping_objects[sequence.seq_hash] = mapping

    def get_or_create_mapping_obj(self,sequence):
        try:
            return self.mapping_objects[Sequence.get_sequence_hash(sequence)]
        except KeyError:
            seq = Sequence(sequence)
            mapping = SequenceMapping(seq, None, None, None, None, None, None)
            seqdf = pd.DataFrame.from_records([seq.to_series()])
            mapdf = pd.DataFrame.from_records([mapping.to_series()])
            seqdf.to_sql('sequences',con=self.con.get_engine(), if_exists='append',index=False)
            mapdf.to_sql('sequence_mapping', con=self.con.get_engine(), if_exists='append',index=False)
            self.mapping_objects[seq.seq_hash] = mapping
            return mapping


    def process_job(self,job_dir : Path,
                    input_name='input.faa', positions_pu_name='output.txt.positions_pu_1.mat',
                    prop_pvals_pu_name='output.txt.prop_pvals_pu_1.mat',
                    prop_scores_pu_name='output.txt.prop_scores_pu_1.mat',
                    prop_types_pu_name='output.txt.prop_types_pu_1.mat',
                    motif_info_name='output.txt.motif_info_1.mat',
                    models_name='output.txt.models_1.mat',
                    notes_name='output.txt.notes_1.mat',
                    features_name='output.txt.feats_1.mat') -> List[MutPred2Output.T]:
        """
        Write the MutPred2 for a given job to the Database

        Required Arguments:
        - job_dir : Path : Directory containing the job's input and output files

        Optional Arguments:
        - input_name : str : Default 'input.faa'
        - positions_pu_name : str : Default 'output.txt.positions_pu_1.mat'
        - prop_pvals_pu_name : str : Default 'output.txt.prop_pvals_pu_1.mat'
        - prop_scores_pu_name : str : Default 'output.txt.prop_scores_pu_1.mat'
        - prop_types_pu_name : str : Default 'output.txt.prop_types_pu_1.mat'

        """
        sequence = str(next(iter(parse(job_dir/input_name, 'fasta').records)).seq)
        mapping = self.get_or_create_mapping_obj(sequence)
        output_df = pd.read_csv(job_dir / 'output.txt')
        mutations = [Mutation(mapping,sub,score) for sub,score in zip(output_df.Substitution.values,
                                                                  output_df['MutPred2 score'].values)]
        positions_pu = loadmat(job_dir/positions_pu_name)['positions_pu']
        positions_pu = np.concatenate((positions_pu,
                                    np.ones((positions_pu.shape[0], 1)) * -1),
                                    axis=1)
        pvals_pu = loadmat(job_dir/prop_pvals_pu_name)['prop_pvals_pu']
        scores_pu = loadmat(job_dir/prop_scores_pu_name)['prop_scores_pu']
        prop_types_pu = loadmat(job_dir/prop_types_pu_name)['prop_types_pu']
        motif_info = loadmat(job_dir/motif_info_name)['motif_info']
        models = loadmat(job_dir / models_name)['models'].ravel()
        notes = loadmat(job_dir / notes_name)['notes']
        feats = loadmat(job_dir / features_name)['feats']
        output_records = []
        motif_records = []
        sequence_feature_sets = []
        substitution_feature_sets = []
        pssm_feature_sets = []
        conservation_feature_sets = []
        homology_feature_sets = []
        structure_feature_sets = []
        function_feature_sets = []
        for mutation,pos,pval,score,types,motifs,note, feat in zip(mutations,positions_pu, pvals_pu,scores_pu,prop_types_pu,motif_info,notes, feats):
            mechanisms = MutPred2Output.read_mechanisms(pos,pval,score,types)
            output_record = MutPred2Output(mutation,mechanisms)
            output_records.append(output_record)
            motif_mech = [m for m in output_record.mechanisms if m.name == "Motifs"][0]
            motif_records += Motif.motifs_from_string(mapping, mutation,motifs, motif_mech.posterior, motif_mech.pvalue)
            sequence_feature_sets.append(mapping, mutation, Features_Sequence(feat[:184]))
            substitution_feature_sets.append(mapping, mutation, Features_Substitution(feat[184:630]))
            pssm_feature_sets.append(mapping, mutation, Features_PSSM(feat[630:799]))
            conservation_feature_sets.append(mapping, mutation, Features_Conservation(feat[799:1036]))
            homology_feature_sets.append(mapping, mutation, Features_Homology(feat[1036:1056]))
            structure_feature_sets.append(mapping, mutation, Features_Structure(feat[1056:1135]))
            function_feature_sets.append(mapping, mutation, Features_Function(feat[1135:]))
        self.write_mutations(output_records)
        self.write_formatted_outputs(output_records)
        self.write_features(sequence_feature_sets, substitution_feature_sets, pssm_feature_sets, conservation_feature_sets, homology_feature_sets, structure_feature_sets, function_feature_sets)

    def write_mutations(self,output_records):
        mutation_df = pd.DataFrame.from_records([o.mutation.to_series() for o in output_records])
        mutation_df.to_sql('mutations',con=self.con.get_engine(),if_exists='append',index=False)

    def write_formatted_outputs(self,output_records):
        output_df = pd.DataFrame.from_records([o.to_series(flatfile_info=False) for o in output_records])
        output_df.to_sql('formatted_output',con=self.con.get_engine(),if_exists='append',index=False)

    def write_features(self,sequence_feature_sets, substitution_feature_sets, pssm_feature_sets, conservation_feature_sets, homology_feature_sets, structure_feature_sets, function_feature_sets):
        pd.DataFrame.from_records([s.to_series() for s in sequence_feature_sets]).to_sql('features_sequence', con=self.con.get_engine(),index=False,if_exists='append')
        pd.DataFrame.from_records([s.to_series() for s in substitution_feature_sets]).to_sql('features_substitution', con=self.con.get_engine(),index=False,if_exists='append')
        pd.DataFrame.from_records([s.to_series() for s in pssm_feature_sets]).to_sql('features_pssm', con=self.con.get_engine(),index=False,if_exists='append')
        pd.DataFrame.from_records([s.to_series() for s in conservation_feature_sets]).to_sql('features_conservation', con=self.con.get_engine(),index=False,if_exists='append')
        pd.DataFrame.from_records([s.to_series() for s in homology_feature_sets]).to_sql('features_homology', con=self.con.get_engine(),index=False,if_exists='append')
        pd.DataFrame.from_records([s.to_series() for s in structure_feature_sets]).to_sql('features_structure', con=self.con.get_engine(),index=False,if_exists='append')
        pd.DataFrame.from_records([s.to_series() for s in function_feature_sets]).to_sql('features_function', con=self.con.get_engine(),index=False,if_exists='append')
