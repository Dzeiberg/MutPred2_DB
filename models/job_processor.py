from models.sequence_mapping import SequenceMapping
from .sequence import Sequence
from .mutation import Mutation
from .output_record import MutPred2Output
from .sql_connection import SQL_Connection
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
            seqdf.to_sql('sequences',con=self.con, if_exists='append',index=False)
            mapdf.to_sql('sequence_mapping', con=self.con, if_exists='append',index=False)
            self.mapping_objects[seq.seq_hash] = mapping
            return mapping


    def process_job(self,job_dir : Path,
                    input_name='input.faa', positions_pu_name='output.txt.positions_pu_1.mat',
                    prop_pvals_pu_name='output.txt.prop_pvals_pu_1.mat',
                    prop_scores_pu_name='output.txt.prop_scores_pu_1.mat',
                    prop_types_pu_name='output.txt.prop_types_pu_1.mat') -> List[MutPred2Output.T]:
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
        output_records = []
        for mutation,pos,pval,score,types in zip(mutations,positions_pu, pvals_pu,scores_pu,prop_types_pu):
            mechanisms = MutPred2Output.read_mechanisms(pos,pval,score,types)
            output_records.append(MutPred2Output(mutation,mechanisms))
        self.write_mutations(output_records)
        self.write_formatted_outputs(output_records)

    def write_mutations(self,output_records):
        mutation_df = pd.DataFrame.from_records([o.mutation.to_series() for o in output_records])
        mutation_df.to_sql('mutations',con=self.con.get_engine(),if_exists='append',index=False)

    def write_formatted_outputs(self,output_records):
        output_df = pd.DataFrame.from_records([o.to_series(flatfile_info=False) for o in output_records])
        output_df = output_df.drop([c for c in output_df.columns if "Motif" in c],axis=1)
        output_df.to_sql('formatted_output',con=self.con.get_engine(),if_exists='append',index=False)