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
import json
import re
import os
import logging


class Processor:
    def __init__(self, con : SQL_Connection):
        self.con = con
        self.init_mapping_objects()
        self.initialize_existing_tuples()
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
        

    def create_status_file(self, job_dir : Path, status : str):
        with open(f"{job_dir}/status.json",'w') as s:
            status_output = json.dumps({ "status": status })
            s.write(status_output)

    def get_job_status(self, job_dir: Path):
        if os.path.isfile(f"{job_dir}/status.json"):
            job_status = json.load(open(f"{job_dir}/status.json"))["status"]
        else:
            job_status = "INCOMPLETE"
        return job_status

    

    def read_input_file(self, job_dir : Path):

        directory = os.listdir(job_dir)

        re_pattern = re.compile('.*.faa$')

        files = sorted([f for f in directory if re_pattern.match(f)])

        if len(files) == 1:
            logging.warning(f'Reading input file: {files[0]}')
            return files[0]
        
        elif len(files)==0:

            logging.warning(f"No input file found ending in '.faa' for job {job_dir}")
            return ""

        else:
            logging.warning(f"Too many input files with '.faa' for job {job_dir}, only one expected")
            return ""
    

    def read_output_file(self, job_dir : Path):

        directory = os.listdir(job_dir)

        re_pattern = re.compile('.*.txt$')

        files = sorted([f for f in directory if re_pattern.match(f)])

        if len(files) == 1:
            return files[0]
        
        elif len(files)==0:
            logging.warning(f"No output file found ending in '.txt' for job {job_dir}")

        else:
            logging.warning(f"Too many output files with '.txt' for job {job_dir}, only one expected")


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
            data = np.concatenate([loadmat(f"{job_dir}/{f}")[key_value].ravel() for f in files if len(loadmat(f"{job_dir}/{f}")[key_value].ravel()) > 0])
        else:
            data = np.concatenate([loadmat(f"{job_dir}/{f}")[key_value] for f in files if len(loadmat(f"{job_dir}/{f}")[key_value]) > 0])

        return data


    def process_job(self, job_dir : Path) -> List[MutPred2Output.T]:
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

        if self.get_job_status() == "COMPLETE":
            pass
        else:
            input_name = self.read_input_file(job_dir=job_dir)
            if input_name != "":
                sequence = str(next(iter(parse(job_dir/input_name, 'fasta').records)).seq)
                mapping = self.get_or_create_mapping_obj(sequence)
                input_pass = True
            else:
                input_pass = False


            output_name = self.read_output_file(job_dir=job_dir)
            if output_name != "":
                output_df = pd.read_csv(job_dir / output_name)
                output_pass = True
            else:
                output_pass = False


            if input_pass and output_pass:
            
                mutations = [Mutation(mapping,sub,score) for sub,score in zip(output_df.Substitution.values,
                                                                        output_df['MutPred2 score'].values)]
                
                positions_pu = self.read_mat_files(job_dir=job_dir, pattern='.*.txt.positions_pu_\d+.mat', key_value='positions_pu') #loadmat(job_dir/positions_pu_name)['positions_pu']
                positions_pu = np.concatenate((positions_pu,
                                            np.ones((positions_pu.shape[0], 1)) * -1),
                                            axis=1)
                pvals_pu        = self.read_mat_files(job_dir=job_dir, pattern='.*.txt.prop_pvals_pu_\d+.mat', key_value='prop_pvals_pu') #loadmat(job_dir/prop_pvals_pu_name)['prop_pvals_pu']
                scores_pu       = self.read_mat_files(job_dir=job_dir, pattern='.*.txt.prop_scores_pu_\d+.mat', key_value='prop_scores_pu') #loadmat(job_dir/prop_scores_pu_name)['prop_scores_pu']
                prop_types_pu   = self.read_mat_files(job_dir=job_dir, pattern='.*.txt.prop_types_pu_\d+.mat', key_value='prop_types_pu') #loadmat(job_dir/prop_types_pu_name)['prop_types_pu']
                motif_info      = self.read_mat_files(job_dir=job_dir, pattern='.*.txt.motif_info_\d+.mat', key_value='motif_info') #loadmat(job_dir/motif_info_name)['motif_info'].ravel()
                models          = self.read_mat_files(job_dir=job_dir, pattern='.*.txt.models_\d+.mat', key_value='models') #loadmat(job_dir / models_name)['models'].ravel()
                notes           = self.read_mat_files(job_dir=job_dir, pattern='.*.txt.notes_\d+.mat', key_value='notes') #loadmat(job_dir / notes_name)['notes']
                feats           = self.read_mat_files(job_dir=job_dir, pattern='.*.txt.feats_\d+.mat', key_value='feats') #loadmat(job_dir / features_name)['feats']

                if len(set([len(pvals_pu),len(scores_pu),len(prop_types_pu),len(motif_info),len(models),len(notes),len(feats)])) == 1:
                    
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
                        motif_records += Motif.motifs_from_string(mapping, mutation,motifs.item(), motif_mech.posterior, motif_mech.p_value)
                        sequence_feature_sets.append(Features_Sequence(mapping, mutation, feat[:184]))
                        substitution_feature_sets.append(Features_Substitution(mapping, mutation, feat[184:630]))
                        pssm_feature_sets.append(Features_PSSM(mapping, mutation, feat[630:799]))
                        conservation_feature_sets.append(Features_Conservation(mapping, mutation, feat[799:1036]))
                        homology_feature_sets.append(Features_Homology(mapping, mutation, feat[1036:1056]))
                        structure_feature_sets.append(Features_Structure(mapping, mutation, feat[1056:1135]))
                        function_feature_sets.append(Features_Function(mapping, mutation, feat[1135:]))
                    
                    self.write_mutations(output_records)
                    self.write_formatted_outputs(output_records)
                    self.write_features(sequence_feature_sets, substitution_feature_sets, pssm_feature_sets, conservation_feature_sets, homology_feature_sets, structure_feature_sets, function_feature_sets)

                    self.create_status_file(job_dir, "COMPLETE")
                
                else:
                    logging.warning(f'{job_dir} not written to database because files are mismatched in size. See concatenated file lengths below.')
                    logging.warning(f"prop_pvals_pu: {len(pvals_pu)}")
                    logging.warning(f"prop_scores_pu: {len(scores_pu)}")
                    logging.warning(f"prop_types_pu: {len(prop_types_pu)}")
                    logging.warning(f"motif_info: {len(motif_info)}")
                    logging.warning(f"models: {len(models)}")
                    logging.warning(f"notes: {len(notes)}")
                    logging.warning(f"features: {len(feats)}")
                    self.create_status_file(job_dir, "INCOMPLETE")

            else:
                self.create_status_file(job_dir, "INCOMPLETE")

    def initialize_existing_tuples(self):
        logging.warning('initializing mutations')
        self.existing_mutations = NamedSet("existing_mutations")
        self.existing_mutations.update(pd.read_sql('select seq_hash,mutation from mutations',con=self.con.get_engine()).apply(tuple,axis=1).values)
        logging.warning('initializing formatted_output')
        self.existing_formatted_outputs = NamedSet("existing_formatted_outputs")
        self.existing_formatted_outputs.update(pd.read_sql('select seq_hash,mutation from formatted_output',con=self.con.get_engine()).apply(tuple,axis=1).values)
        logging.warning('initializing features_sequence')
        self.existing_features_sequence = NamedSet("existing_features_sequence")
        self.existing_features_sequence.update(pd.read_sql('select seq_hash,mutation from features_sequence',con=self.con.get_engine()).apply(tuple,axis=1).values)
        logging.warning('initializing features_substitution')
        self.existing_features_substitution = NamedSet("existing_features_substitution")
        self.existing_features_substitution.update(pd.read_sql('select seq_hash,mutation from features_substitution',con=self.con.get_engine()).apply(tuple,axis=1).values)
        logging.warning('initializing features_pssm')
        self.existing_features_pssm = NamedSet("existing_features_pssm")
        self.existing_features_pssm.update(pd.read_sql('select seq_hash,mutation from features_pssm',con=self.con.get_engine()).apply(tuple,axis=1).values)
        logging.warning('initializing features_conservation')
        self.existing_features_conservation = NamedSet("existing_features_conservation")
        self.existing_features_conservation.update(pd.read_sql('select seq_hash,mutation from features_conservation',con=self.con.get_engine()).apply(tuple,axis=1).values)
        logging.warning('initializing features_homology')
        self.existing_features_homology = NamedSet("existing_features_homology")
        self.existing_features_homology.update(pd.read_sql('select seq_hash,mutation from features_homology',con=self.con.get_engine()).apply(tuple,axis=1).values)
        logging.warning('initializing features_structure')
        self.existing_features_structure = NamedSet("existing_features_structure")
        self.existing_features_structure.update(pd.read_sql('select seq_hash,mutation from features_structure',con=self.con.get_engine()).apply(tuple,axis=1).values)
        logging.warning('initializing features_function')
        self.existing_features_function = NamedSet("existing_features_function")
        self.existing_features_function.update(pd.read_sql('select seq_hash,mutation from features_function',con=self.con.get_engine()).apply(tuple,axis=1).values)


    def filter_existing_and_update(self,df, tupset):
        tuple_series = df.loc[:,['seq_hash','mutation']].apply(tuple,axis=1)
        new_df = df.loc[tuple_series.apply(lambda tup: tup not in tupset)]
        new_len = len(new_df)
        old_len = len(df)
        if old_len != new_len:
            logging.warning(f'skipping {old_len - new_len} records already existing in table {tupset}')
        tupset.update(tuple_series.values)
        return new_df

    def write_mutations(self,output_records):
        mutation_df = pd.DataFrame.from_records([o.mutation.to_series() for o in output_records])
        out = self.filter_existing_and_update(mutation_df, self.existing_mutations)
        out.to_sql('mutations',con=self.con.get_engine(),if_exists='append',index=False)

    def write_formatted_outputs(self,output_records):
        output_df = pd.DataFrame.from_records([o.to_series(flatfile_info=False) for o in output_records])
        out = self.filter_existing_and_update(output_df, self.existing_formatted_outputs)
        out.to_sql('formatted_output',con=self.con.get_engine(),if_exists='append',index=False)

    def write_features(self,sequence_feature_sets, substitution_feature_sets,
                       pssm_feature_sets, conservation_feature_sets, homology_feature_sets,
                       structure_feature_sets, function_feature_sets):
        self.filter_existing_and_update(pd.DataFrame.from_records([s.to_series() for s in sequence_feature_sets]),
                                        self.existing_features_sequence).to_sql('features_sequence',
                                                                                con=self.con.get_engine(),
                                                                                index=False,
                                                                                if_exists='append')
        
        self.filter_existing_and_update(pd.DataFrame.from_records([s.to_series() for s in substitution_feature_sets]),
                                        self.existing_features_substitution).to_sql('features_substitution',
                                                                                    con=self.con.get_engine(),
                                                                                    index=False,
                                                                                    if_exists='append')
        
        self.filter_existing_and_update(pd.DataFrame.from_records([s.to_series() for s in pssm_feature_sets]),
                                        self.existing_features_pssm).to_sql('features_pssm',
                                                                            con=self.con.get_engine(),
                                                                            index=False,
                                                                            if_exists='append')
        
        self.filter_existing_and_update(pd.DataFrame.from_records([s.to_series() for s in conservation_feature_sets]),
                                        self.existing_features_conservation).to_sql('features_conservation',
                                                                                    con=self.con.get_engine(),
                                                                                    index=False,
                                                                                    if_exists='append')
        
        self.filter_existing_and_update(pd.DataFrame.from_records([s.to_series() for s in homology_feature_sets]),
                                        self.existing_features_homology).to_sql('features_homology',
                                                                                con=self.con.get_engine(),
                                                                                index=False,
                                                                                if_exists='append')
        
        self.filter_existing_and_update(pd.DataFrame.from_records([s.to_series() for s in structure_feature_sets]),
                                        self.existing_features_structure).to_sql('features_structure',
                                                                                 con=self.con.get_engine(),
                                                                                 index=False,
                                                                                 if_exists='append')
        
        self.filter_existing_and_update(pd.DataFrame.from_records([s.to_series() for s in function_feature_sets]),
                                        self.existing_features_function).to_sql('features_function',
                                                                                con=self.con.get_engine(),
                                                                                index=False,
                                                                                if_exists='append')

class NamedSet(set):
    def __init__(self,name):
        super().__init__()
        self.name = name
    def __repr__(self):
        return f"set {self.name} "# + super().__repr__()
