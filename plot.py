import streamlit as st
import pandas as pd
import numpy as np
from Bio.PDB.Polypeptide import one_to_index,index_to_three,protein_letters_3to1
from scipy.sparse import csr_array
import numpy as np
# import pymysql.cursors
import matplotlib.pyplot as plt
import seaborn as sns
from sqlalchemy import create_engine

user = st.secrets['DB_USER']
ip = st.secrets['DB_ADDR']
st.write(
    "ip_addr_from_app_secrets: ",ip)

cnx = create_engine(f'mysql+pymysql://{user}@{ip}/MutPred2_DB')

proteins = pd.read_sql_query("SELECT DISTINCT sequence_mapping.ensembl_prot_id from variant INNER JOIN sequence_mapping on \
variant.seq_hash=sequence_mapping.seq_hash", con=cnx)

def on_change():
    # QUERY DATABASE for variants
    query = st.session_state.ensp
    variants = pd.read_sql_query("SELECT * from variant INNER JOIN sequence_mapping on variant.seq_hash=sequence_mapping.seq_hash where sequence_mapping.ensembl_prot_id='{}';".format(query), con=cnx)
    score_mat = csr_array((variants['score'].values,
                    ([one_to_index(r) for r in variants.alternate_aa], variants.position.values)))
    # create inticators for missing values
    mask = csr_array((np.ones(len(score_mat.nonzero()[0])), score_mat.nonzero())).todense() == 0
    # Get gene symbol
    symbol = variants.gene_symbol.unique()[0]
    # Set up the matplotlib figure
    fig,ax = plt.subplots(figsize=(11, 3))
    # draw figure
    # colormap
    cmap = sns.color_palette("vlag", as_cmap=True)
    # draw heatmap
    g= sns.heatmap(score_mat.todense(),vmin=0, vmax=1,cmap=cmap,mask=mask,square=False, ax=ax)
    # set alternate aa y-axis labels
    _ = ax.set_yticks(np.arange(20) + .5,[protein_letters_3to1[index_to_three(i)] for i in np.arange(20)],rotation=0)
    ax.set_title('MutPred2 scores for protein {} ({})'.format(query, symbol))
    g.set_facecolor('xkcd:silver')
    # clear_output(wait=True)
    # display(ax.figure)
    # plt.close(ax.figure)
    return fig


text = st.selectbox('Ensembl Protein ID',
    options=sorted(list(proteins.ensembl_prot_id.values)),
    index=0,on_change=on_change,key='ensp',
)

fig = on_change()
st.pyplot(fig)