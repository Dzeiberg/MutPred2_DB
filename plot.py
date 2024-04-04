import streamlit as st
import pandas as pd
import numpy as np
from Bio.PDB.Polypeptide import one_to_index,index_to_three,protein_letters_3to1
from scipy.sparse import csr_array
import numpy as np
import mysql.connector
import matplotlib.pyplot as plt
import seaborn as sns

user = st.secrets['DB_USER']
db_ipaddr = st.secrets['DB_ADDR']
db_pwd = st.secrets['DB_PWD']

con = mysql.connector.connect(user=user,
                              host=db_ipaddr,
                              password=db_pwd,
                              database='MutPred2_DB')

proteins = pd.read_sql_query("SELECT DISTINCT sequence_mapping.ensembl_prot_id,sequence_mapping.gene_symbol from variant INNER JOIN sequence_mapping on \
variant.seq_hash=sequence_mapping.seq_hash", con=con)

def on_change():
    # QUERY DATABASE for variants
    query = st.session_state.ensp.split("(")[0].strip()
    variants = pd.read_sql_query("SELECT * from variant INNER JOIN sequence_mapping on variant.seq_hash=sequence_mapping.seq_hash where sequence_mapping.ensembl_prot_id='{}';".format(query), con=con)
    variants = variants.loc[:,~variants.columns.duplicated()]
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
    return fig,variants


text = st.selectbox('Ensembl Protein ID',
    options=sorted(list(proteins.apply(lambda r: '{} ({})'.format(r['ensembl_prot_id'],r['gene_symbol']),axis=1).values),key=lambda s: s.split("(")[1][:-1]),
    index=0,on_change=on_change,key='ensp',
)

fig,variants = on_change()
df = variants.loc[:,['reference_aa','position', 'alternate_aa', 'score', 'ensembl_prot_id', "ensembl_nuc_id", 'ensembl_gene_id','gene_symbol','MANE_select']].reset_index(drop=True)
st.pyplot(fig)
st.dataframe(df)
