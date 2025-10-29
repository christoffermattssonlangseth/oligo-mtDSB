#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import pandas as pd
import os
import mygene
mg = mygene.MyGeneInfo()


# In[2]:


metadata = pd.read_csv('/Users/christoffer/work/karolinska/development/data/abc_atlas/metadata/WMB-10X/20241115/cell_metadata.csv')
cluster_to_anno = pd.read_csv('/Users/christoffer/work/karolinska/development/Allen_ABC/data/abc_atlas/metadata/WMB-taxonomy/20231215/views/cluster_to_cluster_annotation_membership_pivoted.csv')
base_dir = '/Users/christoffer/work/karolinska/development/data/abc_atlas/expression_matrices/WMB-10Xv2/20230630/'
files = os.listdir(base_dir)


# In[3]:


files = [
    'WMB-10Xv2-Isocortex-1-raw.h5ad',
    'WMB-10Xv2-OLF-raw.h5ad',
    'WMB-10Xv2-CTXsp-raw.h5ad',
    'WMB-10Xv2-TH-raw.h5ad',
    #'WMB-10Xv2-Isocortex-2-raw.h5ad',
    #'WMB-10Xv2-Isocortex-3-raw.h5ad',
    'WMB-10Xv2-MB-raw.h5ad',
    'WMB-10Xv2-HY-raw.h5ad',
    'WMB-10Xv2-HPF-raw.h5ad'
]


# In[4]:


ad_list = []
for file in files: 
    print(file)
    ad_ = sc.read_h5ad(base_dir+file)
    ad_list.append(ad_)


# In[5]:


ad = sc.concat(
    ad_list,
)


# In[6]:


ad.obs


# In[7]:


del ad_list


# In[8]:


for meta in ['cluster_alias', 'donor_sex', 'dataset_label','x','y']:
    mapping_dict = dict(zip(metadata['cell_barcode'], metadata[meta]))
    ad.obs[meta] = ad.obs['cell_barcode'].map(mapping_dict)


# In[9]:


for meta in ['neurotransmitter', 'class', 'subclass', 'supertype','cluster']:
    mapping_dict = dict(zip(cluster_to_anno['cluster_alias'], cluster_to_anno[meta]))
    ad.obs[meta] = ad.obs['cluster_alias'].map(mapping_dict)


# In[10]:


ad.obs['class'].value_counts()


# In[13]:


# Normalizing to median total counts
sc.pp.normalize_total(ad)
# Logarithmize the data
sc.pp.log1p(ad)


# In[14]:


ad.write('/Users/christoffer/work/karolinska/development/data/abc_atlas/combined_scRNAseq.h5ad')


# In[2]:


ad = sc.read_h5ad('/Users/christoffer/work/karolinska/development/data/abc_atlas/combined_scRNAseq.h5ad')


# In[3]:


ad.var


# In[4]:


ad_sp = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_LLM_anno.h5ad')


# In[5]:


import re
import pandas as pd
from pathlib import Path

def load_ensembl_gene_map(gtf_path):
    """
    Parse an Ensembl GTF and return a dict: {ensembl_gene_id: gene_symbol}.
    Keeps only 'gene' features; strips version suffixes.
    """
    gene_map = {}
    pat_id   = re.compile(r'gene_id "([^"]+)"')
    pat_name = re.compile(r'gene_name "([^"]+)"')
    with open(gtf_path, "r") as fh:
        for line in fh:
            if line.startswith("#"): 
                continue
            # Only keep 'gene' feature lines to avoid huge memory use
            # GTF columns: seqname, source, feature, start, end, score, strand, frame, attributes
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            m_id = pat_id.search(attrs)
            m_nm = pat_name.search(attrs)
            if not (m_id and m_nm):
                continue
            gid = m_id.group(1).split('.')[0]   # strip version, e.g. ENSMUSG... .1 → base
            gnm = m_nm.group(1)
            # keep first occurrence (usually fine). If you want, prefer protein_coding using parts[1]/attrs
            if gid not in gene_map:
                gene_map[gid] = gnm
    return gene_map


# In[6]:


# --- use it ---
# Point to your Mus musculus Ensembl GTF (e.g., Mus_musculus.GRCm39.109.gtf)
gtf_file = "//Users/christoffer/Downloads/Mus_musculus.GRCm39.109.gtf"
gene_map = load_ensembl_gene_map(gtf_file)

# ids_to_map: list/Series/Index of 32,285 Ensembl IDs (with/without version)
def map_ids_to_symbols(ids, gene_map):
    # vectorized mapping via pandas for speed and NaN handling
    s = pd.Series(ids, dtype="string")
    base = s.str.replace(r"\.\d+$", "", regex=True)  # strip version suffixes
    symbols = base.map(gene_map).astype("string")
    return symbols.tolist()

# Example:
# ids_to_map = adata.var_names.tolist()
# symbols = map_ids_to_symbols(ids_to_map, gene_map)
# 


# In[7]:


ids_to_map = ad.var_names.tolist()


# In[8]:


symbols = map_ids_to_symbols(ids_to_map, gene_map)


# In[9]:


ad.var["gene_symbol"] = symbols


# In[10]:


ad.var = ad.var.set_index('gene_symbol')


# In[15]:


ad._inplace_subset_var(~ad.var.index.isna())


# In[18]:


common_genes = ad.var_names.intersection(ad_sp.var_names)


# In[23]:


ad_sp_sub = ad_sp[:, common_genes]


# In[25]:


ad.var_names_make_unique()  # in-place; appends -1, -2, ...


# In[26]:


ad_sp_sub = ad_sp[:, common_genes]
ad_sub = ad[:, common_genes]


# In[28]:


ad_sub.write('/Users/christoffer/work/karolinska/development/data/abc_atlas/combined_scRNAseq_for_tangram_subset_genes.h5ad')


# In[11]:


import pandas as pd

# make sure var_names are string, replace NaN with something safe
ad.var_names = pd.Index(ad.var_names.astype(str))
ad_sp.var_names = pd.Index(ad_sp.var_names.astype(str))



# In[ ]:


# optionally drop "nan" entries
ad = ad[:, ad.var_names != '<NA>'].copy()
ad_sp = ad_sp[:, ad_sp.var_names != "nan"].copy()


# In[ ]:


common_genes = ad.var_names.intersection(ad_sp.var_names)

ad_sub    = ad[:, common_genes].copy()
ad_sp_sub = ad_sp[:, common_genes].copy()


# In[ ]:




