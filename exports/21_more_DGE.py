#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import pandas as pd
import numpy as np


# In[4]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_manual_annotation.h5ad')


# In[12]:


adata_raw = sc.read_h5ad('../data/mtDNA_DSB_5k_raw.h5ad')


# In[15]:


import numpy as np
import pandas as pd
from scipy import sparse

def _report_dups(name, idx):
    dup = pd.Index(idx)
    d = dup[dup.duplicated()].unique().tolist()
    if d:
        print(f"[WARN] {name}: {len(d)} duplicated names e.g. {d[:5]}")
    else:
        print(f"[OK] {name}: all unique")

def attach_counts_intersection_safe(adata, adata_raw, layer_name="counts", collapse_var=False):
    # 1) report duplicates
    _report_dups("adata.obs_names", adata.obs_names)
    _report_dups("adata.var_names", adata.var_names)
    _report_dups("adata_raw.obs_names", adata_raw.obs_names)
    _report_dups("adata_raw.var_names", adata_raw.var_names)

    # 2) handle duplicates
    if collapse_var:
        # collapse duplicated var names by summing columns (good for counts)
        def _collapse(A, names):
            df = pd.DataFrame.sparse.from_spmatrix(A) if sparse.issparse(A) else pd.DataFrame(A)
            df.columns = names
            return df.groupby(level=0, axis=1).sum()
        # collapse both objects on var (genes)
        A1 = _collapse(adata.X, adata.var_names)
        A2 = _collapse(adata_raw.X, adata_raw.var_names)
        # keep obs in original order
        A1.index = adata.obs_names
        A2.index = adata_raw.obs_names
        # replace X and var_names with collapsed
        adata = adata[:, []].copy()
        adata.X = A1.values
        adata.var_names = A1.columns.astype(str)
        adata.obs_names = A1.index.astype(str)

        adata_raw = adata_raw[:, []].copy()
        adata_raw.X = A2.values
        adata_raw.var_names = A2.columns.astype(str)
        adata_raw.obs_names = A2.index.astype(str)
    else:
        # just force uniqueness by appending suffixes
        adata.var_names_make_unique()
        adata.obs_names_make_unique()
        adata_raw.var_names_make_unique()
        adata_raw.obs_names_make_unique()

    # 3) intersect axes
    common_obs = adata.obs_names.intersection(adata_raw.obs_names)
    common_var = adata.var_names.intersection(adata_raw.var_names)

    ad_sub  = adata[common_obs, common_var].copy()
    raw_sub = adata_raw[common_obs, common_var]

    # 4) attach raw into a layer
    Xraw = raw_sub.X
    if sparse.issparse(Xraw):
        Xraw = Xraw.copy()
    else:
        Xraw = np.asarray(Xraw).copy()
    ad_sub.layers[layer_name] = Xraw

    print(f"[attach_counts] kept {ad_sub.n_obs} cells and {ad_sub.n_vars} genes "
          f"(intersected). Raw stored in .layers['{layer_name}'].")
    return ad_sub


# In[16]:


# If you want to simply make names unique (fastest):
adata = attach_counts_intersection_safe(adata, adata_raw, layer_name="counts", collapse_var=False)

# If your var (gene) names are symbols with duplicates and you prefer to SUM duplicates first:
# adata = attach_counts_intersection_safe(adata, adata_raw, layer_name="counts", collapse_var=True)


# In[20]:


adata.layers['counts']


# In[121]:


adata.write('../data/mtDNA_DSB_5k_clustered_manual_annotation_with_raw.h5ad')


# In[105]:


list(adata.obs.cell_class.unique())


# In[106]:


adata_OL = adata[adata.obs.cell_class.isin(['Mature oligodendrocytes','Oligodendrocytes precursor cells'])]


# In[107]:


adata_sub = adata_OL[adata_OL.obs["condition"]=='mtDSB']
adata_sub = adata_sub[adata_sub.obs["age"]=='60']

adata_rest = adata_OL[adata_OL.obs["condition"]!='mtDSB']
adata_rest = adata_rest[adata_rest.obs["age"]=='60']


# In[112]:


import numpy as np
import pandas as pd

# mean expression per gene for each group
mean_interest = np.asarray(adata_sub.layers['counts'].mean(axis=0)).ravel()
mean_rest = np.asarray(adata_rest.layers['counts'].mean(axis=0)).ravel()

diff = mean_interest - mean_rest
fold = (mean_interest + 1e-6) / (mean_rest + 1e-6)

# rough z-score of differences (per gene relative to all genes)
zscore = (diff - diff.mean()) / diff.std()

diff_df = pd.DataFrame({
    "gene": adata.var_names,
    "mean_interest": mean_interest,
    "mean_rest": mean_rest,
    "diff": diff,
    "log2FC": np.log2(fold),
    "zscore": zscore
}).sort_values("zscore", ascending=False)

diff_df.head()


# In[114]:


import numpy as np
import pandas as pd

# mean expression per gene for each group
mean_interest = np.asarray(adata_sub.layers['counts'].mean(axis=0)).ravel()
mean_rest = np.asarray(adata_rest.layers['counts'].mean(axis=0)).ravel()

diff = mean_interest - mean_rest
fold = (mean_interest + 1e-6) / (mean_rest + 1e-6)

# rough z-score of differences (per gene relative to all genes)
zscore = (diff - diff.mean()) / diff.std()

diff_df = pd.DataFrame({
    "gene": adata.var_names,
    "mean_interest": mean_interest,
    "mean_rest": mean_rest,
    "diff": diff,
    "log2FC": np.log2(fold),
    "zscore": zscore
}).sort_values("zscore", ascending=False)


# In[115]:


diff_df_sorted = diff_df.sort_values(by=[ "mean_interest",'log2FC'], ascending=[False, False])
print(diff_df_sorted[diff_df_sorted.log2FC > 0.3].head(60))


# In[111]:


diff_df_sorted = diff_df.sort_values(by=['log2FC',"mean_interest"], ascending=[False, False])
print(diff_df_sorted[diff_df_sorted.log2FC > 0.2].head(60))


# In[120]:


genes = ["Hif1a","Nfe2l2","Akt1","Pik3cd","Hmox1","Txnip","Slc16a1","Slc16a3",
         "Ldha","Ldhb",'Aldoa',"Mfn1","Mfn2","Opa1","Pkm","Sirt2","Ppargc1a"]
print(diff_df_sorted[diff_df_sorted.gene.isin(genes)])


# ## Functional Modules in mtDSB Oligodendrocytes
# 
# ### 1. Oxidative Stress & Detoxification
# - **Mt2** – metallothionein, binds Zn/Cu, scavenges ROS  
# - **Gstp1** – glutathione detox enzyme  
# - **Sqstm1 (p62)** – oxidative stress sensor, links ROS to autophagy  
# - **Nfe2l1** – TF regulating antioxidant genes  
# - **Hspd1, Hspa9** – mitochondrial chaperones for ROS stress  
# ➡️ Evidence for **oxidative stress and mitochondrial redox imbalance**
# 
# ---
# 
# ### 2. Integrated Stress Response (ISR) & UPRmt
# - **Atf4, Jun** – stress-activated TFs driving ISR/UPRmt  
# - **Hspa5 (BiP/GRP78), Hspd1 (HSP60), Hspa9** – ER/mitochondrial chaperones  
# - **Hsph1** – HSP110 family, protein folding/stress tolerance  
# ➡️ Indicates **protein misfolding and translational stress downstream of mitochondrial dysfunction**
# 
# ---
# 
# ### 3. Antigen Presentation & Immune Signaling
# - **B2m, H2-D1, H2-K1** – MHC-I antigen presentation  
# - **Ctss (Cathepsin S)** – lysosomal protease for antigen processing  
# - **Cd40, Cd44, Cd22, Cd3e** – co-stimulatory/immune interaction molecules  
# ➡️ Suggests **stressed OLs present antigen and engage immune surveillance**
# 
# ---
# 
# ### 4. Axonal/Myelin Transport Stress
# - **Kif5a, Kif5b** – kinesin motors for axonal/myelin cargo  
# - **Dync1li1** – dynein light intermediate chain, retrograde transport  
# - **Ptpra, Itgb1, Mpzl1, Sorbs1** – adhesion/cell interaction molecules  
# ➡️ Points to **disturbed axonal transport and OL–axon coupling**
# 
# ---
# 
# ### 5. Astrocytic Reactivity & Glial Crosstalk
# - **Ndrg2, Gfap, Mlc1** – hallmark astrocytic/reactive gliosis genes  
# - **S100a1, Calb2** – calcium-binding proteins in reactive astrocytes  
# ➡️ Reflects **astrocyte activation secondary to OL mtDSB stress**
# 
# ---
# 
# ### 6. Inflammatory Signaling & Cytokines
# - **Nmu (Neuromedin U)** – neuropeptide with pro-inflammatory activity  
# - **Ccl3 (MIP-1α)** – chemokine recruiting myeloid cells  
# - **Irf4** – immune TF controlling cytokine expression  
# ➡️ Suggests **immune-modulatory signaling in the OL/astro niche**
# 
# ---
# 
# ### 7. Metabolic Remodeling
# - **Hadhb** – mitochondrial β-oxidation enzyme  
# - **Parvb, Sorbs1** – cytoskeletal/metabolic adaptors  
# - **Fosl1, Myc** – TFs linked to metabolic reprogramming and proliferation  
# ➡️ Evidence for **shifts in mitochondrial and metabolic regulation**
# 
# ---
# 
# ## ✅ Overall Takeaway
# mtDSB oligodendrocytes show a **coherent multi-pathway stress program**:
# - **Redox stress** → metallothioneins, glutathione enzymes  
# - **Mitochondrial/ER stress** → ISR/UPRmt activation  
# - **Immune presentation** → MHC-I and lysosomal proteases  
# - **Axonal transport disruption** → kinesins/dyneins  
# - **Astro reactivity** → Ndrg2, Gfap upregulation  
# - **Pro-inflammatory signaling** → Nmu, Ccl3  
# ➡️ Together, these changes suggest OLs under mtDNA damage are **alive but stressed**, tipping the microenvironment toward **immune activation and glial crosstalk**.

# In[ ]:




