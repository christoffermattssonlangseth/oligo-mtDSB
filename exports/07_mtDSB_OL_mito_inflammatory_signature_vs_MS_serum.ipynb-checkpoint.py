#!/usr/bin/env python
# coding: utf-8

# # From mtDNA damage to the pre-onset MS signature
# 
# This notebook explores how mitochondrial DNA double-strand breaks (mtDSBs) in oligodendrocytes (OLs) recapitulate molecular signatures observed in **pre-symptomatic multiple sclerosis (MS)** serum proteomics (DOI: 10.1038/s41591-025-04014-w).
# 
# Using Xenium spatial transcriptomics data from the mtDSB-OL model, we profile key stress-response and inflammatory modules:
# 
# - **Integrated stress response (ISR) & UPR:** Atf4, Ddit3, Hspa5, Atf5  
# - **Mitochondrial & redox regulation:** Mtf1, Mt1/2, Sod2, Hif1a  
# - **Metabolic reprogramming:** Hk2, Pfkl, Ldha, Slc16a1/3  
# - **Innate immune sensing:** Tmem173, Ifit1–3, Isg15, Rsad2  
# - **Inflammatory signaling:** Nfkb2, Tnfaip3, Ccl2  
# - **Myelin program integrity:** Mbp, Plp1, Mog
# 
# Together, these data support the concept that **mitochondrial genome instability in OLs initiates a mito-inflammatory transcriptional program**, mirroring early systemic MS-associated changes long before symptomatic onset or axonal injury.

# In[2]:


import scanpy as sc
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib import colors as mcolors
import warnings
warnings.filterwarnings("ignore")


def plot_spatial_compact_fast(
    ad,
    color="leiden_2",
    groupby="sample_id",
    spot_size=8,
    cols=3,
    height=8,
    legend_col_width=1.2,
    palette=None,            # dict {cat:"#hex"} or list
    rasterized=True,         # big speedup for vectors/PDFs
    invert_y=True,           # match Scanpy orientation
    dpi=120,                 # lower dpi → faster
):
    # ----- 0) Preconditions -----
    if "spatial" not in ad.obsm:
        raise ValueError("ad.obsm['spatial'] is required.")

    coords = np.asarray(ad.obsm["spatial"])[:, :2]
    cats = ad.obs[color].astype("category")
    cat_names = cats.cat.categories
    cat_codes = cats.cat.codes.to_numpy()  # -1 for NaN

    # ----- 1) Build shared palette (RGBA array) -----
    if isinstance(palette, dict):
        col_list = [palette[c] for c in cat_names]
    elif isinstance(palette, (list, tuple)):
        if len(palette) < len(cat_names):
            raise ValueError("Palette shorter than number of categories.")
        col_list = list(palette)[:len(cat_names)]
    elif f"{color}_colors" in ad.uns:
        col_list = list(ad.uns[f"{color}_colors"])
        if len(col_list) != len(cat_names):
            raise ValueError(f"{color}_colors length != categories.")
    else:
        base = getattr(getattr(__import__("scanpy").pl.palettes, "default_64", []), "__iter__", None)
        base = list(__import__("scanpy").pl.palettes.default_64
                    if hasattr(__import__("scanpy").pl.palettes, "default_64")
                    else __import__("scanpy").pl.palettes.default_102)
        reps = int(np.ceil(len(cat_names) / len(base)))
        col_list = (base * reps)[:len(cat_names)]

    # store for consistency elsewhere
    ad.uns[f"{color}_colors"] = col_list

    # convert to RGBA float array for fast indexing
    rgba = np.array([mcolors.to_rgba(c) for c in col_list], dtype=float)
    # map codes -> rgba; handle -1 (NaN) as transparent
    colors_arr = np.empty((cat_codes.size, 4), dtype=float)
    colors_arr[cat_codes >= 0] = rgba[cat_codes[cat_codes >= 0]]
    colors_arr[cat_codes < 0] = (0, 0, 0, 0)

    # ----- 2) Precompute group indices (no per-iteration masks) -----
    gvals = ad.obs[groupby].astype(str).to_numpy()
    uniq_groups, gcodes = np.unique(gvals, return_inverse=True)
    # list of index arrays per group
    group_indices = [np.flatnonzero(gcodes == gi) for gi in range(len(uniq_groups))]

    # ----- 3) Figure layout -----
    n = len(uniq_groups)
    rows = int(np.ceil(n / cols))
    panel_w = height * cols * 0.6 / rows
    fig_w = panel_w + legend_col_width

    plt.ioff()  # speed: disable interactive redraws
    fig = plt.figure(figsize=(fig_w, height), dpi=dpi, constrained_layout=False)
    gs = GridSpec(
        rows, cols + 1, figure=fig,
        width_ratios=[1]*cols + [legend_col_width / (fig_w - legend_col_width)],
        wspace=0.02, hspace=0.02
    )

    # ----- 4) Panels (fast scatter) -----
    for i, sid in enumerate(uniq_groups):
        r, c = divmod(i, cols)
        ax = fig.add_subplot(gs[r, c])

        idx = group_indices[i]
        if idx.size:
            xy = coords[idx]
            # note: turn off edgecolor & use rasterized for speed
            sca = ax.scatter(
                xy[:, 0], xy[:, 1],
                c=colors_arr[idx],
                s=spot_size,
                marker='o',
                linewidths=0,
                rasterized=rasterized
            )
        ax.set_title(str(sid), fontsize=9, pad=2)
        ax.set_aspect("equal")
        if invert_y:
            ax.invert_yaxis()
        ax.set_axis_off()

    # blank unused
    for j in range(n, rows * cols):
        r, c = divmod(j, cols)
        fig.add_subplot(gs[r, c]).axis("off")

    # ----- 5) Legend (single, shared) -----
    ax_leg = fig.add_subplot(gs[:, -1])
    ax_leg.axis("off")
    handles = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=col_list[k], markersize=7, label=str(cat))
        for k, cat in enumerate(cat_names)
    ]
    ax_leg.legend(handles=handles, title=color, frameon=False, loc="center left")

    fig.subplots_adjust(left=0.01, right=0.98, top=0.98, bottom=0.02, wspace=0.02, hspace=0.02)
    plt.ion()  # ✅ correct
    plt.show()


# In[3]:


adata = sc.read_h5ad('../data/mtDNA_DSB_5k_clustered_annotation_with_rbd_OL_updated.h5ad')


# In[4]:


list(adata.obs.cell_class.unique())


# In[6]:





# In[7]:


adata_OL


# In[21]:


import scanpy as sc
import pandas as pd

# your AnnData: adata_OL  (obs: cell_types / sample_id / condition, etc.)
panels = {
    "IL3_axis": ["Il3ra","Csf2rb","Jak2","Stat5a","Stat5b","Socs1","Socs3","Ptpn6","Icam1","Vcam1","Ccl2","Ccl5","Cxcl10"],
    "NFkB": ["Nfkb1","Nfkb2","Rela","Relb","Chuk","Ikbkb","Nfkbia","Tnfaip3","Tnf","Il1b","Ccl2","Cxcl10"],
    "ISR_UPR": ["Atf4","Ddit3","Eif2ak3","Ppp1r15a","Asns","Hspa5","Xbp1","Atf6","Herpud1","Dnajb9","Atf5","Hspd1","Hspe1","Lonp1","Clpp"],
    "Mito_Redox": ["Hif1a","Sod2","Cat","Gpx1","Hmox1","Nqo1","Mt1","Mt2","Mtf1","Ppargc1a","Tfam","Dnm1l","Opa1","Mfn1","Mfn2"],
    "Lactate_Glycolysis": ["Hk2","Pfkl","Pkm","Pdk1","Ldha","Ldhb","Slc16a1","Slc16a3"],
    "cGAS_STING_ISG": ["Mb21d1","Tmem173","Irf3","Irf7","Stat1","Isg15","Ifit1","Ifit2","Ifit3","Oasl1","Oasl2","Rsad2","Mx1"],
    "Acute_Complement_Microglia": ["Aif1","Cx3cr1","Trem2","Cst7","Lgals3","C1qa","C1qb","C1qc","C3","C4b","Serpina3n","Lcn2","Saa3"],
    "Myelin_OLcore": ["Mbp","Plp1","Mog","Mag","Cnp","Ugt8a","Gpr17","Sox10","Id2","Id4","Bcl2l1","Bax","Bbc3"]
}

# keep only genes present
panels = {k:[g for g in v if g in adata.var_names] for k,v in panels.items()}

# quick dotplot by condition and cell type (e.g., OL only)
ad_sub = adata[adata.obs.cell_class.str.contains('Oligod')]


# In[24]:


sc.pl.dotplot(
       ad_sub,
       var_names=panels,
       groupby="condition",     # e.g., ["control","mtDSB"]
       standard_scale="var",
       title=name,
        swap_axes=True,
   figsize=(3,14)
   )


# In[ ]:





# In[ ]:




