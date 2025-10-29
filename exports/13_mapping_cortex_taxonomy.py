#!/usr/bin/env python
# coding: utf-8

# This notebook will demonstrate the full functionality of MapMyCells by downoading one mouse and one human "legacy" (read: pre-MapMyCells) taxonomy from brain-map.org, splitting the data defining the taxonomies into test and training sets, and using MapMyCells to map the test set back onto the taxonomy, thus verifying the accuracy of the MapMyCells mapping. We will then Map the human legacy data onto the mouse legacy taxonomy and examine the implied correspondence between human and mouse cell types. Finally, we will map the human legacy data onto the [Yao et al. 2023](https://www.nature.com/articles/s41586-023-06812-z) Whole Mouse Brain taxonomy and perform a similar comparison.
# 
# This notebook will take several hours to run in total (most of that time is spent downloading data and splitting it into test and training sets). It will download and/or generate roughly 40 GB of data.

# First, let's configure out compute system to efficiently use MapMyCells' parallelization scheme.

# In[1]:


import os
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'


# ## Legacy mouse taxonomy
# 
# We will download the [Mouse Whole Cortex and Hippocampus 10x](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x) data from [Yao et al. 2021.](https://www.sciencedirect.com/science/article/pii/S0092867421005018?dgcid=rss_sd_all) We will then split the data into a test set and a training set. We will use the training set to create the data files necessary to use MapMyCells to map data onto this taxonomy. Finally, we will map the test set (which has already been annotated according to the taxonomy) back onto the taxonomy, using the pre-existing annotations in the test dataset to assess the accuracy of the MapMyCells mapping.

# In[2]:


import pathlib
import subprocess


# Create a directory into which we can download ABC Atlas assets as needed.

# In[3]:


abc_cache_dir = pathlib.Path('../data/abc_cache')


# In[4]:


data_dir = pathlib.Path('../data/cortex')
if not data_dir.exists():
    data_dir.mkdir(parents=True)


# ### Download the data
# 
# Define a helper function to download the data.

# In[5]:


def download_file(
        src_url,
        dst_path,
        force_download=False):
    """
    Download a file.

    Parameters
    ----------
    src_url:
        URL of the file to be downloaded
    dst_path:
        path where the file will be saved
    force_download:
        if True and dst_path exists, overwrite
        (otherwise, assume dst_path is the file you want)
    """
    dst_path = pathlib.Path(dst_path)
    if force_download and dst_path.exists():
        dst_path.unlink()

    if not dst_path.is_file():
        args = [
            "wget",
            src_url,
            "-O",
            str(dst_path),
            "-q"
        ]
        p = subprocess.Popen(args)
        exit_status = p.wait()
        if exit_status != 0:
            if metadata_path.exists():
                metadata_path.unlink()
            raise RuntimeError("Failure downloading metadata file")
        print(f"SUCCESSFULLY DOWNLOADED {dst_path}")
    else:
        print(f"{dst_path} ALREADY EXISTS; NO ACTION TAKEN")


# Below we will use the helper function we just defined to download the cell-by-gene matrix and cell annotations associated with the Allen Institute's 2020 [Mouse Whole Cortex and Hippocampus 10x study](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x) as published in [Yao et al. 2021](https://www.biorxiv.org/content/10.1101/2020.03.30.015214v1). This will download a 491 MB CSV file and a 5.3 GB HDF5 file.

# In[6]:


mouse_metadata_path = data_dir / "mouse_cortex_metadata.csv"


# In[7]:


get_ipython().run_cell_magic('time', '', 'download_file(\n    src_url="https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/metadata.csv",\n    dst_path=mouse_metadata_path\n)\nassert mouse_metadata_path.is_file()\n')


# In[8]:


mouse_data_path = data_dir / "mouse_cortex_expression_matrix.hdf5"


# In[9]:


get_ipython().run_cell_magic('time', '', 'download_file(\n    src_url="https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/expression_matrix.hdf5",\n    dst_path=mouse_data_path\n)\nassert mouse_data_path.is_file()\n')


# The metadata file we just downloaded contains many columns describing the individual cells in the cell-by-gene data. For our purposes, we are only insterested in the unique identifier of each cell, along with its taxonomic annotations (`class_label`, `subclass_label`, and `cluster_label`).

# In[11]:


import pandas as pd


# In[12]:


mouse_metadata_df = pd.read_csv(mouse_metadata_path)


# In[13]:


mouse_metadata_df.columns


# In[14]:


mouse_metadata_df = mouse_metadata_df[['sample_name', 'class_label', 'subclass_label', 'cluster_label']].set_index('sample_name')


# The cell-by-gene data also carries gene and sample name annotations.

# In[19]:


import h5py
import numpy as np


# In[20]:


with h5py.File(mouse_data_path, 'r') as src:
    mouse_gene_list = np.array([g.decode('utf-8') for g in src['data/gene'][()]])
    sample_name = np.array([s.decode('utf-8') for s in src['data/samples'][()]])


# Unfortunately, not every sample represented in the cell-by-gene dataset has a corresponding row in the metadata table.

# In[21]:


len(mouse_metadata_df)


# In[22]:


len(sample_name)


# The missing cells represent cells that failes quaity control for one reason or another. When splitting up our data into test and training sets, we must be careful to preserve the relationship between cells in the metadata table and cells in the cell-by-gene dataset.

# ### Map the genes to ENSEMBL IDs
# 
# We will now use the `abc_atlas_access` API to download a sqlite database which we will use to map data both across species (i.e. from mouse to human via orthologous gene relationships) and between different identifiers (gene symbols, NCBI IDs, and ENSEMBL IDs).
# 
# **Note:** at the moment, we are using a wrapper of `AbcProjectCache` defined in this codebase. This is because the gene_mapper_db file has not yet been made part of an official `abc_atlas_access` release.
# 
# This database contains all of the
# 
# - gene symbol to gene identifier
# - NCIB to ENSEMBL
# - cross-species ortholog
# 
# Relationships recognized by MapMyCells. The database is created and queried using [this codebase](https://github.com/AllenInstitute/mmc_gene_mapper). The version you are downloading represents a snapshot of the data published by NCBI and ENSEMBL as of August 4, 2025. Should you want to create your own database, [this Jupyter notebook](https://github.com/AllenInstitute/mmc_gene_mapper/blob/main/notebooks/gene_mapper_demo.ipynb) will show you how.
# 
# The path to this file will be passed in through the `gene_mapping.db_path` configuration parameter when actually mapping the data.

# In[23]:


import cell_type_mapper.test_utils.cache_wrapper as cache_module
import mmc_gene_mapper.mapper.mapper as mapper_module


# In[24]:


abc_cache = cache_module.AbcCacheWrapper.from_cache_dir(abc_cache_dir)


# In[25]:


gene_mapper_db_path = abc_cache.get_data_path(
    directory='mmc-gene-mapper',
    file_name='mmc_gene_mapper.2025-08-04'
)


# Let's instantiate an instance of the `MMCGeneMapper` class defined in [this codebase](https://github.com/AllenInstitute/mmc_gene_mapper).

# In[106]:


gene_mapper = mapper_module.MMCGeneMapper(db_path=gene_mapper_db_path)


# This class will allow us to map the genes from the cell-by-gene dataset into ENSEMBL IDs.

# In[21]:


mouse_genes_as_ensembl = gene_mapper.map_genes(
    gene_list=mouse_gene_list,
    dst_species='Mus musculus',
    dst_authority='ENSEMBL'
)['gene_list']


# In[107]:


for ii in range(10):
    print(f'original gene: {mouse_gene_list[ii]} -> {mouse_genes_as_ensembl[ii]}')


# ### Split the data into test and training sets
# 
# Define a function to help us split the cell-by-gene dataset into a test set and a training set.

# In[26]:


def h5ad_from_subset_of_hdf5(
        hdf5_path,
        row_idx,
        obs,
        gene_list,
        dst_path,
        chunk_size=20000):
    """
    Split off a subset of the cell-by-gene HDF5 file as an h5ad file

    Parameters
    ----------
    hdf5_path:
        path to the cell-by-gene hdf5
    row_idx:
        np.array of integers indicating the rows of the CSV to be
        saved to his h5ad file
    obs:
        obs dataframe for this h5ad file
    gene_list:
        list of gene identifiers for this h5ad file
    dst_path:
        path to the h5ad file being written
    chunk_size:
        number of cells to read/write at a time
    """

    n_cells = len(obs)
    n_genes = len(gene_list)

    last_printed = 0
    r0 = 0
    with h5py.File(hdf5_path, 'r') as src:
        with h5py.File(dst_path, 'w') as dst:
            x_dataset = dst.create_dataset(
                'X',
                shape=(n_cells, n_genes),
                dtype=src['data/counts'].dtype,
                chunks=True,
                compression='gzip',
                compression_opts=4
            )
            x_dataset.attrs.create(
                name='encoding-type',
                data='array'
            )
            x_dataset.attrs.create(
                name='encoding-version',
                data='0.2.0'
            )
            for i0 in range(0, len(row_idx), chunk_size):
                idx_chunk = row_idx[i0:i0+chunk_size]
                chunk = src['data/counts'][:, idx_chunk].transpose()
                x_dataset[r0: r0+chunk.shape[0], :] = chunk
                r0 += chunk.shape[0]
                if i0 >= (last_printed + n_cells//10):
                    last_printed = i0
                    print(f'    wrote {i0+chunk_size:.3e} of {n_cells:.3e} cells')

    # write metadata dataframes to h5ad file
    var = pd.DataFrame(
        [{'gene': gene} for gene in gene_list]
    ).set_index('gene')

    anndata_utils.write_df_to_h5ad(
        h5ad_path=dst_path,
        df_name='var',
        df_value=var)
    anndata_utils.write_df_to_h5ad(
        h5ad_path=dst_path,
        df_name='obs',
        df_value=obs
    )
    print(f'WROTE {dst_path}')
    print(f'contains {n_cells} cells and {n_genes} genes')



# In[27]:


mouse_test_path = data_dir / 'mouse_test.h5ad'
mouse_training_path = data_dir / 'mouse_training.h5ad'


# Create mappings from sample_id to the row index in both the metadata and the cell-by-gene datasets

# In[28]:


metadata_sample_to_idx = {s: ii for ii, s in enumerate(mouse_metadata_df.index.values)}
data_sample_to_idx = {s: ii for ii, s in enumerate(sample_name)}


# Make sure that each cell that is in the metadata table is in the cell-by-gene dataset.

# In[ ]:





# In[29]:


assert len(set(metadata_sample_to_idx.keys())-set(data_sample_to_idx.keys())) == 0


# Creatie a mapping from row in the metadata table to row in the cell-by-gene dataset

# In[30]:


data_idx_to_metadata_idx = dict()
for sample in metadata_sample_to_idx:
    data_idx = data_sample_to_idx[sample]
    metadata_idx = metadata_sample_to_idx[sample]
    data_idx_to_metadata_idx[data_idx] = metadata_idx


# Split the data so that 70% of the cells are in the training set and 30% are in the test set.

# In[31]:


rng = np.random.default_rng(554412)
training_fraction = 0.7
sample_name_list = np.sort(np.array(list(metadata_sample_to_idx)))
n_cells = len(sample_name_list)
training_mask = np.zeros(n_cells, dtype=bool)
chosen_idx = rng.choice(np.arange(n_cells), np.round(training_fraction*n_cells).astype(int), replace=False)
training_mask[chosen_idx] = True
training_samples = sample_name_list[training_mask]
test_samples = sample_name_list[np.logical_not(training_mask)]


# For pedagogical purposes, we are constructing our training set using the ENSEMBL IDs for the genes, however leaving the gene symbols in the test set. This is so we can demonstrate MapMyCells' ability to seamlessly translate from gene symbols to ENSEMBL IDs.
# 
# Generally speaking: if you can identify your genes using stable gene identifiers (i.e. ENSEMBL IDs or NCBI IDs), you shoudl do so. It will minimize ambiguity in the relationship between the genes in the reference dataset and the genes in the unlabeled dataset.

# In[32]:


get_ipython().run_cell_magic('time', '', 'test_data_idx = np.sort(np.array([data_sample_to_idx[s] for s in test_samples]))\ntest_metadata_idx = np.array([data_idx_to_metadata_idx[idx] for idx in test_data_idx])\nnp.testing.assert_array_equal(mouse_metadata_df.index.values[test_metadata_idx], np.array(sample_name)[test_data_idx])\nh5ad_from_subset_of_hdf5(\n    hdf5_path=mouse_data_path,\n    dst_path=mouse_test_path,\n    row_idx=test_data_idx,\n    obs=mouse_metadata_df.iloc[test_metadata_idx],\n    gene_list=mouse_gene_list\n)\n')


# In[29]:


get_ipython().run_cell_magic('time', '', 'training_data_idx = np.sort(np.array([data_sample_to_idx[s] for s in training_samples]))\ntraining_metadata_idx = np.array([data_idx_to_metadata_idx[idx] for idx in training_data_idx])\nnp.testing.assert_array_equal(mouse_metadata_df.index.values[training_metadata_idx], sample_name[training_data_idx])\nh5ad_from_subset_of_hdf5(\n    hdf5_path=mouse_data_path,\n    dst_path=mouse_training_path,\n    row_idx=training_data_idx,\n    obs=mouse_metadata_df.iloc[training_metadata_idx],\n    gene_list=mouse_genes_as_ensembl\n)\n')


# ### Create the data files needed to run MapMyCells
# 
# First, create a `precomputed_stats` file that defines the taxonomy and average gene expression profie of every cell type in the taxonomy.

# In[24]:


from cell_type_mapper.cli.precompute_stats_scrattch import PrecomputationScrattchRunner


# In[25]:


mouse_precomputed_stats_path = data_dir / 'mouse_precomputed_stats.h5'


# In[108]:


# note that the taxonomy is made of classes which are subdivided into subclasses which are subdivided into clusters,
# this is how the 'hierarchy' config parameter is set.

config = {
    "h5ad_path": str(mouse_training_path),
    "output_path": str(mouse_precomputed_stats_path),
    "hierarchy": ["class_label", "subclass_label", "cluster_label"],
    "n_processors": 4,
    "normalization": "raw",
    "clobber": True
}


# In[33]:


get_ipython().run_cell_magic('time', '', 'runner = PrecomputationScrattchRunner(\n    args=[],\n    input_data=config\n)\nrunner.run()\n')


# #### Find reference marker genes
# 
# Now create the reference marker file which defines every marker gene between every pair of taxons in the taxonomy. Generally, this step scales in complexity as`N**2` where `N` is the number of leaf nodes (clusters) in the taxonomy. This taxonomy is small enough that it only takes ~ 2 minutes on 4 cores. For taxonomies with a few thousand clusters, the process will take a few hours instead of a few minutes.

# In[111]:


for level in ('class_label', 'subclass_label', 'cluster_label'):
    unq = set(mouse_metadata_df[level].values)
    print(f'{len(unq)} taxons at level "{level}"')


# In[34]:


from cell_type_mapper.cli.reference_markers import ReferenceMarkerRunner


# In[30]:


mouse_reference_marker_dir = data_dir / "reference_markers/mouse"
if not mouse_reference_marker_dir.exists():
    mouse_reference_marker_dir.mkdir(parents=True)
config = {
    "precomputed_path_list": [str(mouse_precomputed_stats_path)],
    "output_dir": str(mouse_reference_marker_dir),
    "n_processors": 4,
    "max_gb": 10,
    "clobber": True
}


# In[36]:


get_ipython().run_cell_magic('time', '', 'runner = ReferenceMarkerRunner(\n    args=[],\n    input_data=config\n)\nrunner.run()\n')


# #### Find query marker genes
# 
# Now take the reference marker set and select a minimal subset of marker genes to use for the actual mapping. This step can be very memory intensive for complex taxonomies (note: resource usage does not scale with the size of the training data; it scales with the size of the taxonomy). This taxonomy will require only ~ 20 GB of memory. A taxonomy with a few thousand clusters will require ~ 48 GB of memory.

# In[31]:


from cell_type_mapper.cli.query_markers import QueryMarkerRunner


# In[32]:


mouse_query_marker_path = data_dir / "mouse_query_markers.json"
mouse_reference_marker_path = mouse_reference_marker_dir / "reference_markers.h5"
assert mouse_reference_marker_path.is_file()


# In[33]:


config = {
    "reference_marker_path_list": [str(mouse_reference_marker_path)],
    "output_path": str(mouse_query_marker_path),
    "n_processors": 4
}


# In[40]:


get_ipython().run_cell_magic('time', '', 'runner = QueryMarkerRunner(\n    args=[],\n    input_data=config\n)\nrunner.run()\n')


# ### Actually perform the mapping
# 
# Now we will map our test dataset onto the taxonomy defined by our training dataset

# In[34]:


from cell_type_mapper.cli.from_specified_markers import FromSpecifiedMarkersRunner


# In[35]:


mouse_test_mapping_path = data_dir / 'mouse_test_mapping.json'


# In[36]:


config = {
    "precomputed_stats": {
        "path": str(mouse_precomputed_stats_path)
    },
    "query_markers": {
        "serialized_lookup": str(mouse_query_marker_path)
    },
    "query_path": str(mouse_test_path),
    "extended_result_path": str(mouse_test_mapping_path),
    "type_assignment": {
        "n_processors": 4,
        "normalization": "raw"
    },
    "verbose_stdout": True,
    "gene_mapping": {
        "db_path": str(gene_mapper_db_path)
    }
}


# In[44]:


get_ipython().run_cell_magic('time', '', 'runner = FromSpecifiedMarkersRunner(\n    args=[],\n    input_data=config\n)\nrunner.run()\n')


# ### Assess mapping quality
# 
# Because the test dataset was already annotated with class, subclaass, and cluster labeles for this taxonomy, we can use these pre-existing annotations as a ground truth to assess the accuracy of the mapping performed by MapMyCells.

# In[37]:


import matplotlib.pyplot as plt
import json
import cell_type_mapper.utils.anndata_utils as anndata_utils


# Create a dict mapping the cell's unique identifier to its "true" cell type annotations.

# In[38]:


mouse_truth = anndata_utils.read_df_from_h5ad(mouse_test_path, df_name='obs').reset_index().to_dict(orient='records')
mouse_truth = {cell['sample_name']: cell for cell in mouse_truth}


# Load the mapping produced by MapMyCells.

# In[39]:


with open(mouse_test_mapping_path, 'rb') as src:
    mapping = json.load(src)


# At each level of the taxonomy, MapMyCells actually maps each cell 100 times, each time using different random subsample of the available marker genes. The assigned cell type is the cell type that is chosen by the plurality of these 100 "bootstrap iterations" (written documentation of this algorithm can be found [on this page](https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/algorithms/hierarchical_mapping.md)).  This allows MapMyCells to return two quality metrics quantifying it confidence in its mappings:
# 
# - `aggregate_probability`: the fraction of bootstrapping iterations that chose the assigned cell type ([see here](https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/output.md#aggregate_probability))
# - `avg_correlation`: the Pearson's correlation coefficient between the cell and the chosen cell type, averaged over the boostrapping iterations that chose the cell type.
# 
# Either of these metricks can be used as a gatekeeper for the quality of mappings the user is willing to accept. Below, we define a function to plot the [F1 score](https://www.geeksforgeeks.org/machine-learning/f1-score-in-machine-learning/) as a function of cuts in the two quality metrics, as well as the true positive versus false positive rates as we vary those cuts.

# In[40]:


def assess_mapping(
        truth,
        mapping,
        taxonomy_level,
        fontsize=15):
    """
    Generate plots assessing mapping quality

    Parameters
    ----------
    truth:
        a dict mapping cell_label to the correct assignment at each level of the taxonomy
    mapping:
        the 'results' from the mapping run
    taxonomy_level:
        the level at which we are mapping
    """
    truth_arr = np.array(
        [truth[cell['cell_id']][taxonomy_level] for cell in mapping]
    )
    mapping_arr = np.array(
        [cell[taxonomy_level]['assignment'] for cell in mapping]
    )
    prob = np.array(
        [cell[taxonomy_level]['aggregate_probability'] for cell in mapping]
    )
    corr = np.array(
        [cell[taxonomy_level]['avg_correlation'] for cell in mapping]
    )

    is_true = (mapping_arr == truth_arr)

    corr_grid = np.linspace(0.0, 1.0, 20)
    prob_grid = np.linspace(0.0, 1.0, 20)

    n_rows = 1
    n_cols = 2
    fig = plt.figure(figsize=(n_cols*5, n_rows*5))
    axis_list = [fig.add_subplot(n_rows, n_cols, ii+1) for ii in range(n_rows*n_cols)]

    axis_list[0].set_xlabel('metric cut', fontsize=fontsize)
    axis_list[0].set_ylabel('F1', fontsize=fontsize)
    axis_list[1].set_xlabel('False positive rate', fontsize=fontsize)
    axis_list[1].set_ylabel('True positive rate', fontsize=fontsize)

    for metric_grid, metric, metric_name in [(corr_grid, corr, 'metric = Correlation'),
                                             (prob_grid, prob, 'metric = Probability')]:

        tp = np.zeros(metric_grid.shape, dtype=int)
        fp = np.zeros(metric_grid.shape, dtype=int)
        fn = np.zeros(metric_grid.shape, dtype=int)
        for ii, metric_value in enumerate(metric_grid):
            considered_true = (metric >= metric_value)
            tp[ii] = np.logical_and(
                is_true,
                considered_true).sum()
            fp[ii] = np.logical_and(
                np.logical_not(is_true),
                considered_true).sum()
            fn[ii] = np.logical_and(
                is_true,
                np.logical_not(considered_true)).sum()
            fn[ii] += np.logical_not(is_true).sum()

        f1 = tp/(tp+0.5*(fp+fn))
        axis_list[0].plot(metric_grid, f1, label=metric_name)
        axis_list[1].plot(fp/len(mapping_arr), tp/len(mapping_arr), label=metric_name)

    axis_list[0].set_title(taxonomy_level, fontsize=fontsize)
    for axis in axis_list:
        axis.legend(loc=0, fontsize=fontsize)
    fig.tight_layout()


# In[41]:


assess_mapping(truth=mouse_truth, mapping=mapping['results'], taxonomy_level='class_label')


# In[42]:


assess_mapping(truth=mouse_truth, mapping=mapping['results'], taxonomy_level='subclass_label')


# In[43]:


assess_mapping(truth=mouse_truth, mapping=mapping['results'], taxonomy_level='cluster_label')


# As you can see, MapMyCells is more than 95% accurate at the class and subclass level and nearly 80% accurate at the cluster level of this taxonomy.
# 
# It is worth reiterating that this is not one of the taxonomies that MapMyCells was designed against.

# ## Human data
# 
# Let us repeat this process with the Human M1 10x data from brain-map.org

# In[44]:


import re
import tempfile


# In[45]:


human_metadata_path = data_dir / "human_metadata.csv"
human_data_path = data_dir / "human_data.csv"


# User our data download helper function to download the cell-by-gene expression matrix and cell annotation metadata associated with the Allen Institutes [2020 human primary motor cortex study](https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x). This will download 2 CSV files, one that is 23 MB and one that is 7.2 GB in size.

# In[46]:


get_ipython().run_cell_magic('time', '', 'download_file(\n    src_url="https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv",\n    dst_path=human_metadata_path\n)\n')


# In[47]:


get_ipython().run_cell_magic('time', '', 'download_file(\n    src_url="https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv",\n    dst_path=human_data_path\n)\n')


# In[48]:


human_metadata_df = pd.read_csv(human_metadata_path)


# In[49]:


human_metadata_df.columns


# In[50]:


human_metadata_df = human_metadata_df[['sample_name', 'class_label', 'subclass_label', 'cluster_label']].set_index('sample_name')


# In this case, the cell-by-gene data is in a CSV file in which each column is a gene and each row is a cell. Let's confirm that the cell IDs in the metadata table match the cell IDs in the cell-by-gene table

# In[51]:


get_ipython().run_cell_magic('time', '', "sample_name_pattern = re.compile('[^,]+')\ndata_sample_names = []\nwith open(human_data_path, 'r') as src:\n    human_gene_list = src.readline().split(',')[1:]\n    for line in src:\n        match = sample_name_pattern.match(line)\n        data_sample_names.append(line[:match.span()[1]])\n")


# In[52]:


np.array_equal(human_metadata_df.index.values, np.array(data_sample_names))


# ### Map the human genes to ENSEMBL IDs

# In[53]:


human_genes_as_ensembl = gene_mapper.map_genes(
    gene_list=human_gene_list,
    dst_species='human',
    dst_authority='ENSEMBL'
)['gene_list']


# ### Split the human data into test and training sets

# In[54]:


def h5ad_from_subset_of_csv(
        csv_path,
        row_idx,
        obs,
        gene_list,
        dst_path):
    """
    Split off a subset of the human cell-by-gene CSV as an h5ad file

    Parameters
    ----------
    csv_path:
        path to the human cell-by-gene CSV
    row_idx:
        np.array of integers indicating the rows of the CSV to be
        saved to his h5ad file
    obs:
        obs dataframe for this h5ad file
    gene_list:
        list of gene identifiers for this h5ad file
    dst_path:
        path to the h5ad file being written
    """

    row_idx_set = set(row_idx)
    n_cells = len(obs)
    n_genes = len(gene_list)

    (descriptor,
     tmp_path) = tempfile.mkstemp(
         suffix='.csv'
     )
    os.close(descriptor)

    tmp_path = pathlib.Path(
        tmp_path
    )

    try:
        # copy chosen rows of CSV to a new CSV at tmp_path
        dtype = dict()
        with open(csv_path, 'r') as src:
            header = src.readline()

            header_params = header.strip().split(',')
            dtype[header_params[0]] = str
            for hh in header_params[1:]:
                dtype[hh] = int

            with open(tmp_path, 'w') as dst:
                dst.write(header)
                for idx, line in enumerate(src):
                    if idx in row_idx_set:
                        dst.write(line)

        # use pd.read_csv with non-zero chunksize to iterate over
        # the tmp CSV one chunk at a time
        chunk_iterator = pd.read_csv(
            tmp_path,
            chunksize=5000,
            dtype=dtype
        )

        with h5py.File(dst_path, 'w') as dst:
            x_dataset = dst.create_dataset(
                'X',
                shape=(n_cells, n_genes),
                dtype=int,
                chunks=True,
                compression='gzip',
                compression_opts=4
            )
            x_dataset.attrs.create(
                name='encoding-type',
                data='array'
            )
            x_dataset.attrs.create(
                name='encoding-version',
                data='0.2.0'
            )
            r0 = 0
            for chunk in chunk_iterator:
                chunk = chunk.set_index('sample_name')
                chunk = chunk.to_numpy()
                x_dataset[r0:r0+chunk.shape[0], :] = chunk
                r0 += chunk.shape[0]

        # write metadata dataframes to h5ad file
        var = pd.DataFrame(
            [{'gene': gene} for gene in gene_list]
        ).set_index('gene')

        anndata_utils.write_df_to_h5ad(
            h5ad_path=dst_path,
            df_name='var',
            df_value=var)
        anndata_utils.write_df_to_h5ad(
            h5ad_path=dst_path,
            df_name='obs',
            df_value=obs
        )

    finally:
        if tmp_path.exists():
            print(f'deleting {tmp_path}')
            tmp_path.unlink()



# In[55]:


rng = np.random.default_rng(3442213)
training_fraction = 0.7
n_cells = len(data_sample_names)
training_mask = np.zeros(n_cells, dtype=bool)
chosen_idx = rng.choice(np.arange(n_cells), np.round(training_fraction*n_cells).astype(int), replace=False)
training_mask[chosen_idx] = True


# In[56]:


human_training_path = data_dir / 'human_training_data.h5ad'
human_test_path = data_dir / 'human_test_data.h5ad'


# In[65]:


get_ipython().run_cell_magic('time', '', 'training_idx = np.where(training_mask)[0]\nh5ad_from_subset_of_csv(\n    csv_path=human_data_path,\n    row_idx=training_idx,\n    obs=human_metadata_df.iloc[training_idx],\n    dst_path=human_training_path,\n    gene_list=human_genes_as_ensembl\n)\n')


# In[66]:


get_ipython().run_cell_magic('time', '', 'test_idx = np.where(np.logical_not(training_mask))[0]\nh5ad_from_subset_of_csv(\n    csv_path=human_data_path,\n    row_idx=test_idx,\n    obs=human_metadata_df.iloc[test_idx],\n    dst_path=human_test_path,\n    gene_list=human_gene_list\n)\n')


# ### Create the data files for MapMyCells

# In[57]:


human_precomputed_stats_path = data_dir / 'human_precomputed_stats.h5'


# In[58]:


config = {
    "h5ad_path": str(human_training_path),
    "output_path": str(human_precomputed_stats_path),
    "hierarchy": ["class_label", "subclass_label", "cluster_label"],
    "n_processors": 4,
    "normalization": "raw",
    "clobber": True
}


# In[69]:


get_ipython().run_cell_magic('time', '', 'runner = PrecomputationScrattchRunner(\n    args=[],\n    input_data=config\n)\nrunner.run()\n')


# #### Reference markers (human)

# In[59]:


human_reference_marker_dir = data_dir / "reference_markers/human"
if not human_reference_marker_dir.exists():
    human_reference_marker_dir.mkdir(parents=True)
config = {
    "precomputed_path_list": [str(human_precomputed_stats_path)],
    "output_dir": str(human_reference_marker_dir),
    "n_processors": 4,
    "max_gb": 10
}


# In[71]:


get_ipython().run_cell_magic('time', '', 'runner = ReferenceMarkerRunner(\n    args=[],\n    input_data=config\n)\nrunner.run()\n')


# #### Query markers (human)

# In[60]:


human_query_marker_path = data_dir / "human_query_markers.json"
human_reference_marker_path = human_reference_marker_dir / "reference_markers.h5"
assert human_reference_marker_path.is_file()


# In[61]:


config = {
    "reference_marker_path_list": [str(human_reference_marker_path)],
    "output_path": str(human_query_marker_path),
    "n_processors": 4
}


# In[74]:


get_ipython().run_cell_magic('time', '', 'runner = QueryMarkerRunner(\n    args=[],\n    input_data=config\n)\nrunner.run()\n')


# ### Map human test set onto human taxonomy

# In[62]:


human_test_mapping_path = data_dir / 'human_test_mapping.json'


# In[63]:


config = {
    "precomputed_stats": {
        "path": str(human_precomputed_stats_path)
    },
    "query_markers": {
        "serialized_lookup": str(human_query_marker_path)
    },
    "query_path": str(human_test_path),
    "extended_result_path": str(human_test_mapping_path),
    "type_assignment": {
        "n_processors": 4,
        "normalization": "raw"
    },
    "verbose_stdout": True,
    "gene_mapping": {
        "db_path": str(gene_mapper_db_path)
    }
}


# In[77]:


get_ipython().run_cell_magic('time', '', 'runner = FromSpecifiedMarkersRunner(\n    args=[],\n    input_data=config\n)\nrunner.run()\n')


# ### Assess quality of human mapping

# In[112]:


human_truth = anndata_utils.read_df_from_h5ad(human_test_path, df_name='obs').reset_index().to_dict(orient='records')
human_truth = {cell['sample_name']: cell for cell in human_truth}


# In[113]:


with open(human_test_mapping_path, 'rb') as src:
    human_mapping = json.load(src)


# In[114]:


assess_mapping(mapping=human_mapping['results'], truth=human_truth, taxonomy_level='class_label')


# In[115]:


assess_mapping(mapping=human_mapping['results'], truth=human_truth, taxonomy_level='subclass_label')


# In[116]:


assess_mapping(mapping=human_mapping['results'], truth=human_truth, taxonomy_level='cluster_label')


# Again we see that MapMyCells is ~ 95% accurate at the class and subclass level and ~ 80% accurate at the cluster level.

# ## Mapping human data onto mouse taxonomy
# 
# Now let us try mapping our human training data (since there are more cells than in the human test data) onto the taxonomy defined by our mouse training data. By inspecting the results and comparing the mapping to the original human annotation, we can hopefully find some meaningful correspondence between the cell types in the human taxonomy and the cell types in the mouse taxonomy.

# In[69]:


cross_species_mapping_path = data_dir / 'human_to_mouse_mapping.json'


# In[70]:


config = {
    "precomputed_stats": {
        "path": str(mouse_precomputed_stats_path)
    },
    "query_markers": {
        "serialized_lookup": str(mouse_query_marker_path)
    },
    "query_path": str(human_training_path),
    "extended_result_path": str(cross_species_mapping_path),
    "type_assignment": {
        "n_processors": 4,
        "normalization": "raw"
    },
    "verbose_stdout": True,
    "gene_mapping": {
        "db_path": str(gene_mapper_db_path)
    }
}


# In[85]:


get_ipython().run_cell_magic('time', '', 'runner = FromSpecifiedMarkersRunner(\n    args=[],\n    input_data=config\n)\nrunner.run()\n')


# Let's look at how many unique mouse subclasses our data was assigned to and compare that to teh number of unique subclasses in the original human taxonomy.

# In[71]:


with open(cross_species_mapping_path, 'rb') as src:
    cross_species_mapping = json.load(src)


# In[72]:


full_human_truth = human_metadata_df.reset_index().to_dict(orient='records')
full_human_truth = {cell['sample_name']: cell for cell in full_human_truth}


# In[73]:


human_subclasses = sorted(set([cell['subclass_label'] for cell in full_human_truth.values()]))
mouse_subclasses = set()
for cell in cross_species_mapping['results']:
    mouse_subclasses.add(cell['subclass_label']['assignment'])
mouse_subclasses = sorted(mouse_subclasses)


# In[74]:


print(len(human_subclasses))


# In[75]:


print(len(mouse_subclasses))


# ### Visualize cross-species mapping
# 
# Below we will define some helper functions to help us plot a confusion matrix comparing the human subclasses to which the data was originally annotated with the mouse subclasses to which the data was mapped by MapMyCells. We will focus on the subclass level of the taxonomy because there are few enough of them to make the visualization intelligible.

# In[76]:


human_assignments = [full_human_truth[cell['cell_id']]['subclass_label'] for cell in cross_species_mapping['results']]
mouse_assignments = [cell['subclass_label']['assignment'] for cell in cross_species_mapping['results']]


# In[77]:


import matplotlib.cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable


# In[78]:


def create_confusion_matrix(
        row_labels,
        col_labels):
    """
    Create a confusion matrix in which the numerical values are the Jaccard index between the corresponding
    cell type labels.

    Parameters
    ----------
    row_labels:
        an (n_cells, ) array listing the labels of each cell in the cell types that will be along the rows of
        the confusion matrix
    col_abels:
        an (n_cells, ) array listing the labels of each cell in the cell types that will be along the columns
        of the confusion matrix

    Returns
    -------
    A dict
        {'array': the (n_rows, n_cols) confusion matrix (jaccard index)
         'ct': the (n_rows, n_cols) count matrix (keeps track of how many cells are in each grid point)
         'rows': the label of each row (in the order they occur in the array)
         'cols': the label of each column (in the order they occur in the array)
         }
    """
    row_values = sorted(set(row_labels))
    col_values = sorted(set(col_labels))
    nrows = len(row_values)
    ncols = len(col_values)
    arr = np.zeros((nrows, ncols), dtype=float)
    ct_arr = np.zeros((nrows, ncols), dtype=float)
    row_to_idx = {r:ii for ii, r in enumerate(row_values)}
    col_to_idx = {c:ii for ii, c in enumerate(col_values)}
    row_idx = np.array([row_to_idx[r] for r in row_labels])
    col_idx = np.array([col_to_idx[c] for c in col_labels])
    for ir in range(nrows):
        r_set = set(np.where(row_idx==ir)[0])
        for ic in range(ncols):
            c_set = set(np.where(col_idx==ic)[0])
            jj = len(r_set.intersection(c_set))/len(r_set.union(c_set))
            arr[ir, ic] = jj
            ct_arr[ir, ic] = len(r_set.intersection(c_set))
    assert ct_arr.sum() == len(row_labels)
    assert len(row_labels) == len(col_labels)
    return {
        "array": arr,
        "ct": ct_arr,
        "rows": row_values,
        "cols": col_values
    }


# In[79]:


def get_confusion_order(input_arr, ct_arr, rng, n_iter=500, inverse_temp=10.0, order_wgt=0.1):
    """
    Re-order the rows and columns of a confusion matrix so that it appears
    as diagnoal as possible. Use the Metropolis-Hastings algorithm to try
    to optimize a cost function.

    Parameters
    ----------
    input_arr:
        the confusion matrix. Probably the Jaccard index between
        two cel types
    ct_arr:
        an array the same shape as the confusion matrix showing how many
        cells are in each (row, col) pair of input_arr
    rng:
        a numpy random number generator
    n_iter:
        number of random steps to take with the Metropolis-Hastings
        algorithm
    inverse_temp:
        1/Temperature to use when deciding whether or not
        to accept a random step
    order_wgt:
        how important is ordering cell types by population.
        Higher values of order_wgt will privilege ordering cell types
        by n_cells over keeping the confusion matrix diagonal.

    Returns
    -------
    a dict
        {'rows': np.array of indexes of rows as they should be ordered
         'cols': np.array of indexes of cols as they should be ordered
         }
    """
    print("OPTIMIZING ORDER OF ROWS/COLUMNS IN CONFUSION MATRIX")
    def cost_fn(arr, dst, ct, order_wgt):

        val = (arr*dst).sum()
        row_sum = ct.sum(axis=1)
        row_sorted_idx = np.argsort(row_sum)[-1::-1]
        col_sum = ct.sum(axis=0)
        col_sorted_idx = np.argsort(col_sum)[-1::-1]

        desired_row = np.arange(arr.shape[0])
        row_term = ((desired_row-row_sorted_idx)**2).sum()
        desired_col = np.arange(arr.shape[1])
        col_term =((desired_col-col_sorted_idx)**2).sum()

        norm = ct.shape[0]+ct.shape[1]
        col_wgt = ct.shape[0]/norm
        row_wgt = ct.shape[1]/norm

        val += order_wgt*(col_wgt*col_term+row_wgt*row_term)

        return val

    baseline = np.copy(input_arr)
    candidate = np.copy(input_arr)

    candidate_ct = np.copy(ct_arr)
    baseline_ct = np.copy(ct_arr)

    grid = np.meshgrid(
        np.arange(input_arr.shape[0]),
        np.arange(input_arr.shape[1]),
        indexing='ij'
    )
    row_grid = grid[0]
    col_grid = grid[1]
    dst_grid = (row_grid-col_grid)**2
    print(f'raw cost {cost_fn(arr=input_arr, dst=dst_grid, ct=ct_arr, order_wgt=order_wgt):.2e}')

    rows = np.arange(input_arr.shape[0], dtype=int)
    cols = np.arange(input_arr.shape[1], dtype=int)

    # just try to put brightest pixels at top
    #row_max = baseline.max(axis=1)
    #assert row_max.shape == (baseline.shape[0], )
    #sorted_dex = np.argsort(row_max)[-1::-1]

    row_sum = ct_arr.sum(axis=1)
    sorted_dex = np.argsort(row_sum)[-1::-1]

    rows = rows[sorted_dex]
    baseline = baseline[rows, :]
    candidate = candidate[rows, :]
    baseline_ct = baseline_ct[rows, :]
    candidate_ct = candidate_ct[rows, :]
    print(f'cost after row shuffler {cost_fn(arr=baseline, dst=dst_grid, ct=baseline_ct, order_wgt=order_wgt):.2e}')

    #col_max = baseline.max(axis=0)
    #assert col_max.shape == (baseline.shape[1], )
    #sorted_dex = np.argsort(col_max)[-1::-1]

    col_sum = ct_arr.sum(axis=0)
    sorted_dex = np.argsort(col_sum)[-1::-1]

    cols = cols[sorted_dex]

    baseline = baseline[:, cols]
    candidate = candidate[:, cols]
    baseline_ct = baseline_ct[:, cols]
    candidate_ct = candidate_ct[:, cols]

    n_rows = input_arr.shape[0]
    n_cols = input_arr.shape[1]

    best_cost = cost_fn(arr=baseline, dst=dst_grid, ct=baseline_ct, order_wgt=order_wgt)
    actual_best_cost = best_cost
    actual_rows = np.copy(rows)
    actual_cols = np.copy(cols)
    print(f'base cost {best_cost}')
    for i_iter in range(n_iter):
        c0 = None
        c1 = None
        r0 = None
        r1 = None
        row_or_col = rng.integers(0, 2)
        if row_or_col == 0:
            to_swap = rng.choice(np.arange(n_rows), 2, replace=False)
            r0 = baseline[to_swap[0], :]
            r1 = baseline[to_swap[1], :]
            candidate[to_swap[1], :] = r0
            candidate[to_swap[0], :] = r1

            r0 = baseline_ct[to_swap[0], :]
            r1 = baseline_ct[to_swap[1], :]
            candidate_ct[to_swap[1], :] = r0
            candidate_ct[to_swap[0], :] = r1

        else:
            to_swap = rng.choice(np.arange(n_cols), 2, replace=False)
            c0 = baseline[:, to_swap[0]]
            c1 = baseline[:, to_swap[1]]
            candidate[:, to_swap[1]] = c0
            candidate[:, to_swap[0]] = c1

            c0 = baseline_ct[:, to_swap[0]]
            c1 = baseline_ct[:, to_swap[1]]
            candidate_ct[:, to_swap[1]] = c0
            candidate_ct[:, to_swap[0]] = c1

        candidate_cost = cost_fn(arr=candidate, dst=dst_grid, ct=candidate_ct, order_wgt=order_wgt)
        accept = False
        if candidate_cost < best_cost:
            accept = True
        else:
            roll = rng.random()
            delta = inverse_temp*(candidate_cost-best_cost)
            if np.exp(-0.5*(delta)) > roll:
                accept = True

        if accept:
            best_cost = candidate_cost
            if row_or_col == 0:
                baseline[to_swap[0], :] = candidate[to_swap[0], :]
                baseline[to_swap[1], :] = candidate[to_swap[1], :]
                baseline_ct[to_swap[0], :] = candidate_ct[to_swap[0], :]
                baseline_ct[to_swap[1], :] = candidate_ct[to_swap[1], :]

                r0 = rows[to_swap[0]]
                r1 = rows[to_swap[1]]
                rows[to_swap[1]] = r0
                rows[to_swap[0]] = r1
            else:
                baseline[:, to_swap[0]] = candidate[:, to_swap[0]]
                baseline[:, to_swap[1]] = candidate[:, to_swap[1]]

                baseline_ct[:, to_swap[0]] = candidate_ct[:, to_swap[0]]
                baseline_ct[:, to_swap[1]] = candidate_ct[:, to_swap[1]]

                c0 = cols[to_swap[0]]
                c1 = cols[to_swap[1]]
                cols[to_swap[1]] = c0
                cols[to_swap[0]] = c1
            if best_cost < actual_best_cost:
                actual_best_cost = best_cost
                actual_rows = np.copy(rows)
                actual_cols = np.copy(cols)
        else:
            if row_or_col == 0:
                candidate[to_swap[0], :] = baseline[to_swap[0], :]
                candidate[to_swap[1], :] = baseline[to_swap[1], :]

                candidate_ct[to_swap[0], :] = baseline_ct[to_swap[0], :]
                candidate_ct[to_swap[1], :] = baseline_ct[to_swap[1], :]
            else:
                candidate[:, to_swap[0]] = baseline[:, to_swap[0]]
                candidate[:, to_swap[1]] = baseline[:, to_swap[1]]

                candidate_ct[:, to_swap[0]] = baseline_ct[:, to_swap[0]]
                candidate_ct[:, to_swap[1]] = baseline_ct[:, to_swap[1]]

        if i_iter % (n_iter//10) == 0:
            print(f'iteration {i_iter} -- best_cost {actual_best_cost:.2e}')

    print(f'best_cost {actual_best_cost:.2e}')
    test = np.copy(input_arr)
    test = test[actual_rows, :]
    test = test[:, actual_cols]
    test_ct = np.copy(ct_arr)
    test_ct = test_ct[actual_rows, :]
    test_ct = test_ct[:, actual_cols]
    validate = cost_fn(arr=test, dst=dst_grid, ct=test_ct, order_wgt=order_wgt)
    print(f'validating {validate} -- {actual_best_cost}')

    return {'rows': actual_rows, 'cols': actual_cols}


# In[117]:


def plot_confusion_matrix(
        row_assignments,
        col_assignments,
        row_axis_label,
        col_axis_label,
        label_elements=True,
        fontsize=15,
        colorbar_title='Jaccard index',
        col_cut=None,
        row_cut=None,
        order_wgt=0.1,
        ordering_iterations=300000):
    """
    Create a confusion matrix plot colored by the Jaccard index
    to show the correspondence between cell types in two taxonomies.

    Parameters
    ----------
    row_assignments:
        a list of str. One for each cell. The assignment of the cells to the
        cell types listed along the rows of the confusion matrix
    col_assignments:
        as above. The assignment of the cells to the cell types
        listed along the columns of the confusion matrix
    row_axis_label:
        label for the row axis (the vertical axis) of the
        confusion matrix
    col_axis_label:
        label for the column axis (the horizontal axis)
        of the confusion matrix
    label_elements:
        a boolean. If True, label the rows and columns of the confusion
        matrix. If False, do not (probably because there are too many of them)
    color_bar_title:
        title to list along the color bar of the confusion matrix
    col_cut:
        an optional int. If non-None, only show the first {col_cut}
        columns of the confusion matrix
    row_cut:
        an optional int. If non-None, only show the first {row_cut}
        rows of the founction matrix
    order_wgt:
        a float. Prioritize how important it is that cell types be in
        sorted order according to how many cells are assigned to those
        types when ordering rows and columns of the cofnusion matrix.
    ordering_iterations:
        how many random steps to take when trying to find the most
        aesthetically pleasing order for rows and columns of the
        confusion matrix
    """
    confusion = create_confusion_matrix(
        row_labels=row_assignments,
        col_labels=col_assignments
    )
    confusion_matrix = confusion['array']
    row_labels = confusion['rows']
    col_labels = confusion['cols']
    swapping = get_confusion_order(
        input_arr=confusion_matrix,
        ct_arr=confusion['ct'],
        rng=np.random.default_rng(8711121),
        n_iter=ordering_iterations,
        inverse_temp=10.0,
        order_wgt=order_wgt)

    row_labels = np.array(row_labels)[swapping['rows']]
    col_labels = np.array(col_labels)[swapping['cols']]

    fig = plt.figure(figsize=(20, 20))
    to_plot = np.copy(confusion_matrix)
    to_plot = to_plot[swapping['rows'], :]
    to_plot = to_plot[:, swapping['cols']]

    if col_cut is not None:
        col_labels = col_labels[:col_cut]
        to_plot = to_plot[:, :col_cut]
    if row_cut is not None:
        row_labels = row_labels[:row_cut]
        to_plot = to_plot[:row_cut, :]

    print(f'to_plot shape {to_plot.shape}')
    heatmap_axis = fig.add_subplot(2,1,1)
    heatmap_axis.set_title(f'{row_axis_label} vs {col_axis_label}', fontsize=fontsize)
    img = heatmap_axis.imshow(np.ma.masked_array(to_plot, to_plot==0.0), cmap=matplotlib.cm.Reds)

    divider = make_axes_locatable(heatmap_axis)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(
        img,
        ax=heatmap_axis,
        cax=cax,
        label=colorbar_title)

    cbar.ax.tick_params(axis='both', which='both', labelsize=fontsize)
    cbar.ax.set_ylabel(ylabel=colorbar_title, fontdict={'fontsize': fontsize})

    heatmap_axis.set_xlabel(col_axis_label, fontsize=fontsize)
    heatmap_axis.set_ylabel(row_axis_label, fontsize=fontsize)

    if label_elements:
        row_ticks = [ii+0.5 for ii in range(len(row_labels))]
        heatmap_axis.set_yticks(row_ticks, labels=row_labels, va='bottom', ha='right')

        col_ticks = [ii-0.5 for ii in range(len(col_labels))]
        heatmap_axis.set_xticks(col_ticks, labels=col_labels, va='top', ha='left')

        heatmap_axis.tick_params(
            which='both',
            axis='both',
            labelsize=fontsize,
            pad=0,
            size=10
        )
        heatmap_axis.tick_params(
            axis='x',
            labelrotation=-90.0)

    row_hist_axis = fig.add_subplot(2, 2, 3)
    col_hist_axis = fig.add_subplot(2, 2, 4)
    for axis, assignments, labels in [(row_hist_axis, row_assignments, row_labels),
                                      (col_hist_axis, col_assignments, col_labels)]:
        l_to_idx = {l: ii+0.5 for ii, l in enumerate(labels)}
        assn_idx = [l_to_idx[a] for a in assignments if a in l_to_idx]
        _ = axis.hist(assn_idx, bins=np.arange(len(labels)+1))
        if label_elements:
            ticks = [ii for ii in range(len(labels))]
            axis.set_xticks(ticks, labels, va='top', ha='left')
            axis.tick_params(
                which='both',
                axis='both',
                labelsize=fontsize
            )
            axis.tick_params(
                axis='x',
                labelrotation=-90.0
            )
    row_hist_axis.set_xlabel(row_axis_label, fontsize=fontsize)
    col_hist_axis.set_xlabel(col_axis_label, fontsize=fontsize)
    row_hist_axis.set_ylabel('N cells', fontsize=fontsize)
    col_ylim = col_hist_axis.get_ylim()
    row_ylim = row_hist_axis.get_ylim()
    if col_ylim[1] > row_ylim[1]:
        row_hist_axis.set_ylim(col_ylim)
    else:
        col_hist_axis.set_ylim(row_ylim)
    fig.tight_layout()


# Let's plot the confusion matrix comparing the Human M1 subclasses to which the data was originally annotated with the mouse cortex subclasses to which MapMyCells mapped the data.

# In[81]:


plot_confusion_matrix(
    row_assignments=human_assignments,
    col_assignments=mouse_assignments,
    row_axis_label='Human subclasses (expert annotated)',
    col_axis_label='Mouse subclasses (MapMyCells)',
    order_wgt=0.05
)


# The salient features of this visualizaiton are
# - There is a way to arrange the data such that the confusion matrix appears (mostly) diagonal, indicating good correspondence between the human and mouse subclass assignments.
# - Naively comparing the names of the subclasses, tehy appear to represent the same biological reality.

# ## Map Human data to Whole Mouse Brain Taxonomy
# 
# Let's perform the same analysis, this time mapping the human M1 data onto the [Yao et al. 2023](https://www.nature.com/articles/s41586-023-06812-z) Whole Mouse Brain taxonomy.
# 
# ### Download the files defining the Whole Mouse Brain Taxonomy

# In[82]:


get_ipython().run_cell_magic('time', '', "\nwmb_precomputed_stats_path = abc_cache.get_data_path(\n    directory='WMB-taxonomy',\n    file_name='precomputed_stats_ABC_revision_230821'\n)\n")


# In[83]:


get_ipython().run_cell_magic('time', '', "\nwmb_query_marker_path = abc_cache.get_data_path(\n    directory='WMB-taxonomy',\n    file_name='mouse_markers_230821'\n)\n")


# ### Perform the mapping

# In[84]:


human_to_wmb_mapping_path = data_dir / 'human_to_wmb_mapping.json'


# In[85]:


config = {
    "precomputed_stats": {
        "path": str(wmb_precomputed_stats_path)
    },
    "query_markers": {
        "serialized_lookup": str(wmb_query_marker_path)
    },
    "drop_level": "CCN20230722_SUPT",
    "query_path": str(human_training_path),
    "extended_result_path": str(human_to_wmb_mapping_path),
    "type_assignment": {
        "n_processors": 4,
        "normalization": "raw"
    },
    "verbose_stdout": True,
    "gene_mapping": {
        "db_path": str(gene_mapper_db_path)
    }
}


# In[434]:


get_ipython().run_cell_magic('time', '', 'runner = FromSpecifiedMarkersRunner(\n    args=[],\n    input_data=config\n)\nrunner.run()\n')


# ### Visualize confusion matrices

# In[86]:


with open(human_to_wmb_mapping_path, 'rb') as src:
    human_to_wmb = json.load(src)


# In[87]:


wmb_assignment_set = set([cell['CCN20230722_SUBC']['assignment'] for cell in human_to_wmb['results']])
print(len(wmb_assignment_set))


# In[88]:


assert [cell['cell_id'] for cell in human_to_wmb['results']] == [cell['cell_id'] for cell in cross_species_mapping['results']]


# In[89]:


import cell_type_mapper.taxonomy.taxonomy_tree as taxonomy_module


# In[90]:


wmb_taxonomy_tree = taxonomy_module.TaxonomyTree(data=human_to_wmb['taxonomy_tree'])


# In[91]:


wmb_assignments = [
    wmb_taxonomy_tree.label_to_name(
        level='CCN20230722_SUBC',
        label=cell['CCN20230722_SUBC']['assignment']
    )
    for cell in human_to_wmb['results']
]


# In[118]:


# Turn off the row/column labels because there are too many Whole Mouse Brain subclasses
# to reasonaby visualize

plot_confusion_matrix(
    row_assignments=human_assignments,
    col_assignments=wmb_assignments,
    row_axis_label='Human subclass (expert annotation)',
    col_axis_label='WMB subclass (MapMyCells)',
    fontsize=15,
    label_elements=False,
    order_wgt=0.05
)


# Because most of the cells appear to cluster on the far left side of the confusion matrix, let's impose a cut, only visualizing the first 30 Whole Mouse Brain subclasses and turn the axis labeling back on.

# In[101]:


plot_confusion_matrix(
    row_assignments=human_assignments,
    col_assignments=wmb_assignments,
    row_axis_label='Human subclass (expert annotation)',
    col_axis_label='WMB subclass (MapMyCells)',
    fontsize=15,
    label_elements=True,
    col_cut=30,
    order_wgt=0.05
)


# Again, we see that the confusion matrix is mostly diagonal, and the names of the corresponding subclasses appear to indicate good agreement between the two taxonomies.

# ### Cut on quality
# 
# We can use the quality metrics returned by MapMyCells to filter out mappings that are suspect. First, let's look at the distribution of the `aggregate_probability` quality metric at the subclass level.

# In[119]:


mouse_subclass_prob = np.array([cell['CCN20230722_SUBC']['aggregate_probability'] for cell in human_to_wmb['results']])


# In[123]:


fontsize = 15
fig = plt.figure(figsize=(5, 5))
axis = fig.add_subplot(1,1,1)
_ = axis.hist(mouse_subclass_prob, bins=50)
axis.set_xlabel('aggregate probability', fontsize=fontsize)
axis.set_ylabel('N cells', fontsize=fontsize)


# In[96]:


np.quantile(mouse_subclass_prob, (0.1, 0.5, 0.75, 0.8))


# Let's say we want to reject cells with `aggregate_probability < 0.25`

# In[97]:


human_assignments_best = []
wmb_assignments_best = []
human_assignments_worst = []
for cell in human_to_wmb['results']:
    hh = full_human_truth[cell['cell_id']]['subclass_label']
    if cell['CCN20230722_SUBC']['aggregate_probability'] < 0.25:
        human_assignments_worst.append(hh)
        continue
    human_assignments_best.append(hh)
    mm = wmb_taxonomy_tree.label_to_name(
        level='CCN20230722_SUBC',
        label=cell['CCN20230722_SUBC']['assignment']
    )
    wmb_assignments_best.append(mm)


# In[98]:


len(set(wmb_assignments_best))


# In[122]:


print('cells lost to cut on probability\n=======')
print(f'{len(human_assignments_worst)} of {len(human_assignments)} cells dropped by cut\n=====')
full_unq, full_ct = np.unique(human_assignments, return_counts=True)
full_lookup = {u:c for u, c in zip(full_unq, full_ct)}
unq, ct = np.unique(human_assignments_worst, return_counts=True)
for u, c in zip(unq, ct):
    print(f'human subclass {u} -- {c} dropped cells of {full_lookup[u]} total')


# As you can see, we only loose 1/50th of our data to this cut. However, this cut reduces the number of Whole Mouse Brain subclasses in the mapping from > 80 to 47. Clearly, we are dropping Whole Mouse Brain subclasses that have only a few cells (probably spuriously) assigned to them.
# 
# Let's plot the confusion matrix for this new, higher quality set of mappings.

# In[100]:


plot_confusion_matrix(
    row_assignments=human_assignments_best,
    col_assignments=wmb_assignments_best,
    row_axis_label='Human subclass (expert annotation)',
    col_axis_label='WMB subclass (MapMyCells)',
    fontsize=15,
    label_elements=True,
    order_wgt=0.05,
    ordering_iterations=1000000
)


# Now the visualizaiton is much more tractable and we can see that, even mapping across species, MapMyCells reveals a good correspondence between taxonomies.

# ## More questions?
# 
# Hopefully this notebook provided you with enough information to run MapMyCells on arbitrary cell type taxonomies. If you have any other questions, do not hesitate to post them to the MapMyCells channel in the [Allen Brain Map community forum,](https://community.brain-map.org/c/how-to/mapmycells/20) and we will respond to you promptly.

# In[ ]:




