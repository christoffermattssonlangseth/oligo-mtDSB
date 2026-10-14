"""Helper modules for the oligo-mtDSB analysis (maintained).

    sample_classifier  - leave-one-section-out age / condition classifiers (pseudobulk and single cell)
    ol_response        - OL responder calls (ISR score), classifier gene weights, within-section driver genes
    region_time_de     - treatment x region x time-point DE on pseudobulk (pyDESeq2)
    region_time_plots  - figures built from the region_time_de tables
    spatial_utils      - compact spatial plotting (plot_spatial_compact_fast etc.)

The auto-harvested legacy helpers (de_analysis_utils, io_utils, misc_utils,
plotting_utils, qc_utils) live in legacy/utils/.
"""
