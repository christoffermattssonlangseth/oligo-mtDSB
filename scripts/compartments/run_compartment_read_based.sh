#!/usr/bin/env bash
set -euo pipefail

# Activate your env
source /home/christopher/miniconda3/etc/profile.d/conda.sh
conda activate spatial   # <-- change if you use a different env

# Navigate to project directory
cd /home/christopher/projects/oligo-mtDSB

# Run your Python job
python /home/christopher/projects/oligo-mtDSB/run_compartment_read_based.py \
    --data-base /date/gcb/gcb_CML/oligo-mtDSB \
    --adata data/mtDNA_DSB_5k_clustered_manual_annotation.h5ad \
    --results-subdir results/rbd_runs/${1:-test_run} \
    --checkpoints-subdir results/rbd_checkpoints_mtDSB \
    --sample-key sample_id \
    --cell-id-key cell_id \
    --threads 8 \
    --force-rerun 0 \
    --tag monod
