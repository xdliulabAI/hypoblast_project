#!/bin/bash

conda activate pyscenic

chmod +x 3.1_arboreto_with_multiprocessing.py
./3.1_arboreto_with_multiprocessing.py \
    exprMat_filtered_D0Fm-P20NR-P20PR-P15.loom \
    cisTarget_databases/hg38/allTFs_hg38.txt \
    --method grnboost2 \
    --output adj.tsv \
    --num_workers 8 \
    --seed 123

Rscript 3.2_generate_regulon.R
