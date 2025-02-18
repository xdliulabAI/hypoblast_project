#!/bin/bash

export PATH=cellranger-8.0.1:$PATH

cellranger count --id=E1_3 \
                 --transcriptome=refdata-gex-GRCh38-2024-A \
                 --fastqs=E1_3 \
                 --create-bam=true \
                 --sample=E1_3
