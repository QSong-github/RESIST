#!/bin/bash

#SBATCH --job-name=resist_test_share
#SBATCH --output=resist_test_share_%j.log

#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=90gb
#SBATCH --nodes=1

## for step 9, it may need larger memory, 80G or more, if it fails, please resubmit with larger memory
## some R packages may to install first, please install them in your home directory

####
module load R/4.5
cd /blue/qsong1/sen.guo/resist_02/resist_test_share

### Run all scripts in order
Rscript /blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/pipeline1_UMAP_replicate.R
Rscript /blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/pipeline2_DEG_replicate.R
Rscript /blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/pipeline3_enrichment_replicate.R
Rscript /blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/pipeline3_enrichment_barchart_side_by_side_replicate.R
Rscript /blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/pipeline4_ITH_replicate.R
Rscript /blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/pipeline5_EMT_replicate.R
Rscript /blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/pipeline7_miRNA_replicate.R
Rscript /blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/pipeline7_miRNA_plot_replicate.R
Rscript /blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/pipeline8_RBP1_enrichment_replicate.R
Rscript /blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/pipeline8_RBP_plot1_bubble_replicate.R
Rscript /blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/pipeline9_drug_enrichment_replicate.R
Rscript /blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/pipeline9_drug_enrichment_plot_for_web_replicate.R



