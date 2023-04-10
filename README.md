# Single B cell analysis of BND2 population
Writen by Kristen L. Wells

All scripts necessary to replicate the analysis from the manuscript Identification of an anergic BND cell-derived activate B cell population (BND2) in young-onset type 1 diabetes.

## Data processing steps

### Starting from raw fastq files

If you choose to start from the raw fastq files, you can run these through an existing snakemake pipeline. To run this pipeline, first download the raw fastq files from GEO (link) and then update the config file below. The snakemake pipline is built to work on an lsf server.

Otherwise, you can update the config file with your own structure.


#### To use the snakemake pipeline

1. Download and install miniconda3: For Linux
```{bash}
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh bash Miniconda3-latest-Linux-x86_64.sh
```
2. Install Snakemake:
```{bash}
conda install snakemake -c bioconda -c conda-forge
```

3. Update the config file (config.yaml) 
>* RAW_DATA: Location of the raw data 
>* SAMPLES: the list of samples you want to test. This is the name that will be in the output files. The order must be identical to the order of RNA_SAMPLES, ADT_SAMPLES, and VDJ_SAMPLES
>* AGGR_SAMPLES: Which samples should be aggregated together at the end of the pipeline. These should be samples that were split between multiple 10x runs.
>* RNA_SAMPLES: The name of the RNA samples. Not the full fastq name, but the name that is the same for all samples. Likely ends in "GEX"
>* ADT_SAMPLES: *optional* The name of the samples from CITE-seq or hashtagging. If CITE-seq and hashtagging files are separate, include them both separated by a comma. Put samples in the same order as their RNA counterparts. If CITE-seq or hashtagging were not performed, leave this blank.
>* VDJ_SAMPLES: *optional* The name of the samples from VDJ-seq. Put samples in the same order as their RNA counterparts. If VDJ-sequencing was not peformed, leave this blank.
>* RESULTS: Path to the output directory
>* GENOME: Path to the cellranger genome reference
>* ADT_REF: *optional* Path to the ADT-reference. Should be a comma separated file with the columns described in the 10x tutorial: id, name, read, pattern, sequence, and feature_type. The feature_type will be Antibody Capture. The name will be the name in the final output matrix. Leave this blank if CITE-seq or hashtagging were not performed.
>* VDJ_REF: Path to the cellranger VDJ reference. If VDJ sequencing were not performed, leave this blank.
>* MAX_JOBS: The maximum number of jobs that can be submitted by cell ranger at a time
>* LSF_TEMPLATE: Path to an LSF template. One is included in this git repo.
>* CHEMISTRY: *optional* Arguments to the `--chemstiry` flag in cellranger count. If left blank, chemistry will be `auto`. Only use this if the pipeline failed because of a failure to detect the chemistry. Can be filled in for only some samples.

Note: we sequenced the same capture two times. To use the config file as is for the raw data, create the following directory structure

```
- full_raw_data
|--- round_1
|---|--- JH_GEX_S34_L001_R1_001.fastq.gz # JH_round1_GEX_S34_L001_R1_001.fastq.gz
|---|--- JH_GEX_S34_L001_R2_001.fastq.gz # JH_round1_GEX_S34_L001_R2_001.fastq.gz
|---|--- JH_ADT_S7_L003_R1_001.fastq.gz # JH_round1_ADT_S7_L003_R1_001.fastq.gz
|---|--- JH_ADT_S7_L003_R2_001.fastq.gz # JH_round1_ADT_S7_L003_R2_001.fastq.gz
|--- round_2
|---|--- JH_GEX_S6_L002_R1_001.fastq.gz # JH_round2_GEX_S6_L002_R1_001.fastq.gz
|---|--- JH_GEX_S6_L002_R2_001.fastq.gz # JH_round2_GEX_S6_L002_R2_001.fastq.gz
|---|--- JH_ADT_S5_L002_R1_001.fastq.gz # JH_round2_ADT_S5_L002_R1_001.fastq.gz
|---|--- JH_ADT_S5_L002_R2_001.fastq.gz # JH_round2_ADT_S5_L002_R2_001.fastq.gz
```

4. Update snakecharmer.sh to your specific cluster specs. 
>* change the -q argument to the queue you want to use 

5. submit the job using `bsub < snakecharmer.sh`

6. I highly recommend looking at the csv files that are generated and passed to cell ranger to ensure that the correct fastq files have been detected for each sample.

### Starting from the output of cellranger
If you start from the files output from cellranger either by running the above pipeline or by downloading the files from geo (here) you can run just the r scripts. 

All the scripts rely on the `scAnalysisR` package on [github](https://github.com/CUAnschutzBDC/scAnalysisR)

```R
library(devtools)
install_github("CUAnschutzBDC/scAnalysisR")
```

To directly run the scripts, you must first put the filtered outpts (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`) into a folder called `results/healthy_bcells_all/outs/filtered_feature_bc_matrix`)

Once this is done, you can run the scripts in the `healthy_bcells_all` [directory](https://github.com/CUAnschutzBDC/Smith_210825_b_cell_analysis/tree/main/src/scripts/healthy_bcells_all) in numeric order. And then run the scripts in the `healthy_bcells_all_subset` [directory](https://github.com/CUAnschutzBDC/Smith_210825_b_cell_analysis/tree/main/src/scripts/healthy_bcells_all_subset) also in numeric order. This will produce the final UMAP from the figures.

If you only want the final umap, you can just download the metadata, subset to only cells that are "TRUE" for the column "in_final_obj". The UMAP and cluster data are contained with the metadata.



