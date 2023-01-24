---
title: VTAM, Benchmark
author: Emese Meglécz, Aitor González
date: Jan 23, 2023
---

This is the code to reproduce the benchmark analysis of metabarcoding software.

# Quick start

## Install dependencies

Download this repository

~~~
git clone git@github.com:aitgon/vtam_benchmark.git
cd vtam_benchmark
~~~

Prepare output and data directories

~~~
mkdir -p out
mkdir -p "${HOME}"/Software/process
mkdir -p "${HOME}"/Software/public
~~~

Build singularity container

~~~
sudo singularity build out/vtam_benchmark.sif vtam_benchmark.def
~~~

Install and enter snakemake as describe here: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

~~~
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
~~~

## Downloads

Download bat and fish datasets

~~~
mkdir -p ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/62829 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/62831 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/62833 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/62835 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/62837 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/62839 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/62841 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/62843 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/62845 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/62847 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/62849 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/62851 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/64339 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/64340 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/64341 -P ${HOME}/Software/public
wget -c -q -r datadryad.org/stash/downloads/file_stream/64342 -P ${HOME}/Software/public
~~~

Unzip bat datasets

~~~
mkdir -p out/data_bat
unzip -o -j ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/64339 -d  out/data_bat
unzip -o -j ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/64340 -d  out/data_bat
unzip -o -j ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/64341 -d  out/data_bat
unzip -o -j ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/64342 -d  out/data_bat
~~~

Untar fish datasets

~~~
mkdir -p out/data_fish
tar zxvf ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/62829 --to-stdout |gzip -f >out/data_fish/MFZR1_S4_L001_R1_001.fastq.gz
tar zxvf ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/62831 --to-stdout |gzip -f >out/data_fish/MFZR1_S4_L001_R2_001.fastq.gz
tar zxvf ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/62833 --to-stdout |gzip -f >out/data_fish/MFZR2_S5_L001_R1_001.fastq.gz
tar zxvf ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/62835 --to-stdout |gzip -f >out/data_fish/MFZR2_S5_L001_R2_001.fastq.gz
tar zxvf ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/62837 --to-stdout |gzip -f >out/data_fish/MFZR3_S6_L001_R1_001.fastq.gz
tar zxvf ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/62839 --to-stdout |gzip -f >out/data_fish/MFZR3_S6_L001_R2_001.fastq.gz
tar zxvf ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/62841 --to-stdout |gzip -f >out/data_fish/ZFZR1_S1_L001_R1_001.fastq.gz
tar zxvf ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/62843 --to-stdout |gzip -f >out/data_fish/ZFZR1_S1_L001_R2_001.fastq.gz
tar zxvf ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/62845 --to-stdout |gzip -f >out/data_fish/ZFZR2_S2_L001_R1_001.fastq.gz
tar zxvf ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/62847 --to-stdout |gzip -f >out/data_fish/ZFZR2_S2_L001_R2_001.fastq.gz
tar zxvf ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/62849 --to-stdout |gzip -f >out/data_fish/ZFZR3_S3_L001_R1_001.fastq.gz
tar zxvf ${HOME}/Software/public/datadryad.org/stash/downloads/file_stream/62851 --to-stdout |gzip -f >out/data_fish/ZFZR3_S3_L001_R2_001.fastq.gz
~~~

Download shark dataset

~~~
mkdir -p out/data_shark/200309
mkdir -p out/data_shark/200513
~~~

~~~
mkdir -p ${HOME}/Software/public
wget -cr https://osf.io/ymswg/download -P ${HOME}/Software/public
wget -cr https://osf.io/7yf2g/download -P ${HOME}/Software/public
wget -cr https://osf.io/wp6u8/download -P ${HOME}/Software/public
wget -cr https://osf.io/q69de/download -P ${HOME}/Software/public
wget -cr https://osf.io/7pt3x/download -P ${HOME}/Software/public
wget -cr https://osf.io/aqrkb/download -P ${HOME}/Software/public
~~~

Untar shark dataset

~~~
mkdir -p ${HOME}/Software/process/osf.io/ymswg/
tar zxvf ${HOME}/Software/public/osf.io/ymswg/download -C ${HOME}/Software/process/osf.io/ymswg/
mkdir -p ${HOME}/Software/process/osf.io/7yf2g/
tar zxvf ${HOME}/Software/public/osf.io/7yf2g/download -C ${HOME}/Software/process/osf.io/7yf2g/
mkdir -p ${HOME}/Software/process/osf.io/wp6u8/
tar zxvf ${HOME}/Software/public/osf.io/wp6u8/download -C ${HOME}/Software/process/osf.io/wp6u8/
mkdir -p ${HOME}/Software/process/osf.io/q69de/
tar zxvf ${HOME}/Software/public/osf.io/q69de/download -C ${HOME}/Software/process/osf.io/q69de/
mkdir -p ${HOME}/Software/process/osf.io/7pt3x/
tar zxvf ${HOME}/Software/public/osf.io/7pt3x/download -C ${HOME}/Software/process/osf.io/7pt3x/
mkdir -p ${HOME}/Software/process/osf.io/aqrkb/
tar zxvf ${HOME}/Software/public/osf.io/aqrkb/download -C ${HOME}/Software/process/osf.io/aqrkb/
~~~

Pool 2 runs of the shark dataset and gzip files

~~~
mkdir -p out/data_shark
cat ${HOME}/Software/process/osf.io/ymswg/REQPOL-R1_S1_L001_R1_001.fastq ${HOME}/Software/process/osf.io/q69de/REQPOL-R1_S1_L001_R1_001.fastq |gzip -c >out/data_shark/reqpol_R1_fw.fastq.gz
cat ${HOME}/Software/process/osf.io/ymswg/REQPOL-R1_S1_L001_R2_001.fastq ${HOME}/Software/process/osf.io/q69de/REQPOL-R1_S1_L001_R2_001.fastq |gzip -c >out/data_shark/reqpol_R1_rv.fastq.gz
cat ${HOME}/Software/process/osf.io/7yf2g/REQPOL-R2_S2_L001_R1_001.fastq ${HOME}/Software/process/osf.io/7pt3x/REQPOL-R2_S2_L001_R1_001.fastq |gzip -c >out/data_shark/reqpol_R2_fw.fastq.gz
cat ${HOME}/Software/process/osf.io/7yf2g/REQPOL-R2_S2_L001_R2_001.fastq ${HOME}/Software/process/osf.io/7pt3x/REQPOL-R2_S2_L001_R2_001.fastq |gzip -c >out/data_shark/reqpol_R2_rv.fastq.gz
cat ${HOME}/Software/process/osf.io/wp6u8/REQPOL-R3_S3_L001_R1_001.fastq ${HOME}/Software/process/osf.io/aqrkb/REQPOL-R3_S3_L001_R1_001.fastq |gzip -c >out/data_shark/reqpol_R3_fw.fastq.gz
cat ${HOME}/Software/process/osf.io/wp6u8/REQPOL-R3_S3_L001_R2_001.fastq ${HOME}/Software/process/osf.io/aqrkb/REQPOL-R3_S3_L001_R2_001.fastq |gzip -c >out/data_shark/reqpol_R3_rv.fastq.gz
~~~

Download and extract COI blast db and taxonomy file
~~~
mkdir -p out/vtam_db
wget -cr https://osf.io/46f8j/download -P ${HOME}/Software/public
tar zxvf ${HOME}/Software/public/osf.io/46f8j/download -C out/vtam_db
~~~

Download and extract 16S blast db and taxonomy file
~~~
wget -cr https://osf.io/sqkda/download -P ${HOME}/Software/public
tar zxvf ${HOME}/Software/public/osf.io/sqkda/download -C out/vtam_db
~~~


## Prepare data 

VTAM, DALU, OBIbaR - bat, fish, shark

~~~
snakemake -p -c all -s 01snkfl_prep.yml --config data_bat=out/data_bat data_fish=out/data_fish data_shark=out/data_shark process_data_dir="${HOME}"/Software/process public_data_dir="${HOME}"/Software/public min_readcount=0 outdir=out container=out/vtam_benchmark.sif --resources db_bat=1 db_fish=1 db_shark=1 --use-singularity
~~~

Delete necessary files

~~~
rm -r out/dalu_fish/fastq_demultiplexed_trimmed_fw_mfzr
rm -r out/dalu_fish/fastq_demultiplexed_trimmed_fw_zfzr
rm -r out/dalu_fish/fastq_demultiplexed_trimmed_rev_mfzr
rm -r out/dalu_fish/fastq_demultiplexed_trimmed_rev_zfzr
rm -r out/dalu_fish/fastq_demultiplexed_untrimmed_mfzr
rm -r out/dalu_fish/fastq_demultiplexed_untrimmed_zfzr
rm -r out/dalu_shark/fastq_demultiplexed_untrimmed
rm -r out/vtam_fish/merged_mfzr
rm -r out/vtam_fish/merged_zfzr
rm -r out/vtam_shark/merged
~~~


## Analyse 

VTAM - bat, fish, shark

~~~
snakemake -p -c all -s 02snkfl_analysis_vtam.yml --config process_data_dir="{HOME}"/Software/process public_data_dir="{HOME}"/Software/public outdir=out container=out/vtam_benchmark.sif taxonomy_16S=out/vtam_db/rdp_16StrainssetNo18_vtam_dbV4/rdp_taxonomy.tsv blastdbdir_16S=out/vtam_db/rdp_16StrainssetNo18_vtam_dbV4 blastdbname_16S=16S_rdp_db taxonomy=out/vtam_db/COInr_for_vtam_2022_05_06_dbV4/COInr_for_vtam_taxonomy.tsv blastdbdir=out/vtam_db/COInr_for_vtam_2022_05_06_dbV4 blastdbname=COInr_for_vtam --resources db_bat=1 db_fish=1 db_shark=1 --use-singularity
~~~

DALU  - bat, fish, shark

~~~
snakemake -p -c all -s 03snkfl_analysis_dalu.yml --config process_data_dir="{HOME}"/Software/process public_data_dir="{HOME}"/Software/public min_readcount=0 outdir=out container=out/vtam_benchmark.sif taxonomy_16S=out/vtam_db/rdp_16StrainssetNo18_vtam_dbV4/rdp_taxonomy.tsv blastdbdir_16S=out/vtam_db/rdp_16StrainssetNo18_vtam_dbV4 blastdbname_16S=16S_rdp_db taxonomy=out/vtam_db/COInr_for_vtam_2022_05_06_dbV4/COInr_for_vtam_taxonomy.tsv blastdbdir=out/vtam_db/COInr_for_vtam_2022_05_06_dbV4 blastdbname=COInr_for_vtam --resources db_bat=1 db_fish=1 db_shark=1 --use-singularity
~~~


OBIbaR  - bat, fish, shark

~~~
snakemake -p -c all -s 04snkfl_analysis_obibar.yml --config data_bat=out/data_bat data_fish=out/data_fish data_shark=out/data_shark process_data_dir="{HOME}"/Software/process public_data_dir="{HOME}"/Software/public min_readcount=0 outdir=out container=out/vtam_benchmark.sif taxonomy_16S=out/vtam_db/rdp_16StrainssetNo18_vtam_dbV4/rdp_taxonomy.tsv blastdbdir_16S=out/vtam_db/rdp_16StrainssetNo18_vtam_dbV4 blastdbname_16S=16S_rdp_db taxonomy=out/vtam_db/COInr_for_vtam_2022_05_06_dbV4/COInr_for_vtam_taxonomy.tsv blastdbdir=out/vtam_db/COInr_for_vtam_2022_05_06_dbV4 blastdbname=COInr_for_vtam --resources db_bat=1 db_fish=1 db_shark=1 --use-singularity
~~~

Change min_read_count and make plots

min_readcount=0

~~~
snakemake -p -c all -s 05snkfl_min_readcount.yml --config data_bat=out/data_bat data_fish=out/data_fish data_shark=out/data_shark process_data_dir="{HOME}"/Software/process public_data_dir="{HOME}"/Software/public min_readcount=0 outdir=out container=out/vtam_benchmark.sif taxonomy_16S=out/vtam_db/rdp_16StrainssetNo18_vtam_dbV4/rdp_taxonomy.tsv blastdbdir_16S=out/vtam_db/rdp_16StrainssetNo18_vtam_dbV4 blastdbname_16S=16S_rdp_db taxonomy=out/vtam_db/COInr_for_vtam_2022_05_06_dbV4/COInr_for_vtam_taxonomy.tsv blastdbdir=out/vtam_db/COInr_for_vtam_2022_05_06_dbV4 blastdbname=COInr_for_vtam --resources db_bat=1 db_fish=1 db_shark=1 --use-singularity
~~~

**Note** : It is possible to have an error message similar to this:

~~~
...
Job xxx  completed successfully, but some output files are missing. Missing files after 5 seconds. This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait:
...
~~~

In this case, just rerun the same command without deleting the files created by the first execution.

min_readcount=10

~~~
snakemake -p -c all -s 05snkfl_min_readcount.yml --config data_bat=out/data_bat data_fish=out/data_fish data_shark=out/data_shark process_data_dir="{HOME}"/Software/process public_data_dir="{HOME}"/Software/public min_readcount=10 outdir=out container=out/vtam_benchmark.sif taxonomy_16S=out/vtam_db/rdp_16StrainssetNo18_vtam_dbV4/rdp_taxonomy.tsv blastdbdir_16S=out/vtam_db/rdp_16StrainssetNo18_vtam_dbV4 blastdbname_16S=16S_rdp_db taxonomy=out/vtam_db/COInr_for_vtam_2022_05_06_dbV4/COInr_for_vtam_taxonomy.tsv blastdbdir=out/vtam_db/COInr_for_vtam_2022_05_06_dbV4 blastdbname=COInr_for_vtam --resources db_bat=1 db_fish=1 db_shark=1 --use-singularity
~~~

min_readcount=40

~~~
snakemake -p -c all -s 05snkfl_min_readcount.yml --config data_bat=out/data_bat data_fish=out/data_fish data_shark=out/data_shark process_data_dir="{HOME}"/Software/process public_data_dir="{HOME}"/Software/public min_readcount=40 outdir=out container=out/vtam_benchmark.sif taxonomy_16S=out/vtam_db/rdp_16StrainssetNo18_vtam_dbV4/rdp_taxonomy.tsv blastdbdir_16S=out/vtam_db/rdp_16StrainssetNo18_vtam_dbV4 blastdbname_16S=16S_rdp_db taxonomy=out/vtam_db/COInr_for_vtam_2022_05_06_dbV4/COInr_for_vtam_taxonomy.tsv blastdbdir=out/vtam_db/COInr_for_vtam_2022_05_06_dbV4 blastdbname=COInr_for_vtam --resources db_bat=1 db_fish=1 db_shark=1 --use-singularity
~~~

min_readcount=60

~~~
snakemake -p -c all -s 05snkfl_min_readcount.yml --config data_bat=out/data_bat data_fish=out/data_fish data_shark=out/data_shark process_data_dir="{HOME}"/Software/process public_data_dir="{HOME}"/Software/public min_readcount=60 outdir=out container=out/vtam_benchmark.sif taxonomy_16S=out/vtam_db/rdp_16StrainssetNo18_vtam_dbV4/rdp_taxonomy.tsv blastdbdir_16S=out/vtam_db/rdp_16StrainssetNo18_vtam_dbV4 blastdbname_16S=16S_rdp_db taxonomy=out/vtam_db/COInr_for_vtam_2022_05_06_dbV4/COInr_for_vtam_taxonomy.tsv blastdbdir=out/vtam_db/COInr_for_vtam_2022_05_06_dbV4 blastdbname=COInr_for_vtam --resources db_bat=1 db_fish=1 db_shark=1 --use-singularity
~~~

Plot sensitivity and precision

~~~
export MINREAD0_TSV_PATH=out/summary_min_readcount_0/control_counts.tsv
export MINREAD10_TSV_PATH=out/summary_min_readcount_10/control_counts.tsv
export MINREAD40_TSV_PATH=out/summary_min_readcount_40/control_counts.tsv
export MINREAD60_TSV_PATH=out/summary_min_readcount_60/control_counts.tsv
export OUTDIR=out

singularity exec -u out/vtam_benchmark.sif python scripts/plt_sensitivity_precision.py ${MINREAD0_TSV_PATH} ${MINREAD10_TSV_PATH} ${MINREAD40_TSV_PATH} ${MINREAD60_TSV_PATH} ${OUTDIR}
~~~

## Run time and memory usage

Run time and memory useage on a desktop computer with 32Gb RAM, 8 CPU, Intel(TM) Core(TM) i7-7700 CPU @ 3.60GHz

| Pipeline                    | Task                                                         | User (s) | System (s) | Elapsed (hh:mm:ss) | CPU % ((U+S)/E) | maxresident memory (M) |
| --------------------------- | ------------------------------------------------------------ | -------- | ---------- | ------------------ | --------------- | ---------------------- |
| 01snkfl_prep.yml            | Preparation of trimmed, demultiplexed fasta or fastq files   | 18969    | 343        | 01:25:47           | 375%            | 1309                   |
| 02snkfl_analysis_vtam.yml   | VTAM (filter and taxassign)                                  | 1805     | 81         | 00:11:01           | 285%            | 2815                   |
| 03snkfl_analysis_dalu.yml   | DALU (DADA2, LULU, supplementary filtres and taxasign)       | 9666     | 178        | 00:57:53           | 283%            | 1808                   |
| 04snkfl_analysis_obibar.yml | OBIbaR (OBITools3, metabaR, supplementary filtres and taxasign) | 12172    | 208        | 01:01:43           | 334%            | 2689                   |

# Detailed analyses
## Datasets

VTAM was tested on three published datasets of metabarcoding : 

1. The fish dataset (Corse et al., 2017) contained 35 samples of excrement for *Zingel asper*, a fresh-water species, and 5 samples of excrement from *Pomatoschistus microps*, from brackish water. 
2. The bat dataset (Galan et al., 2018) had 357 bat fecal pellet samples from different bat species. 
2. The Shark dataset (Esposito et al., 2022) had 63 shark gut content samples from the Black-Tip Reef Shark (*Carcharhinus melanopterus*).

All three datasets included negative controls (6 for fish, 19 for the bat and 7 for shark datastes) and mock samples (2 for fish, 24 for bat and 2 for shark) with known composition. All samples in all three studies had three PCR replicates. 

The fish samples were amplified by two primer sets (markers: MFZR and ZFZR), amplifying the first 158-182 bases of the COI gene (Meusnier et al., 2008; Zeale, et al., 2011), while from the bat samples a 133-bp minibarcode of the COI (Gillet et al., 2015) was amplified. 

The shark samples wee amplified by the 515FY (Parada et al., 2016) and 806RB (Apprill et al., 2015) primers, which target a 253 bp fragment of the 16S ribosomal (rDNA).

The detailed description of samples and laboratory protocols are found in the original studies (Corse et al., 2017; Galan et al., 2018, Esposito et al., 2022) and the full raw datasets are available in Dryad (https://datadryad.org/stash/dataset/doi:10.5061/dryad.f40v5, https://datadryad.org/stash/dataset/doi:10.5061/dryad.kv02g) and OSF (https://osf.io/3txv7/).


## Preparation of fish data

**VTAM**

Sequences in the fastq files were merged, demultiplexed and trimmed from tags and primers by VTAM using the merge and sortreads commands with default values.

The output fasta files are the input data for the filter step of VTAM.

**DALU**

DADA2 (Callahan et al., 2016) requires demultiplexed and trimmed fastq files pairs with forward and reverse reads in separate files. 
We used the output of sortreads of VTAM to recreate fastq files pairs. The fastq file pairs were then trimmed by cutadapt (Martin, 2011). Fastq reads are not oriented in the original dataset, thus we checked both orientations. 

Cutadapt was run with the following parameters:

- discard-untrimmed 
- no-indels
- minimum-length 160 for MFZR, 140 for ZFZR
- length 170 for MFZR,  150 for ZFZR
- min_overlap: length of the primer
- g TCCACTAATCACAARGATATTGGTAC, AGATATTGGAACWTTATATTTTATTTTTGG for MFZR and ZFZR, respecivelly
- G WACTAATCAATTWCCAAATCCTCC for both markers

**OBIbaR**

OBITools3 (Boyer et al. 2016) can use only tags of identical length for demultiplexing, which is not the case for the fish dataset. We reformatted the output of the sortreads of command of VTAM to feed already demultiplexed, trimmed sequences to OBITools3. 

## Preparation of bat data

**DALU**

The fastq file pairs were trimmed by cutadapt (Martin, 2011) using the following parameters:

- discard-untrimmed 
- no-indels
- minimum-length 100
- length 120
- min_overlap: length of the primer
- g ATTCHACDAAYCAYAARGAYATYGG
- G ACTATAAAARAAAYTATDAYAAADGCRTG

Output files are named in *marker-sample_replicate_fw.fastq* and *marker-sample_replicate_rev.fastq* format. They are the input for DALU pipeline.

**VTAM**

Fastq files were already demultiplexed in the original dataset but primers should be trimmed from the reads.

The above created fastq file pairs were merged by the fastq_mergepairs command of vseach  (Rognes et al., 2016) using default parameters to produce demultiplexed, merged fasta files. These files were the input of the filter command of  VTAM.

**OBIbaR**

The input of vtam (merged fasta files) were formatted for OBITools3 as for the fish dataset.

## Preparation of shark data

**VTAM**

Sequences in the fastq files were merged, demultiplexed and trimmed from tags and primers by VTAM using the merge and sortreads commands with default values.

The output fasta files are the input data for the filter step of VTAM.

**DALU**

DADA2 (Callahan et al., 2016) requires demultiplexed and trimmed fastq files pairs with forward and reverse reads in separate files. 
We used the output of sortreads of VTAM to recreate fastq files pairs. The fastq file pairs were then trimmed by cutadapt (Martin, 2011). 

Cutadapt was run with the following parameters:

- discard-untrimmed 
- no-indels
- minimum-length 200
- length 210
- min_overlap: length of the primer
- g GTGYCAGCMGCCGCGGTAA
- G GGACTACNVGGGTWTCTAAT

**OBIbaR**

OBITools3 (Boyer et al. 2016) can use only tags of identical length for demultiplexing, which is not the case for the shark dataset. We reformatted the output of the sortreads of command of VTAM to feed already demultiplexed, trimmed sequences to OBITools3. 


## Analysis with the VTAM pipeline

### Bat data

**Filter by default**

First bat data were filtered using the *filter* command of VTAM with default parameters.

**Optimize paramaters**

The *known_occurrences.tsv* file was prepared from the output of the filter command using *make_known_occurrences command* of VTAM. This file conyains true positives in mock samples false positives in mock and negative controls to be used in the optimize step.

The optimization of  *lfn_sample_replicate_cutoff* and *pcr_error_var_prop* parametres is independent of the other parameters. Optimal paramaters for them were established first so they could be taken into account when optimizing the other parameters. Then the  *lfn_variant_cutoff* and *lfn_read_count_cutoff* were opimized.

The optimal parameters for the bat dataset were the following:

- lfn_sample_replicate_cutoff: 0.002
- pcr_error_var_prop: 0.2
- lfn_variant_cutoff:  0.006
- lfn_read_count_cutoff: 60

**Filter with optimized paramaters**

The *filter* command was run again using optimized parameters, then sequences were assigned to taxa by the *taxassign* command using the COInr COI reference database (Meglécz, 2022; Meglécz, 2023).

**Note:** The complete pipline with Filter (default) - Optimize - Filter (optimized) was run command by command. The 02snkfl_analysis_vtam.yml is intended to produce the output files rapidly, therefore, we did not inlcude the Filter (default) and optimize steps, only the Filter with optimized paramaters and taxassign.

### Fish data

The MFZR and ZFZR datasets were analyzed separately, using the same procedure as for the bat dataset.

The optimal parameters for both markers were the following:

- lfn_sample_replicate_cutoff: 0.003
- pcr_error_var_prop: 0.1
- lfn_variant_replicate_cutoff:  0.006
- lfn_read_count_cutoff: 10

The two markers were pooled by the *pool* command and assigned to taxa by the *taxassign* command using the COInr COI reference database (Meglécz, 2022; Meglécz, 2023).

**Note:** 
As for the bat dataset, only the Filter with optimized parameters and the taxassign steps are included in the 02snkfl_analysis_vtam.yml pipeline.

### Shark data

The shark data was analysed using the same procedure as the bat and fish datasets with the exception of deactivating the FilterIndel and FilterCodonStop filters, since they can be applied only for coding sequences. 

The optimal parameters for both markers were the following:

- lfn_sample_replicate_cutoff: 0.001
- pcr_error_var_prop: 0.2
- lfn_variant_replicate_cutoff:  0.016
- lfn_read_count_cutoff: 40
- skip_filter_codon_stop: 1
- skip_filter_indel: 1

Sequences were assigned to taxa by the *taxassign* command using the 16S referece dataset of RDP Classifier (RDPClassifier_16S_trainsetNo18) formatted for VTAM (https://osf.io/ns2w8/ )

**Note:** 
As for the bat dataset, only the Filter with optimized parameters and the taxassign steps are included in the 02snkfl_analysis_vtam.yml pipeline.

## Analysis with the DALU pipeline

### DADA2

First DADA2 (v1.12.1, Callahan et al. 2016) was run on all four datasets (bat, fish-MZFR, fish-ZFZR, shark).

Forward and reverse reads were truncated to 170 nucleotides for MFZR and 150 for ZFZR markers in the fish dataset, to 120 in the bat dataset and to 150 in the shark datasets. The accepted merged read length range was 175-190 for MFZR, 154-169 for ZFZR, 124-163 for the bat dataset and 252-254 for the shark dataset.

### LULU

ASV sequences returned by DADA2 were BLASTed against themselves using 80% identity and 80 query sequence coverage, then LULU  (Frøslev et al. 2017) was run using default parameters.

### Replicate pooling

Occurrences were accepted if the ASV is present in at least 2 replicates of the sample. This step is the equivalent of the *FilterMinReplicateNumber* in VTAM.

### FilterIndel,  FilterCodonStop, Read count cutoff

Chimera filtering was already included in DADA2.

FilterIndel,  FilterCodonStop and the LFN_read_count filters of VTAM were reimplemented in a perl script (add_filters.pl). For the shark dataset,  FilterIndel,  FilterCodonStop was not used, since it is a non-coding marker.

This step was run four times each time using a different arbitrary read count cutoff (0, 10, 40, 60) for the LFN_read_count filter.

ASVs are matched to ASVs in the VTAM database and the same sequence IDs are used.

### Pool MFZR and ZFZR markers

Two markers of the fish dataset were pooled by *pool_markers.pl* perl script using the same algorithm as the pool command of VTAM.

### Taxassign

Taxonomic assignment were done using VTAM taxassign command.

## Analysis with the OBIbaR pipeline

### OBITools3

First, OBITools3 (Boyer et al. 2016) was run on the demultiplexed trimmed fasta files using the following commands to dereplicate reads, filter ASVs with low read counts (stricktly less than 10) and denoise (-r 0.1) the dataset:

- obi import
- obi uniq
- obi annotate (-k COUNT -k MERGED_sample)
- obi grep ("sequence['COUNT']>=10")
- obi clean (-s MERGED_sample -r 0.1)
- obi export (--fasta-output)

The three datasest (bat, fish-MZFR, fish-ZFZR) were run separatelly.

### metabaR

Taxonomic assignment were done using VTAM taxassign command as for the other pipelines and the ouput was adapted to metabaR by *adapt_taxa.pl* perl script.

The output of OBITools and the taxonomic assignment were formatted for metabaR using *make_merabaR_input_from_obi.pl* perl script. All occurrences classed as intermediate status were eliminated from the dataset, even if the variant was head or singleton in some other PCRs.

The datasets denoised by OBITools3 were further filtered using the metabaR R package (Zinger et al. 2020). 

The following filters were run from metabaR:
 - Contaslayer and eliminate pcrs (sample-replicate) with total contaminant relative abundance greater than 10%
 - Eliminate motus (ASVs) not assigned to Eukaryotes (for the bat and fish datasets), not assigned to Bacteria or Archaea for the shark dataset
 - Eliminate motus (ASVs) if LTG could not be established or the % identity level was lower than 80%
 - Eliminate pcrs with sequencing depth bellow 500 reads for the bat and fish datasets, and 600 for the shark dataset.
 - Pcrslayer to eliminate replicates too distant from other replicates of the same sample (similar to *FilterRenkonen* in VTAM)
 - Tagjumpslayer (0.005 for bat and fish, 0.01 for shark)

### Replicate pooling

Occurrences were accepted if the ASV is present in at least 2 replicates of the sample. This step is the equivalent of the *FilterMinReplicateNumber* in VTAM.

### FilterIndel,  FilterCodonStop, FilterChimera

FilterIndel,  FilterCodonStop, FilterChimera  and the LFN_read_count filters of VTAM were reimplemented in a perl script (add_filters.pl). The FilterIndel,  FilterCodonStop was not activated for the shark dataset.

The analyses were run four times each time using a different arbitrary read count cutoff (0, 10, 40, 60) for the LFN_read_count filter.

ASVs are matched to ASVs in the VTAM database and the same sequence IDs are used.

### Pool MFZR and ZFZR markers

Two markers of the fish dataset were pooled by *pool_markers.pl* perl script using the same algorithm as the pool command of VTAM.

### Taxassign

Taxonomic assignment were done using VTAM taxassign command.

## Statistical analysis

The previous sections resulted in these main files

- out/summary_min_readcout_XX/dalu_bat_final_taxa.tsv
- out/summary_min_readcout_XX/dalu_fish_final_taxa.tsv
- out/summary_min_readcout_XX/dalu_shark_final_taxa.tsv
- out/summary_min_readcout_XX/obibar_bat_final_taxa.tsv
- out/summary_min_readcout_XX/obibar_fish_final_taxa.tsv
- out/summary_min_readcout_XX/obibar_shark_final_taxa.tsv
- out/summary_min_readcout_XX/vtam_bat_final_taxa.tsv
- out/summary_min_readcout_XX/vtam_fish_final_taxa.tsv
- out/summary_min_readcout_XX/vtam_shark_final_taxa.tsv

Based on these output files, we calculated the Specificity and Sensitivity using the count_false_occ option of the add_filters.pl script.

Alpha and beta diversiy was calculated using the bxplt_diversity_intraspecies.R script.

# References

Apprill, A., McNally, S., Parsons, R., & Weber, L. (2015). Minor revision to V4 region SSU rRNA 806R gene primer greatly increases detection of SAR11 bacterioplankton. *Aquatic Microbial Ecology*, *75*(2), 129–137. https://doi.org/10.3354/ame01753

Boyer, F., Mercier, C., Bonin, A., Bras, Y. L., Taberlet, P., & Coissac, E. (2016). obitools: A unix-inspired software package for DNA metabarcoding. *Molecular Ecology Resources*, *16*(1), 176–182. https://doi.org/10.1111/1755-0998.12428

Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, *13*(7), 581–583. https://doi.org/10.1038/nmeth.3869

Corse, E., Meglécz, E., Archambaud, G., Ardisson, M., Martin, J.-F., Tougard, C., Chappaz, R., & Dubut, V. (2017). A from-benchtop-to-desktop workflow for validating HTS data and for taxonomic identification in diet metabarcoding studies. *Molecular Ecology Resources*, *17*(6), e146–e159. https://doi.org/10.1111/1755-0998.12703

Esposito, A., Sasal, P., Clua, É., Meglécz, E., & Clerissi, C. (2022). Shark Provisioning Influences the Gut Microbiota of the Black-Tip Reef Shark in French Polynesia. Fishes, 7(6), Article 6. https://doi.org/10.3390/fishes7060312

Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Brunbjerg, A. K., Pietroni, C., & Hansen, A. J. (2017). Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. *Nature Communications*, *8*(1), 1–11. https://doi.org/10.1038/s41467-017-01312-x

Galan, M., Pons, J.-B., Tournayre, O., Pierre, E., Leuchtmann, M., Pontier, D., & Charbonnel, N. (2017). *Data from: Metabarcoding for the parallel identification of several hundred predators and their preys: application to bat species diet analysis* (Version 1, p. 1363569679 bytes) [Data set]. Dryad. https://doi.org/10.5061/DRYAD.KV02G

Gillet, F., Tiouchichine, M.-L., Galan, M., Blanc, F., Némoz, M., Aulagnier, S., & Michaux, J. R. (2015). A new method to identify the endangered Pyrenean desman (Galemys pyrenaicus) and to study its diet, using next generation sequencing from faeces. *Mammalian Biology*, *80*(6), 505–509. https://doi.org/10.1016/j.mambio.2015.08.002

Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.Journal*, *17*(1), 10–12. https://doi.org/10.14806/ej.17.1.200

Meglécz, E. (2023). COInr and mkCOInr: Building and customizing a non-redundant barcoding reference database from BOLD and NCBI using a semi-automated pipeline. Molecular Ecology Resources, n/a(n/a). https://doi.org/10.1111/1755-0998.13756

Meglécz, E. (2022). COInr a comprehensive, non-redundant COI database from NCBI-nt and BOLD [Data set]. Zenodo. https://doi.org/10.5281/zenodo.6555985

Meusnier, I., Singer, G. A., Landry, J.-F., Hickey, D. A., Hebert, P. D., & Hajibabaei, M. (2008). A universal DNA mini-barcode for biodiversity analysis. *BMC Genomics*, *9*(1), 214. https://doi.org/10.1186/1471-2164-9-214

Parada, A. E., Needham, D. M., & Fuhrman, J. A. (2016). Every base matters: Assessing small subunit rRNA primers for marine microbiomes with mock communities, time series and global field samples. *Environmental Microbiology*, *18*(5), 1403–1414. https://doi.org/10.1111/1462-2920.13023

Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: A versatile open source tool for metagenomics. *PeerJ*, *4*, e2584. https://doi.org/10.7717/peerj.2584

Zeale, M. R. K., Butlin, R. K., Barker, G. L. A., Lees, D. C., & Jones, G. (2011). Taxon-specific PCR for DNA barcoding arthropod prey in bat faeces. *Molecular Ecology Resources*, *11*(2), 236–244. https://doi.org/10.1111/j.1755-0998.2010.02920.x

Zinger, L., Lionnet, C., Benoiston, A.-S., Donald, J., Mercier, C., & Boyer, F. (2020). metabaR: An R package for the evaluation and improvement of DNA metabarcoding data quality. *BioRxiv*, 2020.08.28.271817. https://doi.org/10.1101/2020.08.28.271817
