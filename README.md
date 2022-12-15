---
title: VTAM, Benchmark
author: Emese Meglecz, Aitor González
date: Oct 12, 2022
---

This is the code to reproduce the benchmark analysis of metabarcoding software presented here: TODO

# Quick start

Run

~~~
snakemake -p -c 15 -s snkfl_all.yml --config process_data_dir="${HOME}"/Software/process public_data_dir="${HOME}"/Software/public min_readcount=0 outdir=out/min_readcount_0 --resources db_fish=1 db_bat=1
~~~

You must repeat the previous command with min read counts 10, 20 and 60.
Then you run this command:

~~~
python scripts/plt_sensitivity_precision.py out_min_readcount_0/summary_min_readcount_0/control_counts.tsv out_min_readcount_10/summary/control_counts.tsv out_min_readcount_20/summary/control_counts.tsv out_min_readcount_60/summary/control_counts.tsv
~~~

# Datasets

VTAM was tested on two published datasets of metabarcoding : 

1. The fish dataset (Corse et al., 2017) contained 35 samples of excrement for *Zingel asper*, a fresh-water species, and 5 samples of excrement from *Pomatoschistus microps*, from brackish water. 

2. The bat dataset (Galan et al., 2018) had 357 bat fecal pellet samples from different bat species. 

Both datasets included negative controls (6 for fish and 19 for the bat dataset) and mock samples (2 for the fish and 24 for bat) with known composition. All samples in both studies had three PCR replicates. The fish samples were amplified by two primer sets (markers: MFZR and ZFZR), amplifying the first 158-182 bases of the COI gene (Meusnier et al., 2008; Zeale, et al., 2011), while from the bat samples a 133-bp minibarcode of the COI (Gillet et al., 2015) was amplified. The detailed description of samples and laboratory protocols are found in the original studies (Corse et al., 2017; Galan et al., 2018) and the full raw datasets are available in Dryad (https://datadryad.org/stash/dataset/doi:10.5061/dryad.f40v5, https://datadryad.org/stash/dataset/doi:10.5061/dryad.kv02g).

# Requirements

Ubuntu packages

~~~
sudo apt-get install r-cran-jpeg
sudo apt-get install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
~~~

R packages

~~~
install.packages('jpeg')
install.packages("ragg")
install.packages("pkgdown")
install.packages("devtools")
install.packages('latticeExtra')
install.packages('dplyr')
install.packages("BiocManager")

BiocManager::install("ShortRead")
BiocManager::install("dada2")
devtools::install_github("tobiasgf/lulu")
BiocManager::install("biomformat")
install.packages("remotes")
remotes::install_github("metabaRfactory/metabaR")
~~~

python

~~~
pip install snakemake
pip install vtam
pip install obitools
~~~

We used the following versions of the different software for the benchmarking analyses:
 - VTAM 0.2.0
 - cutadpat 4.1
 - vsearch 2.21.1
 - DADA2 1.12.1
 - LULU 0.1.0
 - OBITools3 3.0.1b19
 - metabaR 1.0.0

# The ONE command

This commands runs everything from the data download to the figures of the alpha and beta diversity.

All script necessary to run the benchmarking are available at <http://github.com/aitgon/vtam_bechmark>.

~~~
snakemake -p -c all -s snkfl_all.yml --config process_data_dir=out_min_readcount_0/processed public_data_dir=out_min_readcount_0/public min_readcount=0 outdir=out_min_readcount_0 --resources db_fish=1 db_bat=1 -n
~~~

The snakefile "snkfl_all.yml" includes reference to partial snakefiles that are organized like this documents

Retrieval of raw data

~~~
include: 'snkfl_dwnld_bat.yml'
include: 'snkfl_dwnld_fish.yml'
~~~

Conversion of raw data for the bennchmark software

~~~
include: 'snkfl_prep_bat.yml'
include: 'snkfl_prep_fish_vtam.yml'
include: 'snkfl_prep_fish_obibar.yml'
include: 'snkfl_prep_fish_dalu.yml'
~~~

Real analysis of the data

~~~
include: 'snkfl_analyse_obibar.yml'

include: 'snkfl_analyse_dalu.yml'

include: 'snkfl_analyse_fish_vtam_pool.yml'
include: 'snkfl_analyse_fish_vtam_zfzr.yml'
include: 'snkfl_analyse_fish_vtam_mfzr.yml'
include: 'snkfl_analyse_bat_vtam.yml'
~~~

**The following sections describe the essential steps of the benchmarking protocol**

# Data download

The fish and bat fastq files were dowloaded from :
https://datadryad.org/stash/dataset/doi:10.5061/dryad.f40v5, https://datadryad.org/stash/dataset/doi:10.5061/dryad.kv02g

# Preparation of fish data

## VTAM
Sequences in the fastq files were merged, demultiplexed and trimmed from tags and primers by VTAM using the merge and sortreads commands with default values.

The output fasta files are the input data for the filter step of VTAM.

## OBIbaR

OBITools3 (Boyer et al. 2016) can use only tags of identical length for demultiplexing, which is not the case for the fish dataset. We reformatted the output of the sortreads of command of VTAM to feed already demultiplexed, trimmed sequences to OBITools3. 

## DALU

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

# Preparation of bat data

## DALU

The fastq file pairs were trimmed by cutadapt (Martin, 2011) using the following parameters:

- discard-untrimmed 
- no-indels
- minimum-length 100
- length 120
- min_overlap: length of the primer
- g ATTCHACDAAYCAYAARGAYATYGG
- G ACTATAAAARAAAYTATDAYAAADGCRTG

Output files are named in *marker-sample_replicate_fw.fastq* and *marker-sample_replicate_rev.fastq* format. **They are the input for DALU pipeline**

## VTAM

Fastq files were already demultiplexed in the original dataset but primers should be trimmed from the reads.

The above created fastq file pairs were merged by the fastq_mergepairs command of vseach  (Rognes et al., 2016) using default parameters to produce demultiplexed, merged fasta files. These files were the **input of the filter command of  VTAM**.

## OBIbaR

The input of vtam (merged fasta files) were formatted for OBITools3 as for the fish dataset.

# Analysis with the VTAM pipeline

## Bat data

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

The *filter* command was run again using optimized parameters, then sequences were assigned to taxa by the *taxassign* command.

## Fish data

The MFZR and ZFZR datasets were analysed seperetell, using the same procedure as for the bat dataset.

The optimal parameters for both markers were the following:

- lfn_sample_replicate_cutoff: 0.003
- pcr_error_var_prop: 0.1
- lfn_variant_cutoff:  0.006
- lfn_read_count_cutoff: 10

The results of the two markers were pooled by the *pool* command, and sequences were assigned to taxa by the *taxassign* command.

# Analysis with the DALU pipeline

## DADA2

First DADA2 (v1.12.1, Callahan et al. 2016) was run on all three datasest (bat, fish-MZFR, fish-ZFZR).

Forward and reverse reads were truncated to 170 nucleotides for MFZR and 150 for ZFZR markers in the fish dataset and to 120 in the bat dataset. The accepted merged read length range was 175-190 for MFZR, 154-169 for ZFZR and 124-163 for the bat dataset using the following scripts.

## LULU

ASV sequences returned by DADA2 were BLASTed against themselves using 80% identity and 80 query sequence coverage, then LULU  (Frøslev et al. 2017) was using default parameteres.

## Replicate pooling

Occurrences were accepted if the ASV is present in at least 2 replicates of the sample. This step is the equivalent of the *FilterMinReplicateNumber* in VTAM.

## FilterIndel,  FilterCodonStop, Read count cutoff

Chimera filtering was already included in DADA2.

FilterIndel,  FilterCodonStop and the LFN_read_count filters of VTAM were reimplemented in a perl script (add_filters.pl).

The analyses were run four times each time using a different arbitrary read count cutoff (0, 10, 20, 50) for the LFN_read_count filter.

ASVs are matched to ASVs in the VTAM database and the same sequence IDs are used.

## Pool MFZR and ZFZR markers

Two markers of the fish dataset were pooled by *pool_markers.pl* perl script using the same algorithm as the pool command of VTAM.

## Taxassign

Taxonomic assignment were done using VTAM taxassign command.

# Analysis with the OBIbaR pipeline

## OBITools3

First, OBITools3 (Boyer et al. 2016) were run on the demultiplexed trimmed fasta files using the following commands to dereplicate reads, filter ASVs with low read counts (stricktly less than 10) and denoise (-r 0.1) the dataset:

- obi import
- obi uniq
- obi annotate (-k COUNT -k MERGED_sample)
- obi grep ("sequence['COUNT']>=10")
- obi clean (-s MERGED_sample -r 0.1)
- obi export (--fasta-output)

The three datasest (bat, fish-MZFR, fish-ZFZR) were run separatelly.

## metabaR

Taxonomic assignment were done using VTAM taxassign command as for the other pipelines and the ouput was adapted to metabaR by *adapt_taxa.pl* perl script.

The output of OBITools and the taxonomic assignment were formatted for metabaR using *make_merabaR_input_from_obi.pl* perl script. All occurrences classed as intermediate status were eliminated from the dataset, even if the variant was head or singleton in some other PCRs.

The datasets denoised by OBITools3 were further filtered using the metabaR R package (Zinger et al. 2020). 

The following filters were run from metabaR:
 - Contaslayer and eliminate pcrs (sample-replicate) with total contaminant relative abundance greater than 10%
 - Eliminate motus (ASVs) not assigned to Eukaryotes
 - Eliminate motus (ASVs) if LTG could not be established or the % identity level was lower than 80%
 - Eliminate pcrs with sequencing depth bellow 500 reads
 - Pcrslayer to eliminate replicates too distant from other replicates of the same sample (similar to *FilterRenkonen* in VTAM)
 - Tagjumpslayer (0.005)

## Replicate pooling

Occurrences were accepted if the ASV is present in at least 2 replicates of the sample. This step is the equivalent of the *FilterMinReplicateNumber* in VTAM.

## FilterIndel,  FilterCodonStop, FilterChimera

FilterIndel,  FilterCodonStop, FilterChimera  and the LFN_read_count filters of VTAM were reimplemented in a perl script (add_filters.pl).

The analyses were run four times each time using a different arbitrary read count cutoff (0, 10, 20, 50) for the LFN_read_count filter.

ASVs are matched to ASVs in the VTAM database and the same sequence IDs are used.

## Pool MFZR and ZFZR markers

Two markers of the fish dataset were pooled by *pool_markers.pl* perl script using the same algorithm as the pool command of VTAM.

## Taxassign

Taxonomic assignment were done using VTAM taxassign command.

# Statistical analysis

The previous sections resulted in these main files


- out/dalu_bat/dalu_bat_final_taxa.tsv
- out/dalu_fish/dalu_fish_final_taxa.tsv
- out/obibar_bat/obibar_bat_final_taxa.tsv
- out/obibar_fish/obibar_fish_final_taxa.tsv
- out/vtam_bat/vtam_bat_final_taxa.tsv
- out/vtam_fish/vtam_fish_final_taxa.tsv

Based on these output files, we calculated the Specificity and Sensitivity using the count_false_occ option of the add_filters.pl script.

Alpha and beta diversiy was calculated using the bxplt_diversity.R script.

# References

Boyer, F., Mercier, C., Bonin, A., Bras, Y. L., Taberlet, P., & Coissac, E. (2016). obitools: A unix-inspired software package for DNA metabarcoding. *Molecular Ecology Resources*, *16*(1), 176–182. https://doi.org/10.1111/1755-0998.12428

Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, *13*(7), 581–583. https://doi.org/10.1038/nmeth.3869

Corse, E., Meglécz, E., Archambaud, G., Ardisson, M., Martin, J.-F., Tougard, C., Chappaz, R., & Dubut, V. (2017). A from-benchtop-to-desktop workflow for validating HTS data and for taxonomic identification in diet metabarcoding studies. *Molecular Ecology Resources*, *17*(6), e146–e159. https://doi.org/10.1111/1755-0998.12703

Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Brunbjerg, A. K., Pietroni, C., & Hansen, A. J. (2017). Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. *Nature Communications*, *8*(1), 1–11. https://doi.org/10.1038/s41467-017-01312-x

Galan, M., Pons, J.-B., Tournayre, O., Pierre, E., Leuchtmann, M., Pontier, D., & Charbonnel, N. (2017). *Data from: Metabarcoding for the parallel identification of several hundred predators and their preys: application to bat species diet analysis* (Version 1, p. 1363569679 bytes) [Data set]. Dryad. https://doi.org/10.5061/DRYAD.KV02G

Gillet, F., Tiouchichine, M.-L., Galan, M., Blanc, F., Némoz, M., Aulagnier, S., & Michaux, J. R. (2015). A new method to identify the endangered Pyrenean desman (Galemys pyrenaicus) and to study its diet, using next generation sequencing from faeces. *Mammalian Biology*, *80*(6), 505–509. https://doi.org/10.1016/j.mambio.2015.08.002

Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.Journal*, *17*(1), 10–12. https://doi.org/10.14806/ej.17.1.200

Meusnier, I., Singer, G. A., Landry, J.-F., Hickey, D. A., Hebert, P. D., & Hajibabaei, M. (2008). A universal DNA mini-barcode for biodiversity analysis. *BMC Genomics*, *9*(1), 214. https://doi.org/10.1186/1471-2164-9-214

Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: A versatile open source tool for metagenomics. *PeerJ*, *4*, e2584. https://doi.org/10.7717/peerj.2584

Zeale, M. R. K., Butlin, R. K., Barker, G. L. A., Lees, D. C., & Jones, G. (2011). Taxon-specific PCR for DNA barcoding arthropod prey in bat faeces. *Molecular Ecology Resources*, *11*(2), 236–244. https://doi.org/10.1111/j.1755-0998.2010.02920.x

Zinger, L., Lionnet, C., Benoiston, A.-S., Donald, J., Mercier, C., & Boyer, F. (2020). metabaR: An R package for the evaluation and improvement of DNA metabarcoding data quality. *BioRxiv*, 2020.08.28.271817. https://doi.org/10.1101/2020.08.28.271817

git clone git@github.com:aitgon/vtam_benchmark.git
cd vtam_benchmark

Prepare output and data directories
mkdir -p out
mkdir -p "${HOME}"/Software/process
mkdir -p "${HOME}"/Software/public

Build singularity container
sudo /home/gonzalez/Software/singularity_3.5.1/bin/singularity build out/vtam_benchmark.sif vtam_benchmark.def

Install and enter snakemake as describe here: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake

Run benchmark
snakemake -p -c 15 -s snkfl_all.yml --config process_data_dir="${HOME}"/Software/process public_data_dir="${HOME}"/Software/public min_readcount=0 outdir=out/min_readcount_0 container=out/vtam_benchmark.sif --resources db_fish=1 db_bat=1 --rerun-incomplete --use-singularity --singularity-args "\-u"
