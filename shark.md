
## download data

~~~
mkdir -p out/data_shark
~~~

Download

~~~
mkdir -p ${HOME}/Software/public
wget -cr https://osf.io/ymswg/download -P ${HOME}/Software/public
wget -cr https://osf.io/7yf2g/download -P ${HOME}/Software/public
wget -cr https://osf.io/wp6u8/download -P ${HOME}/Software/public
wget -cr https://osf.io/q69de/download -P ${HOME}/Software/public
wget -cr https://osf.io/7pt3x/download -P ${HOME}/Software/public
wget -cr https://osf.io/aqrkb/download -P ${HOME}/Software/public
~~~

Untar

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

Pool 2 runs and gzip files

~~~
mkdir -p out/data_shark
cat ${HOME}/Software/process/osf.io/ymswg/REQPOL-R1_S1_L001_R1_001.fastq ${HOME}/Software/process/osf.io/q69de/REQPOL-R1_S1_L001_R1_001.fastq |gzip -c >out/data_shark/reqpol_R1_fw.fastq.gz
cat ${HOME}/Software/process/osf.io/ymswg/REQPOL-R1_S1_L001_R2_001.fastq ${HOME}/Software/process/osf.io/q69de/REQPOL-R1_S1_L001_R2_001.fastq |gzip -c >out/data_shark/reqpol_R1_rv.fastq.gz
cat ${HOME}/Software/process/osf.io/7yf2g/REQPOL-R2_S2_L001_R1_001.fastq ${HOME}/Software/process/osf.io/7pt3x/REQPOL-R2_S2_L001_R1_001.fastq |gzip -c >out/data_shark/reqpol_R2_fw.fastq.gz
cat ${HOME}/Software/process/osf.io/7yf2g/REQPOL-R2_S2_L001_R2_001.fastq ${HOME}/Software/process/osf.io/7pt3x/REQPOL-R2_S2_L001_R2_001.fastq |gzip -c >out/data_shark/reqpol_R2_rv.fastq.gz
cat ${HOME}/Software/process/osf.io/wp6u8/REQPOL-R3_S3_L001_R1_001.fastq ${HOME}/Software/process/osf.io/aqrkb/REQPOL-R3_S3_L001_R1_001.fastq |gzip -c >out/data_shark/reqpol_R3_fw.fastq.gz
cat ${HOME}/Software/process/osf.io/wp6u8/REQPOL-R3_S3_L001_R2_001.fastq ${HOME}/Software/process/osf.io/aqrkb/REQPOL-R3_S3_L001_R2_001.fastq |gzip -c >out/data_shark/reqpol_R3_rv.fastq.gz
~~~

Prepare data

~~~
conda activate snakemake
snakemake -p -c all -s 01snkfl_prep_shark.yml --config data_shark=out/data_shark process_data_dir="${HOME}"/Software/process public_data_dir="${HOME}"/Software/public min_readcount=0 outdir=out/min_readcount_0 container=out/vtam_benchmark.sif --resources db_shark=1 --use-singularity
~~~

