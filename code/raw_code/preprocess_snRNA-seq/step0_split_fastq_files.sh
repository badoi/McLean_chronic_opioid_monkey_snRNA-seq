#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
#SBATCH --exclude=compute-1-11,compute-1-12,compute-1-35
#SBATCH --time 3-00:00:00
#SBATCH --job-name=grnboost
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=46G
#SBATCH --error=logs/untar_%A.txt
#SBATCH --output=logs/untar_%A.txt

###################
## 1) set up paths 
PROJDIR=/projects/pfenninggroup/singleCell/McLean_chronic_opioid_monkey_snRNA-seq
DATADIR=$PROJDIR/data/raw_data/
TMPDIR=/scratch/${USER}

cd $PROJDIR
mkdir -p $DATADIR/fastq ${DATADIR}

##################################################
## 2) untar the fastq's into the raw fastq folder
cd ${PROJDIR}/data/tidy_data/fenster_round1_analyses/1_input
tar -zxvf 200305_RF7824_R1_R2_fastq.tar.gz
tar -zxvf 200305_RF7824_R3_R4_fastq.tar.gz

mv */*.fastq.gz ${PROJDIR}/data/tidy_data/fenster_round1_analyses/1_input

############################################################
## 3) demultiplex the sequencing reads by sample index in R3
cd ${TMPDIR}
rsync -Pav ${PROJDIR}/data/tidy_data/fenster_round1_analyses/1_input *.fastq.gz .
R1=20200305_pool_RF7824_S1_R1_001.fastq.gz
R2=20200305_pool_RF7824_S1_R2_001.fastq.gz
R3=20200305_pool_RF7824_S1_R3_001.fastq.gz
R4=20200305_pool_RF7824_S1_R4_001.fastq.gz

########################################################################
## 2) paste the R3, R2, R4 together to create combined cell barcode/UMI
## R3 is first, used to demultiplex the reads by sample later on
R5=20200305_pool_RF7824_S1_R5_001.fastq.gz
if [[ !-f ${R5} ]]; then
paste <(zcat $R3) <(zcat $R2) <(zcat $R4) | \
awk '{if (NR%4==1) {
print $1 " " $2
} else if (NR%4==3) {
print "+"
} else {
print $1 $2 $3
}}' | gzip > ${R5}
fi

##########################################################################
## 3) demultiplex the reads based on the w/ first 8bp of new cell barcode
# https://cutadapt.readthedocs.io/en/stable/guide.html#demultiplexing
mamba activate cutadaptenv
cutadapt -e 1 -g ^file:$DATADIR/tables/barcodes.fasta \
-o 20200305_RF7824-{name}.CBUMI.fastq.gz \
-p 20200305_RF7824-{name}.cDNA.fastq.gz \
-j 16 $R5 $R1
mamba deactivate

########################################
## 4) copy back the demultiplexed files
ls 20200305_RF7824*.fastq.gz | \
xargs -n1 -P8 -I% rsync -Pa % $DATADIR/fastq
rsync -Pav 20200305_pool_RF7824* \
${PROJDIR}/data/tidy_data/fenster_round1_analyses/1_input/

