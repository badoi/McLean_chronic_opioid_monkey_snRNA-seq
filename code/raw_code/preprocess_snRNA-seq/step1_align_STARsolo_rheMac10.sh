#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
#SBATCH --time 3-00:00:00
#SBATCH --job-name=grnboost
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=46G
#SBATCH --error=logs/untar_%A.txt
#SBATCH --output=logs/untar_%A.txt
#SBATCH --array=1-16

###################
## 1) set up paths 
PROJDIR=/projects/pfenninggroup/singleCell/McLean_chronic_opioid_monkey_snRNA-seq
DATADIR=$PROJDIR/data/raw_data/
TMPDIR=/scratch/${USER}
SOLODIR=$DATADIR/STARsolo_out

GENOME_DIR=/home/bnphan/resources/genomes/rheMac10
BARCODES=/home/bnphan/resources/cell_ranger_barcodes/inDrops_v3_gel_barcode2.16bp.txt

mkdir -p $TMPDIR $SOLODIR

cd $TMPDIR
NAME=$(awk)

SAMPLE_ID=Fenster_2020

if [[ ! -f "$PROJDIR/data/raw_data/STARsolo_out/${SAMPLE_ID}.Log.final.out" ]]; then
# copy over the fastq files, preserving Run file structure
rsync $DATADIR/fastq/20200305_RF7824-${NAME}.*.fastq.gz $TMPDIR

# find all the files
cDNA_FASTQ=20200305_RF7824-${NAME}.cDNA.fastq.gz
CB_FASTQ=20200305_RF7824-${NAME}.CBUMI.fastq.gz

echo "Aligning samples w/ STARsolo for: ${SAMPLE_ID}."
~/src/STAR-2.7.9a/bin/Linux_x86_64/STAR \
--outFileNamePrefix ${SAMPLE_ID}. \
--readFilesIn $cDNA_FASTQ $CB_FASTQ \
--readFilesCommand zcat \
--genomeDir $GENOME_DIR \
--limitOutSJcollapsed 5000000 \
--runThreadN 12 \
--clipAdapterType CellRanger4 \
--soloType CB_UMI_Simple \
--soloCBwhitelist $BARCODES \
--soloFeatures GeneFull Velocyto \
--soloStrand Forward \
--soloBarcodeReadLength 0 \
--soloCBmatchWLtype 1MM \
--soloCellFilter EmptyDrops_CR \
--soloMultiMappers EM \
--soloCBstart 1 --soloCBlen 16 \
--soloUMIstart 17 --soloUMIlen 6 \
--soloUMIdedup 1MM_CR \
--outSAMtype None

rsync --remove-source-files -Pauv $TMPDIR/${SAMPLE_ID}* $PROJDIR/data/raw_data/STARsolo_out
rm -rf ${SAMPLE_ID}._STARtmp */${SAMPLE_ID}* ${SAMPLE_ID}*
else 
	echo "A completed STARsolo already exists: ${SOLODIR}/${SAMPLE_ID}.Solo.out"
fi
