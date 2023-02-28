#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1,pfen_bigmem,pfen3
#SBATCH --job-name=STARsolo
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --error=logs/align_STARsolo_%A_%a.txt
#SBATCH --output=logs/align_STARsolo_%A_%a.txt
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

SAMPLE_ID=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $1}' \
${PROJDIR}/data/raw_data/tables/metadata.csv )

if [[ ! -f "$PROJDIR/data/raw_data/STARsolo_out/${SAMPLE_ID}.Log.final.out" ]]; then
# copy over the fastq files, preserving Run file structure
cd $TMPDIR
rsync -Paq $DATADIR/fastq/20200305_RF7824-${SAMPLE_ID}.*.fastq.gz $TMPDIR

## find all the files
cDNA_FASTQ=$(ls 20200305_RF7824-${SAMPLE_ID}.cDNA*fastq.gz | tr '[[:space:]]' ',' | sed 's/,$//g')
CB_FASTQ=$(ls 20200305_RF7824-${SAMPLE_ID}.CBUMI*fastq.gz | tr '[[:space:]]' ',' | sed 's/,$//g')

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