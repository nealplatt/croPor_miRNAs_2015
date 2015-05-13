#miRNA prediction in croc small RNA libraries


FASTX_DIR="/lustre/work/apps/fastx_toolkit-0.0.14/bin"
MIRDEEP_BOWTIE="/lustre/work/apps/mirdeep2_0_0_7/essentials/bowtie-1.1.1"
MIRDEEP_BIN="/lustre/work/apps/mirdeep2_0_0_7/"

ADAPTER_SEQ=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

WORK_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015"
RESULTS_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015/results"
DATA_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015/data"
BIN_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015/bin"


#create a home directory and structure
mkdir $WORK_DIR $RESULTS_DIR $DATA_DIR $BIN_DIR
cd $WORK_DIR

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#populate data directory

#add the current croc genome
wget 							\
	--directory-prefix=$DATA_DIR 			\
	-qc						\
	ftp://crocgenomes.org/pub/ICGWG/Genome_drafts/crocodile.current/croc_sub2.assembly.fasta.gz


GENOME=$DATA_DIR/croc_sub2.assembly.fasta
gunzip $GENOME

# and build an index (in the background) w/ bowtie v.1.1.1
cd $DATA_DIR
$MIRDEEP_BOWTIE/bowtie-build $GENOME $GENOME &

cd $WORK_DIR

################################################################################
#get the raw sequence data (this should be on NCBI eventually)
#
# This is to be done manually at this point.  Once data is publically available...then...
#
################################################################################


#remove all execute and write privelages
chmod a-wx $DATA_DIR/*

mkdir $DATA_DIR/seqRaw

#sequences combined (manually) into necessary files by tissue type (only R1 reads)
zcat $DATA_DIR/D10*R1*.gz > $RESULTS_DIR/seqRaw/D10_testis.fastq
zcat $DATA_DIR/D6*R1*.gz > $RESULTS_DIR/seqRaw/D6_testis.fastq
zcat $DATA_DIR/H6*R1*.gz > $RESULTS_DIR/seqRaw/H6_skinJaw.fastq
zcat $DATA_DIR/J4*R1*.gz > $RESULTS_DIR/seqRaw/J4_cloaca.fastq
zcat $DATA_DIR/J6*R1*.gz > $RESULTS_DIR/seqRaw/J6_skinBelly.fastq
zcat $DATA_DIR/K1*R1*.gz > $RESULTS_DIR/seqRaw/K1_smIntestine.fastq
zcat $DATA_DIR/M2*R1*.gz > $RESULTS_DIR/seqRaw/M2_brain.fastq
zcat $DATA_DIR/M6*R1*.gz > $RESULTS_DIR/seqRaw/M6_tongue.fastq
zcat $DATA_DIR/O4*R1*.gz > $RESULTS_DIR/seqRaw/O4_liver.fastq
zcat $DATA_DIR/P6*R1*.gz > $RESULTS_DIR/seqRaw/P6_heart.fastq
zcat $DATA_DIR/P7*R1*.gz > $RESULTS_DIR/seqRaw/P7_stomach.fastq
zcat $DATA_DIR/P8*R1*.gz > $RESULTS_DIR/seqRaw/P8_heart.fastq
zcat $DATA_DIR/Q3*R1*.gz > $RESULTS_DIR/seqRaw/Q3_brain.fastq
zcat $DATA_DIR/D2*R1*.gz > $RESULTS_DIR/seqRaw/D2_testis.fastq
zcat $DATA_DIR/E1*R1*.gz > $RESULTS_DIR/seqRaw/E1_testis.fastq
zcat $DATA_DIR/J2*R1*.gz > $RESULTS_DIR/seqRaw/J2_spleen.fastq


#count raw data
echo "#####################"  >$RESULTS_DIR/seqCountData.tab
echo "###   RAW READS   ###" >>$RESULTS_DIR/seqCountData.tab
echo "#####################" >>$RESULTS_DIR/seqCountData.tab
wc -l $RESULTS_DIR/seqRaw/*fastq | awk '{print $1/4"\t"$2}' >>$RESULTS_DIR/seqCountData.tab
echo "" >>$RESULTS_DIR/seqCountData.tab
echo "" >>$RESULTS_DIR/seqCountData.tab
echo "" >>$RESULTS_DIR/seqCountData.tab


