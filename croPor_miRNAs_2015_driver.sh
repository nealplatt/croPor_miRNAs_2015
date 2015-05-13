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
cd $WORK_DIR
mkdir $RESULTS_DIR $DATA_DIR $BIN_DIR


#populate data directory
#add the current croc genome
wget 							\
	--directory-prefix=$DATA_DIR 			\
	-bqc						\
	ftp://crocgenomes.org/pub/ICGWG/Genome_drafts/crocodile.current/croc_sub2.assembly.fasta.gz

#@@@@@@@@@@@@@
GENOME=$WORK_HOME/data/genomes/croc_sub2.assembly.fasta
