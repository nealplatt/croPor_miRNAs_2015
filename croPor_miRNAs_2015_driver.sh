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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

################################################################################
#now raw data needs to be cleaned and processed


mkdir $RESULTS_DIR/seqQCd


for RAW_READS in $RESULTS_DIR/seqRaw/*.fastq
do

      SAMPLE_ID=$(basename $RAW_READS .fastq)


        $FASTX_DIR/fastx_clipper                                               \
                -v                                                              \
                -a $ADAPTER_SEQ                                                 \
                -i $RAW_READS                                                   \
                | $FASTX_HOME/fastq_quality_filter                              \
                        -Q33                                                    \
                        -q 20                                                   \
                        -p 50							\
                        -o $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clipped.fastq"  

        $FASTX_DIR/fastx_quality_stats                                \
                -Q33                                                   \
                -o $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clipped.stats"      \
                -i $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clipped.fastq"

        $FASTX_DIR/fastx_nucleotide_distribution_graph.sh             \
                -i $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clipped.stats"      \
                -o $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clippedNUC.png"     \
                -t $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clipped"

        $FASTX_DIR/fastq_quality_boxplot_graph.sh                     \
                -i $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clipped.stats"      \
                -o $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clippedBOX.png"     \
                -t $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clipped"
done


#count raw data
echo "#####################" >>$RESULTS_DIR/seqCountData.tab
echo "###    CLIPPED    ###" >>$RESULTS_DIR/seqCountData.tab
echo "#####################" >>$RESULTS_DIR/seqCountData.tab
wc -l $RESULTS_DIR/seqQCd/*clipped.fastq | awk '{print $1/4"\t"$2}' >>$RESULTS_DIR/seqCountData.tab
echo "" >>$RESULTS_DIR/seqCountData.tab
echo "" >>$RESULTS_DIR/seqCountData.tab
echo "" >>$RESULTS_DIR/seqCountData.tab

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#start running miRNA prediction (inital on all samples combined).

mkdir $RESULTS_DIR/initalPredictions
cd $RESULTS_DIR/initalPredictions

#create config file

echo "$RESULTS_DIR/seqQCd/D10_testis_clipped.fastq	D10"  >config.txt
echo "$RESULTS_DIR/seqQCd/D6_testis_clipped.fastq	D06" >>config.txt
echo "$RESULTS_DIR/seqQCd/H6_jawSkin_clipped.fastq	H06" >>config.txt
echo "$RESULTS_DIR/seqQCd/J4_cloaca_clipped.fastq	J04" >>config.txt
echo "$RESULTS_DIR/seqQCd/J6_bellySkin_clipped.fastq	J06" >>config.txt
echo "$RESULTS_DIR/seqQCd/K1_smIntestine_clipped.fastq	K01" >>config.txt
echo "$RESULTS_DIR/seqQCd/M2_brain_clipped.fastq	M02" >>config.txt
echo "$RESULTS_DIR/seqQCd/M6_tongue_clipped.fastq	M06" >>config.txt
echo "$RESULTS_DIR/seqQCd/O4_liver_clipped.fastq	O04" >>config.txt
echo "$RESULTS_DIR/seqQCd/P6_heart_clipped.fastq	P06" >>config.txt
echo "$RESULTS_DIR/seqQCd/P7_stomach_clipped.fastq	P07" >>config.txt
echo "$RESULTS_DIR/seqQCd/P8_heart_clipped.fastq	P08" >>config.txt
echo "$RESULTS_DIR/seqQCd/Q3_brain_clipped.fastq	Q03" >>config.txt
echo "$RESULTS_DIR/seqQCd/D2_testis_clipped.fastq	D02" >>config.txt
echo "$RESULTS_DIR/seqQCd/E1_testis_clipped.fastq	E01" >>config.txt
echo "$RESULTS_DIR/seqQCd/J2_spleen_clipped.fastq	J02" >>config.txt


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#get mirbase sequences
wget 								\
	--directory-prefix=$DATA_DIR 				\
	ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz

wget 								\
	--directory-prefix=$DATA_DIR 				\
	ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz

#these files have white spaces that need to be removed (for miRDeep2)
zcat $DATA_DIR/hairpin.fa.gz | cut -f1 -d" ">$DATA_DIR/hairpin_noSpace.fa
zcat $DATA_DIR/mature.fa.gz | cut -f1 -d" " >$DATA_DIR/mature_noSpace.fa

#and non-canonical nucleotides
fastaparse.pl $DATA_DIR/hairpin_noSpace.fa -b >$DATA_DIR/hairpin_cleaned.fa
fastaparse.pl $DATA_DIR/mature_noSpace.fa -b >$DATA_DIR/mature_cleaned.fa


MIRBASE_MATURE=$DATA_DIR/mature_cleaned.fa
MIRBASE_HAIRPIN=$DATA_DIR/hairpin_cleaned.fa

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#map all small RNA reads to the genome using mapper.pl
$MIRDEEP_BIN/mapper.pl 			\
	config.txt 			\
	-d 				\
	-e 				\
	-h 				\
	-m 				\
	-j 				\
	-l 18 				\
	-v 				\
	-n 				\
	-s config_mapperProcessed.fa 	\
	-t config_mapper.arf 		\
	-p $GENOME  					

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

$MIRDEEP_BIN/miRDeep2.pl 		\
	config_mapperProcessed.fa	\
	$GENOME				\
	config_mapper.arf	 	\
	$MIRBASE_MATURE			\
	none				\
	$MIRBASE_HAIRPIN		\
	-z .initialPred			\
	-P

