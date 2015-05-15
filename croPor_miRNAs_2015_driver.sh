# RNPlatt
# 13 May 2015
#
# Pipeline to identify miRNAs in Crocodylus porusus small RNA libraries


################################################################################
#
# Step 1) setting env variables and create directory structure
#
################################################################################
# Set program/directory paths here

# directory to FASTX toolkit exes
FASTX_DIR="/lustre/work/apps/fastx_toolkit-0.0.14/bin"

# directory to miRDeep2
MIRDEEP_DIR="/lustre/work/apps/mirdeep2_0_0_7"

# directory to bowtie that is included with miRDeep2 
MIRDEEP_BOWTIE="$MIRDEEP_DIR/essentials/bowtie-1.1.1"

# directory to blast
BLAST_DIR=/lustre/work/apps/blast/bin

# directory to bedtools
BEDTOOLS_DIR=/lustre/work/apps/bedtools-2.17.0/bin

# general work directory...think of as project "home"
WORK_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015"

# directory for all manipulated data
RESULTS_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015/results"

# directory for (raw) data that will not be modified
DATA_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015/data"

# bin directory for custom scripts
BIN_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015/bin"

# Illumina sequencing adapter - will be used for clipping in fastx steps
ADAPTER_SEQ=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

# Create subdirectories in project directory 
mkdir $WORK_DIR $RESULTS_DIR $DATA_DIR $BIN_DIR
cd $WORK_DIR


################################################################################
#
# Step 2) Populate directories with available data
#
################################################################################

#Step 2 -1
# Download the croc genome
wget 							\
	--directory-prefix=$DATA_DIR 			\
	-qc						\
	ftp://crocgenomes.org/pub/ICGWG/Genome_drafts/crocodile.current/croc_sub2.assembly.fasta.gz

GENOME=$DATA_DIR/croc_sub2.assembly.fasta

# Go ahead and create a bowtie index for the genome (since thsi takes forever)
gunzip $GENOME

cd $DATA_DIR

$MIRDEEP_BOWTIE/bowtie-build $GENOME $GENOME &

cd $WORK_DIR

#------------------------
#Step 2-2
# Get all fo the sequence data (manually).  This will eventually be available
#   through NCBI.


# Change permissions on this raw data so that it is accidentally modified
chmod a-wx $DATA_DIR/*

#------------------------
#Step 2-3
# All of the seqeunce data needs to be combined into library specific files
#   and placed in a seqRaw directory
 
mkdir $DATA_DIR/seqRaw

# sequences combined (manually) into necessary files by tissue type (only R1 reads)
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


# count the available data (per library)
echo "#####################"  >$RESULTS_DIR/seqCountData.tab
echo "###   RAW READS   ###" >>$RESULTS_DIR/seqCountData.tab
echo "#####################" >>$RESULTS_DIR/seqCountData.tab
wc -l $RESULTS_DIR/seqRaw/*fastq | awk '{print $1/4"\t"$2}' >>$RESULTS_DIR/seqCountData.tab
echo "" >>$RESULTS_DIR/seqCountData.tab
echo "" >>$RESULTS_DIR/seqCountData.tab
echo "" >>$RESULTS_DIR/seqCountData.tab

################################################################################
#
# Step 3) Clean up the raw data and judge quality
#
################################################################################

#Step 3-1
# Create a seqQCd folder to hold all cleaned sequence reads
mkdir $RESULTS_DIR/seqQCd

#------------------------
#Step 3-2
# Using a for loop clip and quality trim each library
for RAW_READS in $RESULTS_DIR/seqRaw/*.fastq
do

      SAMPLE_ID=$(basename $RAW_READS .fastq)
	
	#------------------------
	# Step 3-2a - remove adapters and trim sequences were half the read is
	#   Q20 or less
        $FASTX_DIR/fastx_clipper                                               	\
                -n								\
		-l 10                                                           \
                -a $ADAPTER_SEQ                                                 \
                -i $RAW_READS                                                   \
                | $FASTX_HOME/fastq_quality_filter                              \
                        -Q33                                                    \
                        -q 20                                                   \
                        -p 50							\
                        -o $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clipped.fastq"  

	#------------------------
	# Step 3-2b - caclulate quality stats and output .png figures
        $FASTX_DIR/fastx_quality_stats                                	\
                -Q33                                                   	\
                -o $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clipped.stats"      	\
                -i $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clipped.fastq"

        $FASTX_DIR/fastx_nucleotide_distribution_graph.sh             	\
                -i $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clipped.stats"      	\
                -o $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clippedNUC.png"     	\
                -t $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clipped"

        $FASTX_DIR/fastq_quality_boxplot_graph.sh                     	\
                -i $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clipped.stats"      	\
                -o $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clippedBOX.png"     	\
                -t $RESULTS_DIR/seqQCd/$SAMPLE_ID"_clipped"
done


# count the available data (per library) after QC
echo "#####################" >>$RESULTS_DIR/seqCountData.tab
echo "###    CLIPPED    ###" >>$RESULTS_DIR/seqCountData.tab
echo "#####################" >>$RESULTS_DIR/seqCountData.tab
wc -l $RESULTS_DIR/seqQCd/*clipped.fastq | awk '{print $1/4"\t"$2}' >>$RESULTS_DIR/seqCountData.tab
echo "" >>$RESULTS_DIR/seqCountData.tab
echo "" >>$RESULTS_DIR/seqCountData.tab
echo "" >>$RESULTS_DIR/seqCountData.tab

################################################################################
#
# Step 4) Initial miRNA prediction
#
################################################################################

# Step 4-1 - make folder in results for the inital predictions
mkdir $RESULTS_DIR/initalPredictions
cd $RESULTS_DIR/initalPredictions

#------------------------
#Step 4-2
# miRDeep2 needs a config file to run analysis on combined data
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

#------------------------
#Step 4-3
# Get and process the miRBase miRNAs in a way acceptable to miRDeep2

#get mirbase sequences
wget 								\
	--directory-prefix=$DATA_DIR 				\
	ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz

wget 								\
	--directory-prefix=$DATA_DIR 				\
	ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz

#these files have white spaces that need to be removed (for miRDeep2)
zcat $DATA_DIR/hairpin.fa.gz | cut -f1 -d" ">$RESULTS_DIR/hairpin_noSpace.fa
zcat $DATA_DIR/mature.fa.gz | cut -f1 -d" " >$RESULTS_DIR/mature_noSpace.fa

#and non-canonical nucleotides and convert to SL (for easier parsing)
fastaparse.pl $RESULTS_DIR/hairpin_noSpace.fa -b 	\
	| $FASTX_DIR/fasta_formatter -w 0	\
	 >$RESULTS_DIR/hairpin_cleaned.fa

fastaparse.pl $RESULTS_DIR/mature_noSpace.fa -b 	\
	| $FASTX_DIR/fasta_formatter -w 0	\
	 >$RESULTS_DIR/mature_cleaned.fa

#grep out the chicken miRNAs
grep -A1 "gga" $RESULTS_DIR/mature_cleaned.fa  | grep -v -- "^--$" >$RESULTS_DIR/gga_matureMirnas.fa
grep -A1 "gga" $RESULTS_DIR/hairpin_cleaned.fa | grep -v -- "^--$" >$RESULTS_DIR/gga_hairpinMirnas.fa

# designate the processed miRBase files for future analyses
MIRBASE_MATURE=$RESULTS_DIR/gga_matureMirnas.fa
MIRBASE_HAIRPIN=$RESULTS_DIR/gga_hairpinMirnas.fa

#------------------------
#Step 4-4 
# Begin miRDeep process with mapper.pl

$MIRDEEP_DIR/mapper.pl 			\
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

#------------------------
#Step 4-5 
# Identify conserved and predict novel miRNAs with miRDeep2

$MIRDEEP_DIR/miRDeep2.pl 		\
	config_mapperProcessed.fa	\
	$GENOME				\
	config_mapper.arf	 	\
	$MIRBASE_MATURE			\
	none				\
	$MIRBASE_HAIRPIN		\
	-z .initialPred			\
	-P



################################################################################
#
# Step 5) Process miRDeep2 results
#
################################################################################

#quality filtering done in 5ish steps.

#1) Filter out any miRNAs with mirdeep score less than 1 (awk)

#2) All predicted miRNAs that are similar to tRNAs or rRNAs are removed. (blast)

#3) Keep on the highest scoring in overlapping pairs (perl/bedtools)

#4) Overlap with the known/identified galGal miRNAs, retain galGal designation (bedtools)

#5) create a fasta file to run with quantifier.pl (bedtools)

#...) down the line remove miRNAs that are only expressed in one tissue


MIRDEEP_RESULT_FILE_TSV=result_14_05_2015_t_15_35_12.test.csv
$MIRDEEP_RESULT_FILE_BED=result_14_05_2015_t_15_35_12.test.bed

TR_RNA_DB=$RESULTS_DIR/tRNA_rRNA_2015-05-15.fas

#------------------------
#Step 5-1
# make a file of miRNA hairpins >= 1 scrore 

#Scroll through results file to find the lines containing novel miRNAs
cat -n $MIRDEEP_RESULT_FILE_TSV
#In this case lines 27-462 contain novel miRNAs

#cat only the novel miRNAs and print out the ones score >=1 to a fasta file
cat $MIRDEEP_RESULT_FILE_TSV | head -462 | tail -435 	\
	| awk '{if ($2 >= 1) print">"$1"\n"$18}' 	\
	>$RESULTS_DIR/predictedHairpin.fas

#------------------------
#Step 5-2
# Compare miRNA hairpins to tRNAs and rRNAs via Blast

#dowload rRNA and tRNA dbs - prep for blast
wget 									\
	--directory-prefix=$DATA_DIR					\
	http://gtrnadb.ucsc.edu/download/tRNAs/eukaryotic-tRNAs.fa.gz

#rRNA db done by hand from http://www.arb-silva.de/download/
# make sure to get both the Large and Small subunits.
# Downloaded all vertebrate rRNAs to >$DATA_DIR/eukaryotic-rRNAs.fa.gz

#combine files
zcat $DATA_DIR/eukaryotic-tRNAs.fa.gz $DATA_DIR/eukaryotic-rRNAs.fa.gz \
	>$TR_RNA_DB

#then make a blast db
$BLAST_DIR/makeblastdb -in $TR_RNA_DB -dbtype nucl

$BLAST_DIR/blastn \
	-db $TR_RNA_DB \
	-query $RESULTS_DIR/predictedHairpin.fas \
	-outfmt 6 \
	>$RESULTS_DIR/predictedHairpin_vs_trRNA_blastn.out

#in this case there were no hits to r or tRNAs that need to be removed


#------------------------
#Step 5-3 & 5.4
# Remove overlapping miRNA predictions

#use bedtools to sort the miRNA prediction results (high scoring).  This will be 
# done in two overall steps.  overlap predicted miRNAs, then overlap the galGal
# miRNAs with the merged predictions.  Precedence goes Known -> High Expression

#an initial test shows that there are not any overlaps between known and novel
# miRNAs, 
#
#  cat $MIRDEEP_RESULT_FILE_BED | bedtools sort -i - | bedtools merge -d 1 -nms
# | grep known | grep novel  


# so can proceed with merging only the novel miRNAs
cat $MIRDEEP_RESULT_FILE_BED \
	| grep novel \
	| bedtools sort -i - \
	| bedtools merge -d 1 -i - \
	| bedtools intersect -wao -a - -b $MIRDEEP_RESULT_FILE_BED \
	>$RESULTS_DIR/novel_miRNAoverlap.bed

#use a custom perl script to retain the highest scoring partner from each
# overlapping pair (script is available on github)
cat $RESULTS_DIR/novel_miRNAoverlap.bed \
	| perl $BIN_DIR/removeOverlappingMirnas.pl \
	>$RESULTS_DIR/novel_miRNA_nonoverlap.bed


#------------------------
#Step 5-5
# Create fasta file of HQ miRNAs

#using this information, two files were created (one for mature miRNAs and
# the other for haiprins) that included the known and non-overlapping high qual
# predicted miRNAs.  At this point I am kind of feeling lazy so this was done
# manually.  Will find a CL solution at a later time. File names are listed below

HIGHQUAL_MATURE=$RESULTS_DIR/hq_matureMirna.fas
HIGHQUAL_HAIRPIN=$RESULTS_DIR/hq_hairpinMirna.fas
#these files will be used to quantify miRNA expression levels




################################################################################
#
# Step 6) Quantify miRNA expression levels
#
################################################################################

#re-run miRDeep2 so that it will map to the HQ miRNAs - no new predictions
$MIRDEEP_DIR/miRDeep2.pl \
	config_mapperProcessed.fa \
	$GENOME \
	config_mapper.arf \
	$HIGHQUAL_MATURE \
	none \
	$HIGHQUAL_HAIRPIN \
	-z .secondRndPred
























