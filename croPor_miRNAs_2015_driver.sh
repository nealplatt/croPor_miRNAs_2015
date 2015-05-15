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
MIRDEEP_DIR="/lustre/work/apps/mirdeep2_0_0_7/"

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

#------------------------
#Step 4-5 
# Identify conserved and predict novel miRNAs with miRDeep2

$MIRDEEP_BIN/miRDeep2.pl 		\
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

#quality filtering done in XXXXX steps.

#1) Filter out any miRNAs with mirdeep score less than 1 (awk)

#2) All predicted miRNAs that are similar to tRNAs or rRNAs are removed. (blast)

#3) Keep on the highest scoring in overlapping pairs (perl/bedtools)

#4) Overlap with the known/identified galGal miRNAs, retain galGal designation (bedtools)

#5) create a fasta file to run with quantifier.pl (bedtools)

#6) down the line remove miRNAs that are only expressed in one tissue


MIRDEEP_RESULT_FILE_TSV= #<-insert this here

cat -n $MIRDEEP_RESULT_FILE_TSV
#scrolling through the result.csv file shows that the newly predicted miRNAs
# are on lines 27-462.  Need to use the fasta sequences for these to quality filter

 

cat $MIRDEEP_RESULT_FILE_TSV | head -462 | tail -435 	\
	| awk '{if ($2 >= 1) print">"$1"\n"$18}' 	\
	>$RESULTS_DIR/predictedHairpin.fas


#dowload rRNA and tRNA dbs - prep for blast
wget 					\
	--directory-prefix=$DATA_DIR	\
	http://gtrnadb.ucsc.edu/download/tRNAs/eukaryotic-tRNAs.fa.gz

#rRNA db done by hand from http://www.arb-silva.de/download/
# make sure to get both the Large and Small subunits.
# Downloaded all vertebrate rRNAs to >$DATA_DIR/eukaryotic-rRNAs.fa.gz

zcat $DATA_DIR/eukaryotic-tRNAs.fa.gz $DATA_DIR/eukaryotic-rRNAs.fa.gz \
	>$RESULTS_DIR/tRNA_rRNA_2015-05-15.fas

TR_RNA_DB=$RESULTS_DIR/tRNA_rRNA_2015-05-15.fas



$BLAST_DIR/makeblastdb -in tRNA_rRNA_2015-04-14.fas -dbtype nucl

$BLAST_DIR/blastn \
	-db tRNA_rRNA_2015-04-14.fas \
	-query initialAll_nonOverlap_miRDeep2Predictions.fas \
	-outfmt 6 \
		| sort -k1,1 -k12,12gr -k11,11g -k3,3gr \
		| sort -u -k1,1 --merge >initialAll_vs_trRNA_blastn.out



 bedtools sort -i result_05_05_2015_t_15_10_49.initialALL.bed >result_05_05_2015_t_15_10_49.initialALL_sorted.bed 
merge -d -1 -i result_05_05_2015_t_15_10_49.initialALL_sorted.bed | bedtools intersect -wao -a - -b result_05_05_2015_t_15_10_49.initialALL_sorted.bed >tmp
cat tmp | perl ../bin/removeOverlappingMirnas.pl >initialAll_nonOverlap_miRDeep2Predictions.bed
bedtools getfasta -name -fi ../../genomes/croc_sub2.assembly.fasta  -bed initialAll_nonOverlap_miRDeep2Predictions.bed -fo initialAll_nonOverlap_miRDeep2Predictions.fas


/lustre/work/apps/blast/bin/makeblastdb -in miRBase_v21.fas -dbtype nucl

#go to local computer to create DB and combine.
# delete first table in mirdeep results
remove "NOVEL:"

combine in sql

vertical lookup in excel

create fasta -> cdhitsuite
awk '{print ">"$17"\n"$14}' miRBaseAnnotations_miRDeepPreds_toCDHit.txt >toCDHit_mature.fas
 awk '{print ">"$17"\n"$15}' miRBaseAnnotations_miRDeepPreds_toCDHit.txt >toCDHit_precursor.fas



for SORTED_BED in *_sorted.bed
do

	#create a base name for downstream files
	BASE=$(basename $SORTED_BED _sorted.bed)
	SAMPLE_ID=$(cut -c1,2)

	
	#for each file find the BEST scoring miRNA from those that overlap
	$BEDTOOLS_HOME/bedtools merge 	\
		-d -1 			\
		-i $SORTED_BED 		\
		| $BEDTOOLS_HOME/bedtools intersect 	\
			-wao 				\
			-a - 				\
			-b $SORTED_BED 			\
		| perl removeOverlappingMirnas.pl \
			>$BASE"_nonOverlapMirnas.bed"


	#extract the fasta sequence from the best scoring miRNA hairpin
	$BEDTOOLS_HOME/bedtools getfasta 		\
		-name					\
		-fi $GENOME 				\
		-bed $BASE"_nonOverlapMirnas.bed"	\
		-fo $BASE".fasta"

	#and blast that seqeunce against a tRNA and rRNA db
	$BLAST_HOME/blastn 					\
		-db $TR_RNA_DB		 			\
		-query $BASE".fasta"  				\
		-out $BASE"_hairpins_vs_trRNA_blastn.out" 	\
		-outfmt 6 					\
		| sort -k1,1 -k12,12gr -k11,11g -k3,3gr 	\
		| sort -u -k1,1 --merge >$BASE"_hairpins_vs_trRNA_blastn.out"
		
	#...then against a miRNA db (only output highest scoring seq)
	$BLAST_HOME/blastn 					\
		-db $MIRNA_DB	 				\
		-query $BASE".fasta"  				\
		-out $BASE"_hairpins_vs_miRBase21_blastn.out" 	\
		-outfmt 6 					\
		| sort -k1,1 -k12,12gr -k11,11g -k3,3gr 	\
		| sort -u -k1,1 --merge >$BASE"_hairpins_vs_miRBase21_blastn.out"
	
done
