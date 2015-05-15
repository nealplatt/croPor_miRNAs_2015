FASTX_DIR="/lustre/work/apps/fastx_toolkit-0.0.14/bin"
GENOME=$DATA_DIR/croc_sub2.assembly.fasta
MIRDEEP_DIR="/lustre/work/apps/mirdeep2_0_0_7/"
MIRDEEP_BOWTIE="$MIRDEEP_DIR/essentials/bowtie-1.1.1"
BLAST_DIR=/lustre/work/apps/blast/bin
BEDTOOLS_DIR=/lustre/work/apps/bedtools-2.17.0/bin
WORK_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015"
RESULTS_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015/results"
DATA_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015/data"
BIN_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015/bin"
ADAPTER_SEQ=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
MIRBASE_MATURE=$RESULTS_DIR/gga_matureMirnas.fa
MIRBASE_HAIRPIN=$RESULTS_DIR/gga_hairpinMirnas.fa
MIRDEEP_RESULT_FILE_TSV=result_14_05_2015_t_15_35_12.test.csv
MIRDEEP_RESULT_FILE_BED=result_14_05_2015_t_15_35_12.test.bed
TR_RNA_DB=$RESULTS_DIR/tRNA_rRNA_2015-05-15.fas
HIGHQUAL_MATURE=$RESULTS_DIR/hq_matureMirna.fas
HIGHQUAL_HAIRPIN=$RESULTS_DIR/hq_hairpinMirna.fas


PLAY=/lustre/scratch/roplatt/crocSmallRNA/crocMirnasFullAnalysis/results/play
cd $PLAY
