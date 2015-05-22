FASTX_DIR="/lustre/work/apps/fastx_toolkit-0.0.14/bin"
MIRDEEP_DIR="/lustre/work/apps/mirdeep2_0_0_7"
MIRDEEP_BOWTIE="$MIRDEEP_DIR/essentials/bowtie-1.1.1"
BLAST_DIR=/lustre/work/apps/blast/bin
BEDTOOLS_DIR=/lustre/work/apps/bedtools-2.17.0/bin
WORK_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015"
RESULTS_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015/results"
DATA_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015/data"
BIN_DIR="/lustre/scratch/roplatt/croPor_miRNAs_2015/bin"
ADAPTER_SEQ=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
GENOME=$DATA_DIR/croc_sub2.assembly.fasta
MIRBASE_MATURE=$RESULTS_DIR/gga_matureMirnas.fa
MIRBASE_HAIRPIN=$RESULTS_DIR/gga_hairpinMirnas.fa
MIRDEEP_RESULT_FILE_TSV=$RESULTS_DIR/initalPredictions/result_20_05_2015_t_13_36_32.initialPred.csv
MIRDEEP_RESULT_FILE_BED=$RESULTS_DIR/initalPredictions/result_20_05_2015_t_13_36_32.initialPred.bed
TR_RNA_DB=$RESULTS_DIR/tRNA_rRNA_2015-05-15.fas
MIRDEEP2SCORE=2
HIGHQUAL_MATURE=$RESULTS_DIR/initalPredictions/hq_matureMirna.fas
HIGHQUAL_HAIRPIN=$RESULTS_DIR/initalPredictions/hq_hairpinMirna.fas
QUANT_DIR=$RESULTS_DIR/allQuantified
cd $WORK_DIR
cd $QUANT_DIR
