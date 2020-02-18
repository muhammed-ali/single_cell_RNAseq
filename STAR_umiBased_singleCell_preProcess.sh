#!/bin/bash -l

#SBATCH -J singCel
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --time=3-00:00:00
#SBATCH -p gpu

# Guide: https://github.com/CGATOxford/UMI-tools/blob/master/doc/Single_cell_tutorial.md

proj_id=mmCortex
thread=12
R1=_S1_R1_001.fastq
R2=_S1_R2_001.fastq
trim_R1=_trimmed_R1.fastq
trim_R2=_trimmed_R2.fastq
pattern=CCCCCCCCCCCCNNNNNNNN
ncells=100

raw_fastq=/scratch/users/mali/Tasks/CS_mm_PD_SNCA/raw_fastq
fastqc_raw=/scratch/users/mali/Tasks/CS_mm_PD_SNCA/fastqc_raw
fastqc_trim=/scratch/users/mali/Tasks/CS_mm_PD_SNCA/fastqc_trim
cutadapt_trim=/scratch/users/mali/Tasks/CS_mm_PD_SNCA/cutadapt_trim
umi_output=/scratch/users/mali/Tasks/CS_mm_PD_SNCA/umi_output
gtf=/scratch/users/mali/Indices/Mus_musculus/Ensembl/GRCm38/Annotation/Archives/archive-2015-07-17-14-32-40/Genes/genes_Star_withCHR.gtf

filename=$(echo $1)
echo $filename

cd /scratch/users/mali/Tasks/CS_mm_PD_SNCA

module use /opt/apps/resif/data/production/v1.2/bioinfo/modules/all
module load bio/FastQC/0.11.8-Java-1.8
module load bio/SAMtools
module load lang/Python

module load compiler/GCC
module load lib/zlib


################
#### PART 1 #### 
################

# While you run PART 1, the PART 2 commands must be commented out

# This bart should be run in parallel jobs so that all the paired-fastq files must be processed simultaneously
# call: sbatch launcher_singleCell_preProcess.sh
# from directory: /scratch/users/mali/Tasks/CS_mm_PD_SNCA 

<<'end_long_comment'

# fastqc -o ./fastqc_raw $raw_fastq/$filename$R1 $raw_fastq/$filename$R2
#echo "fastQC done on raw fastq files"

# cutadapt -j "$thread" -m 10 -q 20 -o $cutadapt_trim/$filename$trim_R1 -p $cutadapt_trim/$filename$trim_R2 $raw_fastq/$filename$R1 $raw_fastq/$filename$R2
#echo "adapters are cut from the raw fastq files"

# Not required, nobody does it.
# fastqc -o ./fastqc_trim $cutadapt_trim/$filename$trim_R1 $cutadapt_trim/$filename$trim_R2
#echo "fastqc done on adapter-trimmed fastq files"
#pwd

end_long_comment

################
#### PART 2 #### 
################

# This bart can be run in one go by calling: sbatch singleCell_preProcess.sh MAE1
# from directory: /scratch/users/mali/Tasks/CS_mm_PD_SNCA

# <<'end_long_comment'

# NOTE: Make sure the file naming conventions are compatible, may very from experiment to experiment
# Whle using trimmed file umi_tools "whitelist" threw error as barcodes of some reads are disturbed during trimming.
# Therefore, I and the tutotial I follow use the raw fastq read not the cutadapt trimmed reads.
#cat $cutadapt_trim/MAE?_trimmed_R1.fastq > $cutadapt_trim/"$proj_id"_R1.fastq;
#cat $cutadapt_trim/MAE?_trimmed_R2.fastq > $cutadapt_trim/"$proj_id"_R2.fastq;
cat $raw_fastq/MAE?_S1_R1_001.fastq > $raw_fastq/"$proj_id"_R1.fastq;
cat $raw_fastq/MAE?_S1_R2_001.fastq > $raw_fastq/"$proj_id"_R2.fastq;


## Identify correct cell barcodes
# NOTE: Make sure pattern is right
# NOTE: Make sure you know the cell number, if you do use: --set-cell-number=100, if not, use: --expect-cells=200 \
# NOTE: 5’ end of the read (this can be changed with --3prime) See guide for details
umi_tools whitelist --stdin $raw_fastq/"$proj_id"_R1.fastq --bc-pattern=$pattern --extract-method=string --set-cell-number=$ncells --log2stderr > $umi_output/whitelist.txt


## Extract barcdoes and UMIs and add to read names
#check paired reads
umi_tools extract --bc-pattern=$pattern --stdin $raw_fastq/"$proj_id"_R1.fastq --stdout $raw_fastq/"$proj_id"_R1_extracted.fastq --read2-in $raw_fastq/"$proj_id"_R2.fastq --read2-out=$raw_fastq/"$proj_id"_R2_extracted.fastq --filter-cell-barcode --whitelist=$umi_output/whitelist.txt


mkdir results_STAR
## Map reads 
STAR --runThreadN $thread --genomeDir /scratch/users/mali/Indices/mm9_Star_Index --readFilesIn $raw_fastq/"$proj_id"_R2_extracted.fastq --readFilesCommand cat --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./results_STAR/ 
#echo "Reads alignment by STAR is done"


mkdir results_Counts
## Assigning reads to genes
featureCounts -a $gtf -o ./results_Counts/gene_assigned -R BAM ./results_STAR/Aligned.sortedByCoord.out.bam -T $thread
echo "Reads are assigned to genes by featureCounts"


## SORT and INDEX the bam file
samtools sort ./results_Counts/Aligned.sortedByCoord.out.bam.featureCounts.bam -o ./results_Counts/gene_assigned_sorted.bam
samtools index ./results_Counts/gene_assigned_sorted.bam
echo "BAM file is sorted and indexed"


## Counting molecules # Count UMIs per gene per cell
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I ./results_Counts/gene_assigned_sorted.bam -L ./results_Counts/gene_counts.log -S ./results_Counts/gene_counts.tsv
# output in wide format, cell names as columns
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --wide-format-cell-counts -I ./results_Counts/gene_assigned_sorted.bam -L ./results_Counts/gene_counts_wide.log -S ./results_Counts/gene_counts_wide.tsv
echo "gene counts are retreived by using umi_tools"

# end_long_comment

pwd
