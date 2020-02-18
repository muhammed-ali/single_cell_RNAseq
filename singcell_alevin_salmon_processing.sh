#!/bin/bash -l

#SBATCH -J singCel
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --time=3-00:00:00
#SBATCH -p gpu

# Tutorial: https://salmon.readthedocs.io/en/latest/alevin.html
# Tutorial: https://combine-lab.github.io/alevin-tutorial/2018/running-alevin/

raw_fastq=/scratch/users/mali/Tasks/CS_mm_PD_SNCA/raw_fastq
salmon_index=/scratch/users/mali/Indices/salmon_mm9_index/salmon_index
trx_gene=/scratch/users/mali/Indices/salmon_mm9_index/txp2gene.tsv
R1=_S1_R1_001.fastq
R2=_S1_R2_001.fastq

#mkdir alevin_output
#salmon alevin -l ISR -1 $raw_fastq/mmCortex_R1.fastq -2 $raw_fastq/mmCortex_R2.fastq --dropseq -i $salmon_index -p 12 -o ./alevin_output --tgMap $trx_gene

#mkdir alevin_output2
#salmon alevin -l ISR -1 $raw_fastq/MAE1_S1_R1_001.fastq $raw_fastq/MAE2_S1_R1_001.fastq $raw_fastq/MAE3_S1_R1_001.fastq $raw_fastq/MAE4_S1_R1_001.fastq -2 $raw_fastq/MAE1_S1_R2_001.fastq $raw_fastq/MAE2_S1_R2_001.fastq $raw_fastq/MAE3_S1_R2_001.fastq $raw_fastq/MAE4_S1_R2_001.fastq --dropseq -i $salmon_index -p 12 -o ./alevin_output2 --tgMap $trx_gene --dumpFeatures

cd /scratch/users/mali/Tasks/CS_mm_PD_SNCA/raw_fastq

ids=( `for i in *.fastq; do a=$(echo $i | cut -d "_" -f 1); echo $a ; done | sort -u` )

for j in "${ids[@]}"
do
 cd /scratch/users/mali/Tasks/CS_mm_PD_SNCA/alevin_salmon_processing
 mkdir $j
 echo "salmon alevin -l ISR -1 $raw_fastq/$j$R1 -2 $raw_fastq/$j$R2 --dropseq -i $salmon_index -p 12 -o ./$j --tgMap $trx_gene"
 salmon alevin -l ISR -1 $raw_fastq/$j$R1 -2 $raw_fastq/$j$R2 --dropseq -i $salmon_index -p 12 -o ./$j --tgMap $trx_gene
done

echo "alevin pipeline for single-cell analysis is finished"

pwd
