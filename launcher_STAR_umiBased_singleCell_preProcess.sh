#!/bin/bash -l

cd /scratch/users/mali/Tasks/CS_mm_PD_SNCA/raw_fastq

# get unique list of fastq file names
ids=( `for i in *.fastq; do a=$(echo $i | cut -d "_" -f 1); echo $a ; done | sort -u` )

for j in "${ids[@]}"
do
	echo $j
	#echo "oarsub -l nodes=1/core=8,walltime=90 -p "cputype='xeon-haswell'" -n $a -S "./TopHatAll.sh $j""
	sbatch /scratch/users/mali/Tasks/CS_mm_PD_SNCA/singleCell_preProcess.sh $j
done

