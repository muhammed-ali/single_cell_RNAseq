#!/bin/bash -l

#SBATCH -J MyLongJob
#SBATCH -N 1
#SBATCH --ntasks-per-node=6
#SBATCH --time=3-00:00:00
#SBATCH -p gpu

#genome_fasta_files = /scratch/users/mali/Indices/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa
#gtf_file = /scratch/users/mali/Indices/Mus_musculus/Ensembl/GRCm38/Annotation/Archives/archive-2015-07-17-14-32-40/Genes/genes_Star_withCHR.gtf
#genome_directory = /scratch/users/mali/Indices/mm9_Star_Index/

#STAR --runThreadN 6 --runMode genomeGenerate --genomeDir $genome_directory --genomeFastaFiles $genome_fasta_file --sjdbGTFfile $gtf_file --sjdbOverhang 100

STAR --runThreadN 6 --runMode genomeGenerate --genomeDir /scratch/users/mali/Indices/mm9_Star_Index --genomeFastaFiles /scratch/users/mali/Indices/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile /scratch/users/mali/Indices/Mus_musculus/Ensembl/GRCm38/Annotation/Archives/archive-2015-07-17-14-32-40/Genes/genes_Star_withCHR.gtf --sjdbOverhang 100
pwd
