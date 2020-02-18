#!/bin/bash -l

#SBATCH -J makIndx
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --time=3-00:00:00
#SBATCH -p gpu

# guide: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

cd /scratch/users/mali/Indices/salmon_mm9_index

# We are first going to download the reference transcriptome and genome for salmon index. 
# As an example we are downloading the gencode mouse reference

# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.transcripts.fa.gz
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz


# Installing Salmon
# see Linux_Svvr_cds.sh


# Preparing metadata
# Salmon indexing requires the names of the genome targets, which is extractable by using the grep command:

# grep "^>" <(gunzip -c GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
# sed -i.bak -e 's/>//g' decoys.txt

# Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index. 
# NOTE: the genome targets (decoys) should come after the transcriptome targets in the reference

# cat gencode.vM23.transcripts.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome.fa.gz


# Salmon Indexing
# We have all the ingredients ready for the salmon recipe. We can run salmon indexing step as follows:
# NOTE: --gencode flag is for removing extra metdata in the target header separated by | from the gencode reference. You can skip it if using other references.

salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode
echo "job finished - index is made"
