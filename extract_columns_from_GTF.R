# https://stackoverflow.com/questions/37033238/parse-gtf-file-from-gencode

# cd /scratch/users/mali/Indices/salmon_mm9_index/

# download the gtf:
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz

# remove header lines with comments manually by using vi

library(rtracklayer)

# As I am interested in extracting "transcript_id" and "gene_id", mention them in tags
gtf <- readGFF("/scratch/users/mali/Indices/salmon_mm9_index/gencode.vM23.annotation.gtf", version=2L, tags = c("transcript_id", "gene_id"))

dim(gtf) # 1868204      12

head(gtf, 2)

gtf <- gtf[,c(9:10)]

any(is.na(gtf)) # TRUE

gtf <- na.omit(gtf)
dim(gtf) # 1812869       2

gtf <- unique(gtf)

dim(gtf) # 142351      2

write.table(gtf, file="txp2gene.tsv", sep="\t", row.names=F, col.names=F, quote=F)


#### Getting a gene name mapping file from Ensembl transcript and gene ids:

gtf <- readGFF("/scratch/users/mali/Indices/salmon_mm9_index/gencode.vM23.annotation.gtf", version=2L, tags = c("transcript_id", "gene_id", "gene_name"))
dim(gtf) # 1868204      12
head(gtf, 2)
gtf <- gtf[,c(9:11)]
any(is.na(gtf))
gtf <- na.omit(gtf)
dim(gtf) # 1812869       3
gtf <- unique(gtf)
dim(gtf) # 142351      3
write.table(gtf, file="EnsT_EnsG_GeneName.tsv", sep="\t", row.names=F, quote=F)
