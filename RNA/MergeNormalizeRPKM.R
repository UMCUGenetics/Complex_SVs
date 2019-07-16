###### MERGE.R

### REPLACE XXX with input directory
bam_folder <- "XXX"
exon_file <- "XXX/Homo_sapiens.GRCh37.74_exon_gene_sizes.txt" # path to "Homo_sapiens.GRCh37.74_exon_gene_sizes.txt"
output_folder <- "XXX"


dirs <- list.dirs(path=bam_folder,recursive=F,full.names=F)
dirs <- dirs[grep("^MP.*", dirs)]

nr_cols=0
sampledir = ''
for ( dir in dirs ){
  sample=basename(dir)
  if ( (sample != 'read_counts') && (sample != 'bamMetrics') && (sample != 'logs') && (sample != 'tmp' ) && (sample != 'jobs') && (sample != 'DEanalysis') && (sample != 'variantCalling') ){
    nr_cols = nr_cols+1
    sampledir = dir
    print(sample)
  }
}
 
nr_rows <- nrow(read.table(paste(bam_folder ,sampledir, '/read_counts/', basename(sampledir), '_htseq_counts.txt', sep=""),sep="\t",header=F))-5

output <- matrix(ncol=nr_cols+1, nrow=nr_rows)
col_count=1
samplenames=c('gene')

for ( dir in dirs ) {
  sample=basename(dir)
  if ( (sample != 'read_counts') && (sample != 'bamMetrics') && (sample != 'jobs') && (sample != 'logs') && (sample != 'tmp' ) && (sample != 'DEanalysis') && (sample != 'variantCalling') ) {
    samplenames <- append(samplenames, as.character(sample))
    countsfile=paste(bam_folder,sample,'/read_counts/',sample,'_htseq_counts.txt',sep="")
    counts=read.table(countsfile,sep="\t",header=F)
    colnames(counts) <- c('gene',sample)
    genes<-counts$gene
    
    if ( col_count == 1 ){
      output[,col_count] <- as.character(genes[1:(length(genes)-5)])
      col_count=col_count+1
    }
    
    output[,col_count] <- as.numeric(counts[1:(nrow(counts)-5),2])
    col_count=col_count+1
  }
}

merged_table <- as.data.frame(output)
colnames(merged_table) <- samplenames

outfile1=paste(bam_folder, "read_counts/Merged_BAMS_readCounts_raw.txt", sep = "")
dir.create(paste(bam_folder, "read_counts/",sep = ""), showWarnings = F)
write.table(merged_table,file=outfile1,row.names=F,col.names=T,quote=F,sep="\t")

###### Normalize.R #######
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")
library("DESeq")


#normalize raw read counts
counts=read.delim(outfile1, header=T, row.names=1)

conds <- factor(c(colnames(counts)))
cds <- newCountDataSet( counts, conds )
cds <- estimateSizeFactors( cds )

normalized_counts <- counts(cds, normalized=TRUE)

outfile2= paste(output_folder, "Bam_files_readCounts_normalized.txt", sep = "")
write.table(normalized_counts, file=outfile2, row.names=T,col.names=NA, quote=F,sep="\t")


###### RPKM.R 
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)

exon_gene_sizes <- read.table(file = exon_file,sep="\t",header=F,row.names=1)
raw_read_counts <- read.table(paste(bam_folder, "/read_counts/", list.files(path = paste(input_folder,"/read_counts", sep = ""), pattern = "readCounts_raw.txt"), sep = ""), sep="\t",header=T,row.names=1)

nrsamples <- ncol(raw_read_counts)
nrrows <- nrow(raw_read_counts)

tab <- matrix(data=NA, nrow=nrrows, ncol=nrsamples)

for (j in 1:nrsamples){
  RPKM = rpkm(raw_read_counts[j], exon_gene_sizes, normalized.lib.sizes=F, log=F)
  for(i in 1:nrrows){
    tab[i,j] = RPKM$V2[i]
  }
}

df <- data.frame(tab)
colnames(df) <- colnames(raw_read_counts)
rownames(df) <- rownames(raw_read_counts)

outfile <- paste(bam_folder, "read_counts/Merged_BAMS_readCounts_RPKM.txt", sep = "")
dir.create(paste(bam_folder, "read_counts/",sep = ""), showWarnings = F)
              
write.table(df, file=outfile, row.names=T,col.names=NA, quote=F,sep="	")
