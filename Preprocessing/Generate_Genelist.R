## This script can be used to obtain gene information from Ensembl using biomart (hg19)
## Raw files containing phenotypic information (pLI, RVIS, HPO, DDD etc) are parsed
## This information is annotated with this phenotypic information

## Important: check the filenames at lines 17 to 32!
# The hg19 biomart is sometimes offline and this script will fail. 

args <- commandArgs(trailingOnly=TRUE)
# The required input files should be stored in [input_folder]/Phenotypes/[]
input_folder <- args[1]
# The output files (Genelists) will be written to this folder:
output_folder <- args[2]

library(biomaRt)
library(GenomicRanges)

# HPO file can be downloaded from http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/
HPO_file <- paste(input_folder, "Phenotypes/HPO/20190212_ALL_SOURCES_FREQUENT_FEATURES_genes_to_phenotype.txt", sep = "")

# pLI scores can be downloaded from: ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint
pLI_file <- paste(input_folder, "Phenotypes/pLI/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt", sep = "")

# RVIS scores can be downloaded from http://genic-intolerance.org/data/RVIS_Unpublished_ExACv2_March2017.txt
RVIS_file <- paste(input_folder, "Phenotypes/RVIS/RVIS_Unpublished_ExACv2_March2017.txt", sep = "")

# Pathogenicity scores from Redin et al can be found in Supplemental Table 6 - doi:10.1038/ng.3720
redin_file <- paste(input_folder, "Phenotypes/Redin_Gene_Annotation.txt", sep = "")

# haploinsufficiency indices can be downloaded from: https://decipher.sanger.ac.uk/files/downloads/HI_Predictions_Version3.bed.gz
haploinsufficiency_file <- paste(input_folder, "Phenotypes/HI_Predictions_Version3.bed", sep = "")

# "curated list of genes reported to be associated with developmental disorders", downloaded from: https://decipher.sanger.ac.uk/about#downloads/data
DDG2P_file <- paste(input_folder, "Phenotypes/DDD/DDG2P_19_2_2019.csv", sep = "")

# Biomart will look for genes on these chromosomes:
chromosomes <- c(1:22, "X", "Y", "MT")

# Download gene information using biomaRt 
#ensembl_hg19 <- useMart(host='feb2014.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
ensembl_hg19 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl") # this one is often offline...
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

print("# Obtaining gene information from grch37.ensembl.org using biomaRt")

filters <- listFilters(ensembl_hg19)
attributes <- listAttributes(ensembl_hg19)

All_genes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "hgnc_id", "entrezgene","chromosome_name", "start_position", "end_position", 
                                  "strand"), 
                   filters = c("chromosome_name"), values = chromosomes, mart = ensembl_hg19)

# Retrieve OMIM and new hgnc_symbols from new ensembl mart (not hg19):
All_genes2 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "hgnc_id", "phenotype_description"), 
                    filters = c("ensembl_gene_id"), values = All_genes$ensembl_gene_id, mart = ensembl)
names(All_genes2) <- c("ensembl_gene_id", "hgnc_symbol_new", "hgnc_id", "phenotype_description")

print("# Adding phenotypic information to the genelist")

# Some genes have multiple phenotype_descriptions, which each are shown on a different row. 
# Here the phenotype_descriptions for each gene are paste together to one OMIM description per gene
OMIM <- aggregate(All_genes2$phenotype_description ~ All_genes2$ensembl_gene_id, FUN=function(x){paste(unique(x), collapse = ";")})
names(OMIM) <- c("ensembl_gene_id", "OMIM")

All_genes_merged <- merge(All_genes, OMIM, by = "ensembl_gene_id", all.x = T)

All_genes_merged <- merge(All_genes_merged, All_genes2[!duplicated(All_genes2$ensembl_gene_id),c("ensembl_gene_id", "hgnc_symbol_new")], all.x = T)

# Remove genes that are not located on chromosomes 1-22, X, Y or MT
genes_filtered <- All_genes_merged[which(All_genes_merged$chromosome_name %in% chromosomes),]
All_genes_merged$ensembl_gene_id <- factor(genes_filtered$ensembl_gene_id)

# Read the pLI, RVIS and Happloinsuffiency score files
if(file.exists(pLI_file)){
  pLI <- read.delim(pLI_file,header = T)
  names(pLI)[names(pLI) == "gene"] <- "hgnc_symbol"
  genes_phenotype <- merge(All_genes_merged, pLI[,c("hgnc_symbol","pLI")], by = "hgnc_symbol", all.x = T)
} else {
  print(paste("# Error! pLI file not found: ", pLI_file, sep = ""))
  genes_phenotype$pLI <- NA
}

if(file.exists(RVIS_file)){
  RVIS <- read.delim(RVIS_file, check.names = F, stringsAsFactors = F)
  RVIS_filtered <- RVIS[, c("CCDSr20", "%RVIS[pop_maf_0.05%(any)]")]
  names(RVIS_filtered) <- c("hgnc_symbol", "RVIS")
  
  # RVIS uses the new hgnc_symbols, while this script uses the hg19 hgnc_symbols. Translate the hgnc_symbols to entrezgene ids and merge on entrezgene:
  RVIS_ID <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), 
                   filters = c("hgnc_symbol"), values = as.vector(RVIS_filtered$hgnc_symbol), mart = ensembl)
  RVIS_ID <- RVIS_ID[!duplicated(RVIS_ID$ensembl_gene_id),]
  RVIS_filtered <- merge(RVIS_filtered, RVIS_ID, by = "hgnc_symbol", all.x = T)
  
  # First merge RVIS based on hgnc_symbol and subsequently on ensembl_gene_id:
  RVIS_temp <- merge(genes_phenotype, RVIS_filtered[,c("hgnc_symbol", "RVIS")], by = "hgnc_symbol", all.x = T)
  RVIS_temp <- merge(RVIS_temp, RVIS_filtered[,c("ensembl_gene_id", "RVIS")], by = "ensembl_gene_id", all.x = T)
  RVIS_temp <- RVIS_temp[!duplicated(RVIS_temp$ensembl_gene_id),]
  RVIS_temp$RVIS.x[is.na(RVIS_temp$RVIS.x)] <- 0
  RVIS_temp$RVIS.y[is.na(RVIS_temp$RVIS.y)] <- 0
  RVIS_temp$RVIS <- ifelse(RVIS_temp$RVIS.x == RVIS_temp$RVIS.y, RVIS_temp$RVIS.x, RVIS_temp$RVIS.x+RVIS_temp$RVIS.y)
  genes_phenotype <- merge(genes_phenotype, RVIS_temp[,c("ensembl_gene_id", "RVIS")], by =  "ensembl_gene_id", all.x = T)
  genes_phenotype$RVIS[which(genes_phenotype$RVIS == 0)] <- NA
} else {
  print(paste("# Error! RVIS file not found: ", RVIS_file, sep = ""))
  genes_phenotype$RVIS <- NA
}

if(file.exists(haploinsufficiency_file)){
  HI <- read.delim(haploinsufficiency_file, header = F, skip = 1, stringsAsFactors = F)
  HI_genes <- stringr::str_split_fixed(string = HI$V4, pattern = "\\|", n = 4) # gene names and HI scores are in the fourth column
  HI$hgnc_symbol <- HI_genes[,1]
  HI$HI <- HI_genes[,3]
  HI$HI <- gsub(x = HI$HI, pattern = "%", replacement = "")
  # HI is also based on newer hgnc symbols. Therefore add ensembl_gene_id ids and merge on these:
  HI_ID <- getBM(attributes = c("hgnc_symbol","hgnc_id","ensembl_gene_id"),
                 filters = c("hgnc_symbol"), values = as.vector(HI$hgnc_symbol), mart = ensembl)
  
  HI_new <- merge(HI, HI_ID, by = "hgnc_symbol", all.x = T)
  
  HI_ID_hg19 <- getBM(attributes = c("hgnc_symbol","hgnc_id","ensembl_gene_id"),
                      filters = c("hgnc_symbol"), values = HI_new$hgnc_symbol[is.na(HI_new$ensembl_gene_id)], mart = ensembl_hg19)
  HI_hg19 <- merge(HI, HI_ID_hg19, by = "hgnc_symbol")
  HI_new <- HI_new[!is.na(HI_new$ensembl_gene_id),]
  HI <- rbind(HI_new,HI_hg19)
  genes_phenotype <- merge(genes_phenotype, HI[,c("ensembl_gene_id","HI")], by =  "ensembl_gene_id", all.x = T)
} else {
  print(paste("# Error! HI file not found: ", haploinsufficiency_file, sep = ""))
  genes_phenotype$HI <- NA
}

if(file.exists(DDG2P_file)){
  DDG2P <- read.csv(DDG2P_file, check.names = F)
  names(DDG2P)[names(DDG2P) == "gene symbol"] <- "hgnc_symbol"
  names(DDG2P)[names(DDG2P) == "hgnc id"] <- "hgnc_id"
  names(DDG2P)[names(DDG2P) == "DDD category"] <- "DDG2P"
  
  genes_phenotype <- merge(genes_phenotype, DDG2P[,c("hgnc_id","DDG2P")], by = "hgnc_id", all.x = T)
} else {
  print(paste("# Error! DDG2P file not found: ", DDG2P_file, sep = ""))
  genes_phenotype$DDG2P <- NA
  stop()
}

if(file.exists(HPO_file)){
  # Read the HPO data:
  HPO_raw <- read.delim(HPO_file, skip = 1, header = F)
  names(HPO_raw) <- c("entrezgene", "Entrez_gene_name", "HPO_Term", "HPO_Term_ID")
  HPO_raw$entrezgene <- factor(HPO_raw$entrezgene)
  
  # Each HPO term for each gene is on a seperate row. Merge and count all HPO terms per gene
  all_HPO <- data.frame()
  for(gene in levels(HPO_raw$entrezgene)){
    gene_HPO <- HPO_raw[HPO_raw$entrezgene == gene,]
    Number_HPO_Terms_Gene <- nrow(gene_HPO)
    merged_HPO_term <- paste(gene_HPO$HPO_Term, collapse = ",")
    merged_HPO_term_ID <- paste(gene_HPO$HPO_Term_ID, collapse = ",")
    
    gene <- data.frame(entrezgene = unique(gene_HPO$entrezgene), 
                       Entrez_gene_name = unique(gene_HPO$Entrez_gene_name), 
                       HPO_Terms = merged_HPO_term, 
                       HPO_Term_IDs = merged_HPO_term_ID,
                       Number_HPO_Terms_Gene = Number_HPO_Terms_Gene)
    
    all_HPO <- rbind(all_HPO,gene)
  }
  # Add ensembl_gene_ids to the genes with HPO terms (not all genes have entrezgene ids):
  HPO_ensembl <- getBM(attributes = c("ensembl_gene_id","entrezgene"), 
                       filters = c("entrezgene"), values = all_HPO$entrezgene, mart = ensembl)
  all_HPO2 <- merge(all_HPO, HPO_ensembl, by = "entrezgene", all.x = T)
  genes_phenotype <- merge(genes_phenotype, all_HPO2[,c("HPO_Terms","HPO_Term_IDs","Number_HPO_Terms_Gene","ensembl_gene_id")], by = "ensembl_gene_id", all.x = T)
  genes_phenotype$Number_HPO_Terms_Gene[is.na(genes_phenotype$Number_HPO_Terms_Gene)] <- 0
} else {
  print(paste("# Error! HPO file not found: ", HPO_file, sep = ""))
  genes_phenotype$HPO_Terms <- NA
  genes_phenotype$HPO_Term_IDs <- NA
  genes_phenotype$Number_HPO_Terms_Gene <- NA
}


if(file.exists(redin_file)){
  # Add the pathogenicity scores of Redin et al to the genelist:
  redin <- read.delim(redin_file, header = F, stringsAsFactors = F, check.names = F)
  names(redin) <- c("hgnc_symbol", "redin_score", "redin_source")
  redin$redin_score <- factor(redin$redin_score, levels = c("1", "2", "3","[3-10[", "[10-20[", "20+"))
  
  # bin the scores to numeric values:
  levels(redin$redin_score)[levels(redin$redin_score) == "[10-20["] <- 10
  levels(redin$redin_score)[levels(redin$redin_score) == "[3-10["] <- 4
  levels(redin$redin_score)[levels(redin$redin_score) == "20+"] <- 20
  levels(redin$redin_score) <- as.numeric(levels(redin$redin_score))
  
  genes_phenotype <- merge(genes_phenotype, redin[,c("hgnc_symbol", "redin_score")], by = "hgnc_symbol", all.x = T)
  genes_phenotype$redin_score <- as.vector(genes_phenotype$redin_score)
  genes_phenotype$redin_score[is.na(genes_phenotype$redin_score)] <- 0
  genes_phenotype$redin_score <- factor(genes_phenotype$redin_score)
} else {
  print(paste("# Error! Redin file not found: ", redin_file, sep = ""))
  genes_phenotype$redin_score <- NA
}

print("# Obtaining transcript information from grch37.ensembl.org using biomaRt")

## add transcription start sites of the longest transcript per gene:
# Add refseq_ids:
All_transcripts <- getBM(attributes = c("ensembl_gene_id","refseq_mrna","chromosome_name", "transcription_start_site", "transcript_start", "transcript_end", "transcript_length",
                                        "strand"), 
                         filters = c("ensembl_gene_id"), values = unique(as.vector(genes_phenotype$ensembl_gene_id)), mart = ensembl_hg19)

genes_phenotype <- merge(genes_phenotype, All_transcripts[,c("ensembl_gene_id","refseq_mrna", "transcription_start_site", "transcript_start", "transcript_end", "transcript_length")], 
                         by = "ensembl_gene_id", all.x = T)

# Collapse the refseq_mrna IDs to one line per ensembl_gene_id:
genes_refseq <- genes_phenotype[which(genes_phenotype$refseq_mrna != ""),]

refseq_mrna_overview <- aggregate(genes_refseq$refseq_mrna ~ genes_refseq$ensembl_gene_id, FUN=function(x){paste(unique(x), collapse = ";")})
names(refseq_mrna_overview) <- c("ensembl_gene_id", "refseq_mrnas")

genes_phenotype <- merge(genes_phenotype, refseq_mrna_overview, by = "ensembl_gene_id", all.x = T)

# Add inheritance mode from HPO and DD2GP
print("# Adding information on the mode-of-inheritance to the genelist")

inheritance_mode <- data.frame(hgnc_symbol = genes_phenotype$hgnc_symbol[genes_phenotype$hgnc_symbol != ""], HPO_Terms = genes_phenotype$HPO_Terms[genes_phenotype$hgnc_symbol != ""], stringsAsFactors = F)
inheritance_mode <- inheritance_mode[!duplicated(inheritance_mode$hgnc_symbol),]
inheritance_mode <- merge(inheritance_mode, DDG2P[,c("hgnc_symbol", "allelic requirement")], all.x = T, by = "hgnc_symbol")

inheritance_mode$Inheritance <- ifelse(grepl(pattern = "Autosomal dominant inheritance", x = inheritance_mode$HPO_Terms) == TRUE,
                                 "AD", "ND")
inheritance_mode$Inheritance <- ifelse(grepl(pattern = "X-linked recessive", x = inheritance_mode$HPO_Terms) == TRUE,
                                 "XR", inheritance_mode$Inheritance )
inheritance_mode$Inheritance <- ifelse(grepl(pattern = "X-linked dominant inheritance", x = inheritance_mode$HPO_Terms) == TRUE,
                                 "XD", inheritance_mode$Inheritance )
inheritance_mode$Inheritance <- ifelse(grepl(pattern = "X-linked dominant inheritance", x = inheritance_mode$HPO_Terms) == TRUE & 
                                         grepl(pattern = "X-linked recessive", x = inheritance_mode$HPO_Terms) == TRUE,
                                       "XD+XR", inheritance_mode$Inheritance )

inheritance_mode$Inheritance <- ifelse(grepl(pattern = "Autosomal recessive inheritance", x = inheritance_mode$HPO_Terms) == TRUE,
                                 "AR", inheritance_mode$Inheritance )
inheritance_mode$Inheritance <- ifelse(grepl(pattern = "Autosomal dominant inheritance", x = inheritance_mode$HPO_Terms) == TRUE &
                                   grepl(pattern = "Autosomal recessive inheritance", x = inheritance_mode$HPO_Terms) == TRUE,
                                 "AD+AR", inheritance_mode$Inheritance )
summary(factor(inheritance_mode$Inheritance))
summary(factor(inheritance_mode$`allelic requirement`[inheritance_mode$Inheritance == "ND"]))
# Some inheritance modes are missing in HPO or are different in DD2GP. Add these from DD2GP:
inheritance_mode$Inheritance[inheritance_mode$Inheritance == "AD" & inheritance_mode$`allelic requirement` == "biallelic" & !is.na(inheritance_mode$`allelic requirement`)] <- "AD+AR"
inheritance_mode$Inheritance[inheritance_mode$Inheritance == "AR" & inheritance_mode$`allelic requirement` == "monoallelic" & !is.na(inheritance_mode$`allelic requirement`)] <- "AD+AR"
inheritance_mode$Inheritance[inheritance_mode$Inheritance == "XR" & inheritance_mode$`allelic requirement` == "x-linked dominant" & !is.na(inheritance_mode$`allelic requirement`)] <- "XD+XR"
inheritance_mode$Inheritance[inheritance_mode$Inheritance == "XD" & inheritance_mode$`allelic requirement` == "hemizygous" & !is.na(inheritance_mode$`allelic requirement`)] <- "XD+XR"

inheritance_mode$Inheritance[inheritance_mode$Inheritance == "ND" & inheritance_mode$`allelic requirement` == "monoallelic" & !is.na(inheritance_mode$`allelic requirement`)] <- "AD"
inheritance_mode$Inheritance[inheritance_mode$Inheritance == "ND" & inheritance_mode$`allelic requirement` == "biallelic" & !is.na(inheritance_mode$`allelic requirement`)] <- "AR"
inheritance_mode$Inheritance[inheritance_mode$Inheritance == "ND" & inheritance_mode$`allelic requirement` == "hemizygous" & !is.na(inheritance_mode$`allelic requirement`)] <- "XR"
inheritance_mode$Inheritance[inheritance_mode$Inheritance == "ND" & inheritance_mode$`allelic requirement` %in% c("x-linked dominant", "x-linked over-dominance") & !is.na(inheritance_mode$`allelic requirement`)] <- "XR"

summary(factor(inheritance_mode$Inheritance))

genes_phenotype <- merge(genes_phenotype, inheritance_mode[,c("hgnc_symbol","Inheritance")], by = "hgnc_symbol", all.x = T)

# Order the genelist. Only the first entry of each gene will be included in the final genelist. This will be the longest transcript:
print("# Sorting the genelist")
genes_ordered <- genes_phenotype[order(genes_phenotype$transcript_length, genes_phenotype$refseq_mrna, genes_phenotype$hgnc_symbol, 
                                       genes_phenotype$entrezgene, decreasing = T),]

# Sort the columns:
genes_ordered <- genes_ordered[,c("hgnc_symbol","hgnc_symbol_new","ensembl_gene_id","refseq_mrna","refseq_mrnas","entrezgene",
                                  "chromosome_name","start_position","end_position","strand","transcription_start_site","transcript_start","transcript_end","transcript_length",
                                  "pLI","RVIS","HI","DDG2P","OMIM", "HPO_Terms", "HPO_Term_IDs" , "Number_HPO_Terms_Gene","redin_score", "Inheritance")]

# Remove all the genes that do not have a hgnc_symbol:
all_hgnc_genes <- genes_ordered[which(genes_ordered$hgnc_symbol != ""),]
all_hgnc_genes <- all_hgnc_genes[!duplicated(all_hgnc_genes$ensembl_gene_id),] 
all_hgnc_genes <- all_hgnc_genes[!duplicated(all_hgnc_genes$hgnc_symbol),]
all_hgnc_genes <- all_hgnc_genes[order(all_hgnc_genes$chromosome_name, all_hgnc_genes$start_position),]
print(paste("  Total number of genes with HGNC symbol: ", nrow(all_hgnc_genes), sep = ""))

# Remove all the genes that do not have a refseq_mrna or a hgnc_symbol:
all_protein_coding_genes <- genes_ordered[which(genes_ordered$refseq_mrna != ""),]
all_protein_coding_genes <- all_protein_coding_genes[which(all_protein_coding_genes$hgnc_symbol != ""),]
all_protein_coding_genes <- all_protein_coding_genes[!duplicated(all_protein_coding_genes$ensembl_gene_id),] 
all_protein_coding_genes <- all_protein_coding_genes[!duplicated(all_protein_coding_genes$hgnc_symbol),]
all_protein_coding_genes <- all_protein_coding_genes[order(all_protein_coding_genes$chromosome_name, all_protein_coding_genes$start_position),]
print(paste("  Total number of protein coding genes with Refseq mRNA ID: ", nrow(all_protein_coding_genes), sep = ""))
      
# Select all Refseq transcripts:
all_transcripts <- genes_ordered[which(genes_ordered$refseq_mrna != ""),]
all_transcripts <- all_transcripts[which(all_transcripts$hgnc_symbol != ""),]
all_transcripts <- all_transcripts[!duplicated(all_transcripts$refseq_mrna),] 
all_transcripts <- all_transcripts[order(all_transcripts$chromosome_name, all_transcripts$start_position),]

# Write the genelists to txt files:
print( paste("# Writing: ", output_folder, "Genelist_Protein_Coding.txt", sep = ""))
write.table(x = all_protein_coding_genes, 
            file = paste(output_folder, "Genelist_Protein_Coding.txt", sep = ""),sep = "\t", quote = F, row.names = F)

print( paste("# Writing: ", output_folder, "Genelist_HGNC.txt", sep = ""))
write.table(x = all_hgnc_genes, 
            file = paste(output_folder, "Genelist_HGNC.txt", sep = ""), sep = "\t", quote = F, row.names = F)

print( paste("# Writing: ",output_folder, "Genelist_refseq_transcripts.txt", sep = ""))
write.table(x = all_transcripts, 
            file = paste(output_folder, "Genelist_refseq_transcripts.txt", sep = ""),sep = "\t", quote = F, row.names = F)


# Retrieve all exons coordinates (mainly for plotting):
all_exons <- getBM(attributes = c("ensembl_gene_id", "ensembl_exon_id","chromosome_name","exon_chrom_start", "exon_chrom_end"), 
                   filters = c("ensembl_gene_id"), values = all_protein_coding_genes$ensembl_gene_id, mart = ensembl_hg19)

print( paste("# Writing: ", output_folder, "Exons_protein_coding.txt", sep = ""))
write.table(x = all_exons, 
            file = paste(output_folder, "Exons_protein_coding.txt", sep = ""),sep = "\t", quote = F, row.names = F)

