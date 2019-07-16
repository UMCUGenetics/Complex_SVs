
library(diagram)
library(biomaRt)
library(GenomicRanges)
library(RColorBrewer)
library(ggbio)


#input_folder <- "/hpc/cog_bioinf/cuppen/project_data/Complex_svs/"

print("## Reading input data for breakpoint browser, this may take a while")
print("# Reading genelist")
genelist <- read.delim(paste(input_folder, "Data/Genes/Genelist_Protein_Coding.txt", sep = ""), stringsAsFactors = F)
genes_g <- GRanges(seqnames = genelist$chromosome_name, IRanges(start = genelist$start_position, end= genelist$end_position), strand= genelist$strand, hgnc_symbol = genelist$hgnc_symbol,
                   ensembl_gene_id =genelist$ensembl_gene_id)

# all_exons <- getBM(attributes = c("ensembl_gene_id", "ensembl_exon_id","chromosome_name","exon_chrom_start", "exon_chrom_end"), 
#                    filters = c("ensembl_gene_id"), values = genelist$ensembl_gene_id, mart = ensembl_hg19)

all_exons <- read.delim(paste(input_folder, "Data/Genes/Exons_protein_coding.txt", sep = ""))

# gaps are used by the Plot_ideogram_derivative_chromosomes script (to plot the centromeres)
gaps <- read.delim(paste(input_folder, "Data/Genome/hg19_gaps.txt", sep = ""), stringsAsFactors = F, check.names = F)


SVs <- read.delim(paste(input_folder, "Data/SV/Denovo_SVs.txt", sep = ""), stringsAsFactors = F)

print("# Reading PCHiC")
PCHiC_data_list <- list()
PCHiC_data_list[["NEC"]] <- read.delim(paste(input_folder, "Data/PCHiC/PCHiC_NEC.txt", sep = ""), stringsAsFactors = F)

print("# Reading Enhancers")

enhancer_folder <- paste(input_folder, "Data/Enhancers/", sep = "")
enhancer_list <- list()
for(enhancer_file in list.files(enhancer_folder, pattern = ".bed|.txt")){
  print(enhancer_file)
  
  enhancer_filepath <- paste(enhancer_folder, enhancer_file, sep = "")
  
  enhancers_celltype <- read.delim(enhancer_filepath, header = F, stringsAsFactors = F)
  enhancer_list[[enhancer_file]] <- enhancers_celltype
}
print("# Reading TADs")
# read_TADs <- function(Input_folder = paste(input_folder, "Data/TADs/", sep = "")){
#   
#   TADs <- list()
#   
#   for(cell_type in list.files(Input_folder, full.names = T, pattern = "TAD")){
#     print(cell_type)
#     TAD_cell <- read.delim(cell_type, header = F, stringsAsFactors = F)
#     TADs[[cell_type]] <- TAD_cell
#   }
#   return(TADs)
# }

TADs_raw <- read.delim(paste( paste(input_folder, "Data/TADs/TADs.txt", sep = "")), stringsAsFactors = F, header = T)

#TADs_raw <- read_TADs()
print("## Finished reading breakpoint browser input files")

