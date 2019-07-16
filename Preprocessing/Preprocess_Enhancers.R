
# Add path to input directory:
input_folder <- ""

# Select only the enhancers from the ChromHMM data of Wilderman et al
# ChromHMM Data can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97752
# Raw files are in .bigbed format, which can be reformatted to .bed files using UCSC bigBedToBed 
Wilderman_folder <- paste(input_folder, "/Common_data/ChromHMM/GSE97752_Imputed_Chromatin_State_Segmentation/", sep = "")
output_folder <-  paste(input_folder,"/Common_data/Enhancers/", sep = "")
for(input_file in list.files(Wilderman_folder)){
  print(paste("# Reading: ", Wilderman_folder, input_file, sep = ""))
  chromhmm <- read.delim(paste(Wilderman_folder, input_file, sep = ""), header = F, stringsAsFactors = F)
  sample_name <- unlist(strsplit(input_file, split = "-"))[1]
  summary(factor(chromhmm$V4))
  # there are several enhancer "categories" in the 25-state chromhmm data. We only use EnhA1, EnhA2 and EnhAF.
  # info about 25-state model: https://egg2.wustl.edu/roadmap/web_portal/imputed.html
  chromhmm_enhancers <- chromhmm[which(chromhmm$V4 %in% c("13_EnhA1", "14_EnhA2", "15_EnhAF")),]
  chromhmm_enhancers <- chromhmm_enhancers[,c(1:3)]
  chromhmm_enhancers[,4] <- paste(sample_name, 1:nrow(chromhmm_enhancers), sep = "_")
  chromhmm_enhancers[,5] <- 1
  
  output_file <- paste(output_folder, "Enhancers_Palate_", sample_name, ".bed", sep = "")
  print(output_file)
  write.table(x = chromhmm_enhancers, file = output_file, sep = "\t", quote = F, row.names = F, col.names = F)
}
