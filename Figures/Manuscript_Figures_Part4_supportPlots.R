## This script is used to generate manuscript figures based on breakpoint browser
# Fig 2: HOXD example
# Fig S2: MP3738 SOX3 insertion
# Fig S4 A,B,C,D: Workflow 3D genome datasets > FOXG1
# Fig S6: SOX9 example
# Fig ??: FOXG1 example

#source("https://bioconductor.org/biocLite.R")
#BiocInstaller::biocLite(c("diagram"))

library(reshape2)

## Add path to input directory here:
input_folder <- ""

source( paste(input_folder, "Scripts/Visualization/Breakpoint_browser.R", sep = ""))
source( paste(input_folder, "Scripts/Visualization/read_browser_data.R", sep = ""))
source( paste(input_folder, "Scripts/Visualization/Plot_ideograms_derivative_chromosomes.R", sep = ""))

## Add path to ini file here:
ini_file <- ""

if(file.exists(ini_file)){
  ini <- read.delim(ini_file, header = T, stringsAsFactors = F, comment.char = "#")
  project_folder <- ini[which(ini$ID == "Project_folder"), "Value"]
} else {
  stop("Could not open ini file: ", ini_file)
}

SV_file <- paste(input_folder, ini[which(ini$ID == "SV_file"), "Value"], sep = "")
SVs <- read.delim(SV_file, stringsAsFactors = F)
Genes_SVs <- read.delim(paste(input_folder, "/Results/Integration/MP/Genes_SVs.txt", sep= ""))
# supplental figure S3A

## Fig 2: HOXD
conf <- read.delim(paste(input_folder,"/Scripts/Visualization/Conf_files/HOXD.conf", sep = ""), stringsAsFactors = F)
pdf(file = paste(input_folder,"/Results/Figures_raw/2B_HOXD.pdf", sep = ""), width = 4, height = 3, pointsize = 8, family = "ArialMT")

## add DHS
# check enhancers / expression hox cluster (maybe mouse enhancers)

breakpoint_plot(conf = conf, 
                patient = "UTR22", 
                chromosome = 2,
                xlim_plot = c(176e6,177.5e6), 
                Genes_SVs = "")
# Highlight the fragment translocated to a different derivative chromosome: 
rect(xleft = 176643961, xright = 176751445, ybottom = 0, ytop = 8.5, col = rgb(255,69,0, 50, maxColorValue = 255), border = NA)
# Highlight the inverted fragment
rect(xleft = 176751445, xright = 176857274, ybottom = 0, ytop = 8.5, col = rgb(255,204,0, 50, maxColorValue = 255), border = NA)
# Highlight the HOXD locus
rect(xleft = 176957532, xright = 177055688, ybottom = 0, ytop = 8.5, col = rgb(0,200,0, 50, maxColorValue = 255), border = NA)

# Plot the HOXD enhancers specified by Rodríguez-Carballo et al, Genes & Development (2017). mm10 coordinates mapped to hg19 using liftover
HOXD_enhancers <- read.delim(paste(input_folder,"/Data/Enhancers/Raw/Limb/HOXD_Enhancers_mouse_hg19.bed", sep = ""), header = F, stringsAsFactors = F)
rect(xleft = HOXD_enhancers[,2], xright = HOXD_enhancers[,3], ybottom = 3.7, ytop = 4, col = "#4B0082", border = NA)
text(x = (HOXD_enhancers[,2]+HOXD_enhancers[,3])/2, y = 4.2, labels = as.character(HOXD_enhancers[,4]), cex = 0.8, col = "#4B0082")
mtext(text = "HOXD", side = 2, line = 0, at = 3.85, las = 1, cex = 0.8)

dev.off()



## Fig S3: Insertion superenhancers near SOX3 in individual MP3738
pdf(file = paste(input_folder,"/Manuscript/Fig_S2A_bottom.pdf", sep = ""), width = 5, height = 2.5, pointsize = 10)
SOX3_Conf <- read.delim(paste(input_folder, "Scripts/Visualization/Conf_files/SOX3_insertion.conf", sep = ""), stringsAsFactors = F)
breakpoint_plot(conf = SOX3_Conf, 
                patient = "MP3738", 
                chromosome = "derX",
                xlim_plot = c(139502953-1e6,139673059+1e6), 
                Genes_SVs = "", breakpoints = SVs)
rect(xleft = 139502953, xright = 139673059, ybottom = 4.3, ytop = 9, col = rgb(255,204,0, 50, maxColorValue = 255), border = NA)
axis(side = 2, at = c( as.numeric(SOX3_Conf[SOX3_Conf$Option == "Genes", "y"])-0.2, as.numeric(SOX3_Conf[SOX3_Conf$Option == "Genes", "y"])+1.5), labels = NA, tcl = 0.5, line = 4)
mtext(side = 2, text = "Genes", at = (as.numeric(SOX3_Conf[SOX3_Conf$Option == "Genes", "y"]) + as.numeric(SOX3_Conf[SOX3_Conf$Option == "Genes", "y"])+1.4)/2, line = 4.5, cex = 0.9)
dev.off()

pdf(file = paste(input_folder,"/Manuscript/Fig_S2A_top.pdf", sep = ""), width = 5, height = 2.5, pointsize = 10)
breakpoint_plot(conf = SOX3_Conf, 
                patient = "MP3738", 
                chromosome = 9,
                xlim_plot = c(16489097-1e6,16659203+1e6), 
                Genes_SVs = "")
axis(side = 2, at = c( as.numeric(SOX3_Conf[SOX3_Conf$Option == "Genes", "y"])-0.2, as.numeric(SOX3_Conf[SOX3_Conf$Option == "Genes", "y"])+1.5), labels = NA, tcl = 0.5, line = 4)
mtext(side = 2, text = "Genes", at = (as.numeric(SOX3_Conf[SOX3_Conf$Option == "Genes", "y"]) + as.numeric(SOX3_Conf[SOX3_Conf$Option == "Genes", "y"])+1.4)/2, line = 4.5, cex = 0.9)
rect(xleft = c(16091200, 16171600, 16501000, 16569200, 16852400), xright = c(16106800, 16289800, 16554600, 16629200, 16885400), 
     ybottom = 5.55, ytop = 5.95, border = "black", col = NA, lty = 1, lwd = 0.8)
rect(xleft = 16489097, xright = 16659203, ybottom = 2, ytop = 9, col = rgb(255,204,0, 50, maxColorValue = 255), border = NA)
dev.off()



# Fig S3A
plot_phenomatch_small <- function(Patient,
                            phenomatch_folder,
                            Genes_SVs = Gene_data,
                            SVs = SVs){
  
  phenomatch_file_patient <- paste(phenomatch_folder, "Phenomatch_", Patient, ".txt", sep = "")
  
  phenomatch_overview <- read.delim(phenomatch_file_patient, check.names = F, stringsAsFactors = F)
  
  patient_data <- Genes_SVs[which(Genes_SVs$Patient == Patient),]
  patient_data$gene_symbol <- patient_data$hgnc_symbol
  
  # Only select one row per gene (removes one of the truncated gene fragments)
  patient_data <- patient_data[!duplicated(patient_data$gene_symbol),]
  
  ## Not including SV type in small version of phenomatch plot
  # Phenomatch automatically looks at genes +/- 3Mb from the SVs and the output therefore contains genes we're not looking at.
  # Select the genes +/- 2Mb from the SVs
  
  phenomatch_overview_2Mb <- phenomatch_overview[which(phenomatch_overview$hgnc_symbol %in% as.vector(patient_data$hgnc_symbol)),]
  phenomatch_overview_2Mb$hgnc_symbol <- factor(phenomatch_overview_2Mb$hgnc_symbol, 
                                                levels = phenomatch_overview_2Mb$hgnc_symbol[order(phenomatch_overview_2Mb$phenoMatchScore, decreasing = F)])
  
  phenomatch_overview2 <- merge(phenomatch_overview_2Mb, patient_data[,c("hgnc_symbol", "Phenomatches_high","pLI","RVIS","HI","DDG2P", "OMIM","Inheritance", "Phenotypic_score")], by = "hgnc_symbol")

  phenomatch_overview_melted <- melt(phenomatch_overview2, id.vars = "hgnc_symbol")

  # The score will determine the fill color
  phenomatch_overview_melted$Score <- phenomatch_overview_melted$value

  phenomatch_overview_melted$Score <- as.numeric(phenomatch_overview_melted$Score)
  
  phenomatch_overview_melted$grid <- ifelse(phenomatch_overview_melted$variable == "phenoMatchScore", "phenoMatchScore", "Phenomatch")
  phenomatch_overview_melted$grid <- ifelse(phenomatch_overview_melted$variable == "Phenomatches", "phenoMatchScore", phenomatch_overview_melted$grid)
  phenomatch_overview_melted$grid <- ifelse(phenomatch_overview_melted$variable %in% c("pLI","RVIS","HI","DDG2P", "OMIM", "Inheritance"), "Other", phenomatch_overview_melted$grid)
  phenomatch_overview_melted$grid <- ifelse(phenomatch_overview_melted$variable == "Phenotypic_score", "Phenotypic_score", phenomatch_overview_melted$grid)
  
  phenomatch_overview_melted$grid <- factor(phenomatch_overview_melted$grid, levels = c("Phenomatch", "phenoMatchScore", "Other", "Phenotypic_score"))
  
  
  phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "Phenomatches"] <- as.numeric(phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "Phenomatches"]) / (ncol(phenomatch_overview_2Mb)-2) * 30
  phenomatch_overview_melted$value[phenomatch_overview_melted$variable == "Phenomatches"] <- 
    paste(phenomatch_overview_melted$value[phenomatch_overview_melted$variable == "Phenomatches"], (ncol(phenomatch_overview_2Mb)-2), sep = "/")
  phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "pLI"] <- ifelse(phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "pLI"] > 0.5, 
                                                                                           (phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "pLI"])*30 ,0)
  
  phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "RVIS"] <- ifelse(as.numeric(phenomatch_overview_melted$value[phenomatch_overview_melted$variable == "RVIS"]) < 50, 
                                                                                           (100 - as.numeric(phenomatch_overview_melted$value[phenomatch_overview_melted$variable == "RVIS"])*2)/100*30,0)
  
  phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "HI"] <- ifelse(as.numeric(phenomatch_overview_melted$value[phenomatch_overview_melted$variable == "HI"]) < 50, 
                                                                                            (100 - as.numeric(phenomatch_overview_melted$value[phenomatch_overview_melted$variable == "HI"])*2)/100*30,0)
  
  phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "Phenotypic_score"] <-   as.numeric(phenomatch_overview_melted$value[phenomatch_overview_melted$variable == "Phenotypic_score"]) / 5* 30
    
    ifelse(as.numeric(phenomatch_overview_melted$value[phenomatch_overview_melted$variable == "Phenotypic_score"]) < 50, 
                                                                                          (100 - as.numeric(phenomatch_overview_melted$value[phenomatch_overview_melted$variable == "HI"])*2)/100*30,0)
  
  phenomatch_overview_melted$value[phenomatch_overview_melted$variable == "DDG2P" &  phenomatch_overview_melted$value == "probable"] <- "Prob"
  phenomatch_overview_melted$value[phenomatch_overview_melted$variable == "DDG2P" &  phenomatch_overview_melted$value == "possible"] <- "Pos"
  phenomatch_overview_melted$value[phenomatch_overview_melted$variable == "DDG2P" &  phenomatch_overview_melted$value == "confirmed"] <- "Yes"
  
  phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "DDG2P"] <- 0
  phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "DDG2P" & phenomatch_overview_melted$value == "Prob"] <- 15
  phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "DDG2P" & phenomatch_overview_melted$value == "Yes"] <- 25
  
  phenomatch_overview_melted$value[phenomatch_overview_melted$variable == "OMIM" & phenomatch_overview_melted$value != ""] <- "Yes"
  phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "OMIM"  & phenomatch_overview_melted$value == "Yes"] <- 20
  
  phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "Inheritance" &  phenomatch_overview_melted$value %in% c("AD", "XD")] <- 20
  phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "Inheritance" &  phenomatch_overview_melted$value == "AD+AR"] <- 15
  phenomatch_overview_melted$Score[phenomatch_overview_melted$variable == "Inheritance" &  phenomatch_overview_melted$value == "AR"] <- 0
  
    
  phenomatch_overview_melted$Score <- ifelse(phenomatch_overview_melted$Score > 30, 30, phenomatch_overview_melted$Score)
  
  phenomatch_overview_melted$value[phenomatch_overview_melted$variable != "DDG2P" & 
                                     phenomatch_overview_melted$variable != "Inheritance" & 
                                     phenomatch_overview_melted$variable != "Phenomatches" &
                                     phenomatch_overview_melted$variable != "OMIM"] <- round(as.numeric(phenomatch_overview_melted$value[phenomatch_overview_melted$variable != "DDG2P" & 
                                                                                                                                           phenomatch_overview_melted$variable != "Inheritance" & 
                                                                                                                                           phenomatch_overview_melted$variable != "Phenomatches" &
                                                                                                                                           phenomatch_overview_melted$variable != "OMIM"] ), 1)
  
  names(phenomatch_overview_melted) <- c("Gene", "Phenotype", "phenoMatchScore", "Score", "grid")
  
  output_height <- nrow(phenomatch_overview) / 3 + 2
  
  #print(output_height)
  output_width <- ncol(phenomatch_overview) / 2 + 2
  #print(output_width)
  
  caption <- paste("Showing data for ", length(unique(phenomatch_overview_melted$Gene)), "/", length(unique(patient_data$hgnc_symbol)), " genes +/- 2 Mb from SVs", sep = "")
  
  bottommargin <-  ifelse(length(unique(phenomatch_overview_melted$Gene)) < 50, 20 - length(unique(phenomatch_overview_melted$Gene)) * 0.5, 0.5)
  
  data_table <- ggplot(phenomatch_overview_melted, aes(x = Phenotype, y = Gene, fill = Score, label = phenoMatchScore)) + 
    geom_tile(lwd = 0.4, colour  = "black") +
    geom_text(size = 2, data = phenomatch_overview_melted, aes(label = phenoMatchScore)) + 
    scale_fill_gradient(low ="white", high = "red", limits=c(0,30), na.value = "gray") +
    scale_x_discrete(position = "top") + 
    theme_minimal(base_size = 8) + 
    theme(axis.text.x = element_text(angle = 60, hjust = 0), axis.title.x=element_blank(),
          axis.text.y = element_text(face = "italic"),
          strip.background = element_blank(),strip.text.x = element_blank(),
          legend.position = "none") + 
    #ggtitle(label = Patient) +
    facet_grid(~grid, scales = "free", space = "free") 
  
  print(data_table)
}
pdf(file = paste(input_folder,"Manuscript/Fig_S4A.pdf", sep = ""), width = 4, height = 3, pointsize = 10)

plot_phenomatch_small(Genes_SVs = Genes_SVs, Patient = "35293", phenomatch_folder = paste(input_folder, "Results/Integration/MP/Phenomatch/", sep = ""))
dev.off()

# Fig S3B
conf <- read.delim(paste(input_folder,"Scripts/Visualization/Conf_files/TADs_only.conf", sep = ""), stringsAsFactors = F)

pdf(file = paste(input_folder,"Results/Figures_raw/Fig_S4B.pdf", sep = ""), width = 3.3, height = 2.3, pointsize = 8)
breakpoint_plot(conf = conf, 
                patient = "35293", 
                chromosome = 14,
                xlim_plot = c(28.5e6,31.5e6), 
                Genes_SVs = "")
dev.off()

## Fig S3C
RNA_expression_data <- read.delim(paste(input_folder,"Data/RNA/Schmitt/primary_cohort.gencode.hg19.FPKM.matrix", sep = ""))
RNA_expression_FOXG1 <- RNA_expression_data[RNA_expression_data$gene == "FOXG1",]
RNA_expression_FOXG1_m <- melt(RNA_expression_FOXG1[,c(2, 4:ncol(RNA_expression_FOXG1))])

RNA_expression_FOXG1_m <- RNA_expression_FOXG1_m[order(RNA_expression_FOXG1_m$value, decreasing = T),][1:5,]
RNA_expression_FOXG1_m$variable <- factor(RNA_expression_FOXG1_m$variable, levels = RNA_expression_FOXG1_m$variable[order(RNA_expression_FOXG1_m$value, decreasing = F)])

ggplot(RNA_expression_FOXG1_m, aes(x= variable, y = value, fill = variable)) +
  geom_bar(stat = "identity", 
           col = "black", width = 0.7, lwd = 0.3) + 
  theme_minimal(base_size = 8) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title = element_text(face = "plain", size = 8),
        legend.position="none",
        panel.grid.major.x = element_line(colour="darkgrey", size = 0.2)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip() +
  labs(y = "RNA expression (FPKM)", x = "Cell type")

ggsave(filename = paste(input_folder,"Results/Figures_raw/Fig_S4C_FOXG1_Expression.pdf", sep = ""), height= 45, width= 45, units = "mm")


## Fig S3D: Enhancers
conf <- read.delim(paste(input_folder,"Scripts/Visualization/Conf_files/TADs_Enhancers.conf", sep = ""), stringsAsFactors = F)
pdf(file = paste(input_folder,"Results/Figures_raw/Fig_S4D.pdf", sep = ""), width = 3.5, height = 2, pointsize = 8)

breakpoint_plot(conf = conf, 
                patient = "35293", 
                chromosome = 14,
                xlim_plot = c(28.5e6,31.5e6), 
                Genes_SVs = "")
dev.off()

## Fig S3E: Virtual 4C
conf <- read.delim(paste(input_folder,"Scripts/Visualization/Conf_files/FOXG1_V4C.conf", sep = ""), stringsAsFactors = F)
pdf(file = paste(input_folder,"Results/Figures_raw/Fig_S4F_V4C.pdf", sep = ""), width = 3.5, height = 2, pointsize = 8)

breakpoint_plot(conf = conf, 
                patient = "35293", 
                chromosome = 14,
                xlim_plot = c(28.5e6,31.5e6), 
                Genes_SVs = "")
dev.off()

## Fig S3F: PCHiC + DHS
conf <- read.delim(paste(input_folder,"Scripts/Visualization/Conf_files/PCHiC.conf", sep = ""), stringsAsFactors = F)
pdf(file = paste(input_folder,"Results/Figures_raw/Fig_S4E.pdf", sep = ""), width = 3.5, height = 2, pointsize = 8)

breakpoint_plot(conf = conf, 
                patient = "35293", 
                chromosome = 14,
                xlim_plot = c(28.5e6,31.5e6), 
                Genes_SVs = "")
dev.off()




## Fig S7: SOX9

RNA_expression_data <- read.delim(paste(input_folder,"/Common_data/RNA/Schmitt/primary_cohort.gencode.hg19.FPKM.matrix", sep = ""))
RNA_expression_SOX9 <- RNA_expression_data[RNA_expression_data$gene == "SOX9",]
RNA_expression_SOX9_m <- melt(RNA_expression_SOX9[,c(2, 4:ncol(RNA_expression_SOX9))])

RNA_expression_SOX9_m <- RNA_expression_SOX9_m[order(RNA_expression_SOX9_m$value, decreasing = T),]
RNA_expression_SOX9_m$variable <- factor(RNA_expression_SOX9_m$variable, levels = RNA_expression_SOX9_m$variable[order(RNA_expression_SOX9_m$value, decreasing = F)])

conf <- read.delim(paste(input_folder,"/Visualization/conf/SOX9.conf", sep = ""), stringsAsFactors = F)

pdf(file = paste(input_folder,"/Manuscript/Figures_raw/S6A_MP3760_SOX9_karyotype.pdf", sep = ""), width = 6, height = 3.5, pointsize = 8)
plot_der_karyotype(Patient = "MP3760", breakpoints = SVs, xlim = c(0,150e6))
dev.off()

pdf(file = paste(input_folder,"/Manuscript/Figures_raw/S6B_MP3760_SOX9_regulation.pdf", sep = ""), width = 6, height = 4, pointsize = 9)

breakpoint_plot(conf = conf, patient = "MP3760", chromosome = 17,
                xlim_plot = c(68e6,71e6), Genes_SVs = "")

axis(side = 2, at = c( as.numeric(conf[conf$Option == "Genes", "y"])-0.1, as.numeric(conf[conf$Option == "Genes", "y"])+1.4), labels = NA, tcl = 0.5, line = 4)
mtext(side = 2, text = "Genes", at = (as.numeric(conf[conf$Option == "Genes", "y"]) + as.numeric(conf[conf$Option == "Genes", "y"])+1.4)/2, line = 4.5, cex = 0.9)

## Translocations

mtext(text = "Translocations", side = 2, at = 7.2, las = 1, cex = 0.9)

rect(xleft = 68888853, xright = 69091290, ybottom = -1, ytop = 8.3, col = rgb(255,204,0, 50, maxColorValue = 255), border = NA)
text(x = (68888853+69091290) / 2, y = -0.5, labels = "PRS\nCluster", cex = 0.7)

rect(xleft = 69185496, xright = 69531685, ybottom = -1, ytop = 8.3, col = rgb(255,204,0, 50, maxColorValue = 255), border = NA)
text(x = (69185496+69531685) / 2, y = -0.5, labels = "Distal\nCluster", cex = 0.7)

rect(xleft = 69742161, xright = 70067161, ybottom = -1, ytop = 8.3, col = rgb(255,204,0, 50, maxColorValue = 255), border = NA)
text(x = (69742161+70067161) / 2, y = -0.5, labels = "Proximal\nCluster", cex = 0.7)

# T1 - Benko
#rect(xleft = 68919517, xright = 68956279, ybottom = 7.8, ytop = 8.1)
arrows(x0 = 68919517, x1 = 68956279, y0 = 7.9, angle = 90, code = 3, lwd = 1, length = 0.03, col = "darkred")
text(x = (68919517+68956279) / 2, y = 8.2, labels = "T1", cex = 0.7, col = "darkred")
# T2 - Benko

arrows(x0 = 69007280, x1 = 69091290, y0 = 7.9, angle = 90, code = 3, lwd = 1, length = 0.03, col = "darkred")
#arrows(x0 = (69007280+69091290)/2, y0 = 7.8, y1 = 8, code = 0, lwd = 1.5, length = 0.03)
text(x = (69007280+69091290) / 2, y = 8.2, labels = "T2", cex = 0.7, col = "darkred")

#rect(xleft = 69007280, xright = 69091290, ybottom = 7.8, ytop = 8.1, col = "red")
# T3 - Benko
arrows(x0 = (68888853+68888854)/2, y0 = 7.8, y1 = 8, code = 0, lwd = 1.5, length = 0.03, col = "darkred")
#arrows(x0 = 68888853, x1 = 68888854, y0 = 7.9, angle = 90, code = 3, lwd = 1, length = 0.03)
text(x = (68888853+68888854-50000) / 2, y = 8.2, labels = "T3", cex = 0.7, col = "darkred")
# T4 Jakobsen
arrows(x0 = 68973787, x1 = 68979922, y0 = 7.5, angle = 90, code = 3, lwd = 1, length = 0.03, col = "darkred")
text(x = (68973787+68979922) / 2, y = 7.2, labels = "T4", cex = 0.7, col = "darkred")
# Patient_1 (P1) Fonseca (no cleft palate)
arrows(x0 = 69201539, x1 = 69262086, y0 = 6.9, angle = 90, code = 3, lwd = 1, length = 0.03, col = "blue")
text(x = (69201539+69262086) / 2, y = 6.65, labels = "P1", cex = 0.7, col = "blue")
# Patient_2 Fonseca
arrows(x0 = 69516104, x1 = 69531685, y0 = 7.9, angle = 90, code = 3, lwd = 1, length = 0.03, col = "darkred")
text(x = (69516104+69531685) / 2, y = 8.2, labels = "P2", cex = 0.7, col = "darkred")

# Case 1 Leipoldt
arrows(x0 = (69741445+69741446)/2, y0 = 7.8, y1 = 8, code = 0, lwd = 1.5, length = 0.03, col = "blue")

text(x = (69741445+69741446) / 2, y = 8.2, labels = "C1", cex = 0.7, col = "blue")

# Case 2 Leipoldt 
arrows(x0 = (69327862+69327921)/2, y0 = 7.8, y1 = 8, code = 0, lwd = 1.5, length = 0.03, col = "darkred")
text(x = (69327862+69327862) / 2, y = 8.2, labels = "C2", cex = 0.7, col = "darkred")

# Patient 1 Velagaleti and patient MS Hill Harfe 
arrows(x0 = (69217932), y0 = 7.8, y1 = 8, code = 0, lwd = 1.5, length = 0.03, col = "darkred")
text(x = (69217932), y = 8.2, labels = "V1", cex = 0.7, col = "darkred")

# Patient Refai (no cleft palate)
arrows(x0 = 69306161, x1 = 69341161, y0 = 7.5, angle = 90, code = 3, lwd = 1, length = 0.03, col = "blue")
text(x = (69306161+69341161) / 2, y = 7.2, labels = "R1", cex = 0.7, col = "blue")

# Family F Hill-Harfe (= H1)
arrows(x0 = (69185496), y0 = 7.4, y1 = 7.6, code = 0, lwd = 1.5, length = 0.03, col = "darkred")
text(x = (69185496-30000), y = 7.2, labels = "H1", cex = 0.7, col = "darkred")

# Patient MS (= H2) Hill-Harfe (2005) + Pfeifer (2009)
arrows(x0 = (69216927+69217135) / 2, y0 = 7.4, y1 = 7.6, code = 0, lwd = 1.5, length = 0.03, col = "darkred")
text(x = (69216927+69217135) / 2 + 10000, y = 7.2, labels = "H2", cex = 0.7, col = "darkred")

# 
# ## Deletions (not included in figure)
# mtext(text = "Deletions", side = 2, at = 7.2, las = 1)
# 
# # F1
# rect(xleft = 68663405, xright = 68738405, ybottom = 7.15, ytop = 7.3, col = rgb(128, 0, 0, 200, maxColorValue = 255))
# text(x = (68663405+68738405) / 2, y = 7.4 , labels = "F1", cex = 0.7)
# # F2
# rect(xleft = 68676303, xright = 68676304, ybottom = 6.9, ytop = 7.05, col = rgb(128, 0, 0, 200, maxColorValue = 255))
# text(x = (68676303+68676304) / 2, y = 6.65, labels = "F2", cex = 0.7)
# # SP4
# rect(xleft = 68219155, xright = 68538005, ybottom = 6.9, ytop = 7.05, col = rgb(128, 0, 0, 200, maxColorValue = 255))
# text(x = (68219155+68538005) / 2, y = 6.65, labels = "SP4", cex = 0.7)
# 
# # PRS183
# rect(xleft = 68680882, xright = 68965323, ybottom = 6.9, ytop = 7.05, col = rgb(128, 0, 0, 200, maxColorValue = 255))
# text(x = (68680882+68965323) / 2, y = 6.65, labels = "PRS183", cex = 0.7, col = "darkred")
# 
# rect(xleft = 69704137, xright = 69710167, ybottom = 6.9, ytop = 7.05, col = rgb(128, 0, 0, 200, maxColorValue = 255))
# text(x = (69704137+69710167) / 2, y = 6.65, labels = "PRS116", cex = 0.7, col = "darkred")
# 
# # # SP2
# # rect(xleft = 71641405, xright = 71677405, ybottom = 7, ytop = 7.2)

dev.off()

## Fig S8: FOXG1




