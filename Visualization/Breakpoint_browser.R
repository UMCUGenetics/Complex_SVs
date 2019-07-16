
## REPLACE XXX WITH PATH TO INPUT DIRECTORY
input_folder <- "XXX"

breakpoint_plot <- function(conf = standard.conf,
                            patient,
                            chromosome,
                            xlim_plot,
                            TADs = TADs_raw,
                            breakpoints = SVs,
                            Genes_SVs = ""){
  
  ## Read conf
  # Conf can be a string (path to conf file) or a data.frame. If it is a string, read the conf file:
  if(isSingleString(conf)){
    print(paste("Reading:", conf))
    conf <- read.delim(conf, header = T, stringsAsFactors = F, comment.char = "#" )
  }
  
  ## Determine if normal of derivative chromosomes are plotted:
  derivative_mode <- ifelse(startsWith(x = as.character(chromosome), prefix = "der"), TRUE, FALSE)
  
  patient_data <- breakpoints[which(breakpoints$Patient == patient),]
  
  par(mar=c(0, 6, 0, 0))
  #par(mar=c(1, 5, 2, 1))
  plot.new()
  
  xlim <- xlim_plot
  ylim <- c(as.numeric(conf[conf$Option == "ybottom", "Value"]), 
            as.numeric(conf[conf$Option == "ytop", "Value"]))
  
  plot.window(xlim=xlim, ylim=ylim)
  
  title(main = patient, line = -1)
  
  if(conf[conf$Option == "Ideogram", "Value"] == TRUE){
    print("## Plotting ideogram")
    
    if(conf[conf$Option == "Ideogram","Parameters"] != ""){
      parameters_enhancers <- unlist(strsplit(conf[conf$Option == "Ideogram","Parameters"], split = ";"))
      width_ideogram <- as.numeric(unlist(strsplit(unlist(strsplit(parameters_enhancers[grep("width=", x = parameters_enhancers)], split = "="))[2], split = ",")))
      
    } else {
      width_ideogram <- 0.5
    }
    plot_ideogram(chromosome = chromosome, highlight = c(xlim[1], xlim[2]), plot = "inside", xlim = xlim, 
                  y0 = as.numeric(unlist(strsplit(conf[conf$Option == "Ideogram", "y"], split = ",")))[1],
                  y1 = as.numeric(unlist(strsplit(conf[conf$Option == "Ideogram", "y"], split = ",")))[2],
                  y_chromosome = as.numeric(unlist(strsplit(conf[conf$Option == "Chromosome", "y"], split = ",")))[1],
                  ideogram_width  = width_ideogram,
                  Patient = patient, 
                  SVs = patient_data, 
                  der = derivative_mode)
  }
  
  if(conf[conf$Option == "Genes", "Value"] == TRUE){
    print("## Plotting Genes")
    
    plot_genes(chr = chromosome, genes = genes_g, xlim = xlim, y0 = as.numeric(conf[conf$Option == "Genes", "y"]), 
             derivative = derivative_mode, patient = patient, label = TRUE, y_text = 1.2, breakpoints = patient_data)
  }
  
  if("TADs" %in% conf$Option){
    print("# Plotting TADs")
    
    for(TAD_cell_type in conf[conf$Option == "TADs","Value"]){
      print(TAD_cell_type)
      print(head(TADs))
      TADs_plot <- TADs[which(TADs$TAD_cell == TAD_cell_type),]
      # if(length(grep(pattern = TAD_cell_type, x = names(TADs))) > 0){
      #   TADs_plot <- TADs[[names(TADs)[grep(pattern = TAD_cell_type, x = names(TADs))]]]
      # } else {
      #   if(file.exists(TAD_cell_type)){
      #     TADs_plot <- read.delim(TAD_cell_type)
      #     
      #   } else {
      #     print(paste(TAD_cell_type , "not found!"))
      #   }
      # }
      
      plot_TADs(TADs = TADs_plot, 
                chr = chromosome, 
                Patient = patient, 
                SVs = breakpoints, 
                y0 =  as.numeric(conf[conf$Option == "TADs" & conf$Value == TAD_cell_type,"y"]),
                color =  conf[conf$Option == "TADs" & conf$Value == TAD_cell_type,"Color"],
                label = conf[conf$Option == "TADs" & conf$Value == TAD_cell_type,"Label"], 
                xlim = xlim,
                der = derivative_mode)
    }
    min_y_TADs <- as.numeric(min(conf[conf$Option == "TADs","y"]))
    max_y_TADs <- as.numeric(max(conf[conf$Option == "TADs", "y"]))
    axis(side = 2, line = 4, at = c(min_y_TADs-0.1,max_y_TADs+0.1), labels = NA, tcl = 0.5)
    mtext(side = 2, text = "TADs", at = (min_y_TADs+max_y_TADs)/2, line = 4.5, cex = 0.9)
  }
  
  if("Enhancers" %in% conf$Option){
    print("# Plotting Enhancers")
    
    enhancer_tracks <- conf[conf$Option == "Enhancers",]
    
    for(Enhancer_cell_type in enhancer_tracks[,"Value"]){
      print(Enhancer_cell_type)
      
      if(length(grep(Enhancer_cell_type, x = names(enhancer_list))) > 0){
        Enhancers_plot <- enhancer_list[[grep(Enhancer_cell_type, x = names(enhancer_list))]]
        
        parameters_enhancers <- unlist(strsplit(enhancer_tracks[enhancer_tracks$Value == Enhancer_cell_type,"Parameters"], split = ";"))
        
        bg_color <- unlist(strsplit(unlist(strsplit(parameters_enhancers[grep("background=", x = parameters_enhancers)], split = "="))[2], split = ","))
        color <- as.numeric(unlist(strsplit(conf[conf$Option == "Enhancers" & conf$Value ==Enhancer_cell_type ,"Color"], split = ",")))
        
        plot_enhancers(y0 = as.numeric(conf[conf$Option == "Enhancers" & conf$Value ==Enhancer_cell_type ,"y"]),
                       height = as.numeric(unlist(strsplit(parameters_enhancers[grep("height=", x = parameters_enhancers)], split = "="))[2]),
                       enhancers = Enhancers_plot, 
                       chr = chromosome,
                       fill_color = c(color[1], color[2], color[3]),
                       border_color = c(color[1], color[2], color[3]),
                       border_width = 0.8,
                       label = conf[conf$Option == "Enhancers" & conf$Value ==Enhancer_cell_type ,"Label"], 
                       patient = patient, derivative = derivative_mode, threshold = 0,
                       bg = bg_color, xlim = xlim)
        
        if(exists("y_enhancers") == FALSE){
          y_enhancers <- as.numeric(conf[conf$Option == "Enhancers" & conf$Value == Enhancer_cell_type ,"y"])
        } else {
          y_enhancers <- c(y_enhancers, as.numeric(conf[conf$Option == "Enhancers" & conf$Value ==Enhancer_cell_type ,"y"]))
        }
      } else {
        print(paste("# Data for enhancers: ", Enhancer_cell_type, " NOT found. Skipping. ", sep = ""))
      }
    }
    
    print(as.numeric(y_enhancers))
    axis(side = 2, line = 4, at = c(min(as.numeric(y_enhancers)), max(as.numeric(y_enhancers))+0.25), labels = NA, tcl = 0.5)
    mtext(side = 2, line = 4.5, text = "Enhancers", at = (min(as.numeric(y_enhancers)) + max(as.numeric(y_enhancers))+0.25)/2, cex = 0.9)
  }
  
  if("PCHiC" %in% conf$Option){
    print("# Plotting PCHiC")
    
    PCHiC_cell_types <- conf[conf$Option == "PCHiC", "Value"]
    for(PCHiC_cell_type in PCHiC_cell_types){
      if(derivative_mode != TRUE){
        
        if(is.null(PCHiC_data_list[[PCHiC_cell_type]]) == TRUE){
          print(paste("Reading: ",input_folder,"Common_data/PCHiC/PCHiC_",PCHiC_cell_type,".txt", sep = ""))
          PCHiC_data_list[[PCHiC_cell_type]] <- read.delim(paste(input_folder,"Common_data/PCHiC/PCHiC_",PCHiC_cell_type,".txt", sep = ""), stringsAsFactors = F)
        }
        
        PCHiC_interactions <- PCHiC_data_list[[PCHiC_cell_type]]
        
        PCHiC_color <- as.numeric(unlist(strsplit(x = conf[conf$Option == "PCHiC" & conf$Value == PCHiC_cell_type, "Color"], split = ",")))
        PCHiC_parameters <- unlist(strsplit(conf[conf$Option == "PCHiC","Parameters"], split = ";"))
        PCHiC_bait = unlist(strsplit(PCHiC_parameters[grep("bait=", x = PCHiC_parameters)], split = "="))[2]

        plot_interactions(interactions = PCHiC_interactions, 
                          y0 = as.numeric(conf[conf$Option == "PCHiC" & conf$Value == PCHiC_cell_type, "y"]), 
                          lwd = 0.6,
                          bow = -1.5, xlim = xlim, rgb = c(PCHiC_color[1], PCHiC_color[2], PCHiC_color[3]), 
                          chromosome = chromosome, label = conf[conf$Option == "PCHiC" & conf$Value == PCHiC_cell_type, "Label"], bait = PCHiC_bait)
      }
      
    }
    axis(side = 2, line = 4, at = c(as.numeric(conf[conf$Option == "PCHiC" & conf$Value == PCHiC_cell_type, "y"]),as.numeric(conf[conf$Option == "PCHiC" & conf$Value == PCHiC_cell_type, "y"])-1.5), labels = NA, tcl = 0.5)
    mtext(side = 2, text = "PCHiC", at = (as.numeric(conf[conf$Option == "PCHiC" & conf$Value == PCHiC_cell_type, "y"])+0.5)/2, line = 4.5, cex = 0.9)
  }
  
  if("V4C" %in% conf$Option){
    
    V4C_tracks <- conf[conf$Option == "V4C",]
    print("Plotting V4C")
    for(V4C_track in 1:nrow(V4C_tracks)){
      
      parameters_V4C <- unlist(strsplit(V4C_tracks[V4C_track,"Parameters"], split = ";"))
      
      input_folder = unlist(strsplit(parameters_V4C[grep("folder=", x = parameters_V4C)], split = "="))[2]
      V4C_gene = unlist(strsplit(parameters_V4C[grep("gene=", x = parameters_V4C)], split = "="))[2]
      V4C_resolution = unlist(strsplit(parameters_V4C[grep("resolution=", x = parameters_V4C)], split = "="))[2]
      V4C_color <- V4C_tracks[V4C_track,"Color"]
      
      plot_v4C(gene = V4C_gene, cell_type = V4C_tracks[V4C_track,"Value"], y0 = as.numeric(V4C_tracks[V4C_track ,"y"]),
               bin_size = as.numeric(V4C_resolution), color = V4C_color, V4C_folder = input_folder, label = V4C_tracks[V4C_track,"Label"])
      
      axis(side = 2, line = 4, at = c(as.numeric(V4C_tracks[V4C_track, "y"]),as.numeric(V4C_tracks[V4C_track, "y"])+1.4), labels = NA, tcl = 0.5)
      mtext(side = 2, text = "V4C", at = (as.numeric(V4C_tracks[V4C_track, "y"]) + as.numeric(V4C_tracks[V4C_track, "y"])+1.4)/2, line = 4.5, cex = 0.9)
    }
  }
  
  
  if(conf[conf$Option == "Chromosome", "Value"] == TRUE){
    print("## Plotting chromosome axis")
    parameters_chromosome <- unlist(strsplit(conf[conf$Option == "Chromosome","Parameters"], split = ";"))
    
    plot_chromosome(chr = chromosome, 
                    patient = patient, 
                    y0 = as.numeric(conf[conf$Option == "Chromosome", "y"]), 
                    xlim = xlim, 
                    ylim = c(ylim[1], as.numeric(unlist(strsplit(conf[conf$Option == "Chromosome", "y"], split = ",")))[1]), 
                    derivative = derivative_mode,
                    location =  unlist(strsplit(parameters_chromosome[grep("position=", x = parameters_chromosome)], split = "="))[2],
                    SVs = patient_data)
  }
  
}


plot_ideogram <- function(chromosome,
                          SVs,
                          Patient,
                          highlight = c(10e6, 12e6),
                          plot = "inside",
                          ideogram_width = 0.5, 
                          xlim = "",
                          y0 = 1,
                          y1 = 1.3,
                          y_chromosome = "",
                          der = FALSE,
                          chr_colors = c(255,0,0,30,
                                         50,205,50,30,
                                         0,0,255,30,
                                         255,140,0,30,
                                         128,0,128,30)){
  
  # Karyogram bands and colors are obtained from the biovizBase package:
  data(hg19IdeogramCyto, package = "biovizBase")
  seqlevels(hg19IdeogramCyto) <- gsub("chr", "", seqlevels(hg19IdeogramCyto))
  
  # Breakpoints will plotted only if a patient id and SV data is added:
  if(Patient != ""){
    SVs_patient <- SVs[which(SVs$Patient == Patient),]
  } else {
    SVs_patient <- data.frame()
  }
  
  # The plot width is used to fit the ideogram within the plotwindow 
  plot_width <- (xlim[2] - xlim[1]) * ideogram_width
  
  # Calculate the plot coordinates of the bands (separately for normal and derivative chroms)
  if(der == FALSE){
    SVs_chr <- SVs_patient[which(SVs_patient$chr == chromosome),]
    
    chr_data <- hg19IdeogramCyto[seqnames(hg19IdeogramCyto) == chromosome,]
    
    highlight_g <- GRanges(seqnames = chromosome, IRanges(start = highlight[1], end = highlight[2]))
    
    # The start and end bands of the highlight will be labeled on the ideogram. This label is obtained from the chr_data
    # Determine the band(s) overlapping with the X-coordinates of the plot (the highlight)
    olap_chr_highlight <- findOverlaps(chr_data, highlight_g)
    overlaps <- pintersect(chr_data[queryHits(olap_chr_highlight)], highlight_g[subjectHits(olap_chr_highlight)])
    
    start_label <- as.character(overlaps$name[which(start(overlaps) == min(start(overlaps)))])
    end_label <- as.character(overlaps$name[which(start(overlaps) == max(start(overlaps)))])
    
    bands <- as.data.frame(chr_data)
    bands$band_start <- bands$start
    bands$band_end <- bands$end
    
  } else {
    # Fragments rom different chromosomes will get different colors:
    i <- 0
    for(chr in unique(SVs_patient$chr)){
      SVs_patient$chr_fill[SVs_patient$chr == chr] <- rgb(chr_colors[1+i*4], chr_colors[2+i*4], chr_colors[3+i*4], chr_colors[4+i*4], maxColorValue = 255)
      SVs_patient$chr_border[SVs_patient$chr == chr] <- rgb(chr_colors[1+i*4], chr_colors[2+i*4], chr_colors[3+i*4], 200, maxColorValue = 255)
      i <- i + 1
    }
    
    SVs_chr <- SVs_patient[which(SVs_patient$der_chr == chromosome),]
    
    bands <- data.frame()
    
    for(der_fragment in SVs_chr$der_fragment){
      print(der_fragment)
      fragment_data <- SVs_chr[which(SVs_chr$der_fragment == der_fragment),]
      names(fragment_data) <- paste("SV", names(fragment_data), sep = "_")
      
      # Select the bands on the derivative chromosome:
      der_fragment_g <- GRanges(seqnames = fragment_data$SV_chr, IRanges(start =fragment_data$SV_start, end = fragment_data$SV_end))
      olap_fragment_bands <- findOverlaps(hg19IdeogramCyto, der_fragment_g)
      
      bands_fragment <- cbind(as.data.frame(hg19IdeogramCyto[queryHits(olap_fragment_bands)]), as.data.frame(fragment_data[subjectHits(olap_fragment_bands),]))
      bands_fragment$band_start <- ifelse(bands_fragment$SV_type != "Inversion" & bands_fragment$SV_type != "Flanking_Inversion",
                                          bands_fragment$SV_der_start + (bands_fragment$start - bands_fragment$SV_start),
                                          bands_fragment$SV_der_start + (bands_fragment$SV_end - bands_fragment$end))
      bands_fragment$band_end <- ifelse(bands_fragment$SV_type != "Inversion" & bands_fragment$SV_type != "Flanking_Inversion",
                                        bands_fragment$SV_der_start + (bands_fragment$end - bands_fragment$SV_start),
                                        bands_fragment$SV_der_start + (bands_fragment$SV_end - bands_fragment$start))
      bands_fragment$band_start <- ifelse(bands_fragment$band_start < bands_fragment$SV_der_start, bands_fragment$SV_der_start, bands_fragment$band_start)
      bands_fragment$band_end <- ifelse(bands_fragment$band_end > bands_fragment$SV_der_end, bands_fragment$SV_der_end, bands_fragment$band_end)
      
      if(unique(bands_fragment$SV_type) == "Duplication"){
        
        Duplicated_bands <- bands_fragment
        Duplicated_bands$band_start <- Duplicated_bands$band_start + Duplicated_bands$SV_width
        Duplicated_bands$band_end <- Duplicated_bands$band_end + Duplicated_bands$SV_width
        
        bands_fragment <- rbind(bands_fragment, Duplicated_bands)
        bands_fragment$SV_der_width <- bands_fragment$SV_der_width*2
      }
      bands <- rbind(bands, bands_fragment)
    }
  }
  
  # plot_start and plot_end are used to determine where the band has to be plot (corrected for band coordinates and plot window size)
  bands$plot_start <- bands$band_start / max(bands$band_end) * plot_width + xlim[1]
  bands$plot_end <- bands$band_end / max(bands$band_end) *plot_width + xlim[1]
  
  # Add colors to the bands:
  bands$col <- bands$gieStain
  bands$col  <- factor(bands$col)
  levels(bands$col) <- list( "grey100" = "gneg",
                             "grey100" = "gvar",
                             "brown3" = "stalk",
                             "#C1C1C1" = "gpos25",
                             "#808080" = "gpos50",
                             "#404040" = "gpos75",
                             "#000000" = "gpos100",
                             "brown4" = "acen")
    
    # Plot the bands
    rect(xleft = bands$plot_start, xright = bands$plot_end, ybottom = y0, ytop = y1, col = as.character(bands$col), border = NA)
    
    # Plot the chromosome color
    # Current "normal" chromosomes are just black and the derivative chromosomes have different colors
    if(der == FALSE){
      rect(xleft = min(bands$plot_start), xright = max(bands$plot_end), ybottom = y0, ytop = y1, col = NA, border = "black", lwd = 0.3)
    } else {
      for(der_fragment in unique(bands$SV_der_fragment)){
        bands_fragment <- bands[which(bands$SV_der_fragment == der_fragment),]
        
        # Plot the chromosome colors over the bands
        rect(xleft = min(bands_fragment$plot_start), xright = max(bands_fragment$plot_end), ybottom = y0, ytop = y1, 
             col = unique(as.character(bands_fragment$SV_chr_fill)), border = unique(as.character(bands_fragment$SV_chr_border)), lwd = 1.5)
        
        if(max(bands_fragment$plot_end) != max(bands$plot_end)){
          arrows(x0 = max(bands_fragment$plot_end), y0 = y0-0.3, y1 = y1+0.3,lty = 2, code = 0, col = rgb(0,0,0, 150, maxColorValue = 255), lwd = 1.5)
          
        }
      }
    }
    
    # Plot the highlight on the ideogram
    rect(xleft = highlight[1] / max(bands$band_end) * plot_width + xlim[1], xright = highlight[2] /max(bands$band_end) * plot_width + xlim[1], ybottom = y0-0.2,
         ytop = y1+0.2, col = rgb(200,0,0,100, maxColorValue = 255), border = rgb(200,0,0,200, maxColorValue = 255))
    
    # Plot the breakpoint junctions
    if(nrow(SVs_patient) > 0){
      if(der == FALSE){
        bps_chr <- SVs_chr$end / max(bands$band_end) * plot_width + xlim[1]
      } else {
        bps_chr <-  SVs_chr$der_end / max(bands$band_end) * plot_width + xlim[1]
      }
      bps_chr <- as.numeric(bps_chr[which(bps_chr !=  max(bands$plot_end))])
      arrows(x0 = bps_chr, y0 = y0-0.3, y1 = y1+0.3,lty = 2, code = 0, col = "red", lwd = 0.8)
    }
    
    # This part plot connecting lines (zoom in) between the ideogram and the chromosome axis
    if(y_chromosome != ""){
      arrows(x0 = highlight[1] / max(bands$band_end) * plot_width + xlim[1], x1 = xlim[1],
             y0 = y0 - 0.2, y1 = y_chromosome, lty = 2, code = 0, col = rgb(0,0,0, 80, maxColorValue = 255))
      arrows(x0 = highlight[2] / max(bands$band_end) * plot_width + xlim[1], x1 = xlim[2],
             y0 = y0 - 0.2, y1 = y_chromosome, lty = 2, code = 0, col = rgb(0,0,0, 80, maxColorValue = 255))
    }
    
    # Plot the chromosome number at the y-axis label
    mtext(side = 2, text = chromosome, at = (y0+y1)/2, line = 2, las = 1, font = 2)
}


plot_chromosome <- function(patient = "16784",
                            derivative = FALSE,
                            chr = "6", 
                            y0 = 0, 
                            show_breaks = TRUE,
                            xlim = xlim_plot,
                            ylim = ylim_plot,
                            SVs = breakpoints,
                            location = "top",
                            linewidth = 1.5){
  
  ## This function plots a chromosome axis showing fragments of SVs and a coordinate axis
  patient_data <- SVs[which(SVs$Patient == patient),]

  patient_data$type <- factor(patient_data$type)
  levels(patient_data$type) <- list( "DEL" = "Deletion", "DUP" = "Duplication", "DUP" ="Inverted Duplication", "INV" = "Inversion",
                                        "TRUNC" = "Truncation", "INTRA" = "Normal" , "FLANK" = "Flanking_Inversion", "FLANK" = "Flanking_Normal",
                                        "INS" = "Insertion", "INS" = "Inverted_insertion")
  
  # Set the colors for the different chromosomes:
  colors <- colorRampPalette(brewer.pal(8,"PRGn"))(length(levels(factor(patient_data$chr))))
  chromosome_colors <- data.frame(chr = levels(factor(patient_data$chr)), color = colors)
  
  # Determine the column names for the necessary data:
  fragment_start <- ifelse(derivative == FALSE, "start", "der_start")
  fragment_end <- ifelse(derivative == FALSE, "end", "der_end")
  fragment_column <- ifelse(derivative == FALSE, "Fragment", "der_fragment")
  chr_column <- ifelse(derivative == FALSE, "chr", "der_chr")
  
  # Plot the ticks 
  if(xlim[2]-xlim[1] > 2e6){
    major_ticks <- seq(round(xlim[1]/1e6, 0) * 1e6, round(xlim[2]/1e6, 0) * 1e6, by = 1e6)
    minor_ticks <- seq(round(xlim[1]/1e6, 0) * 1e6, round(xlim[2]/1e6, 0) * 1e6, by = 0.5e6)
  } else if(xlim[2]-xlim[1] > 1e6){
    major_ticks <- seq(round(xlim[1]/1e6, 0) * 1e6, round(xlim[2]/1e6, 0) * 1e6, by = 0.4e6)
    minor_ticks <- seq(round(xlim[1]/1e6, 0) * 1e6, round(xlim[2]/1e6, 0) * 1e6, by = 0.2e6)
  } else if(xlim[2]-xlim[1] > 0.5e6){
    major_ticks <- seq(round(xlim[1]/1e6, 0) * 1e6, round(xlim[2]/1e6, 0) * 1e6, by = 0.2e6)
    minor_ticks <- seq(round(xlim[1]/1e6, 0) * 1e6, round(xlim[2]/1e6, 0) * 1e6, by = 0.1e6)
  }
  
  minor_ticks <- minor_ticks[-which(minor_ticks %in% major_ticks)]
  
  arrows(x0 = major_ticks, code = 0, y0 = y0, y1 =  ifelse(location == "top", y0 + 0.2, y0 - 0.2), lwd = 1.5)
  arrows(x0 = minor_ticks, code = 0, y0 = y0, y1 =  ifelse(location == "top", y0 + 0.15, y0 - 0.15), lwd = 1)
  
  # Tick labels
  text(x = major_ticks, y = ifelse(location == "top", y0 + 0.4, y0 - 0.4), labels = paste(round(major_ticks / 1e6,1), "Mb"), font = 3, cex = 0.8)
  
  
  # Plot a gray line from the left to the right of the plot if there's no patient data (instead of arrows for SV fragments):
  arrows(x0 = xlim[1],
         x1 = xlim[2],
         y0 = y0,
         code = 0,
         lwd = linewidth,
         col = "darkgray")
  
  if(nrow(patient_data[which(patient_data[,chr_column] == chr),]) > 0){
    
    # Only plot arrows for each fragment if breakpoints have been found for the specific patient
    for(fragment in patient_data[,fragment_column][which(patient_data[,chr_column] == chr)]){
      print(fragment)
      
      # Determine the directionaly of the arrow (left to right, or right to left in case of inversion). If the fragment is less than 200kb no arrow will be drawn (just a line)
      if(patient_data[patient_data[,fragment_column] == fragment, "type"] == "INV" & 
         patient_data[patient_data[,fragment_column] == fragment, "width"] > 0.15e6 & 
         (xlim[2]-xlim[1]) > 2e6){
        arrow_code <- 1
      } else if (patient_data[patient_data[,fragment_column] == fragment, "type"] == "INV" & 
           patient_data[patient_data[,fragment_column] == fragment, "width"] > 0.10e6 & 
           (xlim[2]-xlim[1]) > 1e6){
          arrow_code <- 1
      } else if (patient_data[patient_data[,fragment_column] == fragment, "width"] < 0.15e6 & (xlim[2]-xlim[1]) > 2e6){
        arrow_code <- 0
      } else if (patient_data[patient_data[,fragment_column] == fragment, "width"] < 0.10e6 & (xlim[2]-xlim[1]) > 1e6){
        arrow_code <- 0
      } else {
        arrow_code <- 2
      }
      
      # Plot an arrow for each genomic fragment within the plot window
      arrows(x0 = ifelse(patient_data[patient_data[,fragment_column] == fragment, fragment_start] > xlim[1],
                         patient_data[patient_data[,fragment_column] == fragment, fragment_start],
                         xlim[1]-1e6),
             
             x1 = ifelse(patient_data[patient_data[,fragment_column] == fragment, fragment_end] < xlim[2],
                         patient_data[patient_data[,fragment_column] == fragment, fragment_end],
                         xlim[2]+1e6),
             y0 = y0,
             code = arrow_code,
             lwd = linewidth,
             length = 0.08, angle = 50,
             col = as.vector(chromosome_colors$color[chromosome_colors$chr == patient_data$chr[patient_data[,fragment_column] == fragment]]))
      
      # The start (5') of each fragment starts with a vertical line
      arrows(x0 = ifelse(patient_data[patient_data[,fragment_column] == fragment, "type"] == "INV",
                         patient_data[patient_data[,fragment_column] == fragment, fragment_end],
                         patient_data[patient_data[,fragment_column] == fragment, fragment_start]),
             y0 = y0 - 0.2,
             y1 = y0 + 0.2,
             code = 0,
             lwd = linewidth,
             col = as.vector(chromosome_colors$color[chromosome_colors$chr == patient_data$chr[patient_data[,fragment_column] == fragment]]))
      
      # Determine the coordinates of the text label of the fragment:
      # If the start of the fragment is less than the start coordinate of the plot, the fragment label will be placed between the start of the plot and the end of the fragment
      text_start <- ifelse(patient_data[patient_data[,fragment_column] == fragment, fragment_start] < xlim[1],
                           ((xlim[1] + patient_data[patient_data[,fragment_column] == fragment, fragment_end]) / 2),
                           (patient_data[patient_data[,fragment_column] == fragment, fragment_end] + patient_data[patient_data[,fragment_column] == fragment, fragment_start]) / 2)
      
      # If the end of the fragment is more than the end coordinate of the plot, the fragment label will be placed between the start of the fragment and the end of the plot
      text_start <- ifelse(patient_data[patient_data[,fragment_column] == fragment, fragment_end] > xlim[2],
                           ((xlim[2] + patient_data[patient_data[,fragment_column] == fragment, fragment_start]) / 2),
                           text_start)
      
      # If both the start and the end of the fragment are outside the plot limits, the label will be placed in the center of the plot coordinates
      text_start <- ifelse(patient_data[patient_data[,fragment_column] == fragment, fragment_end] > xlim[2] & patient_data[patient_data[,fragment_column] == fragment, fragment_start] < xlim[1],
                           ((xlim[2] +  xlim[1]) / 2),
                           text_start)
      
      # Only plot labels if the fragment is more than 200 kb in length:
      if(patient_data$width[which(patient_data[,fragment_column] == fragment)] > 0.2e6){
        text(x = text_start,
             y = ifelse(location == "top", y0 + 0.6, y0 - 0.6),
             label = patient_data[,fragment_column][which(patient_data[,fragment_column] == fragment)],
             col = "red", cex = 0.8)
        
        text(x = text_start,
             y =  ifelse(location == "top", y0 + 0.8, y0 - 0.8),
             label = patient_data$type[which(patient_data[,fragment_column] == fragment)],
             col = "red", cex = 0.8)
      }
      
      # plot the breakpoint junctions 
      if(show_breaks == TRUE){
        if(nrow(patient_data) > 0){
          arrows(x0 = c(patient_data[patient_data[,fragment_column] == fragment, fragment_start],patient_data[patient_data[,fragment_column] == fragment, fragment_end]),
                 y0 = ylim[1]-0.1,
                 y1 = ylim[2]+0.75,
                 code = 0,
                 col = rgb(200,0,0, 150, maxColorValue = 255), 
                 lty =  2, 
                 lwd = 1)
        }
        
      }
    }
    

    
  } 


 

  # Plot the chromosome number at the y-axis label
  mtext(side = 2, text = ifelse(startsWith(prefix = "der", as.character(chr)), chr, paste("Chr ", chr, sep = "")), at = y0, las = 1, line =2, font = 2)
  
}
    



plot_genes <- function(y0 = 0.5, 
                       genes = genes_g,
                       exons = exons_g,
                       chr, 
                       y_text = 0.9, 
                       derivative = FALSE,
                       patient = "16784",
                       xlim = xlim_plot,
                       gene_info_file = "",
                       label = TRUE, 
                       breakpoints = SVs,
                       exon_height = 0.2){
  
  
  ## streepje begin - eind gen
  ## pijltjes begin - eind gen
  ## exon
  ## gen name
  
  SVs_patient <- breakpoints[which(breakpoints$Patient == patient),]
  
  # fragment_start <- ifelse(derivative == FALSE, "start", "der_start")
  # fragment_end <- ifelse(derivative == FALSE, "end", "der_end")
  # fragment_column <- ifelse(derivative == FALSE, "Fragment", "der_fragment")
  # chr_column <- ifelse(derivative == FALSE, "chr", "der_chr")
  
  if(derivative == FALSE){
    # select the genes on the chromosome within 1mb of the plot window:
    genes_to_plot <- genes[seqnames(genes) == chr,]
    genes_to_plot <- genes_to_plot[start(genes_to_plot) > xlim[1]-0.1e6,]
    genes_to_plot <- genes_to_plot[end(genes_to_plot) < xlim[2]+0.1e6,]
    
    genes_to_plot <- data.frame(genes_to_plot)
    genes_to_plot$plot_start <- ifelse(genes_to_plot$strand == "+", genes_to_plot$start, genes_to_plot$end)
    genes_to_plot$plot_end <- ifelse(genes_to_plot$strand == "+", genes_to_plot$end, genes_to_plot$start)
    
    # genes on + strand will be plotted at a higher y than genes on - strand;
    genes_to_plot$y_gene <- ifelse(genes_to_plot$strand == "+", y0 + 0.5, y0)
    
    exons_to_plot <- all_exons[which(all_exons$ensembl_gene_id %in% genes_to_plot$ensembl_gene_id),]
    exons_to_plot$plot_start <- exons_to_plot$exon_chrom_start
    exons_to_plot$plot_end <- exons_to_plot$exon_chrom_end
    exons_to_plot <- merge(exons_to_plot, genes_to_plot[,c("ensembl_gene_id", "y_gene")])
    
  } else {
    # determine the genes on the derivative chromosome:
    genes_to_plot <- data.frame(stringsAsFactors = F)
    exons_to_plot <- data.frame(stringsAsFactors = F)
    for(fragment in SVs_patient$der_fragment[which(SVs_patient$der_chr == chr)]){
      der_fragment <- SVs_patient[SVs_patient$der_fragment == fragment,]
      fragment_g <- GRanges(seqnames =SVs_patient[SVs_patient$der_fragment == fragment,"chr"], 
                            IRanges(start = SVs_patient[SVs_patient$der_fragment == fragment,"start"],
                                    end = SVs_patient[SVs_patient$der_fragment == fragment,"end"]))
      
      exons_g <- GRanges(seqnames = all_exons$chromosome_name, IRanges(start = all_exons$exon_chrom_start, end = all_exons$exon_chrom_end), ensembl_gene_id = all_exons$ensembl_gene_id)
      
      #print(fragment_g)
      
      olap_genes <- findOverlaps(fragment_g, genes)
      olap_exons <- findOverlaps(fragment_g, exons_g)
      
      overlapping_genes <- genes[subjectHits(olap_genes),]
      overlapping_exons <- exons_g[subjectHits(olap_exons),]
      
      #print(head(overlapping_genes))
      
      if(grepl("nver", x = SVs_patient[SVs_patient$der_fragment == fragment, "type"]) != TRUE){
        gene_starts <- start(overlapping_genes) - SVs_patient[SVs_patient$der_fragment == fragment,"start"] + 
          as.numeric(SVs_patient[SVs_patient$der_fragment == fragment,"der_start"])
        gene_starts <- ifelse(gene_starts < as.numeric(SVs_patient[SVs_patient$der_fragment == fragment,"der_start"]),
                         as.numeric(SVs_patient[SVs_patient$der_fragment == fragment,"der_start"]), gene_starts)
        
        gene_ends <- end(overlapping_genes) - SVs_patient[SVs_patient$der_fragment == fragment,"start"] + 
          as.numeric(SVs_patient[SVs_patient$der_fragment == fragment,"der_start"])
        gene_ends <- ifelse(gene_ends > as.numeric(SVs_patient[SVs_patient$der_fragment == fragment,"der_end"]),
                       as.numeric(SVs_patient[SVs_patient$der_fragment == fragment,"der_end"]), gene_ends)
        
        exon_starts <- start(overlapping_exons) - SVs_patient[SVs_patient$der_fragment == fragment,"start"] + 
          as.numeric(SVs_patient[SVs_patient$der_fragment == fragment,"der_start"])
        exon_starts <- ifelse(exon_starts < as.numeric(SVs_patient[SVs_patient$der_fragment == fragment,"der_start"]),
                              as.numeric(SVs_patient[SVs_patient$der_fragment == fragment,"der_start"]), exon_starts)
        
        exon_ends <- end(overlapping_exons) - SVs_patient[SVs_patient$der_fragment == fragment,"start"] + 
          as.numeric(SVs_patient[SVs_patient$der_fragment == fragment,"der_start"])
        exon_ends <- ifelse(exon_ends > as.numeric(SVs_patient[SVs_patient$der_fragment == fragment,"der_end"]),
                            as.numeric(SVs_patient[SVs_patient$der_fragment == fragment,"der_end"]), exon_ends)
        
      } else {
        gene_starts <- SVs_patient[SVs_patient$der_fragment == fragment,"end"] - start(overlapping_genes) + 
          SVs_patient[SVs_patient$der_fragment == fragment,"der_start"]
        gene_ends <-  SVs_patient[SVs_patient$der_fragment == fragment,"end"] - end(overlapping_genes) + 
          SVs_patient[SVs_patient$der_fragment == fragment,"der_start"]
        
        exon_starts <- SVs_patient[SVs_patient$der_fragment == fragment,"end"] - start(overlapping_exons) + 
          SVs_patient[SVs_patient$der_fragment == fragment,"der_start"]
        
        exon_ends <-  SVs_patient[SVs_patient$der_fragment == fragment,"end"] - end(overlapping_exons) + 
          SVs_patient[SVs_patient$der_fragment == fragment,"der_start"]
      }
      
        genes_to_plot_fragment <- data.frame(hgnc_symbol = overlapping_genes$hgnc_symbol,
                                    ensembl_gene_id = overlapping_genes$ensembl_gene_id,
                                    start = start(overlapping_genes),
                                    end = end(overlapping_genes),
                                    plot_start = gene_starts,
                                    plot_end = gene_ends,
                                    strand = strand(overlapping_genes),
                                    y_gene = ifelse(strand(overlapping_genes) == "+", y0 +0.5, y0))
        
        exons_to_plot_fragment <- data.frame(ensembl_gene_id = overlapping_exons$ensembl_gene_id,
                                             start = start(overlapping_exons),
                                             end = end(overlapping_exons),
                                             plot_start = exon_starts,
                                             plot_end = exon_ends)
        
        genes_to_plot <- rbind(genes_to_plot, genes_to_plot_fragment)
        exons_to_plot <- rbind(exons_to_plot, exons_to_plot_fragment)
        
    }
    
    genes_to_plot <- genes_to_plot[which(genes_to_plot$plot_start > xlim[1]-1e6 & genes_to_plot$plot_start < xlim[2]+1e6),]
    exons_to_plot <- exons_to_plot[which(exons_to_plot$plot_start > xlim[1]-1e6 & exons_to_plot$plot_start < xlim[2]+1e6),]
  }
  
  for(i in 1:nrow(genes_to_plot)){
    gene_to_plot <- genes_to_plot[i,]
    print(gene_to_plot)
    # plot the arrow heads indicating the directionality of the gene (depending if the gene is on + or - strand)
    if(gene_to_plot$plot_end - gene_to_plot$plot_start > 60000){
      arrows(x0 = seq(gene_to_plot$plot_start+30000, gene_to_plot$plot_end-30000, by = 30000), 
             x1 = seq(gene_to_plot$plot_start+30000, gene_to_plot$plot_end-30000, by = 30000)+10000,
             y0 = gene_to_plot$y_gene, code = ifelse(derivative == FALSE, 2, 1),
             col = "gray", length = 0.06, angle = 60)
    } 
    if(gene_to_plot$plot_start - gene_to_plot$plot_end > 60000){
      print("check")
      arrows(x0 = seq(gene_to_plot$plot_start-30000, gene_to_plot$plot_end+30000, by = -30000), 
             x1 = seq(gene_to_plot$plot_start-30000, gene_to_plot$plot_end+30000, by = -30000)-10000,
             y0 = gene_to_plot$y_gene, code = ifelse(derivative == FALSE, 2, 1),
             col = "gray", length = 0.06, angle = 60)
    }
    
    # plot the gene name
    text(x = (gene_to_plot$plot_start + gene_to_plot$plot_end) / 2,
         y = y0 + y_text,
         labels = gene_to_plot$hgnc_symbol, srt = 60, cex = 0.7,
         col = "darkblue")
    
    # plot the gene body (line from gene start to gene end)
    arrows(x0 = gene_to_plot$plot_start,
           x1 = gene_to_plot$plot_end ,
           y0 = gene_to_plot$y_gene,
           code = 0,
           length = 0.08, col = "darkblue")
    if(nrow(exons_to_plot) > 0){
      rect(xleft = exons_to_plot$plot_start[as.character(exons_to_plot$ensembl_gene_id) == as.character(gene_to_plot$ensembl_gene_id)],
           exons_to_plot$plot_start[as.character(exons_to_plot$ensembl_gene_id) == as.character(gene_to_plot$ensembl_gene_id)],
           ybottom = gene_to_plot$y_gene-0.2, ytop = gene_to_plot$y_gene+0.2, col = "darkblue", border = "darkblue", lwd = 0.6)
    }
  }
}

plot_TADs <- function(TADs,
                      Patient,
                      chr,
                      xlim,
                      color = "darkblue",
                      SVs = breakpoints,
                      y0 = 2,
                      label = "TADs",
                      der = FALSE){
  
  # input data needs to contain chr - start - end as first 3 columns
  TADs[,1] <- gsub("chr", "", TADs[,1])
  TADs$TAD_ID <- 1:nrow(TADs)
  names(TADs)[1:3] <- c("chr", "start", "end")
  
  TADs_g <- GRanges(seqnames = TADs[,1], IRanges(start = TADs[,2], end = TADs[,3]))
  
  SVs_patient <- SVs[which(SVs$Patient == Patient),]
  names(SVs_patient) <- paste("SV_", names(SVs_patient), sep = "")
  
  breakpoints_g <- GRanges(seqnames = SVs_patient$SV_chr, IRanges(start = SVs_patient$SV_end, end = SVs_patient$SV_end+1))
  
  if(der == FALSE){
    
    TADs_plot <- TADs[which(TADs[,1] == chr),]
    TADs_plot$TAD_start <- TADs_plot$start
    TADs_plot$TAD_end <- TADs_plot$end
    
    TADs_plot$code <- 3
    TADs_plot$Lwd <- 1
    TADs_plot$Lty <- 1
    
    # TADs overlapping with breakpoint will be highlighted with a red shade
    olap_TADs_BP <- findOverlaps(TADs_g, breakpoints_g)
    disrupted_TADs <- TADs[queryHits(olap_TADs_BP),]
    
    TADs_plot$Highlight <- ifelse(TADs_plot$TAD_ID %in% disrupted_TADs$TAD_ID, "yes", "no")
    
  } else {
    SVs_patient_chr <- SVs_patient[which(SVs_patient$SV_der_chr == chr),]
    
    fragments_g <- GRanges(seqnames = SVs_patient_chr$SV_chr, IRanges(start = SVs_patient_chr$SV_start, end = SVs_patient_chr$SV_end))
    olap_TADs_chr <- findOverlaps(TADs_g, fragments_g)
    TADs_plot <- cbind(TADs[queryHits(olap_TADs_chr),], SVs_patient_chr[subjectHits(olap_TADs_chr),])
    
    for(i in 1:nrow(TADs_plot)){
      TADs_plot$TAD_start[i] <- ifelse(grepl(pattern = "nver", TADs_plot$SV_type[i]) != TRUE, 
                                       TADs_plot$SV_der_start[i] + (TADs_plot$start[i] - TADs_plot$SV_start[i]),
                                       TADs_plot$SV_der_start[i] + (TADs_plot$SV_end[i] - TADs_plot$end[i]))
      TADs_plot$TAD_end[i] <- ifelse(grepl(pattern = "nver", TADs_plot$SV_type[i]) != TRUE, 
                                     TADs_plot$SV_der_start[i] + (TADs_plot$end[i] - TADs_plot$SV_start[i]),
                                     TADs_plot$SV_der_start[i] + (TADs_plot$SV_end[i] - TADs_plot$start[i]))
      
      # TADs overlapping with breakpoint will be highlighted with a red shade
      TADs_plot$Highlight[i] <- ifelse(TADs_plot$TAD_start[i] < TADs_plot$SV_der_start[i] | TADs_plot$TAD_end[i] > TADs_plot$SV_der_end[i], "yes", "no")
      
      if(TADs_plot$TAD_end[i] > TADs_plot$SV_der_end[i] & TADs_plot$TAD_start[i] < TADs_plot$SV_der_start[i]){
        #print(TADs_plot[i,])
        TADs_plot$code[i] <- 0
      } else if (TADs_plot$TAD_start[i] < TADs_plot$SV_der_start[i]) {
        TADs_plot$code[i] <- 2
      } else if (TADs_plot$TAD_end[i] > TADs_plot$SV_der_end[i]) {
        TADs_plot$code[i] <- 1
      } else {
        TADs_plot$code[i] <- 3
      }
      
      TADs_plot$TAD_start[i] <- ifelse(TADs_plot$TAD_start[i] < TADs_plot$SV_der_start[i], TADs_plot$SV_der_start[i], TADs_plot$TAD_start[i] )
      TADs_plot$TAD_end[i] <- ifelse(TADs_plot$TAD_end[i] > TADs_plot$SV_der_end[i], TADs_plot$SV_der_end[i], TADs_plot$TAD_end[i] )
      
      TADs_plot$Lty[i] <- ifelse(TADs_plot$Highlight[i] == "yes", 2,1)
      TADs_plot$Lwd[i] <- ifelse(TADs_plot$Highlight[i] == "yes", 0.6,1.5)
      
    }
  }
  
  TADs_plot <- TADs_plot[which(TADs_plot$TAD_start > (xlim[1]-4e6)),]
  
  TADs_plot <- TADs_plot[which(TADs_plot$TAD_start < (xlim[2]+4e6)),]
  
  #print(head(TADs_plot))
  #print(TADs_plot[,c("TAD_start","TAD_end")])
  for(i in 1:nrow(TADs_plot)){
    #print(i)
    arrows(x0 = TADs_plot$TAD_start[i], x1 = TADs_plot$TAD_end[i],code = as.numeric(TADs_plot$code)[i], y0 = y0, length = 0.05, angle = 90, col = color, lwd = TADs_plot$Lwd[i], lty = TADs_plot$Lty[i])
    
  }
  if(nrow(TADs_plot[which(TADs_plot$Highlight == "yes"), ]) > 0){
    rect(xleft = TADs_plot[which(TADs_plot$Highlight == "yes"), "TAD_start"], 
         xright = TADs_plot[which(TADs_plot$Highlight == "yes"), "TAD_end"], ybottom = y0-0.1, ytop = y0+0.1, border = NA, col = rgb(255,0,0,50, maxColorValue = 255))
  }
  
  mtext(side = 2, text = label, at = y0, las = 1, line = 0.2, cex = 0.8)
  
}



plot_enhancers <- function(y0 = 2,
                           height = 0.5,
                           enhancers,
                           chr,
                           xlim,
                           
                           score = 5,
                           bg = NULL,
                           fill_color = c(0,0,128),
                           border_color = c(0,0,128),
                           border_width = 0.8,
                           label = "",
                           patient = "33008",
                           derivative = FALSE,
                           threshold = 1,
                           breakpoints = SVs){
  
  enhancers[,1] <- gsub(pattern = "chr", replacement = "", enhancers[,1])
  
  if(is.null(bg) == FALSE){
    rect(xleft = xlim[1]-1e6, xright = xlim[2]+1e6, ybottom = y0, ytop = y0+height, 
         col = rgb(as.numeric(bg[1]),as.numeric(bg[2]),as.numeric(bg[3]), as.numeric(bg[4]), maxColorValue = 255), border = NA)
  }
  
  if(derivative == FALSE){
    
    enhancers_to_plot <- enhancers[enhancers[,1] == chr & enhancers[,score] > threshold,]
    
    rect(xleft = enhancers_to_plot[, 2],
         xright = enhancers_to_plot[, 3],
         ybottom = y0,
         ytop = y0 + height, col = rgb(fill_color[1], fill_color[2], fill_color[3], maxColorValue = 255), border = rgb(border_color[1], border_color[2], border_color[3], maxColorValue = 255),
         lwd = border_width)
  } else {
    patient_data <- breakpoints[breakpoints$Patient == patient,]
    for(fragment in patient_data[patient_data$der_chr == chr,"der_fragment"]){
      #print(fragment)
      
      fragment_g <- GRanges(seqnames = patient_data[patient_data$der_fragment == fragment,"chr"], 
                            IRanges(start = patient_data[patient_data$der_fragment == fragment,"start"],
                                    end = patient_data[patient_data$der_fragment == fragment,"end"]))
      
      enhancer_g <- GRanges(seqnames = enhancers[,1],
                            IRanges(start = enhancers[,2], end = enhancers[,3]),
                            score = enhancers[,score])
      
      olap <- findOverlaps(fragment_g, enhancer_g)
      
      overlapping_enhancers <- enhancer_g[subjectHits(olap),]
      
      if(patient_data[patient_data$der_fragment == fragment, "type"] != "Inversion"){
        starts <- start(overlapping_enhancers) - patient_data[patient_data$der_fragment == fragment,"start"] + 
          patient_data[patient_data$der_fragment == fragment,"der_start"]
        
        ends <- end(overlapping_enhancers) - patient_data[patient_data$der_fragment == fragment,"start"] + 
          patient_data[patient_data$der_fragment == fragment,"der_start"]
        
        enhancers_to_plot <- data.frame(der_start = starts,
                                        der_end = ends,
                                        score = overlapping_enhancers$score)
        
        enhancers_to_plot <- enhancers_to_plot[enhancers_to_plot$der_start > xlim[1]-1e6,]
        enhancers_to_plot <- enhancers_to_plot[enhancers_to_plot$der_end < xlim[2]+1e6,]
        enhancers_to_plot <- enhancers_to_plot[enhancers_to_plot$score > threshold,]
      } else {
        starts <- patient_data[patient_data$der_fragment == fragment,"end"] - start(overlapping_enhancers) + 
          patient_data[patient_data$der_fragment == fragment,"der_start"]
        
        ends <- patient_data[patient_data$der_fragment == fragment,"end"] - end(overlapping_enhancers) + 
          patient_data[patient_data$der_fragment == fragment,"der_start"]
        
        enhancers_to_plot <- data.frame(der_start = starts,
                                        der_end = ends,
                                        score = overlapping_enhancers$score)
        
        enhancers_to_plot <- enhancers_to_plot[enhancers_to_plot$der_start > xlim[1]-1e6,]
        enhancers_to_plot <- enhancers_to_plot[enhancers_to_plot$der_end < xlim[2]+1e6,]
        
        enhancers_to_plot <- enhancers_to_plot[enhancers_to_plot$score > threshold,]
      }
      
      
      if(nrow(enhancers_to_plot) > 0){
        rect(xleft = enhancers_to_plot$der_start,
             xright = enhancers_to_plot$der_end,
             ybottom = y0,
             ytop = y0 + height, col = rgb(fill_color[1], fill_color[2], fill_color[3], maxColorValue = 255), border = rgb(border_color[1], border_color[2], border_color[3], maxColorValue = 255),
             lwd = border_width)  
      }
    }
  }
  mtext(text = label, side = 2, line = 0, at = y0 + height /2, las = 1, cex = 0.8)
  #axis(side = 2, at = c(y0, y0 + height), labels = NA, col = "darkgray")
}



plot_interactions <- function(interactions = PCHiC_ESC,
                              genes = genes_g,
                              label = "PCHiC",
                              y0 = 4,
                              chromosome = chromosome,
                              bait_chr = "bait_chr",
                              bait_start = "bait_start",
                              bait_end = "bait_end",
                              bait_ID = interactions$bait_ID,
                              PIR_chr = "PIR_chr",
                              PIR_start = "PIR_start",
                              PIR_end = "PIR_end",
                              PIR_ID = "PIR_ID",
                              score = "Score",
                              score_threshold = 5,
                              lwd = 0.5,
                              rgb = c(255,0,0),
                              bow = 1,
                              bait = NULL,
                              xlim = xlim){
  
  #interactions = PCHiC_ESC
  interactions[,bait_chr] <- gsub("chr", "", interactions[,bait_chr])
  interactions[,PIR_chr] <- gsub("chr", "", interactions[,PIR_chr])
  
  # Select the interactions on the specific chromosome:
  interactions_to_plot <- interactions[which(interactions[,bait_chr] == chromosome),]
  
  # Only plot the interactions that are +- 5 Mb around the plot region (xlim)
  interactions_to_plot <- interactions_to_plot[which(interactions_to_plot[,bait_start] > (xlim[1] - 5e6)),]
  interactions_to_plot <- interactions_to_plot[which(interactions_to_plot[,bait_start] < (xlim[2] + 5e6)),]
  
  # Remove the interactions with a score lower than the score threshold
  interactions_to_plot <- interactions_to_plot[which(interactions_to_plot[,score] > score_threshold),]
  head(interactions_to_plot)
  
  ## This only selects one specific bait.
  if(is.null(bait) == FALSE){
    
    if(bait %in% genes$hgnc_symbol == TRUE){
      bait_info <- genes[genes$hgnc_symbol == bait]
      
      interactions_from_bait <- interactions_to_plot[which(interactions_to_plot[,bait_chr] == as.vector(seqnames(bait_info)) &
                                                           interactions_to_plot[,bait_start] > start(bait_info)-10000 & 
                                                           interactions_to_plot[,bait_end] < end(bait_info)+10000),]
    }
    if(nrow(interactions_from_bait) > 0){
      print(paste("Only plotting PCHiC interaction from bait: ", bait, sep = ""))
      interactions_to_plot <- interactions_from_bait
    } else {
      print(paste("# Warning: no PCHiC interactions found for bait: ", bait, sep = ""))
    }
  } 
  
  #print(head(interactions_to_plot))
  
  coords <- data.frame(x1 = (interactions_to_plot[,bait_start] + interactions_to_plot[,bait_end])/ 2, y1 = y0,
                       x2 = (interactions_to_plot[,PIR_start] + interactions_to_plot[,PIR_end]) / 2, y2 = y0)
  coords$score <- interactions_to_plot[,score]
  head(coords)
  
  for(i in 1:nrow(coords)){
    
    curve <- ifelse(coords[i,1] < coords[i,3], -0.0000005, 0.0000005)
    
    curvedarrow(from = c(coords[i,1], coords[i,2]),
                to =c(coords[i,3], coords[i,4]),
                curve = curve * bow,
                lwd = lwd, 
                #arr.pos = 0.5, arr.type = "simple", arr.length = 0.2,
                arr.length= 0,
                lcol = rgb(rgb[1],rgb[2],rgb[3], 200, coords[i, "score"] / max(coords$score) * 255, maxColorValue = 255))
  }
  
  mtext(side = 2, text = ifelse(is.null(bait) == TRUE, label, paste(label, "\n(",bait, ")", sep = "")), at = y0+bow/2, las = 1, line =0, cex = 0.8)
}


# This function is used to generate a virtual 4C profile based on the polygon function
plot_v4C <- function(gene, genes = Genes_SVs, cell_type, y0 = 5, height = 1.5, bin_size = 40000, color = "65,105,225", 
                     V4C_folder, 
                     label = "NPC"){
  
  # The V4C file should contain the cell type in the first column (first column will be row.names)
  V4C_file <- paste(V4C_folder, gene, ".txt", sep = "")
  if(file.exists(V4C_file)){
    print(paste("Reading: ", V4C_file, sep = ""))
    V4C_data <- read.delim(V4C_file, header = F, stringsAsFactors = F, row.names = 1)
    V4C_data <- V4C_data[cell_type,]

    bins <- seq(from = bin_size , to = ((ncol(V4C_data)) * bin_size), by = bin_size)
    names(V4C_data) <- bins
    
    transcription_start_site <- Genes_SVs[Genes_SVs$hgnc_symbol == gene, "transcription_start_site"]
    viewpoint_start <- max(bins[bins < transcription_start_site])
    viewpoint_end <- min(bins[bins > transcription_start_site])
    
    color <- as.numeric(unlist(strsplit(color, split = ",")))
    
    polygon(x = bins, y = y0 + V4C_data / max(V4C_data) * height, 
            col = rgb(color[1], color[2], color[3], 175, maxColorValue = 255), 
            border = rgb(color[1], color[2], color[3], 255, maxColorValue = 255), lwd = 1)
    
    # Plot a box around the viewpoint
    rect(xleft = viewpoint_start, xright = viewpoint_end, ybottom = y0-0.1, ytop = y0+height+0.1, lty = 2, lwd = 0.5)
    
    axis(side = 2, at = c(y0, y0 + height), labels = c(0, round(max(V4C_data),1)), col = "white", col.ticks = "black", cex.axis = 0.8, las = 1)
    mtext(side = 2, text = label, at = y0+height/2, las = 1, line = 0.2, cex = 0.8)
  } else {
    print(paste("! Cannot find V4C file: ", V4C_folder, gene, ".txt", sep = ""))
  }
}



define_plot_region <- function(patient,
                               SVs,
                               max_width_panel = 7e6,
                               derivative = FALSE){
  
  # First select the SVs in the patient
  
  print(paste("Determining the plot windows for patient ", patient, sep = ""))
  
  SVs_patient <- SVs[which(SVs$Patient == patient),]
  
  all_fragments <- data.frame()
  
  if(derivative == FALSE){
    SVs_patient$plot_chr = SVs_patient$chr
    SVs_patient$plot_start = SVs_patient$start
    SVs_patient$plot_end = SVs_patient$end
    SVs_patient$plot_fragment = SVs_patient$Fragment
  } else {
    SVs_patient$plot_chr = SVs_patient$der_chr
    SVs_patient$plot_start = SVs_patient$der_start
    SVs_patient$plot_end = SVs_patient$der_end
    SVs_patient$plot_fragment = SVs_patient$der_fragment
    SVs_patient <- SVs_patient[which(SVs_patient$type != "Deletion"),]
  }
  
  for(chromosome in unique(SVs_patient$plot_chr)){
    #print(chromosome)
    
    SVs_chr <- SVs_patient[which(SVs_patient$plot_chr == chromosome),]
    
    SVs_chr <- SVs_chr[order(SVs_chr$plot_fragment),]
    
    # Now add 2 Mb in before and after the fragments:
    for(fragment in 1:length(SVs_chr$plot_fragment)){
      #print(fragment)
      
      fragment_data <- SVs_chr[fragment,]
      #print(fragment_data)
      
      # The first fragment is (should be) a flanking region. Only select the breakpoint and not the start of the chromosome:
      if(fragment == 1){
        
        if(derivative == TRUE && fragment_data$type == "Deletion"){
          fragment_plot_start <-  fragment_data$plot_start - 1e6
          fragment_plot_end <- fragment_data$plot_end + 1e6
          fragment_info <- data.frame(chr = chromosome, panel_start = fragment_plot_start, panel_end = fragment_plot_end, stringsAsFactors = F)
        } else {
          fragment_plot_start <- fragment_data$plot_end - 2e6
          fragment_plot_end <- fragment_data$plot_end + 2e6
          fragment_info <- data.frame(chr = chromosome, panel_start = fragment_plot_start, panel_end = fragment_plot_end, stringsAsFactors = F)
        }
      } else if(fragment == length(SVs_chr$Fragment)) {
        # The last fragment should also be a flanking region. Only select the breakpoint (= start of the flank) and not the end of the chromosome:
        fragment_plot_start <- fragment_data$plot_start -2e6
        fragment_plot_end <- fragment_data$plot_start + 2e6
        fragment_info <- data.frame(chr = chromosome, panel_start = fragment_plot_start, panel_end = fragment_plot_end, stringsAsFactors = F)
        
        # Select the fragments in between the flanking regions:
      } else {
        
        # Some fragments in between breakpoints or inversions can be very long. We only want to show +- 2MB of the breakpoints. 
        if(fragment_data$type == "Normal" | fragment_data$type == "Inversion"){
          if(fragment_data$width > 4e6){
            subfragment1 <- data.frame(chr = chromosome,
                                       panel_start = fragment_data$plot_start - 2e6,
                                       panel_end = fragment_data$plot_start + 2e6,
                                       stringsAsFactors = F)
            subfragment2 <- data.frame(chr = chromosome,
                                       panel_start = fragment_data$plot_end - 2e6,
                                       panel_end = fragment_data$plot_end + 2e6,
                                       stringsAsFactors = F)
            fragment_info <- rbind(subfragment1, subfragment2)
            
          } else {
            fragment_plot_start <-  fragment_data$plot_start -2e6
            fragment_plot_end <- fragment_data$plot_end + 2e6
            fragment_info <- data.frame(chr = chromosome, panel_start = fragment_plot_start, panel_end = fragment_plot_end, stringsAsFactors = F)
          }
          
          
        } else {
          fragment_plot_start <-  fragment_data$plot_start -2e6
          fragment_plot_end <- fragment_data$plot_end + 2e6
          fragment_info <- data.frame(chr = chromosome, panel_start = fragment_plot_start, panel_end = fragment_plot_end, stringsAsFactors = F)
        }
      }
      all_fragments <- rbind(all_fragments, fragment_info)
    }
  }
  # Now we need to merge all overlapping fragments in one panel for plotting:
  all_fragments_g <- GRanges(seqnames = all_fragments$chr, IRanges(start = all_fragments$panel_start, end = all_fragments$panel_end))
  plot_region <- as.data.frame(reduce(all_fragments_g))
  
  # Some merged fragments may be too long to plot in just one panel. Split the fragments that are longer than the max_width_panel
  filtered_regions <- data.frame()
  
  for(regions in 1:nrow(plot_region)){
    region <- plot_region[regions,]
    #print(region)
    
    if(region$width > max_width_panel){
      # Determine into how many panels the region has to be split:
      number_subregions <- ceiling(region$width / max_width_panel)
      
      # Each split panel will have the same size:
      size_subregions <- region$width / number_subregions
      
      # Calculate the start and ends of each split panel:
      final_region <- data.frame()
      for(i in 0:(number_subregions-1)){
        if(i == 0){
          start_subregion <- region$start + i * size_subregions
          end_subregion <- region$start + (i+1) * size_subregions + 0.5e6
        } else if (i == (number_subregions-1)){
          start_subregion <- region$start + i * size_subregions - 0.5e6
          end_subregion <- region$start + (i+1) * size_subregions + 0.5e6
        } else {
          start_subregion <- region$start + i * size_subregions - 0.5e6
          end_subregion <- region$start + (i+1) * size_subregions
        }
        subregion <- data.frame(seqnames = region$seqnames, start = start_subregion, end = end_subregion, width = size_subregions, strand = "*", stringsAsFactors = F)
        final_region <- rbind(final_region, subregion)
      }
    } else (
      final_region <- region
    )
    filtered_regions <- rbind(filtered_regions, final_region)
  }
  filtered_regions$order <- gsub(pattern = "der|chr", "", filtered_regions$seqnames)
  filtered_regions <- filtered_regions[order(filtered_regions$order, filtered_regions$start),]
  
  print("Regions to plot:")
  print(filtered_regions)
  
  plot_input <- NULL
  for(i in 1:nrow(filtered_regions)){
    
    panel <- filtered_regions[i,1:3]
    
    panel$seqnames <- as.character(panel$seqnames)
    
    plot_input <- c(plot_input, as.character(panel))
    
  }
  
  return(plot_input)
}

autoplot <- function(patient, derivative_mode = FALSE){
  
  plot_regions <- define_plot_region(patient = patient, derivative = derivative_mode)
  
  number_of_panels <- length(plot_regions)/3
  number_of_pages <- ceiling(number_of_panels/4)
  
  if(number_of_pages > 1){
    for(i in 0:(number_of_pages-1)){
      #print(i)
      max <- 12+i*12
      #print(max)
      if(max > (number_of_panels*3)){
        range <- (1 + i*12) : (number_of_panels*3)
      } else {
        range <- (1 + i*12) : (12 + i*12)
      }
      #print(range)
      plot_regions_page <- plot_regions[range]
      print(plot_regions_page)
      
      print(paste("Plotting page ", i+1, " of ",number_of_pages, sep = ""))
      
      patient_overview_plot(panels = plot_regions_page, patient = patient, derivative = derivative_mode)
    }
  }
}

