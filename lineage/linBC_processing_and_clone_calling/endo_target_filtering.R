# Description ####
# Functions for UMI-based filtering of lineage barcodes on endogenous target genes

# Written by Bastiaan Spanjaard, wrapped up here by Nora Fresmann

suppressPackageStartupMessages(library(ggplot2))
require(igraph)
suppressPackageStartupMessages(library(reshape2))
require(stringdist)
require(plyr)
require(data.table)
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(plyr))
#suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(DescTools))



# Wrapper for running the UMI filtering functions below

UMI_filtering_wrap <- function(input_dat_name, input_dat_wd, dat_name, fraction, gene_names, seur_dat, scars_from_trans = F, read_cutoffs){

    # Bring together gene names and fractions
    fraction <- fraction
    names(fraction) <- gene_names
    
    
    if(scars_from_trans == T){
    
        all.scars.g1 <- read.table(paste0(input_dat_wd,"/",input_dat_name),
                             stringsAsFactors = F, sep = "\t")
        
    }else{
        # In case the full dataset that was obtained from the whitelist extraction approach are used
        all.scars.g1 <- read.table(paste0(input_dat_wd,"/",input_dat_name),
                             stringsAsFactors = F, sep = "\t")

        all.scars.g1$V11 <- trimws(all.scars.g1$V1)
        all.scars.g1$Reads <- sapply(all.scars.g1$V11,
                                       function(x) unlist(strsplit(x, " "))[1])
        all.scars.g1$Barcode <- sapply(all.scars.g1$V11,
                                       function(x) unlist(strsplit(x, " "))[2])
        all.scars.g1$V1 <- NULL
        all.scars.g1$V11 <- NULL

        colnames(all.scars.g1) <- c("UMI", "Gene", "CIGAR", "Sequence","Reads","Barcode")
        all.scars.g1$Reads <- as.numeric(all.scars.g1$Reads)
    }
  
    full_dataset <- all.scars.g1

    full_dataset <- as.data.table(full_dataset[, c("Reads", "Barcode", "UMI", "Gene", "CIGAR", "Sequence")])

    scars <- InitiateFiltering(full_dataset, cell_number = ncol(dat))
    scars <- ReadFilter(scars, read_cutoffs = read_cutoffs)
    scars <- CountUMIs(scars)
    scars <- FilterUMIs(scars, fraction = fraction)

   return(scars)                          
}

           
                                   
# Function for adding probabilities for scars in the data by comparing to the frequency the scar is found in in bulk sequencing data (from larvae).
# If probabilities of certain genes are not available, they will automatically be assigned a low probability to avoid being filtered out.

compare_scars <- function(input_dat_name, dat_name, gene_names, scar_probabilities_in, seur_dat, primer_set, min.presence, ...){

    all.scars.g1 <- input_dat_name
    all.scars.g1 <- all.scars.g1[all.scars.g1$Barcode %in% seur_dat$Barcode,]
    
    print(paste0('Starting with ', nrow(all.scars.g1),' entries in scar file.'))
    
    for(i in 1:length(gene_names)){
    
        gene_it <- gene_names[i]
        
        print(paste0('Comparing gene ', gene_it))
        
        Z1.scars <- f_scars_in_trnscrptm[f_scars_in_trnscrptm$Gene == gene_it,]
        
        
        if(nrow(Z1.scars) == 0){
        
            print(paste0(gene_it, " not found! Moving to next gene."))
            next
        }
        
        # Count presence of scars ####
        unique.scars <- data.frame(table(Z1.scars$Sequence))
        colnames(unique.scars) <- c("Sequence", "Freq.Z1")

        # This comes after merging scars from several datasets - maybe it is not needed here
        unique.scars[is.na(unique.scars)] <- 0
        unique.scars$Presence <- apply(as.data.frame(unique.scars[, -1]), 1,
                                       function(x) sum(x >= min.presence))

        all.CIGARs <- unique(Z1.scars[, c("Sequence", "CIGAR")])
        all.CIGARs <- all.CIGARs[!duplicated(all.CIGARs$Sequence), ]
        unique.scars <- merge(unique.scars, all.CIGARs)
        unique.scars$Sequence <- as.character(unique.scars$Sequence)
        
    
        # Compare presence with probabilities ####
        
         if(primer_set == "01"){                              
                                       
            if(gene_it == "actb1"){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 5, 75)

            }else if(gene_it == "actb2"){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 1, 73)

            }else if(gene_it == "rpl39"){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 14, 75)

            }else if(gene_it == "rpl18a"){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 6, 75)
            }else if(gene_it == "cfl1"){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 1, 75)

                unique.scars$p <- 0.0001
                unique.scars$Embryos <- 0
                unique.scars$p[unique.scars$Sequence.short == wildtype.clip] <- max(unique.scars$p, na.rm = T) + 0.1
                unique.scars <- unique.scars[order(-unique.scars$Presence, -unique.scars$p), ]
                unique.scars$Scar <- paste(0:(nrow(unique.scars)-1), unique.scars$CIGAR, sep = ":")

                col_order <- c('Sequence.short','Sequence','Freq.Z1','Presence','CIGAR','p','Embryos','Scar')
                unique.scars_2 <- unique.scars[, col_order]
            } else if(gene_it == "dsRedRecCas"){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 1, 116)

                unique.scars$p <- 0.0001
                unique.scars$Embryos <- 0
                unique.scars$p[unique.scars$Sequence.short == wildtype.clip] <- max(unique.scars$p, na.rm = T) + 0.1
                unique.scars <- unique.scars[order(-unique.scars$Presence, -unique.scars$p), ]
                unique.scars$Scar <- paste(0:(nrow(unique.scars)-1), unique.scars$CIGAR, sep = ":")

                col_order <- c('Sequence.short','Sequence','Freq.Z1','Presence','CIGAR','p','Embryos','Scar')
                unique.scars_2 <- unique.scars[, col_order]
            }


        }else if(primer_set == "02"){
            
            if(gene_it %in% c('actb1','actb2','rpl39')){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 1, 75)

            }else if(gene_it %in% c('cfl1','cirbpb','ube2e1')){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 1, 75)

                unique.scars$p <- 0.0001
                unique.scars$Embryos <- 0
                unique.scars$p[unique.scars$Sequence.short == wildtype.clip] <- max(unique.scars$p, na.rm = T) + 0.1
                unique.scars <- unique.scars[order(-unique.scars$Presence, -unique.scars$p), ]
                unique.scars$Scar <- paste(0:(nrow(unique.scars)-1), unique.scars$CIGAR, sep = ":")

                col_order <- c('Sequence.short','Sequence','Freq.Z1','Presence','CIGAR','p','Embryos','Scar')
                unique.scars_2 <- unique.scars[, col_order]
        
            }
            
        }else{
            print("Please supply a valid name for the primer set used! At the moment options are Nora or Nina.")
            break   
        }
            
            
        if(!gene_it %in% c('cfl1','cirbpb','ube2e1','dsRedRecCas')){
            # Load scar probabilities. Because bulk scar sequencing is done differently than
            # single-cell scar sequencing, some single cell scars cannot be assigned a
            # probability because they cannot be observed in bulk sequencing. We filter 
            # these out.

            scar.probabilities <- scar_probabilities_in[scar_probabilities_in$Name == gene_it,]

            if(primer_set == "Nora"){      

                if(gene_it == "actb1"){
                    scar.probabilities$Sequence <- substr(scar.probabilities$Sequence,1,71)

                }else if(gene_it == "actb2"){
                    scar.probabilities$Sequence <- substr(scar.probabilities$Sequence,3,75)

                }else if(gene_it == "rpl39"){
                    scar.probabilities$Sequence <- substr(scar.probabilities$Sequence,1,62)

                }else if(gene_it == "rpl18a"){
                    scar.probabilities$Sequence <- substr(scar.probabilities$Sequence,1,75)

                }
            }


            duplicated_scars <- scar.probabilities[duplicated(scar.probabilities$Sequence) == T,]
            duplicated_scars <- scar.probabilities[scar.probabilities$Sequence %in% duplicated_scars$Sequence,]
            duplicated_scars <- duplicated_scars[order(duplicated_scars$Sequence),]
            length(unique(duplicated_scars$CIGAR))
            length(unique(duplicated_scars$Sequence))
            # Even before clipping the sequences in the scar-probabilities a little bit, there are duplicated sequences! This will lead to the creation of more than one row for a unique scar, when datasets are merged below! I'll therefore remove duplicated sequences, keeping the one that has the highest number of "Embryos".

            scar.probabilities <- scar.probabilities[order(scar.probabilities$Embryos, decreasing = T),]
            scar.probabilities <- scar.probabilities[order(scar.probabilities$Sequence, decreasing = T),]
            scar.probabilities <- scar.probabilities[duplicated(scar.probabilities$Sequence) == F,]

            # MERGE
            unique.scars_2 <- merge(unique.scars, scar.probabilities[, c("Sequence", "p", "Embryos")],
                                  by.x = "Sequence.short", by.y = "Sequence", all.x = T)
            unique.scars_2$Embryos[is.na(unique.scars_2$Embryos)] <- 0
            unique.scars_2$p[is.na(unique.scars_2$p)] <- min(unique.scars_2$p, na.rm = T)/100
            unique.scars_2$p[unique.scars_2$Sequence.short == wildtype.clip] <- max(unique.scars_2$p, na.rm = T) + 0.1
            unique.scars_2 <- unique.scars_2[order(-unique.scars_2$Presence, -unique.scars_2$p), ]
            unique.scars_2$Scar <- paste(0:(nrow(unique.scars_2)-1), unique.scars_2$CIGAR, sep = ":")
            
        }
        
        # Count total and unique scars per library ####
        unique.scars.nowt <- unique.scars_2[-1, ]
        scars.Z1 <- unique.scars.nowt$Scar[unique.scars.nowt$Freq.Z1 >= min.presence]

        sum(scars.Z1 %in% unique.scars.nowt$Scar[unique.scars.nowt$Presence == 1])

        # Write output ####
        Z1.scars.compared <- 
          merge(Z1.scars, unique.scars_2[, c("Sequence", "Presence", "p", "Embryos", "Scar")], by = "Sequence", all.x = T)
         
        write.csv(Z1.scars.compared, file = paste0("Z1_scars_compared_",dat_name,"_",gene_it,".csv"),
                  quote = F, row.names = F)
                
    }
    
}







# Functions for UMI-based filtering of lineage barcodes (scars) on endogenous targets extracted with the whitelist-approach in the lineage barcode extraction pipeline
InitiateFiltering <- function(full_dataset, cell_number){
      # Initiates a 'scars' object for filtering. This object includes number of transcriptome
      # cells and the initial scars object, together with a count of how many cells have at 
      # least one scar on a specific gene in absolute numbers (scars$Statistics) and
      # percentages (scars$Perc_statistics). Finally, it counts the number of sequences per 
      # cell per gene (scars$Stats_per_cell)

      scars <- list(Parameters = list(WT_cells = cell_number),
                    Input = full_dataset)
      scars$Statistics <- data.frame(table(unique(full_dataset[, c("Barcode", "Gene")])$Gene))
      colnames(scars$Statistics)[2] <- "Cells_after_whitelisting"
      scars$Perc_statistics <- scars$Statistics
      scars$Perc_statistics$Cells_after_whitelisting <- 100 * scars$Statistics$Cells_after_whitelisting/cell_number

      # filtering_results <- full_dataset %>%
      #   group_by(Barcode, Gene) %>%
      #   summarise(Whitelisted_sequences = n_distinct(Sequence))
      filtering_results <- full_dataset[, .(Whitelisted_sequences = uniqueN(Sequence)), by = .(Barcode, Gene)]
      scars$Stats_per_cell <- filtering_results

      return(scars)
}



ReadFilter <- function(scars, read_cutoffs){
      # Keep one sequence per umi+gene, does this lose us any barcodes? What's the rationale here?
      # Afterwards, removes any sequences with read_cutoff or less reads; used for one-time
      # sequencing error removal.

      # 131976 in full_datasetf before removing scars with only 1 read, 132196 when also taking
      # barcodes into account in Nina1. I am keeping the barcode in here for now.

      full_dataset <- scars$Input
      # full_dataset <- full_dataset %>%
      #   group_by(Gene, UMI) %>% 
      #   mutate(Max_reads = max(Reads)) %>% ungroup()
      full_dataset <- full_dataset[, ':=' (Max_reads = max(Reads)), by = .(Gene, UMI, Barcode)]

      full_datasetf <- full_dataset

      for(gene_num in 1:length(read_cutoffs)){
          gene_it <- names(read_cutoffs)[gene_num]
          read_cutoff <- read_cutoffs[gene_num]
          
          print(gene_it)
          print(read_cutoff)
          full_datasetf <- full_datasetf[full_datasetf$Reads == full_datasetf$Max_reads, ]
          full_datasetf_i <- full_datasetf[full_datasetf$Gene == gene_it & full_datasetf$Reads > read_cutoff, ]
          
          if(gene_num == 1){
              full_datasetf_all <- full_datasetf_i
          }else{
              full_datasetf_all <- rbind(full_datasetf_all,full_datasetf_i)
          }
          
       }
      # Add a minimal read filter just in case any gene was missed
      full_datasetf <- full_datasetf_all[full_datasetf_all$Reads > 1, ]

      scars$Read_filtered <- full_datasetf
      scars$Parameters$Read_cutoff <- read_cutoff
      return(scars)
}



CountUMIs <- function(scars){
      # Collate scars with the same sequence in the same cell on the same gene. Count
      # UMIs per scar and count scars per gene+cell.

      cell_number <- scars$Parameters$WT_cells
      full_datasetf <- scars$Read_filtered

      full_datasetfr <- full_datasetf[, .(UMIs = .N), by = .(Barcode, Gene, CIGAR, Sequence)]

      # full_datasetfr <- full_datasetf[, -(which(colnames(full_datasetf) %in% c("Reads", "Max_reads")))] %>%
      #   group_by(Barcode, Gene, CIGAR, Sequence) %>%
      #   summarise(UMIs = n())
      full_datasetfr <- full_datasetfr[order(-full_datasetfr$UMIs), ]
      scars$UMI_counted <- full_datasetfr

      # filtering_results <- full_datasetfr %>%
      #   group_by(Barcode, Gene) %>%
      #   summarise(Max_read_sequences = n_distinct(Sequence))
      filtering_results <- full_datasetfr[, .(Max_read_sequences = uniqueN(Sequence)), by = .(Barcode, Gene)]

      scars$Stats_per_cell <- merge(scars$Stats_per_cell, filtering_results)

      scars$Statistics <- merge(scars$Statistics,
                                data.frame(table(filtering_results$Gene)))
      colnames(scars$Statistics)[3] <- "Cells_after_readfiltering"

      scars$Perc_statistics <- scars$Statistics
      scars$Perc_statistics[, -1] <- 100 * scars$Statistics[, -1]/cell_number

      return(scars)
}



FilterUMIs <- function(scars, fraction){
      # For each sequence in each cell, calculate how big of a fraction of the total transcripts
      # of that gene in that cell this sequence contains. Calculate the cumulative sums of UMI
      # fractions and keep only high-UMI scars that contain 80% (or other, indicated in the fraction)
      # of the UMIs for that gene+barcode.
      # After that, count scars on endogenous genes. If a gene has >2 scars in a cell, that gene
      # is removed. If a cell has more than one gene with >2 scars, the cell is flagged as doublet
      # and removed altogether.
      # Filtering results are kept in scars$Filtering_output.
      # Statistics are updated with the following fields: after_UMIfilter, indicating the fraction-thresholded
      # scars, after_Ambigene_removal, indicating the scars after removal of genes with too many
      # scars, and after_doublet_removal, indicating the scars after removal of doublet cells.

      scars$Parameters$UMI_fraction <- fraction
      fraction_dt <- data.table(Gene = names(fraction),
                                Fraction = fraction)

      full_datasetfr <- scars$UMI_counted
      # full_dataset_UMIf <- full_datasetfr %>%
      #   group_by(Barcode, Gene) %>%
      #   mutate(UMI_fraction = UMIs/sum(UMIs)) %>% 
      #   mutate(UMI_cfraction = cumsum(UMI_fraction)) %>% 
      #   mutate(UMI_cf_min = min(UMI_cfraction[UMI_cfraction > fraction])) %>% ungroup()
      full_dataset_UMIf <- 
        full_datasetfr[, ':=' (UMI_fraction = UMIs/sum(UMIs)), by = .(Barcode, Gene)]
      full_dataset_UMIf <- full_dataset_UMIf[, ':=' (UMI_cfraction = cumsum(UMI_fraction)), by = .(Barcode, Gene)]
      full_dataset_UMIf <- merge(full_dataset_UMIf, fraction_dt, by = "Gene")

      full_dataset_UMIf <- full_dataset_UMIf[, ':=' (UMI_cf_min = min(UMI_cfraction[UMI_cfraction > Fraction])), by = .(Barcode, Gene)]

      full_dataset_UMIff <- full_dataset_UMIf[full_dataset_UMIf$UMI_cfraction <= full_dataset_UMIf$UMI_cf_min, ]
      scars$UMI_filtered <- full_dataset_UMIff

      ffiltering_results <- full_dataset_UMIff[, .(UMI_filtering_sequences = uniqueN(Sequence)), by = .(Barcode, Gene)]

      # ffiltering_results <- full_dataset_UMIff %>%
      #   group_by(Barcode, Gene) %>%
      #   summarise(UMI_filtering_sequences = n_distinct(Sequence))

      # Remove genes that still have more than two scars
      ffiltering_results$Over <- !(ffiltering_results$Gene %in% c("RFP", "RFP1", "RFP3")) & (ffiltering_results$UMI_filtering_sequences > 2)
      # Remove cells that look like doublets (more than two scars in more than one endogenous gene)
      # doublets <- ffiltering_results %>%
      #   group_by(Barcode) %>%
      #   summarise(Doublet = (sum(Over) > 1))
      doublets <- ffiltering_results[, .(Doublet = sum(Over) > 1), by = .(Barcode)]
      ffiltering_results <- merge(ffiltering_results, doublets, all = T)
      scars$Stats_per_cell <- merge(scars$Stats_per_cell, ffiltering_results, all = T,by = c("Gene", "Barcode"))
      scars$Output_stats_per_cell <- scars$Stats_per_cell[!scars$Stats_per_cell$Doublet & !scars$Stats_per_cell$Over, ]


      scars$Doublets <- doublets

      results <- merge(full_dataset_UMIff, ffiltering_results[!ffiltering_results$Over, c("Barcode", "Gene")],
                       by = c("Barcode", "Gene"))
      results <- results[results$Barcode %in% doublets$Barcode[!doublets$Doublet], ]
      scars$Filtering_output <- results

      # scars$Statistics <- merge(scars$Statistics,
      #                           data.frame(table(filtering_results$Gene)))
      scars$Statistics <- merge(scars$Statistics,
                                data.frame(table(ffiltering_results$Gene)), all = T)
      colnames(scars$Statistics)[4] <- "Cells_after_UMIfilter"
      scars$Statistics <- merge(scars$Statistics,
                                data.frame(table(ffiltering_results$Gene[!ffiltering_results$Over])), all = T)
      colnames(scars$Statistics)[5] <- "Cells_after_Ambigene_removal"
      scars$Statistics <- merge(scars$Statistics,
                                data.frame(table(ffiltering_results$Gene[!ffiltering_results$Over & !ffiltering_results$Doublet])), all = T)
      colnames(scars$Statistics)[6] <- "Cells_after_doublet_removal"
      scars$Statistics[is.na(scars$Statistics)] <- 0

      scars$Perc_statistics <- scars$Statistics
      scars$Perc_statistics[, -1] <- 100 * scars$Statistics[, -1]/scars$Parameters$WT_cells

      return(scars)
}






                
                