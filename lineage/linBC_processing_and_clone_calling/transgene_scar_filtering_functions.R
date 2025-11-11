# Description ####
# Functions for filtering of transgenic target lineage barcodes

# Written mainly by Bastiaan Spanjaard (merged and altered by Nora Fresmann)

suppressPackageStartupMessages(library(ggplot2))
require(igraph)
suppressPackageStartupMessages(library(reshape2))
require(stringdist)
require(plyr)
require(data.table)
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(DescTools))
suppressPackageStartupMessages(library(pheatmap))
                       

## Function for filtering lineage barcode reads on RFP or other transgenic targets

filter_RFP_scars <- function(unfilt_scar_file_1, unfilt_scar_file_2,  all_barcode_UMIs_file, all_cells){

    scars.unfiltered <- unfilt_scar_file_1
    all.scars.g1 <- unfilt_scar_file_2
    all.barcode.UMIs <- all_barcode_UMIs_file
    cells <- all_cells
    
    # Filter scars per cell ####
    # Go through all cells one by one. Within each cell, calculate the Hamming
    # distances between the sequences; if these are one or two and the read 
    # difference is high (log2(read1/read2) > log2.cutoff), flag the scar for 
    # removal.
    for(c in 1:length(cells)){
      cell <- cells[c]
      cell.scars <- scars.unfiltered[scars.unfiltered$Cell == cell, ]
      if(nrow(cell.scars) < 2){next}
      cell.scars <- cell.scars[order(-cell.scars$Reads), ]
      n.scars <- nrow(cell.scars)
      cell.seqdist <- stringdistmatrix(cell.scars$Sequence, method = "hamming")

      # Identify all sequences that have to be removed.
      for(k in 1:length(cell.seqdist)){
        if(cell.seqdist[k] >2){next}
        else{
          # Get scar indices. Calculate log2 ratio between two reads. If too high,
          # flag lowest scar (Keep = F).
          scar.indices <- get.dist.index(k, n.scars)
          high.scar <- min(scar.indices)
          low.scar <- max(scar.indices)
          if(log2(cell.scars$Reads[high.scar]/cell.scars$Reads[low.scar]) > log2.cutoff){
            scars.unfiltered$Keep[scars.unfiltered$Scar.id == cell.scars$Scar.id[low.scar]] <-
              F
          }
        }
      }
    }


    # Remove scars that have HD 1 to scars with much higher reads.
    scars.filter.1 <- scars.unfiltered[scars.unfiltered$Keep, ]



    # Go throught the barcodes and calculate the distances between the scars again. 
    # If the HD is one and the read difference is low, mark the scars as suspect. 
    # If the HD is two and the read difference is high, mark the scars as suspect.
    for(c in 1:length(cells)){
      cell <- cells[c]
      cell.scars <- scars.filter.1[scars.filter.1$Cell == cell, ]
      if(nrow(cell.scars) < 2){next}
      cell.scars <- cell.scars[order(-cell.scars$Reads), ]
      n.scars <- nrow(cell.scars)
      cell.seqdist <- stringdistmatrix(cell.scars$Sequence, method = "hamming")

      for(k in 1:length(cell.seqdist)){
        if(cell.seqdist[k] >2){next}
        else if(cell.seqdist[k] == 1){
          # Get scar indices. Flag both (paste(Pair, other
          # scar)).
          scar.indices <- get.dist.index(k, n.scars)
          high.scar <- min(scar.indices)
          low.scar <- max(scar.indices)
          # If the scars are legitimate scars that are known to be close in HD,
          # skip this pair.
          # if((cell.scars$Sequence[high.scar] %in% allowed.close.scars) &
          #    (cell.scars$Sequence[low.scar] %in% allowed.close.scars)){next}
          scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]] <-
            paste(scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]],
                  scars.filter.1$Scar.id[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]],
                  sep = ".")
          scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]] <-
            paste(scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]],
                  scars.filter.1$Scar.id[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]],
                  sep = ".")
          }
        else if(cell.seqdist[k] == 2){
          # Get scar indices. Calculate log2 ratio between two reads. If too high,
          # flag both (paste(Pair, other scar)).
          scar.indices <- get.dist.index(k, n.scars)
          high.scar <- min(scar.indices)
          low.scar <- max(scar.indices)
          if(log2(cell.scars$Reads[high.scar]/cell.scars$Reads[low.scar]) > log2.cutoff){
            scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]] <-
              paste(scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]],
                    scars.filter.1$Scar.id[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]],
                    sep = ".ratio.")
            scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]] <-
              paste(scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]],
                    scars.filter.1$Scar.id[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]],
                    sep = ".ratio.")
          }
        }
      }
    }


    # Assess and remove suspect scars ####
    # First determine which scars may need to be filtered out - occurring more than
    # once but less than their possible parent scar in the same cell.
    scars.assess <- scars.filter.1[scars.filter.1$Pair != "With", ]
    scars.filter.2 <- scars.filter.1[, 1:7]
    scars.assess$Incidence <-
      sapply(scars.assess$Sequence,
             function(x) sum(scars.filter.1$Sequence == x))
    scars.assess$Pair.incidence <-
      sapply(scars.assess$Pair,
             function(x){
               pairs <- unlist(strsplit(x, "[.]"))[-1]
               return(min(scars.assess$Incidence[scars.assess$Scar.id %in% pairs]))}
      )
    scars.assess$Min.pair.incidence <-
      apply(scars.assess[, c("Incidence", "Pair.incidence")], 1, min)
    scars.assess.2 <- scars.assess[scars.assess$Min.pair.incidence > 1, ]

    # Determine whether scars also occur in cells without the parent sequence 
    # (criterion 1)
    seq.freq <- data.frame(table(scars.assess.2$Sequence))
    colnames(seq.freq)[1] <- "Sequence"
    scars.assess.2 <- merge(scars.assess.2, seq.freq)
    scars.assess.2$Crit.1 <- 
      !(scars.assess.2$Min.pair.incidence == scars.assess.2$Freq)

    # Determine whether scars occur with more than one UMI (criterion 2)
    scars.assess.2 <- merge(scars.assess.2, all.barcode.UMIs)
    scars.assess.2$Crit.2 <- (scars.assess.2$UMIs > 1)

    # Determine whether scars have at least one UMI unrelated to the UMIs of the
    # parent sequence (criterion 3).
    # Select scars to test - the ones that have the lowest incidence in the full
    # dataset.
    scars.assess.2.low.inc <- 
      scars.assess.2[scars.assess.2$Incidence <= scars.assess.2$Pair.incidence, ]
    # If criteria 1 and 2 are both false, the scars are removed automatically
    scars.assess.2.out <- 
      scars.assess.2.low.inc[!scars.assess.2.low.inc$Crit.1 & 
                               !scars.assess.2.low.inc$Crit.2, ]
    # Else if one of criteria 1 and 2 are true, the scars are maybe kept, if
    # they have at least one UMI that is at least 2 HD away from its suspected
    # parent scar UMIs.
    scars.assess.2.maybe <- 
      scars.assess.2.low.inc[scars.assess.2.low.inc$Crit.1 | 
                               scars.assess.2.low.inc$Crit.2, ]

    if(nrow(scars.assess.2.maybe) == 0){

        scars.output <- scars.filter.1[, c(1:6, 8)]

    }else{

        # Loop over all scars to calculate the number of UMIs they have that are more
        # than 1 HD away from UMIs of their 'parent' scar
        scars.assess.2.maybe$Crit.3 <- NA

        for(c.s in 1:nrow(scars.assess.2.maybe)){
              #print(paste(c.s, "out of", nrow(scars.assess.2.maybe)))
              c.sequence <- scars.assess.2.maybe$Sequence[c.s]
              c.cell <- scars.assess.2.maybe$Cell[c.s]
              c.scar <- scars.assess.2.maybe$Scar.id[c.s]
              c.parent.sequence <- 
                scars.assess$Sequence[grepl(paste(c.scar, "\\.", sep = ""), scars.assess$Pair) |
                                          grepl(paste(c.scar, "$", sep = ""), scars.assess$Pair)]

              validate.seq <- 
                all.scars.g1[all.scars.g1$Cell == c.cell &
                               all.scars.g1$Sequence %in% c.sequence, ]
              validate.parent <- 
                all.scars.g1[all.scars.g1$Cell == c.cell &
                               all.scars.g1$Sequence %in% c.parent.sequence, ]
              validate.seq$Min.HD <-
                sapply(validate.seq$UMI,
                       function(x) min(stringdistmatrix(x, validate.parent$UMI, 
                                                        method = "hamming")))
              scars.assess.2.maybe$Crit.3[c.s] <- sum(validate.seq$Min.HD > 1)
        }

        scars.assess.2.maybe$Out <-ifelse(scars.assess.2.maybe$Crit.3 ==0 , T, F)
        # Remove sequencing error scars from the dataset
        sequencing.error.scars <- c(scars.assess.2.out$Scar.id, 
            scars.assess.2.maybe$Scar.id[scars.assess.2.maybe$Out])

        scars.output <- scars.filter.1[!(scars.filter.1$Scar.id %in% sequencing.error.scars),
                                       c(1:6, 8)]

    }


    return(scars.output)

}


                                   
# Function for adding probabilities for scars in the data by comparing to the frequency the scar is found in in bulk sequencing data (from larvae).
# If probabilities of certain genes are not available, they will automatically be assigned a low probability to avoid being filtered out.

compare_scars_transgene <- function(input_dat_name, dat_name,  min.presence, gene_it, scar_probabilities_in){

            
        Z1.scars <- input_dat_name
        
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

        
    
        # Compare presence with probabilities ####
        # For RFP1
        unique.scars$Sequence.short <- 
          substr(unique.scars$Sequence, 25, 75)

        
        # Load scar probabilities. Because bulk scar sequencing is done differently than
        # single-cell scar sequencing, some single cell scars cannot be assigned a
        # probability because they cannot be observed in bulk sequencing. We filter 
        # these out.

        scar.probabilities <- scar_probabilities_in
                                       
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


         