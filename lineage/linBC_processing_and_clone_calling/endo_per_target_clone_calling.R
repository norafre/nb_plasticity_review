# DESCRIPTION
# Functions for grouping cells according to lineage barcode similarity using sequencing data for individual endogenous CRISPR-Cas9 lineage tracing targets:
# Consists of the following functions, carried out for one target gene at a time:
# 1. function 'getclones':
#    Lineage barcode reads from the extraction pipeline are shortened to 30 bp sequence IDs.
#    Cells that have two distinct lineage barcode sequence IDs (i.e. two distinct alleles) are selected.
#    sequence IDs that only appear once are removed. sequence IDs that appear together across many cells are grouped into sequence ID-combinations.
#    sequence ID pairs with sequence IDs that are ambiguous due to co-occurrence with multiple other sequence IDs are removed.
#    Group cells carrying a valid sequence ID pair into 'clones' according to their sequence ID pair.
#    Add cells that only have one sequence ID, if this sequence ID was defined as a clone-defining sequence ID in the previous analysis steps.
#    The output of this function represents the most confident grouping of cells according to one target and is most suited to be used
#    as input for clone calling across all target genes.
#
# 2. function 'addwtclones':
#    Some cells may not have been edited by Cas9 at all and carry two uncut alleles of the target gene.
#    These 'wildtype' cells are not useful for clone calling. But it is useful to identify cells that likely did not incur
#    any edits, to check whether a lack of recovered clones is due to low editing efficiency (i.e. many cells with uncut alleles).
#    Data for all cells is loaded again and lineage barcode reads are shortened into 30 bp sequence IDs.
#    Keep only cell-barcode-sequence IDs combinations that are in the top 70 % of UMI-counts across all cell-barcode-sequence IDs combinations.
#    This ensures that only cells with multiple uncut sequence ID UMIs are assigned as true wildtype cells.
#    Merge the output with the output of 'getclones' to get a full list of cells assigned to a sequence ID-based clone or wildtype.
#
# 3. function 'addonescarclones':
#    Similar to 'wildtype' alleles, cells may have the same edit on both alleles of a gene or one allele may co-occur with another,
#    which isn't detectable. These scenarios are more unlikely, but we can still check if there are any cells that follow this pattern.
#    Data for all cells is loaded again and lineage barcode reads are shortened into 30 bp sequence IDs.
#    Keep only cells with a single sequence ID that is not part of a sequence ID pair used to define clones in 'getclones'.
#    Only keep sequence IDs that appear in at least three cells.
#    Make sure these sequence IDs don't co-occur with other sequence IDs in other cells.
#    Group cells based on these individual sequence IDs and merge them into the results produced by 'getclones'.
#    
# 4. function 'getclones_uncert':
#    In 'getclones' sequence ID pairs with sequence IDs that are ambiguous due to co-occurrence with multiple other sequence IDs are removed.
#    Here, a more lenient cut-off is applied allowing sequence IDs to appear in more distinct sequence ID pairs.
#    The output is added to the output of getclones.
#
# 5. function 'getvalidseqids':
#    Extract all cells with two distinct sequence IDs and assign them their sequence ID pair.
#    This performs in the same way as initial sequence ID selection and grouping into pairs performed by 'getclones'
#    sequence ID pairs are not filtered further, but written out into a new output table.
#    The sequence ID-pairs can then be cross-referenced with an anchor dataset, e.g. with clones called on another set of cells derived
#    from the same biological sample.



# Written by Nora Fresmann 2022-2024

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
require(plyr)
require(data.table)
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(DescTools))
suppressPackageStartupMessages(library(pheatmap))


# Wild-type seq-IDs for each target gene 
wt_seq_ids <- c("rpl39" = "TGGGCTTGTAAACCACTGTTTACTTCAGCA", 
                "actb1" = "CTTTTAGTCATTCCAGAAGCGTTTACCACT", 
                "cfl1" = "CGTATGAAGTTTAGATGGGGAGAGCGATAT",
                "rpl18a" = "TTTGCAGAAGTGTTTCGGCTGGCTTTAAGA",
                "actb2" = "TGCCAGACATTTGGTGGGGCCAACCTGTAC",
                "rpsa" = "TTTTACGTCATACCAATTCACTTAGTTTCT",
                "rps8a" = "TTAAAGAAAACACCCCAACTCATGTCATTG"
                )

# Delimiters for shortening lineage barcode reads into 30 bp seq-IDs for each target gene
seq_id_delim <- list("actb1" = c(36,65), 
                     "actb2" = c(40,69),
                     "rpl39" = c(38,67), 
                     "cfl1" = c(46,75),
                     "rpl18a" = c(44,73)
                    )



# Main function to call clones on a per-endogenous-target basis
getclones <- function(dat_wd, dat_name, wt_seq_ids, seq_id_delim, scar_file_name, heatmap_after_initial_filt = F, heatmap_after_final_filt = F, cooc_fraction_cutoff, seur_dat, tums_list, genes_list) {  
    
    # create plot output directory if missing
    if (!dir.exists(paste0(getwd(),"/pics"))) {
        dir.create(paste0(getwd(),"/pics"), recursive = TRUE)
    }

  for (tum_it in tums_list) {
    for (gene_it in genes_list) {
      
      print(paste0('Now working on tumour ',tum_it,' and gene ', gene_it, '.'))  
        
      # Step 1: Data preparation
      scar.input_full <- read.csv(paste0(dat_wd, "/", scar_file_name, "_", gene_it, ".csv"), stringsAsFactors = F)
      scar.input <- process_scar_input(scar.input_full, tum_it, gene_it)
      
      if (nrow(scar.input) == 0) next
      
      scar.input <- remove_duplicates(scar.input, gene_it, seq_id_delim)
      
      # Step 2: Selection of cells with two distinct alleles
      BCs_2scars <- scar.input$Barcode[duplicated(scar.input$Barcode)]
      scar.input.2scars <- scar.input[scar.input$Barcode %in% BCs_2scars,]
      
      # Identify sequences that appear in only one cell and iteratively remove them
      Scar_clone_overview <- table(scar.input.2scars$Barcode, scar.input.2scars$seq_id)
      seqs_one_cell <- names(which(colSums(Scar_clone_overview) == 1))
      
      for (i in 1:20) {
        scar.input.2scars <- scar.input.2scars[!scar.input.2scars$seq_id %in% seqs_one_cell, ]
        BCs_2scars <- scar.input.2scars$Barcode[duplicated(scar.input.2scars$Barcode)]
        scar.input.2scars <- scar.input[scar.input$Barcode %in% BCs_2scars, ]
        Scar_clone_overview <- table(scar.input.2scars$Barcode, scar.input.2scars$seq_id)
        seqs_one_cell <- names(which(colSums(Scar_clone_overview) == 1))
        if (length(seqs_one_cell) == 0) break
      }
      if (i == 20) next
      
      # Plot heatmap after initial filtering if desired
      if (heatmap_after_initial_filt & length(unique(scar.input.2scars$seq_id)) > 2 & length(unique(scar.input.2scars$Barcode)) > 2) {
        plot_heatmap(Scar_clone_overview, paste0(dat_wd, "/pics/scarsintwocells_afterInitialFilt_heatmap_", dat_name, "_tum_", tum_it, "_", gene_it, ".png"))
      }
      
      # Step 3: Definition of scar pairs
      scar.input.2scars <- scar.input.2scars[order(scar.input.2scars$seq_id, scar.input.2scars$Barcode), ]
      cell_scar_pairs <- get_cell_scar_pairs(scar.input.2scars)
      cell_scar_pairs_filt <- filter_cell_scar_pairs(cell_scar_pairs)
      
      if (nrow(cell_scar_pairs_filt) == 0) next
      
      # Step 4: Calculating scar co-occurrence for all scar combinations
      scar.input.2scars.filt <- scar.input.2scars[scar.input.2scars$Barcode %in% cell_scar_pairs_filt$Barcode,]
      Scar_clone_overview <- as.data.frame.matrix(table(scar.input.2scars.filt$Barcode, scar.input.2scars.filt$seq_id))
      cell_scar_pairs_filt <- calculate_scar_cooccurrence(cell_scar_pairs_filt, Scar_clone_overview)
      
      # Step 5: Removing scar-pairs with scars that are ambiguous due to co-occurrence with multiple other scars
      cell_scar_pairs_filt_filt <- cell_scar_pairs_filt[cell_scar_pairs_filt$combo_fraction_1 > cooc_fraction_cutoff | cell_scar_pairs_filt$combo_fraction_2 > cooc_fraction_cutoff,]
      
      if (nrow(cell_scar_pairs_filt_filt) == 0) next
      
      scar.input.2scars.filt.filt <- scar.input.2scars[scar.input.2scars$Barcode %in% cell_scar_pairs_filt_filt$Barcode,]
      
      # Plot heatmap after final filtering if desired
      if (heatmap_after_final_filt & length(unique(scar.input.2scars$seq_id)) > 2 & length(unique(scar.input.2scars$Barcode)) > 2) {
        plot_heatmap(Scar_clone_overview, paste0(dat_wd, "/pics/scarsintwocells_afterFinalFilt_seqcombs_over_", (cooc_fraction_cutoff * 10), "_", dat_name, "_tum_", tum_it, "_", gene_it, ".png"))
      }
      
      # Step 6: Classify each scar in a pair based on whether it was likely generated first
      cell_scar_pairs_filt_filt <- classify_scar_pairs(cell_scar_pairs_filt_filt, wt_seq_ids[gene_it])
      
      if (nrow(cell_scar_pairs_filt_filt) == 0) next
      
      # Step 7: Merge relevant data and write to file
      merge_and_write_data(cell_scar_pairs_filt_filt, dat_wd, dat_name, tum_it, gene_it, cooc_fraction_cutoff, scar.input.2scars.filt.filt)
      
      # Step 8: Add cells that only have one scar in the data, when this scar was defined as a clone-defining scar in the previous analysis steps
      scar_dat <- read.delim(paste0(dat_wd, "/final_clones_allinfo_per_cell_scar_", dat_name, "_tum_", tum_it, "_", gene_it, "_cutoff_", cooc_fraction_cutoff, ".txt"), stringsAsFactors = F, sep = ',', row.names = 1, header = T)
        
      scar.input.1scar <- add_single_scar_cells(scar.input, scar_dat)  
        
      if (nrow(scar.input.1scar) != 0) {
        print(paste0('Adding ',nrow(scar.input.1scar),' cells with a single scar to clones.'))
        scar_dat_ext <- rbind(scar_dat, scar.input.1scar)
        write.table(scar_dat_ext, paste0(dat_wd, "/final_clones_allinfo_per_cell_scar_", dat_name, "_tum_", tum_it, "_", gene_it, "_cutoff_", cooc_fraction_cutoff, "_singleAlleleExtended.txt"), quote = F, sep = ",")
      }
      
    }
  }
}


addwtclones <- function(dat_wd, dat_name, wt_seq_ids, seq_id_delim, cooc_fraction_cutoff, seur_dat, tums_list, genes_list) {

    for (tum_it in tums_list) {
        for (gene_it in genes_list) {

            print(paste0(tum_it, ",", gene_it))
            
            # Get wildtype seq_id from wt_seq_ids.
            wt_seq_id <- wt_seq_ids[gene_it] 
            
            # Load the information on clones produced so far, if there is any
            file_patterns <- c(paste0("_cutoff_", cooc_fraction_cutoff, "_singleAlleleExtended.txt"), paste0("_cutoff_", cooc_fraction_cutoff, ".txt"))
            
            scars_1 <- load_scar_data(dat_wd, dat_name, tum_it, gene_it, cooc_fraction_cutoff, file_patterns)
            if (is.null(scars_1)) {
                print(paste0("tum. ", tum_it, ", found no scar-based clones for ", gene_it))
                scars_1 <- data.frame(Barcode = c(NA, NA), Empty = c(NA, NA))
            }

            # Load and process scar input data
            scar.input_full <- read.csv(paste0(dat_wd, "/", scar_file_name, "_", gene_it, ".csv"), stringsAsFactors = F)
            scar.input_full <- process_scar_input(scar.input_full, tum_it, gene_it)

            if (nrow(scar.input_full) == 0) next

            # Remove duplicates and get sequence IDs
            scar.input <- remove_duplicates(scar.input_full, gene_it, seq_id_delim)
            scar.input$seq_id <- substr(scar.input$Sequence, seq_id_delim[[gene_it]][1], seq_id_delim[[gene_it]][2])
            wt_seq_id <- wt_seq_ids[gene_it]

            # Get cells that are most likely WT only by selecting those with a high number of UMIs for WT but no other seq-IDs
            BCs_1scar <- as.data.frame(table(scar.input$Barcode))$Var1[as.data.frame(table(scar.input$Barcode))$Freq == 1]
            scar.input.1scar <- scar.input[scar.input$Barcode %in% BCs_1scar,]
            scar.input.1scar <- scar.input.1scar[!scar.input.1scar$Barcode %in% scars_1$Barcode,]
            scar.input.1scar <- scar.input.1scar[scar.input.1scar$seq_id == wt_seq_id,]

            if (nrow(scar.input.1scar) == 0) next

            scar.input.1scar <- scar.input.1scar[scar.input.1scar$UMIs > quantile(scar.input.1scar$UMIs, 0.3),]

            if (nrow(scar.input.1scar) == 0) next

            if (nrow(scar.input.1scar) != 0) {
                # Fill in additional columns
                scar.input.1scar$comb_seq_id <- paste0(wt_seq_id,wt_seq_id)
                scar.input.1scar$seq_id_1 <- scar.input.1scar$seq_id
                scar.input.1scar$seq_id_2 <- NA
                scar.input.1scar$cloneID <- paste0(gene_it, "_WT")
                scar.input.1scar$clone_def_scar <- "WT"
                scar.input.1scar$parent_scar <- "WT"

                if (TRUE %in% names(table(scar.input.1scar$Barcode %in% scars_1$Barcode))) {
                    print("WRONG table merger 1")
                } else if (FALSE %in% table(colnames(scar.input.1scar) == colnames(scars_1))) {
                    print("WRONG table merger 2")
                }

                if (ncol(scars_1) == ncol(scar.input.1scar)) {
                    scar_dat_ext <- rbind(scars_1, scar.input.1scar)
                } else {
                    print('Previous warning about ncol inequality can be ignored. Since there were no true lineage barcode based clones, only the newly found WT cells will contribute to the output.')
                    scar_dat_ext <- scar.input.1scar
                }

                print("add WT seqs")

                write.table(scar_dat_ext, paste0(dat_wd, "/final_clones_allinfo_per_cell_scar_", dat_name, "_tum_", tum_it, "_", gene_it, "_cutoff_", cooc_fraction_cutoff, "_with_wt.txt"), quote = F, sep = ",")
            }
        }
    }
}

addonescarclones <- function(dat_wd, dat_name, wt_seq_ids, seq_id_delim, cooc_fraction_cutoff, seur_dat, tums_list, genes_list) {

  for (tum_it in tums_list) {
    for (gene_it in genes_list) {

      # Step 1: Data preparation
      scar.input_full <- read.csv(paste0(dat_wd, "/", scar_file_name, "_", gene_it, ".csv"), stringsAsFactors = F)
      scar.input <- process_scar_input(scar.input_full, tum_it, gene_it)
      if (nrow(scar.input) == 0) next
      scar.input <- remove_duplicates(scar.input, gene_it, seq_id_delim)
      scar.input$seq_id <- substr(scar.input$Sequence, seq_id_delim[[gene_it]][1], seq_id_delim[[gene_it]][2])
      wt_seq_id <- wt_seq_ids[gene_it]

      # Load the information on clones produced so far, if there is any
      file_patterns <- c(paste0("_cutoff_", cooc_fraction_cutoff, "_with_wt.txt"), paste0("_cutoff_", cooc_fraction_cutoff, "_singleAlleleExtended.txt"), paste0("_cutoff_", cooc_fraction_cutoff, ".txt"))

      scars_1 <- load_scar_data(dat_wd, dat_name, tum_it, gene_it, cooc_fraction_cutoff, file_patterns)
        
      if (is.null(scars_1)) {
        print(paste0("tum. ", tum_it, ", found no previously classified cells for gene ", gene_it))
        scars_1 <- data.frame(Barcode = c(NA, NA), Empty = c(NA, NA))
        empty_ticker <- 'yes'
      } else {
        empty_ticker <- 'no'
      }

      # Select cell barcodes that only have one scar
      BCs_1scar <- as.data.frame(table(scar.input$Barcode))$Var1[as.data.frame(table(scar.input$Barcode))$Freq == 1]
      scar.input.1scar <- scar.input[scar.input$Barcode %in% BCs_1scar,]

      if (nrow(scar.input.1scar) == 0) {
        print("No additional clone-level single scars were found")
        next
      }

      # Select cell barcodes that have a scar that has not been assigned to a double-scar-based clone ID so far.
      scar.input.1scar <- scar.input.1scar[!scar.input.1scar$seq_id %in% c(scars_1$seq_id_1, scars_1$seq_id_2),]

      if (nrow(scar.input.1scar) == 0) {
        print("No additional clone-level single scars were found")
        next
      }

      # Only keep scars that are found in at least 3 cells
      scars_three_cells <- as.data.frame(table(scar.input.1scar$seq_id))
      scars_three_cells <- unique(scars_three_cells$Var1[scars_three_cells$Freq >= 3])
      scar.input.1scar <- scar.input.1scar[scar.input.1scar$seq_id %in% scars_three_cells,]

      # Check that this scar isn't found in presence of other scars too often
      BCs_2scars <- scar.input$Barcode[duplicated(scar.input$Barcode)]
      scar.input.2scars <- scar.input[scar.input$Barcode %in% BCs_2scars,]

      for (i in 1:length(scars_three_cells)) {
        scar_i <- scars_three_cells[i]
        numb_2scarclones <- nrow(scar.input.2scars[scar.input.2scars$seq_id == scar_i,])
        if (numb_2scarclones == 0) next
        ratio <- nrow(scar.input.1scar[scar.input.1scar$seq_id == scar_i,]) / numb_2scarclones
        if (ratio < 3) {
          scar.input.1scar <- scar.input.1scar[scar.input.1scar$seq_id != scar_i,]
        } else {
          next
        }
      }

      if (nrow(scar.input.1scar) == 0) {
        print("No additional clone-level single scars were found")
        next
      } else {
        print(paste0('adding single scar clones for ', gene_it))
        scar.input.1scar$comb_seq_id <- NA
        scar.input.1scar$seq_id_1 <- scar.input.1scar$seq_id
        scar.input.1scar$seq_id_2 <- NA
        scar.input.1scar$cloneID <- paste0(gene_it, "_", scar.input.1scar$CIGAR)
        scar.input.1scar$clone_def_scar <- paste0(gene_it, "_", scar.input.1scar$CIGAR)
        scar.input.1scar$parent_scar <- NA

        if (TRUE %in% names(table(scar.input.1scar$Barcode %in% scars_1$Barcode))) {
          print("WRONG table merger 1")
        } else if (FALSE %in% table(colnames(scar.input.1scar) == colnames(scars_1))) {
          print("WRONG table merger 2")
        }

          
        if (empty_ticker == 'yes') {
          print('Previous warning about ncol inequality can be ignored. Since there were no lineage barcode based clones, only the newly found single allele cells will contribute to the output.')
          print('No two-allele based clones were previously found, but there are some single-allele clones')
          scar_dat_ext <- scar.input.1scar        
        } else if (ncol(scars_1) == ncol(scar.input.1scar)) {
          scar_dat_ext <- rbind(scars_1, scar.input.1scar)
        } else {
          print('Error - files could not be merged!')
        }

        print("Adding cells with only one scar")
        write.table(scar_dat_ext, paste0(dat_wd, "/final_clones_allinfo_per_cell_scar_", dat_name, "_tum_", tum_it, "_", gene_it, "_cutoff_", cooc_fraction_cutoff, "_SingleScarExt.txt"), quote = F, sep = ",")
      }
    }
  }
}


getclones_uncert <- function(dat_wd, dat_name, wt_seq_ids, seq_id_delim, heatmap_after_initial_filt = F, heatmap_after_final_filt = F, cooc_fraction_cutoff, cooc_fraction_cutoff_uncert, seur_dat, tums_list, genes_list) {

  for (tum_it in tums_list) {
    for (gene_it in genes_list) {
        
        print(paste0('Now working on tumour ',tum_it,' and gene ',gene_it))

        # Load the information on clones produced so far, if there is any
      file_patterns <- c(paste0("_cutoff_", cooc_fraction_cutoff, "_SingleScarExt.txt"), paste0("_cutoff_", cooc_fraction_cutoff, "_with_wt.txt"), paste0("_cutoff_", cooc_fraction_cutoff, "_singleAlleleExtended.txt"), paste0("_cutoff_", cooc_fraction_cutoff, ".txt"))

      scars_1 <- load_scar_data(dat_wd, dat_name, tum_it, gene_it, cooc_fraction_cutoff, file_patterns)
        
      if (is.null(scars_1)) {
        print(paste0("tum. ", tum_it, ", found no previously classified cells for gene ", gene_it))
        scars_1 <- data.frame(Barcode = c(NA, NA), Empty = c(NA, NA))
        empty_ticker <- 'yes'
      } else {
        empty_ticker <- 'no'
      }

        
      # Step 1: Data preparation
      scar.input_full <- read.csv(paste0(dat_wd, "/", scar_file_name, "_", gene_it, ".csv"), stringsAsFactors = F)
      scar.input <- process_scar_input(scar.input_full, tum_it, gene_it)
      if (nrow(scar.input) == 0) next
      scar.input <- remove_duplicates(scar.input, gene_it, seq_id_delim)
      wt_seq_id <- wt_seq_ids[gene_it]

      # Remove cell barcodes and seq_ids that were already assigned to a clone
      scar.input <- scar.input[!scar.input$Barcode %in% scars_1$Barcode,]
      scar.input <- scar.input[!scar.input$seq_id %in% scars_1$seq_id,]

      # First select cell barcodes that have been assigned two scars
      BCs_2scars <- scar.input$Barcode[duplicated(scar.input$Barcode)]
      scar.input.2scars <- scar.input[scar.input$Barcode %in% BCs_2scars,]

      # Remove sequences that only appear in one cell iteratively
      Scar_clone_overview <- table(scar.input.2scars$Barcode, scar.input.2scars$seq_id)
      seqs_one_cell <- colnames(t(as.matrix(colSums(Scar_clone_overview)))[, t(as.matrix(colSums(Scar_clone_overview)))[1, ] == 1])

      for (i in 1:20) {
        scar.input.2scars <- scar.input.2scars[!scar.input.2scars$seq_id %in% seqs_one_cell,]
        BCs_2scars <- scar.input.2scars$Barcode[duplicated(scar.input.2scars$Barcode)]
        scar.input.2scars <- scar.input[scar.input$Barcode %in% BCs_2scars,]
        Scar_clone_overview <- table(scar.input.2scars$Barcode, scar.input.2scars$seq_id)
        seqs_one_cell <- names(which(colSums(Scar_clone_overview) == 1))
        if (identical(seqs_one_cell, character(0))) break
      }

      if (i == 20) {
          print('No additional allele-pairs could be found.')
          next
      }
        
        
      if (heatmap_after_initial_filt & length(unique(scar.input.2scars$seq_id)) > 2 & length(unique(scar.input.2scars$Barcode)) > 2) {
        Scar_clone_overview <- as.data.frame.matrix(table(scar.input.2scars$Barcode, scar.input.2scars$seq_id))
        scar_freqs <- sort(colSums(Scar_clone_overview), decreasing = T)
        ploto <- pheatmap(Scar_clone_overview, show_rownames = FALSE, show_colnames = FALSE, breaks = seq(-1, +1, length = 101), cluster_cols = F)
        ggsave(plot = ploto, filename = paste0(dat_wd, "/pics/scarsintwocells_heatmap_", dat_name, "_tum_", tum_it, "_", gene_it, ".png"), width = 6, height = 6, dpi = 300)
        rm(ploto)
      }

      # Define clone-scar-pairs
      
      
      scar.input.2scars <- scar.input.2scars[order(scar.input.2scars$seq_id),]
      scar.input.2scars <- scar.input.2scars[order(scar.input.2scars$Barcode),]
      cell_scar_pairs <- get_cell_scar_pairs(scar.input.2scars)
      cell_scar_pairs$UMI_sum <- cell_scar_pairs$UMIs_1 + cell_scar_pairs$UMIs_2
      cell_scar_pairs$comb_seq_id <- paste0(cell_scar_pairs$seq_id_1, cell_scar_pairs$seq_id_2)
      comb_seq_ids <- as.data.frame(table(cell_scar_pairs$comb_seq_id))
      comb_seq_ids_over2 <- comb_seq_ids$Var1[comb_seq_ids$Freq > 2]
      cell_scar_pairs_filt <- cell_scar_pairs[cell_scar_pairs$comb_seq_id %in% comb_seq_ids_over2,]

      if (nrow(cell_scar_pairs_filt) == 0) next
      
        
      scar.input.2scars.filt <- scar.input.2scars[scar.input.2scars$Barcode %in% cell_scar_pairs_filt$Barcode,]
      Scar_clone_overview <- as.data.frame.matrix(table(scar.input.2scars.filt$Barcode, scar.input.2scars.filt$seq_id))

      for (i in 1:length(unique(cell_scar_pairs_filt$comb_seq_id))) {
        comb_seq_id <- unique(cell_scar_pairs_filt$comb_seq_id)[i]
        scar_1 <- unique(cell_scar_pairs_filt$seq_id_1[cell_scar_pairs_filt$comb_seq_id == comb_seq_id])
        scar_2 <- unique(cell_scar_pairs_filt$seq_id_2[cell_scar_pairs_filt$comb_seq_id == comb_seq_id])
        subs <- Scar_clone_overview[, c(scar_1, scar_2)]
        subs <- subs[rowSums(subs) != 0,]
        scar_freq_1 <- sum(subs[, 1] == 1)
        scar_freq_2 <- sum(subs[, 2] == 1)
        subs$common <- rowSums(subs)
        shared_freq <- sum(subs$common == 2)
        frac_scar_1 <- shared_freq / scar_freq_1
        frac_scar_2 <- shared_freq / scar_freq_2
        if (i == 1) {
          cell_scar_pairs_filt$combo_fraction_1 <- 0
          cell_scar_pairs_filt$combo_fraction_2 <- 0
        }
        cell_scar_pairs_filt$combo_fraction_1[cell_scar_pairs_filt$comb_seq_id == comb_seq_id] <- frac_scar_1
        cell_scar_pairs_filt$combo_fraction_2[cell_scar_pairs_filt$comb_seq_id == comb_seq_id] <- frac_scar_2
      }

      cell_scar_pairs_filt_filt <- cell_scar_pairs_filt[cell_scar_pairs_filt$combo_fraction_1 > cooc_fraction_cutoff_uncert | cell_scar_pairs_filt$combo_fraction_2 > cooc_fraction_cutoff_uncert,]
      if (nrow(cell_scar_pairs_filt_filt) == 0) next

      scar.input.2scars.filt.filt <- scar.input.2scars[scar.input.2scars$Barcode %in% cell_scar_pairs_filt_filt$Barcode,]
      if (heatmap_after_final_filt & length(unique(scar.input.2scars$seq_id)) > 2 & length(unique(scar.input.2scars$Barcode)) > 2) {
        Scar_clone_overview <- as.data.frame.matrix(table(scar.input.2scars.filt.filt$Barcode, scar.input.2scars.filt.filt$seq_id))
        scar_freqs <- sort(colSums(Scar_clone_overview), decreasing = T)
        ploto <- pheatmap(Scar_clone_overview, show_rownames = FALSE, show_colnames = FALSE, breaks = seq(-1, +1, length = 101))
        ggsave(plot = ploto, filename = paste0(dat_wd, "/pics/scarsintwocells_seqcombs_over_", (cooc_fraction_cutoff_uncert * 10), "_", dat_name, "_tum_", tum_it, "_", gene_it, ".png"), width = 6, height = 6, dpi = 300)
      }

      # Define the "child" scar that actually defines the clone
      cell_scar_pairs_filt_filt <- classify_scar_pairs(cell_scar_pairs_filt_filt, wt_seq_id)
      if (nrow(cell_scar_pairs_filt_filt) == 0) next

      clone_seqs <- as.data.frame(unique(cell_scar_pairs_filt_filt$comb_seq_id))
      colnames(clone_seqs)[1] <- "clone_seq"
      clone_seqs$cloneID <- paste0(gene_it, "_", seq_len(nrow(clone_seqs)))
      cell_scar_pairs_filt_filt <- merge(cell_scar_pairs_filt_filt, clone_seqs, by.x = "comb_seq_id", by.y = "clone_seq", all.x = T)
      write.table(cell_scar_pairs_filt_filt, paste0(dat_wd, "/final_clones_allinfo_", dat_name, "_tum_", tum_it, "_", gene_it, "_cutoff_", cooc_fraction_cutoff_uncert, ".txt"), quote = F)

      cell_scar_pairs_filt_filt_sub <- cell_scar_pairs_filt_filt[c("comb_seq_id", "Barcode", "seq_id_1", "seq_id_2", "cloneID", "clone_def_scar", "parent_scar")]
      all_info <- merge(scar.input.2scars.filt.filt, cell_scar_pairs_filt_filt_sub, by = "Barcode", all.x = T)

      if (nrow(all_info) == 0) {
        print("No two-scar clones identified, even with these looser filter settings")
        next
      } else {
        print(paste0('adding loosely filtered two-scar clones clones for ', gene_it))
        if (TRUE %in% names(table(all_info$Barcode %in% scars_1$Barcode))) {
          print("WRONG table merger 1")
        } else if (FALSE %in% table(colnames(all_info) == colnames(scars_1))) {
          print("WRONG table merger 2")
        }
          
        if (empty_ticker == 'yes') {
          print('Previous warning about ncol inequality can be ignored. Since there were no previously found clones, only the newly found WT cells will contribute to the output.')
          scar_dat_ext <- all_info
        }else if (ncol(scars_1) == ncol(all_info)) {
          scar_dat_ext <- rbind(scars_1, all_info)
        }else{
          print('Error - check the code')
        }
            
        write.table(scar_dat_ext, paste0(dat_wd, "/final_clones_allinfo_per_cell_scar_", dat_name, "_tum_", tum_it, "_", gene_it, "_withUncert.txt"), quote = F, sep = ",")
      }
    }
  }
}



# Extract all endogenous target allele combinations

getvalidseqids <- function(dat_wd, dat_name, wt_seq_ids, seq_id_delim, scar_file_name, heatmap_after_initial_filt = F, seur_dat, tums_list, genes_list){

# Loop through samples    
    for(z in 1:length(tums_list)){

        tum_it <- tums_list[z]

# Loop through genes
        for(y in 1:length(genes_list)){
    
            gene_it <- genes_list[y]        

# STEP 1: Data preparation           
            # Read in filtered scar data

            scar.input_full <- read.csv(paste0(dat_wd,"/",scar_file_name,"_",gene_it,".csv"), stringsAsFactors = F)

            scar.input_full <- merge(scar.input_full, seur_dat[,c('Barcode','sample_inf')], by = "Barcode")
            scar_freqs <- as.data.frame(table(scar.input_full$Sequence))
            scar.input_full$Presence <- 0
            scar.input_full$Presence[scar.input_full$Sequence %in% scar_freqs[scar_freqs$Freq > 1,]$Var1] <- 1
            
            print(paste0('entries in scar input object for gene ',gene_it,": ",nrow(scar.input_full)))
            
            scar.input_full$UMIs <- as.numeric(scar.input_full$UMIs)
            
            # Subset for cells from a single sample of interest, if per-sample analysis is desired
            if(tum_it == 'all'){
                scar.input <- scar.input_full
            }else{
                scar.input <- scar.input_full[scar.input_full$sample_inf == tum_it,]   
            }

            if(nrow(scar.input) == 0){
                next
            }



            # Get seq_ids, which are short stretches of sequence around the expected scar site. These will be used as identifiers for each sequence.
            # The length and positioning are defined by the input variable seq_id_delim and need to correspond to the supplied wildtype IDs in wt_seq_ids.
            scar.input$seq_id <- substr(scar.input$Sequence, seq_id_delim[[gene_it]][1],seq_id_delim[[gene_it]][2])
            
            # Get wildtype seq_id from wt_seq_ids.
            wt_seq_id <- wt_seq_ids[gene_it]                   

            # Include hard-clipped reads by giving them a different seq_id, unless they have only one UMI.
            # Then remove entries that still have an empty seq_id slot
            scar.input$seq_id[scar.input$seq_id == "" & scar.input$UMIs > 1] <- substr(scar.input$Sequence[scar.input$seq_id == "" & scar.input$UMIs > 1], 1,30)
            scar.input <- scar.input[scar.input$seq_id != "",]
            

            # Remove cell barcode seq_id duplicates
            scar.input <- scar.input[order(scar.input$UMIs, decreasing = T),]
            scar.input$BC_seq_id <- paste0(scar.input$Barcode, "_", scar.input$seq_id)
            scar.input <- scar.input[!duplicated(scar.input$BC_seq_id),]
            scar.input$BC_seq_id  <- NULL


# STEP 2: Selection of cells with two distinct alleles
            # First select cell barcodes that have been assigned two scars (wildtype sequences are allowed!)
            BCs_2scars <- scar.input$Barcode[duplicated(scar.input$Barcode) == T]
            scar.input.2scars <- scar.input[scar.input$Barcode %in% BCs_2scars,]
            
            # Remove sequences that only appear in one cell. Since some cells that previously had two different alleles, might now only have one, remove cells that have been left with only a single scar.
            # Iterate through these two steps until there are no more scars that can only be found in a single cell. If this does not resolve within 20 rounds of iteration, send a notification and continue to next gene.

            Scar_clone_overview <- table(scar.input.2scars$Barcode, scar.input.2scars$seq_id)

            seqs_one_cell <- as.data.frame(t(as.matrix(colSums(Scar_clone_overview))))
            seqs_one_cell <- colnames(seqs_one_cell[,seqs_one_cell[1,] == 1])

            for(i in 1:20){
                scar.input.2scars <- scar.input.2scars[!scar.input.2scars$seq_id %in% seqs_one_cell,]

                BCs_2scars <- scar.input.2scars$Barcode[duplicated(scar.input.2scars$Barcode) == T]
                scar.input.2scars <- scar.input[scar.input$Barcode %in% BCs_2scars,]
                nrow(scar.input.2scars) == 2*length(BCs_2scars)

                Scar_clone_overview <- table(scar.input.2scars$Barcode, scar.input.2scars$seq_id)
                seqs_one_cell <- colSums(Scar_clone_overview)
                seqs_one_cell <- names(which(seqs_one_cell == 1))

                if(identical(seqs_one_cell,character(0))  == T) {
                  break
                }
            }
            # MOVE TO NEXT gene, IF this converges to 20
            if(i == 20){
                print(paste0('Isolation of cells with two alleles that are also found in other cells did not converge. Skipping clone calling for gene ',gene_it,'.'))
                next
            }

            # Hierarchically cluster cells and make a heatmap, if desired.
            if(heatmap_after_initial_filt == T & length(unique(scar.input.2scars$seq_id)) > 2 & length(unique(scar.input.2scars$Barcode)) > 2){
                # Make a heatmap
                Scar_clone_overview <- as.data.frame.matrix(table(scar.input.2scars$Barcode, scar.input.2scars$seq_id))
                scar_freqs <- sort(colSums(Scar_clone_overview), decreasing = T)

                ploto <-  pheatmap(Scar_clone_overview,
                  show_rownames = FALSE, show_colnames = FALSE,
                  breaks = seq(-1, +1, length = 101),
                  cluster_cols = F
                )

                ggsave(plot = ploto, filename = paste0(dat_wd,"/pics/scarsintwocells_heatmap_",dat_name,"_tum_",tum_it,"_",gene_it,".png"), width = 7, height = 7, dpi = 300)
                rm(ploto)
            }


# STEP 3: Definition of scar pairs

            # Split the data into two dataframes. One contains the first allele of a cell and the second one contains the second allele of the cell.
            # We do this by sorting for seq_id first and then for the barcode to make sure that all cells with the same seq_id combination will have the two seq_ids appear in the same order in the dataframe.
            scar.input.2scars <- scar.input.2scars[order(scar.input.2scars$seq_id),]
            scar.input.2scars <- scar.input.2scars[order(scar.input.2scars$Barcode),]
            cell_scar_pairs <- scar.input.2scars[!duplicated(scar.input.2scars$Barcode),c("Barcode","seq_id","Sequence","CIGAR","UMIs","UMI_fraction")]
            colnames(cell_scar_pairs) <- c("Barcode","seq_id_1","Sequence_1","CIGAR_1","UMIs_1","UMI_fraction_1")

            cell_scar_pairs_2 <- scar.input.2scars[duplicated(scar.input.2scars$Barcode),c("Barcode","seq_id","Sequence","CIGAR","UMIs","UMI_fraction")]
            colnames(cell_scar_pairs_2) <- c("Barcode","seq_id_2","Sequence_2","CIGAR_2","UMIs_2","UMI_fraction_2")
            
            # Merge these two dataframes, so that each cell has a row with the two seq_ids in the columns seq_id_1 and seq_id_2
            cell_scar_pairs <- merge(cell_scar_pairs, cell_scar_pairs_2, by = "Barcode")
            cell_scar_pairs$UMI_sum <- cell_scar_pairs$UMIs_1 + cell_scar_pairs$UMIs_2

            # Remove combinations of seq_ids that only occur once or twice. If no seq_id combination is found more than once, send a note and move to next gene.
            cell_scar_pairs$comb_seq_id <- paste0(cell_scar_pairs$seq_id_1, cell_scar_pairs$seq_id_2) # Combine both seq_ids in a joint column
            comb_seq_ids <- as.data.frame(table(cell_scar_pairs$comb_seq_id))
            comb_seq_ids_over2 <- comb_seq_ids[comb_seq_ids$Freq > 2,]$Var1
            cell_scar_pairs_filt <- cell_scar_pairs[cell_scar_pairs$comb_seq_id %in% comb_seq_ids_over2,]

            if(nrow(cell_scar_pairs_filt) == 0){
                print(paste0("Sample ",tum_it,", gene ",gene_it,": There were no combinations that occur in more than 2 cells. Moving to the next gene"))
                next
            }

            # Now we have all observed sequence-ID-combinations that seem valid at first glance. Instead of filtering the sed-IDs that co-occur consistently and exclusively, we now write out all of them. 
            
            
            # Subset the filtered original scar data for the cells that passed the last filter
            scar.input.2scars.filt.filt <- scar.input.2scars[scar.input.2scars$Barcode %in% cell_scar_pairs_filt$Barcode,]


# STEP 4: Merge relevant data and write to file.

            # Add information on the comined seq_ids defining clones to the scar data
            clone_seqs <- as.data.frame(unique(cell_scar_pairs_filt$comb_seq_id))
            colnames(clone_seqs)[1] <- "clone_seq"

            clone_seqs$cloneID <- c(1:length(unique(clone_seqs$clone_seq)))
            clone_seqs$cloneID <- paste0(gene_it,"_",clone_seqs$cloneID)

            print(paste0('Isolated ',length(unique(clone_seqs$cloneID)), ' clone IDs for gene ',gene_it,'.'))

            
            cell_scar_pairs_filt <- merge(cell_scar_pairs_filt, clone_seqs, by.x = "comb_seq_id", by.y = "clone_seq", all.x = T)
            
            # Write extensive output to file
            write.table(cell_scar_pairs_filt, paste0(dat_wd,"/cells_with_combSeqIDInfo_",dat_name,"_tum_",tum_it,"_",gene_it,".txt"), quote = F)

            # Add metadata information and reduce individual scar-cell information to most relevant.
            cell_scar_pairs_filt_filt_sub <- cell_scar_pairs_filt[c("comb_seq_id","Barcode","seq_id_1","seq_id_2", "cloneID")]

            all_info <- merge(scar.input.2scars.filt.filt, cell_scar_pairs_filt_filt_sub, by = "Barcode", all.x = T)
            
            # Write analysis output to file
            write.table(all_info, paste0(dat_wd,"/cells_with_combSeqIDInfo_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",gene_it,".txt"), quote = F, sep = ",")
            
            
            # Add WILDTYPE CELLS

# Reload original data into variable scar.input. Subset for cells from a single sample of interest, if per-sample analysis is desired
            if(tum_it == 'all'){
                scar.input <- scar.input_full
            }else{
                scar.input <- scar.input_full[scar.input_full$sample_inf == tum_it,]   
            }

            if(nrow(scar.input) == 0){
                next
            }


            # Get seq_ids, which are short stretches of sequence around the expected scar site. These will be used as identifiers for each sequence.
            # The length and positioning are defined by the input variable seq_id_delim and need to correspond to the supplied wildtype IDs in wt_seq_ids.
            scar.input$seq_id <- substr(scar.input$Sequence, seq_id_delim[[gene_it]][1],seq_id_delim[[gene_it]][2])
            
            # Get wildtype seq_id from wt_seq_ids.
            wt_seq_id <- wt_seq_ids[gene_it]                   

            # Include hard-clipped reads by giving them a different seq_id, unless they have only one UMI.
            # Then remove entries that still have an empty seq_id slot
            scar.input$seq_id[scar.input$seq_id == "" & scar.input$UMIs > 1] <- substr(scar.input$Sequence[scar.input$seq_id == "" & scar.input$UMIs > 1], 1,30)
            scar.input <- scar.input[scar.input$seq_id != "",]
            

            # Remove cell barcode seq_id duplicates
            scar.input <- scar.input[order(scar.input$UMIs, decreasing = T),]
            scar.input$BC_seq_id <- paste0(scar.input$Barcode, "_", scar.input$seq_id)
            scar.input <- scar.input[!duplicated(scar.input$BC_seq_id),]
            scar.input$BC_seq_id  <- NULL

            # Get cells that are most likely WT only
            BCs_1scar <- as.data.frame(table(scar.input$Barcode))$Var1[as.data.frame(table(scar.input$Barcode))$Freq == 1]
            scar.input.1scar <- scar.input[scar.input$Barcode %in% BCs_1scar,]
            scar.input.1scar <- scar.input.1scar[!scar.input.1scar$Barcode %in% cell_scar_pairs_filt$Barcode,]
            scar.input.1scar <- scar.input.1scar[scar.input.1scar$seq_id == wt_seq_id,]


            if (nrow(scar.input.1scar) == 0) next

            scar.input.1scar <- scar.input.1scar[scar.input.1scar$UMIs > quantile(scar.input.1scar$UMIs, 0.3),]

            if (nrow(scar.input.1scar) == 0) next

            if (nrow(scar.input.1scar) != 0) {
                # Fill in additional columns
                scar.input.1scar$comb_seq_id <- paste0(wt_seq_id,wt_seq_id)
                scar.input.1scar$seq_id_1 <- scar.input.1scar$seq_id
                scar.input.1scar$seq_id_2 <- scar.input.1scar$seq_id
                scar.input.1scar$cloneID <- paste0(gene_it, "_WT")

                if (TRUE %in% names(table(scar.input.1scar$Barcode %in% all_info$Barcode))) {
                    print("WRONG table merger 1")
                } else if (FALSE %in% table(colnames(scar.input.1scar) == colnames(all_info))) {
                    print("WRONG table merger 2")
                }

                if (ncol(all_info) == ncol(scar.input.1scar)) {
                    scar_dat_ext <- rbind(all_info, scar.input.1scar)
                } else {
                    print('Previous warning about ncol inequality can be ignored. Since there were no scar based clones, only the newly found WT cells will contribute to the output.')
                    scar_dat_ext <- scar.input.1scar
                }

                print("add WT seqs")

                write.table(scar_dat_ext, paste0(dat_wd,"/cells_with_combSeqIDInfo_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",gene_it,"_withWT.txt"), quote = F, sep = ",")

            }
            
            rm(all_info)
            rm(scar_dat_ext)
            

        }

    }
}





# Helper functions to modularize the code

# Helper function to process scar input
process_scar_input <- function(scar.input_full, tum_it, gene_it) {
  scar.input_full <- merge(scar.input_full, seur_dat[, c('Barcode', 'sample_inf')], by = "Barcode")
  scar_freqs <- as.data.frame(table(scar.input_full$Sequence))
  scar.input_full$Presence <- ifelse(scar.input_full$Sequence %in% scar_freqs[scar_freqs$Freq > 1, ]$Var1, 1, 0)
  scar.input_full$UMIs <- as.numeric(scar.input_full$UMIs)
  if (tum_it != 'all') {
    scar.input_full <- scar.input_full[scar.input_full$sample_inf == tum_it, ]
  }
  print(paste0('entries in scar input object for gene ', gene_it, ": ", nrow(scar.input_full)))
  return(scar.input_full)
}

# Helper function to shorten the sequence to a seq_id and to remove duplicates
remove_duplicates <- function(scar.input, gene_it, seq_id_delim) {
  scar.input$seq_id <- substr(scar.input$Sequence, seq_id_delim[[gene_it]][1], seq_id_delim[[gene_it]][2])
  scar.input$seq_id[scar.input$seq_id == "" & scar.input$UMIs > 1] <- substr(scar.input$Sequence[scar.input$seq_id == "" & scar.input$UMIs > 1], 1, 30)
  scar.input <- scar.input[scar.input$seq_id != "", ]
  scar.input <- scar.input[order(scar.input$UMIs, decreasing = TRUE), ]
  scar.input$BC_seq_id <- paste0(scar.input$Barcode, "_", scar.input$seq_id)
  scar.input <- scar.input[!duplicated(scar.input$BC_seq_id), ]
  scar.input$BC_seq_id <- NULL
  return(scar.input)
}

# Function to get cell-scar pairs for a single gene
get_cell_scar_pairs <- function(scar.input.2scars) {
  cell_scar_pairs <- scar.input.2scars[!duplicated(scar.input.2scars$Barcode), c("Barcode", "seq_id", "Sequence", "CIGAR", "UMIs", "UMI_fraction")]
  colnames(cell_scar_pairs) <- c("Barcode", "seq_id_1", "Sequence_1", "CIGAR_1", "UMIs_1", "UMI_fraction_1")
  cell_scar_pairs_2 <- scar.input.2scars[duplicated(scar.input.2scars$Barcode), c("Barcode", "seq_id", "Sequence", "CIGAR", "UMIs", "UMI_fraction")]
  colnames(cell_scar_pairs_2) <- c("Barcode", "seq_id_2", "Sequence_2", "CIGAR_2", "UMIs_2", "UMI_fraction_2")
  merge(cell_scar_pairs, cell_scar_pairs_2, by = "Barcode")
}


# Function to filter cell-scar pairs based on their frequency for a single gene
filter_cell_scar_pairs <- function(cell_scar_pairs) {
  cell_scar_pairs$UMI_sum <- cell_scar_pairs$UMIs_1 + cell_scar_pairs$UMIs_2
  cell_scar_pairs$comb_seq_id <- paste0(cell_scar_pairs$seq_id_1, cell_scar_pairs$seq_id_2)
  comb_seq_ids <- as.data.frame(table(cell_scar_pairs$comb_seq_id))
  comb_seq_ids_over2 <- comb_seq_ids$Var1[comb_seq_ids$Freq > 2]
  cell_scar_pairs[cell_scar_pairs$comb_seq_id %in% comb_seq_ids_over2,]
}


# Function to calculate scar co-occurrence for a single gene
calculate_scar_cooccurrence <- function(cell_scar_pairs_filt, Scar_clone_overview) {
    
  for (i in 1:length(unique(cell_scar_pairs_filt$comb_seq_id))) {
    comb_seq_id <- unique(cell_scar_pairs_filt$comb_seq_id)[i]
    scar_1 <- unique(cell_scar_pairs_filt$seq_id_1[cell_scar_pairs_filt$comb_seq_id == comb_seq_id])
    scar_2 <- unique(cell_scar_pairs_filt$seq_id_2[cell_scar_pairs_filt$comb_seq_id == comb_seq_id])
    subs <- Scar_clone_overview[, c(scar_1, scar_2)]
    subs <- subs[rowSums(subs) != 0,]
    scar_freq_1 <- sum(subs[, 1] == 1)
    scar_freq_2 <- sum(subs[, 2] == 1)
    subs$common <- rowSums(subs)
    shared_freq <- sum(subs$common == 2)
    frac_scar_1 <- shared_freq / scar_freq_1
    frac_scar_2 <- shared_freq / scar_freq_2
    if (i == 1) {
      cell_scar_pairs_filt$combo_fraction_1 <- 0
      cell_scar_pairs_filt$combo_fraction_2 <- 0
    }
    cell_scar_pairs_filt$combo_fraction_1[cell_scar_pairs_filt$comb_seq_id == comb_seq_id] <- frac_scar_1
    cell_scar_pairs_filt$combo_fraction_2[cell_scar_pairs_filt$comb_seq_id == comb_seq_id] <- frac_scar_2
  }
  cell_scar_pairs_filt
}


# Function to classify scar pairs for a single gene
classify_scar_pairs <- function(cell_scar_pairs_filt_filt, wt_seq_id) {
  cell_scar_pairs_filt_filt$clone_def_scar <- cell_scar_pairs_filt_filt$seq_id_1
  cell_scar_pairs_filt_filt$clone_def_scar[cell_scar_pairs_filt_filt$combo_fraction_2 > cell_scar_pairs_filt_filt$combo_fraction_1] <- cell_scar_pairs_filt_filt$seq_id_2[cell_scar_pairs_filt_filt$combo_fraction_2 > cell_scar_pairs_filt_filt$combo_fraction_1]
  cell_scar_pairs_filt_filt$clone_def_scar[abs(cell_scar_pairs_filt_filt$combo_fraction_1 - cell_scar_pairs_filt_filt$combo_fraction_2) < 0.05] <- "both"
  cell_scar_pairs_filt_filt <- cell_scar_pairs_filt_filt[cell_scar_pairs_filt_filt$clone_def_scar != wt_seq_id,]
  cell_scar_pairs_filt_filt$parent_scar <- "both"
  cell_scar_pairs_filt_filt$parent_scar[cell_scar_pairs_filt_filt$clone_def_scar == cell_scar_pairs_filt_filt$seq_id_1] <- cell_scar_pairs_filt_filt$seq_id_2[cell_scar_pairs_filt_filt$clone_def_scar == cell_scar_pairs_filt_filt$seq_id_1]
  cell_scar_pairs_filt_filt$parent_scar[cell_scar_pairs_filt_filt$clone_def_scar == cell_scar_pairs_filt_filt$seq_id_2] <- cell_scar_pairs_filt_filt$seq_id_1[cell_scar_pairs_filt_filt$clone_def_scar == cell_scar_pairs_filt_filt$seq_id_2]
  cell_scar_pairs_filt_filt
}


# Function to plot heatmap
plot_heatmap <- function(Scar_clone_overview, filename) {
  ploto <- pheatmap(Scar_clone_overview, show_rownames = FALSE, show_colnames = FALSE, breaks = seq(-1, +1, length = 101))
  ggsave(plot = ploto, filename = filename, width = 7, height = 7, dpi = 300)
}

# Function to merge and write data for single targets
merge_and_write_data <- function(cell_scar_pairs_filt_filt, dat_wd, dat_name, tum_it, gene_it, cooc_fraction_cutoff, scar.input.2scars.filt.filt) {
  clone_seqs <- data.frame(clone_seq = unique(cell_scar_pairs_filt_filt$comb_seq_id))
  clone_seqs$cloneID <- paste0(gene_it, "_", seq_len(nrow(clone_seqs)))
  cell_scar_pairs_filt_filt <- merge(cell_scar_pairs_filt_filt, clone_seqs, by.x = "comb_seq_id", by.y = "clone_seq", all.x = TRUE)
  write.table(cell_scar_pairs_filt_filt, paste0(dat_wd, "/final_clones_allinfo_", dat_name, "_tum_", tum_it, "_", gene_it, "_cutoff_", cooc_fraction_cutoff, ".txt"), quote = F)
  cell_scar_pairs_filt_filt_sub <- cell_scar_pairs_filt_filt[c("comb_seq_id", "Barcode", "seq_id_1", "seq_id_2", "cloneID", "clone_def_scar", "parent_scar")]
  all_info <- merge(scar.input.2scars.filt.filt, cell_scar_pairs_filt_filt_sub, by = "Barcode", all.x = TRUE)
  write.table(all_info, paste0(dat_wd, "/final_clones_allinfo_per_cell_scar_", dat_name, "_tum_", tum_it, "_", gene_it, "_cutoff_", cooc_fraction_cutoff, ".txt"), quote = F, sep = ",")
}


# Function to add single scar cells for a single target gene
add_single_scar_cells <- function(scar.input, scar_dat) {
  BCs_1scar <- as.data.frame(table(scar.input$Barcode))$Var1[as.data.frame(table(scar.input$Barcode))$Freq == 1]
  scar.input.1scar <- scar.input[scar.input$Barcode %in% BCs_1scar,]
  scar.input.1scar <- scar.input.1scar[scar.input.1scar$seq_id %in% scar_dat$clone_def_scar,]
  scar.input.1scar <- scar.input.1scar[scar.input.1scar$UMIs > 1,]
  if (nrow(scar.input.1scar) != 0) {
    scar.input.1scar$comb_seq_id <- NA
    scar.input.1scar$seq_id_1 <- scar.input.1scar$seq_id
    scar.input.1scar$seq_id_2 <- NA
    clone_def <- scar_dat[!duplicated(scar_dat$cloneID), ]
    clone_def$seq_id <- clone_def$clone_def_scar
    clone_def <- clone_def[, c("seq_id", "cloneID", "clone_def_scar", "parent_scar")]
    scar.input.1scar <- left_join(scar.input.1scar, clone_def, by = "seq_id")
    scar.input.1scar
  } else {
    data.frame()
  }
}


# Load scar data
load_scar_data <- function(dat_wd, dat_name, tum_it, targ, cooc_fraction_cutoff, file_patterns) {
    
    for (pattern in file_patterns) {
        file_path <- paste0(dat_wd, "/final_clones_allinfo_per_cell_scar_", dat_name, "_tum_", tum_it, "_", targ, pattern)
        if (file.exists(file_path)) {
            return(read.delim(file_path, stringsAsFactors = F, sep = ','))
        }
    }
    print(paste0(tum_it, " ,target ", targ, " not found"))
    return(NULL)
}







# Main function to reload output for single-target-clones
load_single_target_clones <- function(dat_wd, dat_name, tum_it, genes_list, cooc_fraction_cutoff, use_uncert = TRUE) {
    file_list <- vector("list", length = length(genes_list))
    j <- 1
    
    for (targ in genes_list) {
        scars_i <- load_scar_data_single_target(dat_wd, dat_name, tum_it, targ, cooc_fraction_cutoff, use_uncert)
        
        if (!is.null(scars_i)) {
            file_list[[j]] <- scars_i
            print(paste0("Loaded clones for gene ", targ))
            j <- j + 1
        }
    }
    
    print(paste0("Highest scar file number: ", (j - 1)))
    return(file_list[1:(j - 1)])  # Return only the filled elements of the list
}

# Example usage
# dat_wd <- "path_to_data"
# dat_name <- "data_name"
# tum_it <- "tumor_identifier"
# genes_list <- c("gene1", "gene2", "gene3")
# cooc_fraction_cutoff <- 0.5
# use_uncert <- TRUE  # or FALSE
# loaded_clones <- load_single_target_clones(dat_wd, dat_name, tum_it, genes_list, cooc_fraction_cutoff, use_uncert)


# Helper function to load scar data for a single target
load_scar_data_single_target <- function(dat_wd, dat_name, tum_it, targ, cooc_fraction_cutoff, use_uncert) {
    file_patterns <- if (use_uncert) {
        c(
            "_withUncert.txt", 
            paste0("_cutoff_", cooc_fraction_cutoff, "_SingleScarExt.txt"), 
            paste0("_cutoff_", cooc_fraction_cutoff, "_with_wt.txt"), 
            paste0("_cutoff_", cooc_fraction_cutoff, "_singleAlleleExtended.txt"),
            paste0("_cutoff_", cooc_fraction_cutoff, ".txt")
            )
    } else {
        c(
            paste0("_cutoff_", cooc_fraction_cutoff, "_SingleScarExt.txt"), 
            paste0("_cutoff_", cooc_fraction_cutoff, "_with_wt.txt"), 
            paste0("_cutoff_", cooc_fraction_cutoff, "_singleAlleleExtended.txt"),
            paste0("_cutoff_", cooc_fraction_cutoff, ".txt")
        )
    }
    
    for (pattern in file_patterns) {
        file_path <- paste0(dat_wd, "/final_clones_allinfo_per_cell_scar_", dat_name, "_tum_", tum_it, "_", targ, pattern)
        if (file.exists(file_path)) {
            print(paste0('Loaded file pattern for target ',targ,': ',pattern))
            return(read.delim(file_path, stringsAsFactors = FALSE, sep = ','))
        }
    }
    
    print(paste0(tum_it, " ,target ", targ, " not found"))
    return(NULL)
}






# Main function to reload output for comb_seq_ids
load_single_target_comb_seq_ids <- function(dat_wd, dat_name, tum_it, genes_list) {
    file_list <- vector("list", length = length(genes_list))
    j <- 1
    
    for (targ in genes_list) {
        scars_i <- load_comb_seq_ids_single_target(dat_wd, dat_name, tum_it, targ)
        
        if (!is.null(scars_i)) {
            file_list[[j]] <- scars_i
            print(paste0("Loaded comb_seq_ids for gene ", targ))
            j <- j + 1
        }
    }
    
    print(paste0("Highest comb_seq_id file number: ", (j - 1)))
    return(file_list[1:(j - 1)])  # Return only the filled elements of the list
}

# Example usage
# dat_wd <- "path_to_data"
# dat_name <- "data_name"
# tum_it <- "tumor_identifier"
# genes_list <- c("gene1", "gene2", "gene3")


# Helper function to load scar data for a single target
load_comb_seq_ids_single_target <- function(dat_wd, dat_name, tum_it, targ) {

        file_path <- paste0(dat_wd, "/cells_with_combSeqIDInfo_", dat_name, "_tum_", tum_it, "_", targ,".txt")
        if (file.exists(file_path)) {
            
            scar_file <- read.delim(file_path, stringsAsFactors = FALSE, sep = ' ')
            scar_file$Gene <- targ
            return(scar_file)
        }

    print(paste0(tum_it, " ,target ", targ, " not found"))
    return(NULL)
}

