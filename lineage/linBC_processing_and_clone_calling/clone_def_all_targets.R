# DESCRIPTION
# Functions for grouping cells into clones using lineage barcode information from multiple target genes.
# Takes as input the cells grouped according to sequence ID-pairs on individual endogenous target genes.
# And cells that passed filtering for transgenic targets.
# All alleles (sequence ID-pairs or transgenic barcodes) are hereafter called seq-IDs.
# Function clone_def:
#    For each sample in tums_list, cluster seq-ID values by shared barcodes, merge/prune clusters by overlap rules,
#    to form clone-defining seq-ID clusters, assign each cell barcode to exactly one clone,
#    derive clone-defining seq-ID clusters.
#  Key inputs:   
#    Paths/names: dat_wd, dat_name, cell_clone_object_file_name.
#    Modes: tum_def_only (only relevant if calling clones on samples from multiple individuals at once,
#           keep scars ≥90% specific to one sample), create_intermediate_plots.
#    Thresholds (defaults values in brackets): merger_frac (default 0.9), overlap_threshold (0.8),
#               secondary_overlap_threshold (0.6), ambiguous_cutoff (0.3).
#    Tweak “resolution” of clustering: allow_large_parent_IDs/child_to_parent_fraction_cutoff
#    Create plots, mainly heatmaps of intermediate steps of cell and seq-ID grouping: create_intermediate_plots
#  Pipeline per sample:
#    Load & subset. Optionally keep only “tumor-defining” scars (≥90% in one sample).
#    Initial clustering. Build Jaccard overlap of cell barcode sets per seq-ID.
#    Build an adjacency matrix with threshold 0.3 and take connected components (igraph) -> preliminary cluster labels.
#    Merge seq-ID-clusters based on barcode overlap:
#      - Keep only clusters that overlap ≥ merger_frac (0.9) with at least one other; drop the rest from the merge graph.
#      - Ambiguity pruning: if a cluster overlaps > overlap_threshold (0.8) with ≥2 other clusters that do not overlap
#                           each other (< secondary_overlap_threshold (0.6)), mark it ambiguous and remove it.
#      - If allow_large_parent_IDs=TRUE, a size-ratio check (child_to_parent_fraction_cutoff) can rescue parent/child seq-ID
#        structures.
#      - Build a graph of remaining pairs (edges where overlap ≥ merger_frac).
#        Connected components define merged “common_cluster” (mapped to the minimum cluster ID).
#    Recheck overlaps between common_clusters, using ambiguous_cutoff (=0.3).
#    Remove clusters or choose the “child” when both directions exceed the cutoff but one is higher.
#    Remove any cellbarcodes that still appear in multiple common_clusters.
#    Remaining common_clusters are clones -> name them.
#    Merge information on cell barcode-clone assignments into original cell-lineage barcode dataframe.
#    If requested, final plots will be made: Heatmap of all cells×seq-ID with gene and clone annotations.
#                                            Heatmap of all cells×clone highlighting kept cells.




# Written by Nora Fresmann 2022-2024

suppressPackageStartupMessages(library(dichromat))
suppressPackageStartupMessages(library(ggplot2))
require(igraph)
suppressPackageStartupMessages(library(reshape2))
require(stringdist)
library(dplyr)
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
suppressPackageStartupMessages(library(reshape2))



# color palettes to use in the function
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788", "#525252", "#969696", "#252525", "#C6DBEF", "#FC9272", "#C7E9C0", "#DBDBDB", "#FDEE02", "#FFAB00", "#4D4DBF", "#BF4DBF", "#999999")


colorblindfriend24 = c( '#003D30',
                        '#005745',
                        '#00735C',
                        '#009175',
                        '#00AF8E',
                        '#00CBA7',
                        '#00EBC1',
                        '#86FFDE',
                        '#00306F',
                        '#00489E',
                        '#005FCC',
                        '#0079FA',
                        '#009FFA',
                        '#00C2F9',
                        '#00E5F8',
                        '#7CFFFA',
                        '#5F0914',
                        '#86081C',
                        '#B20725',
                        '#DE0D2E',
                        '#FF4235',
                        '#FF8735',
                        '#FFB935',
                        '#FFE239'
                        )



# function for grouping cells into clones
clone_def <- function(dat_wd, dat_name, tums_list, cell_clone_object_file_name, tum_def_only = F, create_intermediate_plots = F, merger_frac = 0.9, overlap_threshold = 0.8, secondary_overlap_threshold = 0.6, ambiguous_cutoff = 0.3, allow_large_parent_IDs = FALSE, child_to_parent_fraction_cutoff = 0.2){

    # create plot output directory if missing
    if (!dir.exists(paste0(dat_wd,"/pics"))) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    
    if(allow_large_parent_IDs == F){
        child_to_parent_fraction_cutoff <- 0
    }
    
    if(allow_large_parent_IDs == T) {
            resolution <- 'lowResolution'
    }else{
            resolution <- 'highResolution'
    }

    for(z in 1:length(tums_list)){

        tum_it <- tums_list[z]
        
        print(paste0('Working on sample ', tum_it))
        
        cell_clone_object <- read.delim(paste0(dat_wd, '/', cell_clone_object_file_name,dat_name,'_tum_',tum_it,'.csv'), stringsAsFactors = F, sep = ',', row.names = 1)

        cell_clone_object$sample_inf <- cell_clone_object$fish_all
        
        if(tum_it != 'all'){
            cell_clone_object <- cell_clone_object[cell_clone_object$sample_inf == tum_it,] 
        }
        
        
        if(nrow(cell_clone_object) == 0){
        
            print(paste0('No cells with lineage info found for sample ',tum_it,'. Moving on to the next sample'))
            next
            
        }
        
        if(tum_def_only == TRUE){
            # Get fraction of cells with a certain comb_seq_id that are part of a certain sample

            tum_scar <- as.data.frame(table(cell_clone_object$sample_inf, cell_clone_object$comb_seq_id))

            colnames(tum_scar) <- c('sample_inf','comb_seq_id','freq')

            for(i in c(1:length(unique(cell_clone_object$comb_seq_id)))){

                scar_i <- unique(cell_clone_object$comb_seq_id)[i]

                dat_i <- tum_scar[tum_scar$comb_seq_id == scar_i,]

                dat_i$frac <- dat_i$freq/(sum(dat_i$freq))


                if(i == 1){
                    tum_scar_frac <- dat_i
                }else{
                    tum_scar_frac <- rbind(tum_scar_frac, dat_i)
                }
            }
            
            # Get fish-defining scars and keep only those
            tum_scar_frac <- tum_scar_frac[tum_scar_frac$frac > 0.9,]

            cell_clone_object <- cell_clone_object[cell_clone_object$comb_seq_id %in% tum_scar_frac$comb_seq_id,]

            cell_clone_object$comb_seq_id <- as.character(cell_clone_object$comb_seq_id)
            
            # Remove cells that don't belong to an allowed fish-comb_seq_id combination
            tum_scar_frac$fish_comb_seq_id <- paste0(tum_scar_frac$sample_inf,'_',tum_scar_frac$comb_seq_id)

            
            # Only keep barcodes with allowed fish-comb_seq_id-combos
            cell_clone_object$samp_comb_seq_id <- paste0(cell_clone_object$sample_inf,'_',cell_clone_object$comb_seq_id)
            cell_clone_object <- cell_clone_object[cell_clone_object$samp_comb_seq_id %in% tum_scar_frac$fish_comb_seq_id,]

        }

# GROUP COMB_SEQ_IDS (seq-ID-pairs) THAT CO-OCCUR

        # A network-based clustering approach
        # This approach will group comb_seq_id values based on the overlap of their associated Barcode values, even when the number of rows for each comb_seq_id is imbalanced. Nodes within the same cluster belong to the same group, while nodes in different clusters represent distinct groups.

        # Assuming your data frame is called cell_clone_object
        result <- cell_clone_object %>%
          group_by(comb_seq_id) %>%
          mutate(barcode_list = list(unique(Barcode))) %>%
          ungroup() %>%
          distinct(comb_seq_id, barcode_list) %>%
          arrange(comb_seq_id)

        # Create a matrix to store the overlap fractions
        overlap_matrix <- matrix(NA, nrow = nrow(result), ncol = nrow(result))
        colnames(overlap_matrix) <- result$comb_seq_id
        rownames(overlap_matrix) <- result$comb_seq_id

        for (i in 1:(nrow(result)-1)) {
          for (j in (i+1):nrow(result)) {
            set1 <- result$barcode_list[[i]]
            set2 <- result$barcode_list[[j]]
            overlap_fraction <- length(intersect(set1, set2)) / length(union(set1, set2))
            overlap_matrix[i, j] <- overlap_fraction
            overlap_matrix[j, i] <- overlap_fraction
          }
        }

        # Create an adjacency matrix based on the overlap matrix
        threshold <- 0.3  # Adjust this threshold as needed

        adj_matrix <- as.matrix(overlap_matrix >= threshold)

        # Create an igraph graph from the adjacency matrix
        graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

        # Find connected components in the graph
        clusters <- components(graph)

        # Extract the cluster membership for each comb_seq_id
        cluster_df <- data.frame(comb_seq_id = result$comb_seq_id, cluster = clusters$membership)


         # Merge the original dataframe with the cluster assignments based on comb_seq_id
        merged_df <- left_join(cell_clone_object, cluster_df, by = "comb_seq_id")

        
# Some plots to check that the process works correctly
                
        if(create_intermediate_plots == T &
           length(unique(cell_clone_object$Barcode)) >= 3 &
           length(unique(cell_clone_object$comb_seq_id)) >= 2){

            # Plot heatmap of clones vs. comb_seq_IDs with comb_seq_IDs colored by cluster:

            # Get matrix of cell barcode vs seq_id with integration id
            allls_cbs <- as.data.frame.matrix(table(cell_clone_object$Barcode, cell_clone_object$comb_seq_id))
            table(rowSums(allls_cbs))
            allls_cbs$Barcode <- rownames(allls_cbs)
            allls_cbs <- allls_cbs[!duplicated(allls_cbs$Barcode),]
            dim(allls_cbs)
            table(duplicated(allls_cbs$Barcode))
            allls_cbs$Barcode <- NULL

            # Create a dataframe for the comb_seq_id clusters
            scar_info_clust <- cluster_df
            scar_info_clust$cluster <- as.character(scar_info_clust$cluster)
            scar_info_clust$comb_seq_id <- NULL

            table(unique(colnames(allls_cbs)) == rownames(scar_info_clust))

            cols_clus <- colorRampPalette(tol21rainbow)(length(unique(scar_info_clust$cluster)))
            names(cols_clus) <- unique(scar_info_clust$cluster)

            ann_colors = list(
                cluster = cols_clus
            )

            ploto <- pheatmap(allls_cbs,
              show_rownames = FALSE, show_colnames = FALSE,
              breaks = seq(-1, +1, length = 101),
              annotation_col = scar_info_clust,
              annotation_colors = ann_colors,
              cluster_cols = T,
              main = 'combseqIDs_colored_by_cluster'

            )

            ggsave(plot = ploto, filename = paste0(dat_wd,"/pics/cloneDefPlots_",dat_name,"_",tum_it,"_Heatmap_allCells_initialCombSeqIDclustering.png"), width = 8 , height = 10, dpi = 500)
            rm(ploto)

        }

        
        if(create_intermediate_plots == T){
            
            # Heatmap to show overlap between barcodes. This is the overlap matrix that is being fed into the graph-based clustering

            # Convert overlap matrix to a format suitable for ggplot
            melted_overlap <- melt(overlap_matrix, na.rm = TRUE)

            # Plot heatmap
            ploto <- ggplot(melted_overlap, aes(Var1, Var2, fill = value)) +
                          geom_tile() +
                          scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                               midpoint = 0.5, limit = c(0, 1), space = "Lab", 
                                               name = "Overlap") +
                          theme(axis.text.x=element_blank(),
                                axis.ticks.x=element_blank(),
                                axis.text.y=element_blank(),
                                axis.ticks.y=element_blank()) +
                          coord_fixed() +
                          labs(title = "Comb_seq_ID-Barcode Overlap Heatmap", x = "comb_seq_id", y = "comb_seq_id")

            ggsave(plot = ploto, filename = paste0(dat_wd,"/pics/cloneDefPlots_",dat_name,"_",tum_it,"_combSeqID_overlap_matrix_forGraphBasedClustering.png"), width = 8 , height = 6.5, dpi = 500)
            rm(ploto)
            

            
            # Bar Plots of Cluster Size (Cluster Membership) and number of comb_seq_ids per cloneID

            # Get cluster size
            cluster_sizes <- merged_df %>%
              group_by(cluster) %>%
              summarise(cell_count = n())

            # Plot bar plot
            ploto <- ggplot(cluster_sizes, aes(x = as.factor(cluster), y = cell_count)) +
                  geom_bar(stat = "identity", fill = "steelblue") +
                  labs(title = "Number of Cells per Cluster", x = "Cluster ID", y = "Cell Count") +
                  theme_bw() +
                  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


            ggsave(plot = ploto, filename = paste0(dat_wd,"/pics/cloneDefPlots_",dat_name,"_",tum_it,"_CellCountPer_CombSeqIDCluster.png"), width =9 , height = 5, dpi = 500)
            rm(ploto)

            # Count the number of unique comb_seq_id per cluster
            comb_seq_count_per_cluster <- merged_df %>%
              group_by(cluster) %>%
              summarise(comb_seq_id_count = n_distinct(comb_seq_id))

            # Plot bar plot for comb_seq_id count per cluster
            ploto <- ggplot(comb_seq_count_per_cluster, aes(x = as.factor(cluster), y = comb_seq_id_count)) +
                  geom_bar(stat = "identity", fill = "darkgreen") +
                  labs(title = "Number of comb_seq_ids per Cluster", 
                       x = "Cluster ID", 
                       y = "comb_seq_id Count") +
                  theme_bw() +
                  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

            ggsave(plot = ploto, filename = paste0(dat_wd,"/pics/cloneDefPlots_",dat_name,"_",tum_it,"_CombSeqIDCountPer_CombSeqIDCluster.png"), width = 9 , height = 5, dpi = 500)
            rm(ploto)
        }
        
        
# MERGE CLUSTERS BASED ON BARCODE OVERLAP

        # Create a matrix to store the overlap fractions
        n_clusters <- max(merged_df$cluster, na.rm = TRUE)
        overlap_matrix <- matrix(NA, nrow = n_clusters, ncol = n_clusters)
        colnames(overlap_matrix) <- 1:n_clusters
        rownames(overlap_matrix) <- 1:n_clusters

        # Calculate overlap fractions
        for (cluster_i in 1:n_clusters) {
          for (cluster_j in 1:n_clusters) {
            # Filter the data for comb_seq_id in cluster_i
            cluster_i_barcodes <- merged_df %>%
              filter(cluster == cluster_i) %>%
              pull(Barcode)

            # Filter the data for comb_seq_id in cluster_j
            cluster_j_barcodes <- merged_df %>%
              filter(cluster == cluster_j) %>%
              pull(Barcode)

            # Calculate the overlap fraction for both of the pair separately and keep the maximum
            overlap_fraction_1 <- length(intersect(cluster_i_barcodes, cluster_j_barcodes)) / length(cluster_i_barcodes)
            overlap_fraction_2 <- length(intersect(cluster_i_barcodes, cluster_j_barcodes)) / length(cluster_j_barcodes)

            overlap_fraction <- max(c(overlap_fraction_1, overlap_fraction_2))

            overlap_matrix[cluster_i, cluster_j] <- overlap_fraction
          }
        }
        
        # Set diagonal values to NA
        diag(overlap_matrix) <- 0

        # Find row and column indices with at least one overlap of > merger_frac
        rows_to_keep <- apply(overlap_matrix, 1, function(row) any(row >= merger_frac))
        cols_to_keep <- apply(overlap_matrix, 2, function(col) any(col >= merger_frac))

        # Subset the overlap matrix to keep only the relevant rows and columns
        filtered_overlap_matrix <- overlap_matrix[rows_to_keep, cols_to_keep]


         if(nrow(filtered_overlap_matrix) != 0){
         # Identify ambiguous clusters that have significant overlap with multiple other clusters.
            # If these other clusters are independent of eachother, the overarching cluster is considered to be ambiguous and is removed.

#            overlap_threshold <- 0.8
#            secondary_overlap_threshold <- 0.6

            # List to store the clusters with significant overlaps
            clusters_with_high_overlap <- list()

            for (i in 1:nrow(filtered_overlap_matrix)) {
              overlapping_clusters <- which(filtered_overlap_matrix[i, ] > overlap_threshold)

              # Only consider clusters with overlap > 0.8 with at least 2 other clusters
              if (length(overlapping_clusters) > 1) {
                clusters_with_high_overlap[[rownames(filtered_overlap_matrix)[i]]] <- overlapping_clusters
              }
            }

            # Step 2: Identify ambiguous clusters
            ambiguous_clusters <- c()

            for (cluster in names(clusters_with_high_overlap)) {
              overlapping_clusters <- clusters_with_high_overlap[[cluster]]

              # Check pairwise overlaps between the clusters that overlap with the current cluster
              is_ambiguous <- FALSE

              for (i in 1:(length(overlapping_clusters) - 1)) {
                for (j in (i + 1):length(overlapping_clusters)) {
                  cluster1 <- overlapping_clusters[i]
                  cluster2 <- overlapping_clusters[j]

                  # If the overlap between two of the overlapping clusters is less than 0.6, the current cluster is ambiguous
                  if (filtered_overlap_matrix[cluster1, cluster2] < secondary_overlap_threshold) {
                    is_ambiguous <- TRUE
                    break
                  }
                }
                if (is_ambiguous) break
              }

              if(is_ambiguous == TRUE & allow_large_parent_IDs == TRUE){
                  
                  ratio_child_parent <- nrow(merged_df[merged_df$cluster %in% as.numeric(names(overlapping_clusters)),]) / nrow(merged_df[merged_df$cluster %in% as.numeric(cluster),])
                  
                  if(ratio_child_parent < child_to_parent_fraction_cutoff) {
                      filtered_overlap_matrix[cluster1, cluster2] <- 1
                      filtered_overlap_matrix[cluster2, cluster1] <- 1
                  }else{
                      ambiguous_clusters <- unique(c(ambiguous_clusters, cluster))
                  }
              }else if (is_ambiguous) {
                  ambiguous_clusters <- unique(c(ambiguous_clusters, cluster))
              }
            }

            # Remove the ambiguous clusters from the matrix and from the data in general!
            filtered_overlap_matrix <- filtered_overlap_matrix[!rownames(filtered_overlap_matrix) %in% ambiguous_clusters, !colnames(filtered_overlap_matrix) %in% ambiguous_clusters]

            ambiguous_comb_seq_ids <- cluster_df$comb_seq_id[cluster_df$cluster %in% ambiguous_clusters]
            cell_clone_object_filt <- cell_clone_object %>%
                      filter(!(comb_seq_id %in% ambiguous_comb_seq_ids)) # Remove from original dataframe


            if(!is.null(nrow(filtered_overlap_matrix))){
                # Now remove clusters that no longer have significant overlap with any other clusters from the matrix
                # Step 1: Identify rows and columns that have values greater than overlap_threshold
                rows_to_keep <- apply(filtered_overlap_matrix, 1, function(row) any(row > overlap_threshold))
                cols_to_keep <- apply(filtered_overlap_matrix, 2, function(col) any(col > overlap_threshold))

                # Step 2: Filter the matrix to include only rows and columns that have at least one value > 0.8
                filtered_overlap_matrix <- filtered_overlap_matrix[rows_to_keep, cols_to_keep]
            }  
      
         }else{
         
             cell_clone_object_filt <- cell_clone_object 
             
         }      
                              
                              
        # There may be nothing to merge. In that case move on. Otherwise, merge cloneIDs that are in the same cluster.
                              
        if(is.null(nrow(filtered_overlap_matrix))){
            # Fill the common_cluster column with the old cluster numbers
            cluster_df$common_cluster <- cluster_df$cluster
        }else if(nrow(filtered_overlap_matrix) == 0){
            # Fill the common_cluster column with the old cluster numbers
            cluster_df$common_cluster <- cluster_df$cluster
        }else{
            # First fill the common_cluster column with the old cluster numbers
            cluster_df$common_cluster <- cluster_df$cluster 

            # Find the pairs of clusters with a value over merger_frac in the filtered_overlap_matrix
            pairs_to_merge_indices <- which(filtered_overlap_matrix >= merger_frac, arr.ind = TRUE)
            pairs_to_merge <- rownames(filtered_overlap_matrix)[pairs_to_merge_indices[, "row"]]
            pairs_to_merge_clusters <- colnames(filtered_overlap_matrix)[pairs_to_merge_indices[, "col"]]

            # Combine pairs into a single edge list
            edges <- data.frame(cluster1 = as.numeric(pairs_to_merge), 
                                cluster2 = as.numeric(pairs_to_merge_clusters))


            # Use igraph to find connected components
            graph <- graph_from_data_frame(edges, directed = FALSE)
            components <- clusters(graph)

            # Map each cluster to the lowest cluster in its component
            component_clusters <- split(as.numeric(names(components$membership)), components$membership)

            # Create a mapping of each cluster to the lowest cluster in its component
            cluster_mapping <- unlist(lapply(component_clusters, function(x) {
                rep(min(x), length(x))
            }))

            # Assign names to the mapping for direct lookup
            names(cluster_mapping) <- unlist(component_clusters)

            # Update only the clusters that are part of the merging
            to_merge <- as.numeric(names(cluster_mapping))
            cluster_df$common_cluster[cluster_df$cluster %in% to_merge] <- 
            cluster_mapping[as.character(cluster_df$cluster[cluster_df$cluster %in% to_merge])]
        
        }

        cell_clone_object_ext_full <- left_join(cell_clone_object_filt, cluster_df, by = 'comb_seq_id')
        cell_clone_object_ext_full$comm_clust_BC <- paste0(cell_clone_object_ext_full$Barcode,'_',cell_clone_object_ext_full$common_cluster)

        cell_clone_object_ext <- cell_clone_object_ext_full[!duplicated(cell_clone_object_ext_full$comm_clust_BC),]         
                                  
                                  
# Remove common clusters and specific barcodes to make sure that each barcode is assigned to a single clone:

        #    We calculate the overlap fractions for all pairs of common_cluster values by iterating through each pair, filtering the data for the two clusters, and calculating the overlap in Barcode values.

        #    We store the overlap fractions in the overlap_data data frame.

        #    We determine which common_cluster values to remove based on the specified criteria: if either of the overlap fractions is above 0.25, the common_cluster is marked for removal.

        #    We determine which Barcode values to remove based on both overlap fractions being <= 0.25 and belonging to clusters marked for removal.

        #    Finally, we remove rows corresponding to the common clusters and barcodes to remove, resulting in the filtered_cluster_df.


        # Calculate the overlap fractions for all pairs of common_cluster values
        common_clusters <- unique(cell_clone_object_ext$common_cluster)
        overlap_data <- data.frame()

        clust_combs <- expand.grid(unique(cell_clone_object_ext$common_cluster),unique(cell_clone_object_ext$common_cluster))
        clust_combs <- clust_combs[clust_combs$Var1 != clust_combs$Var2,]

        for (i in 1:nrow(clust_combs)) {

            cluster1 <- clust_combs$Var1[i]
            cluster2 <- clust_combs$Var2[i]

            # Filter the data for the two common_cluster values
            data1 <- cell_clone_object_ext %>% filter(common_cluster == cluster1)
            data2 <- cell_clone_object_ext %>% filter(common_cluster == cluster2)

            # Calculate the overlap in Barcode values
            barcode_overlap <- intersect(data1$Barcode, data2$Barcode)
            total_barcode_count1 <- length(unique(data1$Barcode))
            total_barcode_count2 <- length(unique(data2$Barcode))

            # Calculate overlap fractions
            overlap_fraction1 <- length(barcode_overlap) / total_barcode_count1
            overlap_fraction2 <- length(barcode_overlap) / total_barcode_count2

            # Store the results
            overlap_data <- rbind(overlap_data, data.frame(
              common_cluster1 = cluster1,
              common_cluster2 = cluster2,
              overlap_fraction1 = overlap_fraction1,
              overlap_fraction2 = overlap_fraction2
            ))
        }

        overlap_data <- overlap_data[!overlap_data$common_cluster1 == overlap_data$common_cluster2,]
        overlap_data <- overlap_data[!is.na(overlap_data$common_cluster1) & !is.na(overlap_data$common_cluster2),]

        # Determine which common_cluster values to remove based on the criteria
        clusters_to_remove <- c()

        for (i in 1:length(common_clusters)) {
          cluster1 <- common_clusters[i]

          # Check if either of the overlap fractions is above 0.3.
          # If so, check whether the other overlap-fraction is also above 0.3.
          # If so keep the cluster with a higher overlap fraction (i.e. the one that is likely the 'child' cluster)
          # If not, remove the cluster with overlap-fraction above 0.3

        if (any(overlap_data$common_cluster1 == cluster1 & overlap_data$overlap_fraction1 > ambiguous_cutoff)) {
           
            overlap_data_sub <- overlap_data[overlap_data$common_cluster1 == cluster1 & overlap_data$overlap_fraction1 > ambiguous_cutoff,]
            
            for(row in 1:nrow(overlap_data_sub)){
            
                if(overlap_data_sub[row,'overlap_fraction2'] > ambiguous_cutoff){
                
                    if(overlap_data_sub$overlap_fraction1[row] > overlap_data_sub$overlap_fraction2[row]){
                        clusters_to_remove <- c(clusters_to_remove, cluster2)
                    }else if(overlap_data_sub$overlap_fraction1[row] < overlap_data_sub$overlap_fraction2[row]){
                        clusters_to_remove <- c(clusters_to_remove, cluster1)
                    }else if(overlap_data_sub$overlap_fraction1[row] == overlap_data_sub$overlap_fraction2[row]){
                        clusters_to_remove <- c(clusters_to_remove, cluster1, cluster2)
                    }
                
                }else{
                
                    clusters_to_remove <- c(clusters_to_remove, cluster1)
                
                }   
            }
          }
        }                            


        # Remove rows corresponding to the common_clusters
        filtered_cell_clone_object_ext <- cell_clone_object_ext %>%
          filter(!(common_cluster %in% clusters_to_remove))

        # Remove barcodes that are still duplicated, i.e. associated with multiple clusters entirely
        bcs_to_rem <- as.data.frame(table(filtered_cell_clone_object_ext$Barcode, filtered_cell_clone_object_ext$common_cluster))
        bcs_to_rem <- bcs_to_rem[bcs_to_rem$Freq > 0,]
        bcs_to_rem <- bcs_to_rem$Var1[duplicated(bcs_to_rem$Var1)]

        filtered_cell_clone_object_ext <- filtered_cell_clone_object_ext[!filtered_cell_clone_object_ext$Barcode %in% bcs_to_rem,]

                              
        write.csv(filtered_cell_clone_object_ext, paste0(dat_wd,'/clone_defs_',dat_name,'_',tum_it,'.csv'), quote = F)
        
                              
        ## Name clones per sample and write to file
            
        cell_clone_object_ext_final <- cell_clone_object_ext_full[cell_clone_object_ext_full$common_cluster %in% filtered_cell_clone_object_ext$common_cluster,]
        cell_clone_object_ext_final <- cell_clone_object_ext_final[cell_clone_object_ext_final$Barcode %in% cell_clone_object_ext_final$Barcode,]  
 
        cluster_seq_ids <- as.data.frame(table(cell_clone_object_ext_final$common_cluster, cell_clone_object_ext_final$comb_seq_id))
        cluster_seq_ids <- cluster_seq_ids[cluster_seq_ids$Freq != 0,]
        cluster_seq_ids <- cluster_seq_ids[order(cluster_seq_ids$Var1),]
        cluster_seq_ids$Freq <- NULL

        colnames(cluster_seq_ids) <- c('common_cluster','comb_seq_id')

        cluster_seq_ids$common_cluster <- as.numeric(as.character(cluster_seq_ids$common_cluster))
        cluster_seq_ids$comb_seq_id <- as.character(cluster_seq_ids$comb_seq_id)

        #filtered_cell_clone_object_ext$sample_inf <- 'all'
                              
        fish_clust_ids <- as.data.frame(table(cell_clone_object_ext_final$common_cluster, cell_clone_object_ext_final$sample_inf))
        fish_clust_ids <- fish_clust_ids[fish_clust_ids$Freq != 0,]
        fish_clust_ids <- fish_clust_ids[order(fish_clust_ids$Var1),]
        fish_clust_ids$Freq <- NULL

        colnames(fish_clust_ids) <- c('common_cluster','sample_inf')
                              
        fish_clust_ids$common_cluster <- as.numeric(as.character(fish_clust_ids$common_cluster))
        fish_clust_ids$sample_inf <- as.character(fish_clust_ids$sample_inf)
        fish_clust_ids$fish_clone <- paste0(fish_clust_ids$sample_inf, '_',fish_clust_ids$common_cluster)

        if(tum_def_only == FALSE){
            
            fish_clust_ids <- fish_clust_ids[!duplicated(fish_clust_ids$common_cluster),]
            
        }
                              
        clone_defs <- left_join(cluster_seq_ids, fish_clust_ids, by = 'common_cluster')
        clone_defs$common_cluster <- as.numeric(clone_defs$common_cluster)

            
            
            
   # Merge with original lineage barcode table and save this
        filtered_cell_clone_object_ext_full <- left_join(cell_clone_object, clone_defs[,c('comb_seq_id','common_cluster','fish_clone')], by = 'comb_seq_id')
        
        write.csv(filtered_cell_clone_object_ext_full, paste0(dat_wd,'/original_scar_data_with_clone_defs_final_',dat_name,'_',tum_it,'_',resolution,'Clones_resCutoff_',child_to_parent_fraction_cutoff,'.csv'), quote = F)



    # Make sure each cell was only assigned to one clone
        # Get observed combinations of barcodes and cloneIDs (not taking into account NA-clone entries right now!)
        BC_clone_assign <- as.data.frame(table(filtered_cell_clone_object_ext_full$Barcode[!is.na(filtered_cell_clone_object_ext_full$fish_clone)], filtered_cell_clone_object_ext_full$fish_clone[!is.na(filtered_cell_clone_object_ext_full$fish_clone)]))
        BC_clone_assign <- BC_clone_assign[BC_clone_assign$Freq > 0,]
        BC_clone_assign_ques <- BC_clone_assign

        table(duplicated(BC_clone_assign_ques$Var1)) # if there are TRUE incidents, there are cells that have been assigned to multiple clones

        # Keep only the duplicated ones in the filter table
        ques <- BC_clone_assign_ques$Var1[duplicated(BC_clone_assign_ques$Var1)]
        BC_clone_assign_ques <- BC_clone_assign_ques[BC_clone_assign_ques$Var1 %in% ques,]
        #BC_clone_assign$bc_sample <- paste0(BC_clone_assign$Var1,'_',pt_per_cell$Var2)

        # Rescue the barcodes that have been assigned to a tumor with multiple lines of evidence, i.e. that have many code_combinations pointing to one tumor and few to another.
        # Only keep the truely ambiguous barcode-tumor combinations in the table.
        BC_clone_assign_ques <- BC_clone_assign_ques %>%
          group_by(Var1) %>%
          filter(Freq <= (sum(Freq)+1) - Freq) %>%
          ungroup()

        # Remove the questionable barcodes from the table entirely
        dim(filtered_cell_clone_object_ext_full)
        filtered_cell_clone_object_ext_full <- filtered_cell_clone_object_ext_full[!filtered_cell_clone_object_ext_full$Barcode %in% BC_clone_assign_ques$Var1,]
        dim(filtered_cell_clone_object_ext_full)



        # Merge the inferred clone assignments back to the full dataframe and save this
        BC_clone_assign$Freq <- NULL
        colnames(BC_clone_assign) <- c('Barcode','inferred_clone')

        filtered_cell_clone_object_inf <- left_join(filtered_cell_clone_object_ext_full, BC_clone_assign, by = 'Barcode')

        write.csv(filtered_cell_clone_object_inf, paste0(dat_wd,'/filtered_scar_data_with_clone_defs_final_',dat_name,'_',tum_it,'_',resolution,'Clones_resCutoff_',child_to_parent_fraction_cutoff,'.csv'), quote = F)


    # Identify and extract all clone-defining combinations of comb_seq_ids to use those later on for allograft cell assignment

        # Extract barcode comb_seq_id combos
        bc_seqs <- as.data.frame(table(filtered_cell_clone_object_inf$comb_seq_id, filtered_cell_clone_object_inf$Barcode))
        bc_seqs <- bc_seqs[bc_seqs$Freq > 0,]    
        bc_seqs$Freq <- NULL
        colnames(bc_seqs) <- c('comb_seq_id', 'Barcode')

        table(is.na(bc_seqs$comb_seq_id)) # should all be FALSE

        # Add all combinations of comb_seq_IDs that are seen in the clones

        result <- bc_seqs %>%
          group_by(Barcode) %>%
          summarise(code_combinations = sapply(strsplit(paste(sort(unique(comb_seq_id)), collapse = ","), ","), function(x) {
            unlist(sapply(1:length(x), function(size) combn(x, size, simplify = TRUE, FUN = function(y) paste(y, collapse = ""))))
          })) %>%
          unnest(cols = code_combinations)



        # Add information on the associated clone to the combinations of comb_seq_ids
        filtered_cell_clone_object_inf_min <- filtered_cell_clone_object_inf[!duplicated(filtered_cell_clone_object_inf$Barcode),]
        filtered_cell_clone_object_inf_min <- filtered_cell_clone_object_inf_min[,c('Barcode','fish_clone','inferred_clone')]

        comb_seq_id_combos <- left_join(result, filtered_cell_clone_object_inf_min, by = 'Barcode')

        # Of the sequence_ID combinations in cells that could not be assigned to a clone, only keep the ones that are seen in 10 or more cells

        clone_cells <- comb_seq_id_combos[!is.na(comb_seq_id_combos$inferred_clone),]
        none_clone_cells <- comb_seq_id_combos[is.na(comb_seq_id_combos$inferred_clone),]

        combo_freqs <- as.data.frame(table(none_clone_cells$code_combinations))
        combo_freqs <- combo_freqs$Var1[combo_freqs$Freq >= 10]
        none_clone_cells <- none_clone_cells[none_clone_cells$code_combinations %in% combo_freqs,]

        comb_seq_id_combos <- rbind(clone_cells, none_clone_cells)

        # Get an inferred-clone column without NAs. Replace by 'none'. 'None' will be considered as the background sequence-ID information.
        comb_seq_id_combos$inferred_clone_NAless <- comb_seq_id_combos$inferred_clone
        comb_seq_id_combos$inferred_clone_NAless <- as.character(comb_seq_id_combos$inferred_clone_NAless)
        comb_seq_id_combos$inferred_clone_NAless[is.na(comb_seq_id_combos$inferred_clone_NAless)] <- 'none'



        ## Get fractions of comb_seq_ID combinations per clone

        tum_scar <- as.data.frame(table(comb_seq_id_combos$inferred_clone_NAless, comb_seq_id_combos$code_combinations))

        colnames(tum_scar) <- c('inferred_clone_NAless','code_combinations','freq')

        for(i in c(1:length(unique(comb_seq_id_combos$code_combinations)))){

            scar_i <- unique(comb_seq_id_combos$code_combinations)[i]

            dat_i <- tum_scar[tum_scar$code_combinations == scar_i,]


            dat_i$frac <- dat_i$freq/(sum(dat_i$freq))


            if(i == 1){

                tum_scar_frac <- dat_i

            }else{

                tum_scar_frac <- rbind(tum_scar_frac, dat_i)

            }


        }

        # Only keep sequence combinations that can be assigned to a single cloneID with high certainty
        # And only keep cells that have the combination and were assigned to the correct fish!
        tum_scar_frac <- tum_scar_frac[tum_scar_frac$frac > 0.9,]
        tum_scar_frac$inferred_clone_code_combination <- paste0(tum_scar_frac$inferred_clone_NAless,'_', tum_scar_frac$code_combinations)
        comb_seq_id_combos$inferred_clone_code_combination <- paste0(comb_seq_id_combos$inferred_clone_NAless,'_', comb_seq_id_combos$code_combinations)

        # Label comb_seq_id combinations that are clone defining 
        comb_seq_id_combos$clone_defining <- FALSE
        comb_seq_id_combos$clone_defining[comb_seq_id_combos$inferred_clone_code_combination %in% tum_scar_frac$inferred_clone_code_combination] <- TRUE
        comb_seq_id_combos$clone_defining[comb_seq_id_combos$inferred_clone_NAless == 'none'] <- FALSE


        write.csv(comb_seq_id_combos, paste0(dat_wd,'/all_combseqid_combinations_with_cloneAssignments_',dat_name,'_',tum_it,'_',resolution,'Clones_resCutoff_',child_to_parent_fraction_cutoff,'.csv'), quote = F)



# Plot final clone assignments for all cells
            # First plot the comb_seq_ids and their assignment to the common clusters
            filtered_cell_clone_object_ext_full <- left_join(cell_clone_object, clone_defs[,c('comb_seq_id','common_cluster','fish_clone')], by = 'comb_seq_id')

            # Plot heatmap of clones vs. comb_seq_IDs with comb_seq_IDs colored by cluster and highlighting the cells that made it through filtering:

            # Get matrix of cell barcode vs seq_id with integration id
            allls_cbs <- as.data.frame.matrix(table(filtered_cell_clone_object_ext_full$Barcode, filtered_cell_clone_object_ext_full$comb_seq_id))
            table(rowSums(allls_cbs))
            allls_cbs$Barcode <- rownames(allls_cbs)
            allls_cbs <- allls_cbs[!duplicated(allls_cbs$Barcode),]
            dim(allls_cbs)
            table(duplicated(allls_cbs$Barcode))
            allls_cbs$Barcode <- NULL


            # Create a dataframe for the comb_seq_id clusters

            comb_seq_gene_info <- filtered_cell_clone_object_inf[,c('comb_seq_id','Gene')]
            comb_seq_gene_info <- comb_seq_gene_info[!duplicated(comb_seq_gene_info$comb_seq_id),]
            clone_defs_gene <- left_join(clone_defs,comb_seq_gene_info, by = 'comb_seq_id')
            
            scar_info_clust <- clone_defs_gene
            rownames(scar_info_clust) <- scar_info_clust$comb_seq_id
            scar_info_clust$fish_clone <- as.character(scar_info_clust$fish_clone)
            scar_info_clust$sample_inf <- NULL
            scar_info_clust$comb_seq_id <- NULL
            scar_info_clust$common_cluster <- NULL

            filteredout_ids <- filtered_cell_clone_object_ext_full[!duplicated(filtered_cell_clone_object_ext_full$comb_seq_id),c('comb_seq_id','Gene')]
            filteredout_ids$fish_clone <- 'none'
            filteredout_ids <- filteredout_ids[!filteredout_ids$comb_seq_id %in% rownames(scar_info_clust),]
            rownames(filteredout_ids) <- filteredout_ids$comb_seq_id
            filteredout_ids$comb_seq_id <- NULL

            scar_info_clust_all <- rbind(scar_info_clust, filteredout_ids)
            scar_info_clust_all$fish_clone <- as.factor(scar_info_clust_all$fish_clone)


            scar_info_clust_all <- rbind(scar_info_clust, filteredout_ids)
            scar_info_clust_all$fish_clone <- as.factor(scar_info_clust_all$fish_clone)  

            gene_pal <- colorblindfriend24[c(1,3,8,10,13,15,18,21,24)]

            # Create a dataframe for the barcodes
            cell_info <- as.data.frame(unique(rownames(allls_cbs)))
            colnames(cell_info) <- 'Barcode'
            rownames(cell_info) <- cell_info$Barcode
            cell_info$passed_filter <- 'no'
            cell_info$passed_filter[cell_info$Barcode %in% filtered_cell_clone_object_inf$Barcode[!is.na(filtered_cell_clone_object_inf$inferred_clone)]] <- 'yes'
            cell_info$Barcode <- NULL


            table(unique(colnames(allls_cbs)) %in% rownames(scar_info_clust_all))
            table(unique(rownames(allls_cbs)) %in% rownames(cell_info))
            

            if(length(unique(scar_info_clust_all$fish_clone))-1 > length(tol21rainbow)){
                cols_clus <- c(colorRampPalette(tol21rainbow)(length(unique(scar_info_clust_all$fish_clone))-1), 'light grey')
                names(cols_clus) <- unique(scar_info_clust_all$fish_clone)
                cols_gene <- colorRampPalette(gene_pal)(length(unique(scar_info_clust_all$Gene)))
                names(cols_gene) <- unique(scar_info_clust_all$Gene)
            }else{
                cols_clus <- c(tol21rainbow[1:(length(unique(scar_info_clust_all$fish_clone))-1)], 'light grey')
                names(cols_clus) <- unique(scar_info_clust_all$fish_clone)
                cols_gene <- gene_pal[1:length(unique(scar_info_clust_all$Gene))]
                names(cols_gene) <- unique(scar_info_clust_all$Gene)
            }


            cols_bcs <- c('blue','orange')[1:length(unique(cell_info$passed_filter))]
            names(cols_bcs) <- unique(cell_info$passed_filter)


            ann_colors = list(
                fish_clone = cols_clus,
                passed_filter = cols_bcs,
                Gene = cols_gene

            )

            ploto <- pheatmap(allls_cbs,
              show_rownames = FALSE, show_colnames = FALSE,
              breaks = seq(-1, +1, length = 101),
              annotation_row = cell_info,
              annotation_col = scar_info_clust_all[c('Gene','fish_clone')],
              annotation_colors = ann_colors,
              cluster_cols = T,
              main = 'combseqIDs_colored_by_cluster'

            )

            ggsave(plot = ploto, filename = paste0(dat_wd,"/pics/cloneDefPlots_",dat_name,"_",tum_it,"_",resolution,"_Heatmap_allCells_finalCloneAssignments.png"), width = 8 , height = 10, dpi = 500)
            rm(ploto)

                          
        print(paste0('Now plotting final clones for sample ',tum_it))                 
        # Then plot the cloneIDs for all cells
        filtered_cell_clone_object_inf$inferred_clone_NAless <- as.character(filtered_cell_clone_object_inf$inferred_clone)
        filtered_cell_clone_object_inf$inferred_clone_NAless[is.na(filtered_cell_clone_object_inf$inferred_clone_NAless)] <- 'none'

        # Get matrix of cell barcode vs seq_id with integration id
            allls_cbs <- as.data.frame.matrix(table(filtered_cell_clone_object_inf$Barcode, filtered_cell_clone_object_inf$inferred_clone_NAless))
            table(rowSums(allls_cbs))
            allls_cbs$Barcode <- rownames(allls_cbs)
            allls_cbs <- allls_cbs[!duplicated(allls_cbs$Barcode),]
            dim(allls_cbs)
            table(duplicated(allls_cbs$Barcode))
            allls_cbs$Barcode <- NULL

            scar_info_clust <- as.data.frame(sort(unique(filtered_cell_clone_object_inf$inferred_clone_NAless)))
            colnames(scar_info_clust) <- 'inferred_clone'
            rownames(scar_info_clust) <- scar_info_clust$inferred_clone


            # Create a dataframe for the barcodes
            cell_info <- as.data.frame(unique(rownames(allls_cbs)))
            colnames(cell_info) <- 'Barcode'
            rownames(cell_info) <- cell_info$Barcode
            cell_info$passed_filter <- 'no'
            cell_info$passed_filter[cell_info$Barcode %in% filtered_cell_clone_object_inf$Barcode[!is.na(filtered_cell_clone_object_inf$inferred_clone)]] <- 'yes'
            cell_info$Barcode <- NULL


            table(unique(colnames(allls_cbs)) %in% rownames(scar_info_clust))
            table(unique(rownames(allls_cbs)) %in% rownames(cell_info))


            if(length(unique(scar_info_clust$fish_clone))-1 > length(tol21rainbow)){
                cols_clus <- c(colorRampPalette(tol21rainbow)(length(unique(scar_info_clust$inferred_clone))-1), 'light grey')
                names(cols_clus) <- unique(scar_info_clust$inferred_clone)
            }else{
                cols_clus <- c(tol21rainbow[1:(length(unique(scar_info_clust$inferred_clone))-1)], 'light grey')
                names(cols_clus) <- unique(scar_info_clust$inferred_clone)

            }


            cols_bcs <- c('blue','orange')[1:length(unique(cell_info$passed_filter))]
            if(length(cols_bcs) == 2){
                
                names(cols_bcs) <- c('yes','no')
            
            }else{
                
                names(cols_bcs) <- unique(cell_info$passed_filter)

            }


            ann_colors = list(
                inferred_clone = cols_clus,
                passed_filter = cols_bcs

            )


            ploto <- pheatmap(allls_cbs,
              show_rownames = FALSE, show_colnames = FALSE,
              breaks = seq(-1, +1, length = 101),
              annotation_row = cell_info,
              annotation_col = scar_info_clust,
              annotation_colors = ann_colors,
              cluster_cols = T,
              main = 'cloneIDs_colored_by_cluster'

            )
            
            ggsave(plot = ploto, filename = paste0(dat_wd,"/pics/cloneDefPlots_",dat_name,"_",tum_it,"_",resolution,"Clones_resCutoff_",child_to_parent_fraction_cutoff,"_Heatmap_finalCloneIDs.png"), width = 8 , height = 10, dpi = 500)
            rm(ploto)

            
    }


}
