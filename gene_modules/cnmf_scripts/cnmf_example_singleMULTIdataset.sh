conda activate cnmf_env

# Set variables for blas control to limit number of threads used
export MKL_NUM_THREADS=40
export NUMEXPR_NUM_THREADS=40
export OMP_NUM_THREADS=40

# Prepare - pass file names, locations and ks to work with and number of iterations to perform
cnmf prepare --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes -c ./multi_bAllos_2_allos_nrblstm_counts.tsv -k 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 --n-iter 200 --total-workers 20 --seed 14 --genes-file multi_bAllos_2_allos_nrblstm_varGenes.txt

# Factorize - parallelize by using 20 workers
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 0 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 1 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 2 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 3 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 4 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 5 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 6 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 7 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 8 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 9 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 10 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 11 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 12 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 13 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 14 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 15 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 16 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 17 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 18 --total-workers 20 &
cnmf factorize --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --worker-index 19 --total-workers 20

# Combine results from individual iterations
cnmf combine --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes

# Print k-selection plot
cnmf k_selection_plot --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes

# Generate consensus
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 3 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 4 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 5 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 6 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 7 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 8 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 9 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 10 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 11 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 12 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 13 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 14 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 15 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 16 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 17 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 18 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 19 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 20 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 21 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 22 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 23 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 24 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 25 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 26 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 27 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 28 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 29 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 30 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 31 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 32 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 33 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 34 --local-density-threshold 0.2 --show-clustering
cnmf consensus --output-dir ../nmf_data/out --name multi_bAllos_2_allos_nrblstm_seurVarGenes --components 35 --local-density-threshold 0.2 --show-clustering

