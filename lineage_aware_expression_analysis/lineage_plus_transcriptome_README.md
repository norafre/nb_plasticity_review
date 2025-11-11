# Analysis of transcriptional profiles across clones and transplantation time points
Processed single cell transcriptomes are brought together with gene module scores and cell clonal assignments. This is done for high-resolution primary tumor clones [here](./merge_transcriptome_modScores_clones.ipynb) and used to:  
**[on primary tumor data](./clone_transcriptome_analysis.ipynb):**
- Assess overall inter-clonal expression differences.
- Perform differential module expression analysis between clones.
- Calculate intra-clone module expression variance.
  
Processed single cell transcriptomes are brought together with gene module scores and cell clonal assignments. This is done for low-resolution clones for cells from all sample types [here](./merge_transcriptome_modScores_clones_allTPs.ipynb) and used to:  
**[on data from primary tumors and allograft samples](./clone_transcriptome_analysis_allTPs.ipynb):**
- Assess overall expression differences across sample types in general and within clones.
- Perform differential module expression analysis between cells from different sample types (e.g. primary tumor vs. late allograft tumor).
- Perform differential module expression analysis between clones within allograft tumors.A
- Calculate intra-clone module expression variance.
