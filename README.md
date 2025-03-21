# Marine Food Web Analysis

## Content Description

### Input Files
- **F_KL-77.csv** | Matrix containing 297 marine families across 42 samples (spanning 124k years) from the SO-201-2-77KL Bering Sea sediment core. The abundances represent family clade counts classified using Kraken2 with a confidence threshold of 0.2.
  - Dataset Source (Stella Z. Buchwald, 2024):
    - [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/view/PRJEB66201)
    - [Zenodo](https://zenodo.org/records/10064386)
    - The resulting Kraken2 reports from step2 shotgun sequencing were merged using this command:
        ```bash
          awk  'BEGIN{FS=OFS="\t"} {print FILENAME, $0}' *report | awk 'BEGIN{FS=OFS="\t"} {gsub(/^[ \t]+/, "", $7)}1' > KL-77_nt0.2.txt
        ```
    - the merged Kraken2 report file was then used to filter the families clade counts, resulting in the given F_KL-77.csv file
- **families.csv** | Mapping of 297 families to 16 taxonomic groups.
- **Grant_age_RSL_col_EF.txt** | Relative sea level reconstruction data.
- **ngrip.csv** | NGRIP temperature proxy data.
- **trophic_levels.csv** | Mapping of 16 taxonomic groups to 4 trophic levels.
- **db_allresults_final.csv** | GloBI interaction results (queried 14.11.2024) for the 297 marine families.
- **IP25_77KL.csv** | IP25 values for the 77KL sediment core.
- **interpol_IP25_SSTs_stack_average_praetorius_NAreplaced.txt** | IP25 values for the 12KL sediment core.

### Scripts
- **1_inferESABO.py** | Python script to infer the ESABO network.
- **2_inferNets.R** | R script to infer six network types (Spearman, ecoCopula, Spiec-Easi, CCREPE, SPARCC, Propr).
- **3_consensus_network.py** | Python script for consensus network construction, robustness testing, and energy flow analysis.
- **4_fam_abunds.py** | Script for compositional analysis and generating abundance plots.
- **5_proxie_corr.py** | Script for correlating IP25 proxies from cores 77KL and 12KL with NGRIP.
- **6_globi_interactions.py** | Script for comparing the consensus network with GloBI to identify direct and indirect interactions.
- **7_module_maturity_plots.py** | Script for analyzing module flow, comparing bottom-up, top-down, and intraguild flows across interglacial and glacial clusters.

## Setup Workflow

### Clone Repository
To get started, clone this repository into your working directory:
```bash
git clone https://github.com/v-dinkel/MarineFoodWeb-sedaDNA.git
```

### Environment Setup
Follow the setup instructions in the [FoodWeb_gLV repository](https://github.com/v-dinkel/FoodWeb_gLV).

#### 1. Python Environment Setup (Skip Snakemake Installation)
- Install **scikit-learn (1.6.1)** in the conda environment:
  ```bash
  conda activate FoodWeb_gLV
  conda install scikit-learn
  ```
- Open `config.txt` and update `workdir` to match your local folder path.

#### 2. R Environment Setup
- Open `2_inferNetworks.R` and update `workdir` to match your local folder path.

## Analysis Workflow

### 1. Infer ESABO Network (Python)
- Run `1_inferESABO.py` in a Python editor (e.g., Spyder).
- Outputs: `output/KL-77_esabo.csv`

### 2. Infer Additional Networks (R)
- Run `2_inferNetworks.R` in RStudio.
- Outputs: Six network files, relative abundance matrix, and species correlations with NGRIP (saved in the `output` folder).

### 3. Build Consensus Network and Analyze
- Run `3_consensus_network.py`.
- Performs:
  - Consensus network construction from seven networks.
  - Positive/negative NGRIP correlation linkage analysis.
  - Network coverage statistics.
  - Network export to **Gephi-compatible (.gml) format**.
  - Computation of network statistics (modularity, composition, trophic levels, etc.).
  - Robustness analysis (set `runRobustness = False` to skip).
  - Comparison with Spiec-Easi network.
  - Energy flow metrics computation (relative ascendency of LCC and modules).
- Outputs:
  - Consensus network as a **.csv adjacency matrix**.
  - **Gephi files**, additional results in `/plots` and `/supplementary_information`.

### 4. Plot Family Stratigraphy
- Run `4_fam_abunds.py`.
- Output: Stratigraphy plots saved in `/plots`.

### 5. Correlate IP25 and NGRIP
- Run `5_proxie_corr.py`.
- Compares IP25 values with NGRIP.
- Outputs correlation coefficients and plots in `/plots`.

### 6. Identify Direct and Indirect Interactions (GloBI)
- Run `6_globi_interactions.py`.
- Requires a pre-generated GloBI query results file (included in the repository).
- Matches consensus network edges with GloBI interactions.
- Computes randomized null model for interaction occurrences.
- Compares results with Spiec-Easi.
- Outputs interaction graph (`outputs/`) and null model benchmarks (`plots/`).

### 7. Compute Flow Metrics for Module ig1
- Run `7_module_maturity_plots.py`.
- Analyzes **bottom-up, top-down, and intraguild flows**.
- Saves plots in the `/plots` directory.
