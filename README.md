# Content Description
## Input
- <b>F_KL-77.csv</b> | Matrix of 297 marine families and 42 samples (124kyrs) from the SO-201-2-77KL Bering Sea sediment core. The abundances are the family clade counts given Kraken2 taxonomic classificaiton with confidence threshold of 0.2.
  - The origin of the dataset (Stella Z Buchwald 2024, https://academic.oup.com/ismej/article/18/1/wrad006/7513109):
    - https://www.ebi.ac.uk/ena/browser/view/PRJEB66201
    - https://zenodo.org/records/10064386 
- <b>families.csv</b> | Mapping of 297 families to 16 taxonomic groups
- <b>Grant_age_RSL_col_EF.txt</b> | Relative sea level reconstruction
- <b>ngrip.csv</b> | NGRIP proxy for temperature
- <b>trophic_levels.csv</b> | Mapping of 16 taxonomic groups to 4 trophic levels
## Scripts
- <b>inferNets.R</b> | R script which infers 6 network types (Spearman, ecoCopula, Spiec-Easi, CCREPE, SPARCC, Propr)
- <b>inferESABO.py</b> | python implementation to infer the ESABO network
- <b>consensus_network.py</b> | Main analysis python file to build consensus network and perform robustnes and energy flow analysis
- <b>fam_abunds.py</b> | Compositional analysis & plot genearion
- <b>compare_BN_HT.py</b> | Comparison script of the consensus network with the spiec-easi network being sparsificatied by association strength 
## Output
- <F_KL-77_rel.csv> | the relative abundance matrix 

## Setup Workflow
- Setup your Python and R environment as described in this project: https://github.com/v-dinkel/FoodWeb_gLV
    - 1. Setting Up the Python Environment (skip install Snakemake)
        - Open the config.txt and change workdir to match your folder path
    - 2. Setting Up the R Environment. The corresponding R file in this project is 2_inferNetworks.R
        - Open the 2_inferNetworks.R and change workdir to match your folder path
    
- Clone this repository into your working directory
    - git clone https://github.com/v-dinkel/MarineFoodWeb-sedaDNA

## Analysis Workflow
- infer ESABO network
    - open inferESABO.py in your editor e.g. Sypder run it
    - output: the ESABO network

- infer 6 other networks:
    - 
    
