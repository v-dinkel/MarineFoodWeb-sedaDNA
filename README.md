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
- <b>db_allresutls_final.csv</> | GloBI interactions results for the 297 marine families
- <b>IP25_77KL.csv</b> | IP25 values for the 77KL core
- <b>interpol_IP25_SSTs_stack_average_praetorius_NAreplaced.txt</b> | IP25 values for the 12KL core
## Scripts
- <b>1_inferESABO.py</b> | python implementation to infer the ESABO network
- <b>2_inferNets.R</b> | R script which infers 6 network types (Spearman, ecoCopula, Spiec-Easi, CCREPE, SPARCC, Propr)
- <b>3_consensus_network.py</b> | Main analysis python file to build consensus network and perform robustnes and energy flow analysis
- <b>4_fam_abunds.py</b> | Compositional analysis & plot genearion
- <b>5_proxie_corr.py</b> | Correlates IP25 proxies from core 77KL and 12KL with NGRIP
- <b>6_globi_interactions.py</b> | Compares the consensus network (and spiec easi) with GloBI to identify direct and indirect interactions.
- <b>7_module_maturity_plots.py</b> | Flow analysis of a module, specifically bottom-up, top-down and intraguild flows compared between two module states (interglacial and glacial cluster)


## Setup Workflow
- Setup your Python and R environment as described in this project: https://github.com/v-dinkel/FoodWeb_gLV
    - 1. Setting Up the Python Environment (skip install Snakemake)
        - install scikit-learn (1.6.1) in the conda environment:
            - open terminal and activate the environment
            - conda activate FoodWeb_gLV
            - conda install scikit-learn
        - Open the config.txt and change workdir to match your folder path
    - 2. Setting Up the R Environment. The corresponding R file in this project is 2_inferNetworks.R
        - Open the 2_inferNetworks.R and change workdir to match your folder path
    
- Clone this repository into your working directory
    - git clone https://github.com/v-dinkel/MarineFoodWeb-sedaDNA

## Analysis Workflow
- infer ESABO network with Python:
    - open and run 1_inferESABO.py in your Python editor e.g. Sypder
    - writes the ESABO network in output/KL-77_esabo.csv

- infer 6 other networks using R:
    - open and run 2_inferNetworks.R in RStudio (change the workdir if you haven't already).
    - writes 6 network files, relative abundance matrix and species correlations with ngrip to output folder 
    
- build consensus network and run analysis
    - open and run 3_consensus network.py
    - builds the consensus network from given 7 networks. Base_method defines the base method for the consensus network.
    - computes the positive/negative ngrip correlation linkage between nodes of a network. E.g. how likely an ngrip positively correlated note is linked to another one.
    - computes the percentage, how many edges of in the consensus network are covered by a single network
    - exports the consensus network to gephi compatible format (.gml)
    - computes network statistics such as modularity, module composition, mean trophic levels etc.
    - computes the module robustness by random extinction. Set runRobustness = False to bypass this step as it can be time consuming.
    - compares the modularity and clustering (positive/negative ngrip linkage) of the consensus network with spiec easi network of similiar size. plo
    - computes energy flow metrics such as relative ascendency of the LCC and three modules. Plots them against NGRIP and relative sea level.
    - outputs are the consensus network as .csv adjacency matrix and a gephi file. Other generated outputs are in /plots and /supplementary_information folders

- plot families stratigram
    - open and run 4_fam_abunds.py
    - it saves the stratigram into /plots

- correlate ip25 and ngrip from cores kl12 and kl77
    - open 5_proxie_corr.py and run it
    - loads and compares IP25 values with ngrip
    - prints correlation coefficients and saves the plot to /plots
    
- identify direct and indirect interactions using GloBI
    - open 6_globi_interactions.py and run it
    - this script requires a file with query results of the globi database. in this repository this file is provided.
    - matches the consensus network edges with globi interaction results
    - identifies direct and indirect interactions
    - computes randomized null model of randomly occuring interactions
    - compares with spiec easi results
    - saves interaction graph in outputs and null model benchmark to plots
    
- get flow metrics of the module ig1
    - open 7_module_maturity_plots.py and run it
    - computes the flow metrics of module ig1 and creates plots
    - divides the analysis into bottom-up, top-down and intraguild flows
    - saves the plots