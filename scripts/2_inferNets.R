#################
## 2_inferNets.R ##
#################

# Some installations which might need to be handled manually:
#install.packages("remotes")
#install.packages("ragg")
#install.packages("devtools")
#install.packages("BiocManager")
#remotes::install_version("Matrix", version = "1.6")
#install.packages("pulsar")
#devtools::install_github("zdk123/SpiecEasi")
#install.packages("sn")
#devtools::install_github("GraceYoon/SPRING")
#remotes::install_version("Hmisc", version = "5.1")
#devtools::install_github("stefpeschel/NetComi", repos= c("https://cloud.r-project.org/", BiocManager::repositories()))
#install.packages("RcppGSL") #this is system requirement
#install.packages("ecoCopula")
#devtools::install_github("tpq/propr")


library(analogue)
library(psych)
library(NetCoMi)
library(tidyverse)
library(ecoCopula)

# Install the rstudioapi package if not installed
if (!require(here)) install.packages("here")
library(here)

# Set the working directory to your root of the repository
workdir <- "/home/viktor/project_migration/MarineFoodWeb-sedaDNA/"

# Load KL-77 Dataset
M <- read.csv(here(workdir, "input/F_KL-77.csv"), header = TRUE)
rownames(M) <- M$X
M$X <- NULL

# Calculate relative abundances (Total Sum Scaling)
M_rel <-t(apply(M,1, function(x) x/sum(x)))
write.table(M_rel, here(workdir, "/output/F_KL-77_rel.csv"), row.names = TRUE, col.names = TRUE, sep=";")

normmethod <- "none"

# Spiec-Easi
net_spieceasi <- netConstruct(M_rel,
                              dataType = "counts",
                              measure = "spieceasi",
                              filtSamp = "none",
                              filtSampPar = "none",
                              normMethod = normmethod,
                              sparsMethod = 'threshold',
                              weighted = FALSE,
                              thresh = 0.4,
                              verbose = 3,
                              seed = 100
                              )

A_spieceasi <- net_spieceasi$adjaMat1
diag(A_spieceasi) <- 0
sum(A_spieceasi)
write.table(A_spieceasi, here(workdir, "output/KL-77_spieceasi.csv"), row.names = TRUE, col.names = TRUE, sep=";")

# Weighted Spiec-Easi (for comparison)
w_net_spieceasi <- netConstruct(M_rel,
                                dataType = "counts",
                                measure = "spieceasi",
                                filtSamp = "none",
                                filtSampPar = "none",
                                normMethod = normmethod,
                                sparsMethod = 'threshold',
                                weighted = TRUE,
                                thresh = 0.4,
                                verbose = 3,
                                seed = 100)
w_A_spieceasi <- w_net_spieceasi$adjaMat1
diag(w_A_spieceasi) <- 0
sum(w_A_spieceasi)
write.table(w_A_spieceasi, here(workdir, "output/KL-77_spieceasi_weighted.csv"), row.names = TRUE, col.names = TRUE, sep=";")
write.table(w_net_spieceasi$edgelist1, here(workdir, "output/KL-77_spieceasi_weighted_edgelist.csv"), row.names = TRUE, col.names = TRUE, sep=";")

# CCREPE
net_ccrepe <- netConstruct(M_rel,
                             dataType = "counts",
                             measure = "ccrepe",
                             filtSamp = "none",
                             filtSampPar = "none",
                             normMethod = normmethod,
                             sparsMethod = 'threshold',
                             thresh = 0.57,
                             weighted = FALSE,
                             verbose = 3,
                             seed = 100)
A_ccrepe <- net_ccrepe$adjaMat1
diag(A_ccrepe) <- 0
sum(A_ccrepe)
write.table(A_ccrepe, here(workdir, "output/KL-77_ccrepe.csv"), row.names = TRUE, col.names = TRUE, sep=";")

# PROPR
net_propr <- netConstruct(M_rel,
                            dataType = "counts",
                            measure = "propr",
                            filtSamp = "none",
                            filtSampPar = "none",
                            normMethod = normmethod,
                            sparsMethod = 'threshold',
                            weighted = FALSE,
                            thresh = 0.482,
                            verbose = 3,
                            seed = 100)
A_propr <- net_propr$adjaMat1
diag(A_propr) <- 0
sum(A_propr)
write.table(A_propr, here(workdir, "output/KL-77_propr.csv"), row.names = TRUE, col.names = TRUE, sep=";")

# Sparcc
# Using M_rel returns an empty network. Hence, the raw counts are used with intrinsic "TSS" normalization
net_sparcc <- netConstruct(M,
                             dataType = "counts",
                             measure = "sparcc",
                             filtSamp = "none",
                             filtSampPar = "none",
                             normMethod = "TSS",
                             weighted = FALSE,
                             sparsMethod = 't-test',
                             alpha = 0.25,
                             verbose = 3,
                             seed = 100)

A_sparcc <- net_sparcc$adjaMat1
diag(A_sparcc) <- 0
sum(A_sparcc)
write.table(A_sparcc, here(workdir, "output/KL-77_sparcc.csv"), row.names = TRUE, col.names = TRUE, sep=";")

# spearman
net_spearman <- netConstruct(M_rel,
                               dataType = "counts",
                               measure = "spearman",
                               filtSamp = "none",
                               filtSampPar = "none",
                               normMethod = normmethod,
                               weighted = FALSE,
                               sparsMethod = 'threshold',
                               thresh = 0.52,
                               verbose = 3,
                               seed = 100)

A_spearman <- net_spearman$adjaMat1
diag(A_spearman) <- 0
sum(A_spearman)
write.table(s_A_spearman, here(workdir, "output/KL-77_spearman.csv"), row.names = TRUE, col.names = TRUE, sep=";")

# NGRIP
ngrip <- read.csv(here(workdir, "input/ngrip.csv"), header = TRUE, sep = "\t")
colnames(ngrip)[1]  <- "age"
colnames(ngrip)[2]  <- "temp"
ngrip$age <- ngrip$age * 1000

# GET VALUES FROM NGRIP CLOSEST TO AGES IN SAMPLES
# KL-77
indeces <- c()
for (val in as.integer(rownames(M_rel))){
  indeces <- append(indeces, which.min(abs(ngrip$age - val)))
}
temps <- ngrip[indeces,]$temp
mytemps <- (temps - mean(temps))/sd(temps) # centered & standardized
exp_temps <-  ngrip[indeces,]
exp_temps$temp <- mytemps
write.table(exp_temps, here(workdir, "output/KL-77_ngrip.csv"), sep=";", dec=".")

corrs <- corr.test(x = M_rel, y = mytemps, method ="spearman", adjust="holm")
corrs$indeces <- which(corrs$p %in% corrs$p[corrs$p <= 0.05])
corrs$corr_names <- rownames(corrs$p)[corrs$indeces]
corrs$coeffs <- corrs$r[corrs$indeces]
corrs$ps <- corrs$p[corrs$indeces]
export_corrs = data.frame()
export_corrs <- rbind(corrs$corr_names, corrs$coeffs)
export_ps = data.frame()
export_ps <- rbind(corrs$corr_names, corrs$ps)
write.table(export_corrs, here(workdir, "output/KL-77_temp_spearman_corrs.csv"), sep=";", dec=".")
write.table(export_ps, here(workdir, "output/KL-77_temp_spearman_corrs_p.csv"), sep=";", dec=".")

# ECOCOPULA
ngrip = data.frame(ngrip = mytemps)
ngrip$ngrip <- 1 # set all ngrip datapoints to 1 to remove any temperature input. comment this out if you want to keep environmental factor in
row.names(ngrip) <- row.names(M)
mod <- stackedsdm(M, ~., data = ngrip, ncores=4)
ecocopula = cgr(mod, lambda = 0.5) #0.5 # 0.788
ecocopula_A = ecocopula$best_graph$graph
diag(ecocopula_A) <- 0
sum(ecocopula_A)
write.csv(ecocopula_A , here(workdir, "output/KL-77_ecocopula.csv"), row.names = TRUE, col.names = TRUE, sep=";")

# if you computed ecocopula using unmodified ngrip, then I suggest you include it in the saved filename:
#write.csv(ecocopula_A , here(workdir, "output/KL-77_ecocopula_ngrip.csv"), row.names = TRUE, col.names = TRUE, sep=";")
