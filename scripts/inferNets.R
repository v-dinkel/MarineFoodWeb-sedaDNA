#################
## inferNets.R ##
#################

#install.packages(ecoCopula)
library(analogue)
library(psych)
library(NetCoMi)
library(tidyverse)
library(ecoCopula)

# Load KL-77 Dataset
s_M <- read.csv("C:/Users/vdinkel/Desktop/Manuscript/input/F_KL-77.csv",header = TRUE)
rownames(s_M) <- s_M$X
s_M$X <- NULL

# Calculate relative abundances (Total Sum Scaling)
s_M_rel <-t(apply(s_M,1, function(x) x/sum(x)))
write.table(s_M_rel, "C:/Users/vdinkel/Desktop/Manuscript/output/F_KL-77_rel.csv", row.names = TRUE, col.names = TRUE, sep=";")

normmethod <- "none"

# Spiec-Easi
s_net_spieceasi <- netConstruct(s_M_rel,
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

s_A_spieceasi <- s_net_spieceasi$adjaMat1
diag(s_A_spieceasi) <- 0
sum(s_A_spieceasi)
write.table(s_A_spieceasi, "C:/Users/vdinkel/Desktop/Manuscript/output/KL-77_spieceasi.csv", row.names = TRUE, col.names = TRUE, sep=";")

# Weighted Spiec-Easi (for comparison)
s_net_spieceasi <- netConstruct(s_M_rel,
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
s_A_spieceasi <- s_net_spieceasi$adjaMat1
diag(s_A_spieceasi) <- 0
sum(s_A_spieceasi)
write.table(s_A_spieceasi, "C:/Users/vdinkel/Desktop/Manuscript/output/KL-77_spieceasi_weighted.csv", row.names = TRUE, col.names = TRUE, sep=";")
write.table(s_net_spieceasi$edgelist1, "C:/Users/vdinkel/Desktop/Manuscript/output/KL-77_spieceasi_weighted_edgelist.csv", row.names = TRUE, col.names = TRUE, sep=";")


# CCREPE
s_net_ccrepe <- netConstruct(s_M_rel,
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
s_A_ccrepe <- s_net_ccrepe$adjaMat1
diag(s_A_ccrepe) <- 0
sum(s_A_ccrepe)
write.table(s_A_ccrepe, "C:/Users/vdinkel/Desktop/Manuscript/output/KL-77_ccrepe.csv", row.names = TRUE, col.names = TRUE, sep=";")

# PROPR
s_net_propr <- netConstruct(s_M_rel,
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
s_A_propr <- s_net_propr$adjaMat1
diag(s_A_propr) <- 0
sum(s_A_propr)
write.table(s_A_propr, "C:/Users/vdinkel/Desktop/Manuscript/output/KL-77_propr.csv", row.names = TRUE, col.names = TRUE, sep=";")

# Sparcc
# Using s_M_rel returns an empty network. Hence, the raw counts are used with intrinsic "TSS" normalization
s_net_sparcc <- netConstruct(s_M,
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

s_A_sparcc <- s_net_sparcc$adjaMat1
diag(s_A_sparcc) <- 0
sum(s_A_sparcc)
write.table(s_A_sparcc, "C:/Users/vdinkel/Desktop/Manuscript/output/KL-77_sparcc.csv", row.names = TRUE, col.names = TRUE, sep=";")

# spearman
s_net_spearman <- netConstruct(s_M_rel,
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

s_A_spearman <- s_net_spearman$adjaMat1
diag(s_A_spearman) <- 0
sum(s_A_spearman)
write.table(s_A_spearman, "C:/Users/vdinkel/Desktop/Manuscript/output/KL-77_spearman.csv", row.names = TRUE, col.names = TRUE, sep=";")

# NGRIP
ngrip <- read.csv("C:/Users/vdinkel/Desktop/Manuscript/input/ngrip.csv", header = TRUE, sep = "\t")
colnames(ngrip)[1]  <- "age"
colnames(ngrip)[2]  <- "temp"
ngrip$age <- ngrip$age * 1000

# GET VALUES FROM NGRIP CLOSEST TO AGES IN SAMPLES
# KL-77
s_indeces <- c()
for (val in as.integer(rownames(s_M_rel))){
  s_indeces <- append(s_indeces, which.min(abs(ngrip$age - val)))
}
s_temps <- ngrip[s_indeces,]$temp
s_mytemps <- (s_temps - mean(s_temps))/sd(s_temps) # centered & standardized
s_exp_temps <-  ngrip[s_indeces,]
s_exp_temps$temp <- s_mytemps
write.table(s_exp_temps, "C:/Users/vdinkel/Desktop/Manuscript/output/KL-77_ngrip.csv", sep=";", dec=".")

corrs <- corr.test(x = s_M_rel, y = s_mytemps, method ="spearman", adjust="holm")
corrs$indeces <- which(corrs$p %in% corrs$p[corrs$p <= 0.05]) # p.adj
corrs$corr_names <- rownames(corrs$p)[corrs$indeces]
corrs$coeffs <- corrs$r[corrs$indeces]
corrs$ps <- corrs$p[corrs$indeces]
export_corrs = data.frame()
export_corrs <- rbind(corrs$corr_names, corrs$coeffs)
export_ps = data.frame()
export_ps <- rbind(corrs$corr_names, corrs$ps)
write.table(export_corrs, "C:/Users/vdinkel/Desktop/Manuscript/output/KL-77_temp_spearman_corrs.csv", sep=";", dec=".")
write.table(export_ps, "C:/Users/vdinkel/Desktop/Manuscript/output/KL-77_temp_spearman_corrs_p.csv", sep=";", dec=".")

# ECOCOPULA
s_ngrip = data.frame(ngrip = s_mytemps)
#s_ngrip$ngrip <- 1 # set all ngrip datapoijnts to 1 to remove any temperature input
row.names(s_ngrip) <- row.names(s_M)
s_mod <- stackedsdm(s_M, ~., data = s_ngrip, ncores=4)
s_ecocopula = cgr(s_mod, lambda = 0.5) #0.5 # 0.788
s_ecocopula_A = s_ecocopula$best_graph$graph
diag(s_ecocopula_A) <- 0
sum(s_ecocopula_A)
#write.csv(s_ecocopula_A , "C:/Users/vdinkel/Desktop/Manuscript/output/KL-77_ecocopula.csv", row.names = TRUE, col.names = TRUE, sep=";")
write.csv(s_ecocopula_A , "C:/Users/vdinkel/Desktop/Manuscript/output/KL-77_ecocopula_ngrip.csv", row.names = TRUE, col.names = TRUE, sep=";")
