#############################
# Parent script
### calls:
###### generate_data.R, three_groups_nimble.R, pickrell_data.rds, 
######     competitors_methods.R, evaluation_metrics.R
#############################
# # save arguments from bash 
# bash_args = commandArgs(trailingOnly = T)
# # test if there is at least one argument: if not, return an error
# if (length(bash_args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# } else if (length(bash_args)==1) {
#   # default output file
#   bash_args[2] = "out.txt"
# }
# #############################
# if(!exists("bash_args")){           # for debugging
#   bash_args <- c(127, ".....TEMP....")
# }
# #############################
# my_seed1 <- as.numeric(bash_args[1])/127


library(edgeR)
library(seqgendiff)
library(nimble)
library(tictoc)
library(limma)
library(qvalue)
library(MASS)
library(dplyr)
library(DESeq2)
tic()

my_seed1 <- 1
set.seed(my_seed1*127)

# !!! SAVE NAMES !!! 
save_location <- "./"
save_name <- paste0("simulation_", my_seed1, "_", Sys.Date())


# data shape
num_individuals_RNA  <- 100
num_individuals_GWAS <- 1000
num_genes        <- 250
num_beneficial   <- 10 
GWAS_effect      <- c(rep(0.5, 7), rep(0.1, 3), 
                      rep(0.5, 7), rep(0.1, 3), rep(0, 230))
RNA_effect       <- c(rep(1.05, 3), rep(1.4, 7), 
                      rep(1.05, 3), rep(1.4, 7), rep(0, 230))

# data generation:
source("generate_data.R") 
  # This script loads Montgomery_and_Pickrell.rds which are published data 
  #   Pickrell et. al. (2010) and Montgomery et. al. (2010) and can be accessed
  #   here:   https://bowtie-bio.sourceforge.net/recount/

#############################
# Three_groups
models_for_mcmc  <- c("combined", "RNA_only", "GWAS_only")
priors_for_mcmc <- c("piMOM", "local")
niter           <- 12500
nburnin         <- 2500
thin            <- 2
num_samples     <- (niter-nburnin)/thin

mcmc_save_name <- paste0(save_location, "samples_", save_name)
 
source("TG_nimble_model.R")
# results in post_probs_null_TG

#############################
# competitors
source("competitors_methods.R")

#############################
# Evaluate
source("evaluation_metrics.R")
three_groups_combined_piMOM <- 
  evaluation_function(post_probs_null_TG$combined_piMOM_samples, 
                      groups^2, method = "TG_combined_piMOM", my_seed=my_seed1)
three_groups_RNA_piMOM      <- 
  evaluation_function(post_probs_null_TG$RNA_only_piMOM_samples,
                      groups^2, method = "TG_RNA_piMOM", my_seed=my_seed1)
three_groups_GWAS_piMOM      <- 
  evaluation_function(post_probs_null_TG$GWAS_only_piMOM_samples, 
                      groups^2, method = "TG_GWAS_piMOM", my_seed=my_seed1)
three_groups_combined_local <- 
  evaluation_function(post_probs_null_TG$combined_local_samples,
                      groups^2, method = "TG_combined_local", my_seed=my_seed1)
three_groups_RNA_local      <- 
  evaluation_function(post_probs_null_TG$RNA_only_local_samples, 
                      groups^2, method = "TG_RNA_local", my_seed=my_seed1)
three_groups_GWAS_local      <- 
  evaluation_function(post_probs_null_TG$GWAS_only_local_samples, 
                      groups^2, method = "TG_GWAS_local", my_seed=my_seed1)
GWAS_only                   <- 
  evaluation_function(post_probs_null_GWAS, groups^2, 
                      method = "GWAS", my_seed=my_seed1)
edgeR_only                  <- 
  evaluation_function(post_probs_null_edgeR, groups^2, 
                      method = "edgeR", my_seed=my_seed1)
voom_only                   <- 
  evaluation_function(post_probs_null_voom, groups^2, 
                      method = "voom", my_seed=my_seed1)
DEseq2_only                 <- 
  evaluation_function(post_probs_null_DEseq2, groups^2,
                      method = "DEseq2", my_seed=my_seed1)
GWAS_edgeR                  <- 
  evaluation_function(post_probs_null_GWAS_edgeR, groups^2,
                      method = "GWAS_edgeR", my_seed=my_seed1)
GWAS_voom                   <- 
  evaluation_function(post_probs_null_GWAS_voom, groups^2,
                      method = "GWAS_voom", my_seed=my_seed1)
GWAS_DEseq2                 <- 
  evaluation_function(post_probs_null_GWAS_DEseq2, groups^2,
                      method = "GWAS_DEseq2", my_seed=my_seed1)

myresults <- rbind(three_groups_combined_piMOM, three_groups_RNA_piMOM, 
                   three_groups_GWAS_piMOM, 
                   three_groups_combined_local, three_groups_RNA_local, 
                   three_groups_GWAS_local,
                   GWAS_only, edgeR_only, voom_only, DEseq2_only, 
                   GWAS_edgeR, GWAS_voom, GWAS_DEseq2)

saveRDS(myresults, 
        file=paste0(save_location, "results_", save_name, ".rds"))

myresults_null <- 
  cbind(TG_comb_piMOM = post_probs_null_TG$combined_piMOM_samples, 
        TG_comb_local = post_probs_null_TG$combined_local_samples, 
        TG_RNA_poMOM = post_probs_null_TG$RNA_only_piMOM_samples,
        TG_RNA_local = post_probs_null_TG$RNA_only_local_samples,
        TG_GWAS_piMOM = post_probs_null_TG$GWAS_only_piMOM_samples,
        TG_GWAS_local = post_probs_null_TG$GWAS_only_local_samples,
        GWAS = post_probs_null_GWAS, 
        edgeR = post_probs_null_edgeR,
        voom = post_probs_null_voom,
        DEseq2 = post_probs_null_DEseq2,
        GWAS_edgeR = post_probs_null_GWAS_edgeR,
        GWAS_voom = post_probs_null_GWAS_voom,
        GWAS_DEseq2 = post_probs_null_GWAS_DEseq2)
myresults_ben <- 
  cbind(TG_comb_piMOM = post_probs_ben_TG$combined_piMOM_samples,
        TG_comb_local = post_probs_ben_TG$combined_local_samples,
        TG_RNA_poMOM = post_probs_ben_TG$RNA_only_piMOM_samples,
        TG_RNA_local = post_probs_ben_TG$RNA_only_local_samples,
        TG_GWAS_piMOM = post_probs_ben_TG$GWAS_only_piMOM_samples,
        TG_GWAS_local = post_probs_ben_TG$GWAS_only_local_samples)

saveRDS(myresults_null, 
        file=paste0(save_location, "results_null_", save_name, ".rds"))
saveRDS(myresults_ben, 
        file=paste0(save_location, "results_ben_", save_name, ".rds"))
toc()