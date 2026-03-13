# three_groups_simulations
Houses the simulation scripts for the three groups GWAS and RNA-seq model (publication forthcoming). 

The code for the real data analysis is the same as in the TG_nimble_model.R except for two noticable changes. First, the real data includes addition of other covariates. Second, the script used for the real data anlysis does not loop over different models, each model was run on its own script. Due to data privacy concerns we cannot share the data and thus sharing the script for that analysis is not useful. 


### **WARNING**: these simulations take a long time to run. It takes around 65 minutes to run on an M2pro macbook laptop with the `parent_script.R` options as they are. 


## To Run: 
These scripts were written to be called from a bash script so that simulations could be run in parallel with different random seeds. I have commented out the relevant lines at the beginning.

Below I will show you how to run a single simulation. Feel free to contact me if you want to learn how I ran multiple but I'm sure that you can figure out how to run several after learning how to run one. 

## To run a single simulation:
1. The following files ought to be in the same directory:
   a) `parent_script.R` - calls the relevant scripts, is set up for parallel simulations with bash scripting
   b) `generate_data.R` - generates synthetic data for simulation studies
   c) `Montgomery_and_Pickrell.rds` - RNA-seq data from Pickrell et. al. (2010) and Montgomery et. al. (2010) [accessible here](https://bowtie-bio.sourceforge.net/recount/)
   d) `TG_nimble_model.R` - three-groups model in nimble, this will loop over models and priors and save results after running
   e) `competitors_methods.R` - code to run conventional GWAS, edgeR, limma-voom, and DESeq2 with both the Fisher and Cauchy p-value combinations 
   f) `evaluation_metrics.R` - code to compute the logarithmic score and brier score from the output of `TG_nimble_model.R` and `competitors_methods.R`
3. There are several packages that are required: `edgeR`, `seqgendiff`, `nimble`, `tictoc`, `limma`, `qvalue`, `MASS`, `dplyr`, and `DESeq2`. 
4. Open `parent_script.R` and set your working directory to the directory with the aforementioned files. 
5. Set the `save_location` and `save_name` variables. 
6. Adjust the models to be run and MCMC iterations, etc.
7. Run `parent_script.R`. This will:
   - Generate synthetic data using the `generate_data.R` script. GWAS data are generated from a logistic regression model. RNA-seq data are generated through randomly selecting a subset of a real RNA-seq dataset and adding signal to the specified genes (this is done using the binomial thinning of Gerard (2020)).
   - Run `TG_nimble_model.R` which performs MCMC sampling of the three groups method on the same generated data for each combination of model and gene-effect prior you have specified.
   - Run the competitors methods on that same generated data with `competitors_methods.R`.
   - Compute evaluation metrics (like the logarithmic score, etc.) for all runs using `evaluation_metrics.R`.
   - Save evaluation metrics in a table and posterior probabilities of group membership in `.rds` files.


## Simualation with missingness
The real data has 34 genes that are only measured in the GWAS experiment. To run simulations under this setting we simply delete the RNA-seq information (replace with `NA`'s) for a handful of genes. There is code in `parent_script.R` to do this that is commented out.
