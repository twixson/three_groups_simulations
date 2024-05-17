# three_groups_simulations
Houses the simulation scripts for the three groups GWAS and RNA-seq model (publication forthcoming). 


### **WARNING**: these simulations take a long time to run. It takes X hours to run on my laptop with the `parent_script.R` options as they are.


## To Run: 
These scripts were written to be called from a bash script so that simulations could be run in parallel with different random seeds. I have commented out the relevant lines at the beginning.

Below I will show you how to run a single simulation. Feel free to contact me if you want to learn how I ran multiple but I'm sure that you can figure out how to run several after learning how to run one. 

## To run a single simulation:
1. The following files ought to be in the same directory: `parent_script.R` `generate_data.R`, `Montgomery_and_Pickrell.rds`, `TG_nimble_model.R`, and `competitors_methods.R`.
2. There are several packages that are required: `edgeR`, `seqgendiff`, `nimble`, `tictoc`, `limma`, `qvalue`, `MASS`, `dplyr`, and `DESeq2`. 
3. Open `parent_script.R` and set your working directory to the directory with the aforementioned files. 
4. Set the `save_location` and `save_name` variables. 
5. Adjust the models to be run and MCMC iterations, etc.
6. Run `parent_script.R`. This will:
   - Generate synthetic data using the `generate_data.R` script. GWAS data are generated from a logistic regression model. RNA-seq data are generated through randomly selecting a subset of a real RNA-seq dataset and adding signal to the specified genes (this is done using the binomial thinning of Gerard (2020)).
   - Run `TG_nimble_model.R` which performs MCMC sampling of the three groups method on the same generated data for each combination of model and gene-effect prior you have specified.
   - Run the competitors methods on that same generated data with `competitors_methods.R`.
   - Compute evaluation metrics (like the logarithmic score, etc.) for all runs using `evaluation_metrics.R`.
   - Save evaluation metrics in a table and posterior probabilities of group membership in `.rds` files.

