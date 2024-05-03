library(tidyverse)
library(cmdstanr)
library(rstan)

# Load the data
data <- read_csv("data.csv")

# Pivot wider into A x R matrix
data <- data %>% pivot_wider(id_cols     = c(animal_idx, group_idx),
                             values_from = counts,
                             names_from  = regions,
                             names_sort  = T)

# Set up model data
stan_data <- list("A"                = nrow(data),
                  "R"                = 50,
                  "G"                = 2,
                  "group_idx"        = data$group_idx,
                  "y"                = as.matrix(data[,-c(1:2)])) #  drop first two indexing columns

# Compile the model
mod <- cmdstan_model("models/model_hs.stan")

# Sample the model
fit <- mod$sample(data            = stan_data,
                  chains          = 4,
                  parallel_chains = 4,
                  iter_warmup     = 4000,
                  iter_sampling   = 4000,
                  seed            = 2017) # reproducibility
# Save fit
fit <- read_stan_csv(fit$output_files()) # Save output from cmdstanr in a way that preserves param layout.
saveRDS(fit, "fits/fit_hs.rds")
