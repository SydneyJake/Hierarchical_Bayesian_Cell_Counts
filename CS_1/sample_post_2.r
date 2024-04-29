library(tidyverse)
library(rstan)
library(cmdstanr)

# read the data
data <- read_csv("data.csv", col_types = cols(BrainRegion = col_factor(),
                                              rat_ID      = col_factor(),
                                              group       = col_factor()))

# counts are summed if there are multiple, same applies to the exposure
data <- data %>%
          group_by(BrainRegion, group, rat_ID) %>%
          summarise(y=sum(Cfos_Counts), area = sum( (1-VOID_PERCENTAGE/100) * (0.94*0.67) )) %>%
          ungroup()

# Data for the Stan model
stan_data = list(R          = 23,
                 G          = 4,
                 N          = length(data$y),
                 region_idx = as.numeric(data$BrainRegion),
                 group_idx  = as.numeric(data$group),
                 E          = log(data$area),
                 y          = data$y)

# Compile the model
mod <- cmdstan_model("models/model_poiss.stan")

# fit the model
fit <- mod$sample(data            = stan_data,
                  chains          = 4,
                  parallel_chains = 4,
                  iter_warmup     = 4000,
                  iter_sampling   = 4000)

fit <- read_stan_csv(fit$output_files())

# Save
saveRDS(fit, "fit_poiss.rds")
