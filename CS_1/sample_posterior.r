library(tidyverse)
library(rstan)
library(cmdstanr)

data <- read_csv("data.csv", show_col_types = FALSE)

# counts are summed if there are multiple, same applies to the exposure
data <- data %>%
          group_by(BrainRegion, group, rat_ID) %>%
          summarise(y=sum(Cfos_Counts), area = sum( (1-VOID_PERCENTAGE/100) * (0.94*0.67) )) %>%
          ungroup()

data <- data %>%
  mutate(BrainRegion = factor(BrainRegion), # alphabetical
         group       = factor(group)) # alphabetical

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
                  iter_sampling   = 4000,
                  seed            = 3511)

fit <- read_stan_csv(fit$output_files())

# Save
saveRDS(fit, "fits/fit_poiss.rds")
