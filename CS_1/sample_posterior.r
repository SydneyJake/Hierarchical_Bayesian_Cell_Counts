library(tidyverse)
library(rstan)
library(cmdstanr)

data <- read_csv("Count_Data_Ben.csv", show_col_types = FALSE)

data <- data %>% group_by(BrainRegion, group, rat_ID) %>%
    summarise(y=sum(Cfos_Counts), area = sum( (1-VOID_PERCENTAGE/100) * (0.94*0.67) )) %>% ungroup()

data <- data %>% mutate(brain_region_idx=as.numeric(factor(BrainRegion)))
data <- data %>% mutate(rat_idx = as.numeric(factor(rat_ID)))
data <- data %>% mutate(group_idx = as.numeric(factor(group)))



check_missing <- data %>% group_by(rat_idx, brain_region_idx) %>% summarise(n=n())

n <- length(unique(data$rat_idx))
d <- length(unique(data$brain_region_idx))

# Create a matrix that contains 1 where observations are available, 0 otherwise.
missing_mtx <- matrix(0, ncol = length(unique(data$brain_region_idx)), nrow=length(unique(data$rat_idx)))
for(i in 1:length(check_missing$rat_idx)){
    missing_mtx[check_missing$rat_idx[i], check_missing$brain_region_idx[i]] <- 1
}

contains_missing <- vector(mode = "integer", length = n)
pad_boundary     <- vector(mode = "integer", length = n)
obs_idxs         <- matrix(0, ncol = length(unique(data$brain_region_idx)), nrow=length(unique(data$rat_idx)))

for(i in 1:n){
    if(length(which(missing_mtx[i,]==0)) > 0){
        contains_missing[i] <-1
    }
    else
        contains_missing[i] <- 0

    avail_idxs <- which(missing_mtx[i,]!=0)
    pad_boundary[i] <- length(avail_idxs)

    for(j in 1:length(avail_idxs)){
        obs_idxs[i,j] <- avail_idxs[j]
    }
}

stan_data = list(A = n,
                 R = 23,
                 G = 4,
                 N = length(data$y),
                 animal_idx = data$rat_idx,
                 region_idx = data$brain_region_idx,
                 group_idx  = data$group_idx,
                 animal_group_idx = (data %>% group_by(rat_idx, group_idx) %>% summarise())$group_idx,
                 E          = log(data$area),
                 y          = data$y,
                 contains_missing = contains_missing,
                 obs_idxs         = obs_idxs,
                 pad_boundary     = pad_boundary,
                 nu               = 1.5
                )

mod <- cmdstan_model("model_poiss_t.stan")

fit <- mod$sample(data            = stan_data,
                  chains          = 4,
                  parallel_chains = 4,
                  iter_warmup     = 4000,
                  iter_sampling   = 4000)

fit <- read_stan_csv(fit$output_files())

# Save
saveRDS(fit, "fit_t.rds")
