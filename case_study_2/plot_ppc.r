library(bayesplot)
library(tidyverse)
library(rstan)

# Load the data
data <- read_csv("data.csv")

# Pivot wider into A x R matrix
data <- data %>% pivot_wider(id_cols     = c(animal_idx, group_idx),
                             values_from = counts,
                             names_from  = regions,
                             names_sort  = T)


# Load the fit
fit <- readRDS("fits/fit_hs.rds")
y_rep <- extract(fit, "y_rep")$y_rep
dim(y_rep) <- c(dim(y_rep)[1], prod(dim(y_rep)[-1])) # flatten 3rd dim

# plotting
color_scheme_set("purple")
bayesplot_theme_update(text= element_text(size = 10))

# Proportion of zeroes
p <- ppc_stat(c(as.matrix(data[,-c(1:2)])), y_rep, stat=function(y) mean(y==0)) +
        theme_classic() +
        coord_cartesian(xlim = c(0, 0.1)) +
        xlab("proportion of zeros") +
        ylab("frequency")

ggsave(plot=p, "ppc_zeros.png", dpi=600, units="in", width=5.2, height=2)

# Standard deviation
p <- ppc_stat(c(as.matrix(data[,-c(1:2)])), y_rep, stat=sd) +
        theme_classic() +
        xlab("standard deviation") +
        ylab("frequency")

ggsave(plot=p, "ppc_sd.png", dpi=600, units="in", width=5.2, height=2)
