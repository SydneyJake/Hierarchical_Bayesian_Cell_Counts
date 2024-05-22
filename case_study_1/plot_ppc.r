library(bayesplot)
library(tidyverse)
library(rstan)

# Load the data
data <- read_csv("data.csv", show_col_types = FALSE)

# counts are summed if there are multiple, same applies to the exposure
data <- data %>%
          group_by(BrainRegion, group, rat_ID) %>%
          summarise(y=sum(Cfos_Counts), area = sum( (1-VOID_PERCENTAGE/100) * (0.94*0.67) )) %>%
          ungroup()

# Load the fit
fit <- readRDS("fits/fit_poiss.rds")
y_rep <- extract(fit, "y_rep")$y_rep

# plotting
color_scheme_set("purple")
bayesplot_theme_update(text= element_text(size = 10))

# Proportion of zeroes
p <- ppc_stat(data$y, y_rep, stat=function(y) mean(y==0), binwidth=0.01) +
        theme_classic() +
        coord_cartesian(xlim = c(0, 0.1)) +
        xlab("proportion of zeros") +
        ylab("frequency")

ggsave(plot=p, "ppc_zeros.png", dpi=600, units="in", width=5.2, height=2)

# Standard deviation
p <- ppc_stat(data$y, y_rep, stat=sd) +
        theme_classic() +
        xlab("standard deviation") +
        ylab("frequency")

ggsave(plot=p, "ppc_sd.png", dpi=600, units="in", width=5.2, height=2)
