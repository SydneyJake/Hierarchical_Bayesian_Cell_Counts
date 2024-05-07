library(cmdstanr)
library(tidyverse)
library(rstan)

# What model to fit
mod <- cmdstan_model("models/model_hs.stan")

# Tau values
taus <- seq(from=log(1.001), to=log(1.2), by=0.001)
N <- length(taus)

# Storing the results
gamma_means <- matrix(0, nrow=N, ncol=6)
theta_means <- matrix(0, nrow=N, ncol=2)
theta_sd    <- matrix(0, nrow=N, ncol=2)

# Use real data from a single brain region. This contain zeros.
data <- list("y" = c(770, 820, 712, 541, 0, 0))

# Computation
for(i in 1:N){
  stan_data <- list("y"          = data$y,
                    "N"          = 6,
                    "tau"        = taus[i])

  fit <- mod$sample(data            = stan_data,
                    chains          = 4,
                    parallel_chains = 4,
                    iter_warmup     = 1000,
                    iter_sampling   = 1000)
  # Save fit
  fit <- read_stan_csv(fit$output_files()) # Save output from cmdstanr in a way that preserves param layout.

  gamma_means[i,] <- colMeans(exp(extract(fit, "gamma")$gamma))
  theta_means[i,] <- mean(exp(extract(fit, "theta")$theta))
  theta_sd[i,]    <- sd(exp(extract(fit, "theta")$theta))
}

# Gather results for plotting
gamma_df <- data.frame("mean"       = c(gamma_means),
                       "animal_idx" = rep(c("1":"6"), each=N),
                       "tau"        = rep(taus, 6),
                       "alpha"      = rep(seq(1, 0.17, length.out=N), 6))

theta_df <- data.frame("mean" = c(theta_means),
                       "sd"   = c(theta_sd),
                       "tau"  = rep(taus, 6))

# Plot
p <- ggplot() +
        # theta
        geom_ribbon(data=theta_df, aes(x=tau,y=mean, ymin=mean-sd, ymax=mean+sd), alpha=0.33) +
        geom_line(data=theta_df, aes(x=tau,y=mean), color="black", linewidth=0.5) +
        # gamma
        geom_line(data=gamma_df, aes(x=tau, y=mean, group=animal_idx, alpha=alpha), color="indianred") +
        # data
        geom_point(aes(x=tail(taus, 1), y=data$y), fill=alpha("indianred", 0.17), shape=21, size=0.5) +
        # theme
        theme_classic() +
        ylab(expression(exp(gamma))) +
        xlab(expression(exp(tau))) +
        theme(legend.position="none",
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8))

ggsave(plot = p, filename = "hs_pooling.png",  dpi=600, units = "cm", width = 3.7, height=3.7*0.9)
