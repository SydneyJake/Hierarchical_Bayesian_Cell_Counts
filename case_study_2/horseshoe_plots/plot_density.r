library(tidyverse)

# Horseshoe density using either normal or Cauchy up ot a constant
get_density <- function(gamma, gamma_tilde, kappa, y, norm_param){
    if(norm_param==T){
        dens <- dpois(y, exp(gamma), log = T) +
                dnorm(gamma_tilde, 0, 1, log = T) +
                dnorm(kappa, 0, 1, log = T)
        return(exp(dens))
    }
    else{
        dens <- dpois(y, exp(gamma), log = T) +
                dnorm(gamma_tilde, 0, 1, log = T) +
                dcauchy(kappa, 0, 1, log = T)
        return(exp(dens))
    }
}

data  <- c(770, 820, 712, 541, 0, 0) # Example data from TMd region
theta <- log(mean(data)) # Assume the population mean is the log of the mean of y
tau   <- log(1.05) # Assume the population standard devation is log(1.05)

# Set the grid
gamma_tilde <- seq(-10, 5, by = 0.02)
kappas      <- seq(log(0.01), log(2000), by = 0.02) # needs to exponentiated
param_grid  <- expand.grid("gamma_tilde"=gamma_tilde, "kappa"=kappas)

param_grid$gamma <- param_grid$gamma_tilde*(tau * exp(param_grid$kappa)) + theta

# Now calculate the four densities
param_grid <- param_grid %>%
                  mutate(hs_0  = get_density(gamma, gamma_tilde, exp(kappa), 0, F),
                         hs_y  = get_density(gamma, gamma_tilde, exp(kappa), round(mean(data)), F),
                         mhs_0 = get_density(gamma, gamma_tilde, exp(kappa), 0, T),
                         mhs_y = get_density(gamma, gamma_tilde, exp(kappa), round(mean(data)), T))

# plot the results
theme_set(theme_classic())
theme_update(legend.position="none")

p <- ggplot(data=param_grid, aes(x=gamma_tilde, y=kappa, fill=hs_0)) +
        geom_tile() +
        scale_fill_viridis_c()  +
        ylab(expression(log(kappa))) +
        xlab(expression(tilde(gamma)))

ggsave(plot = p, filename = "Horseshoe_zero.png", dpi=600, units="in", width=5.2/2, height=5.2/2)

p <- ggplot(data=param_grid, aes(x=gamma_tilde, y=kappa, fill=hs_y)) +
        geom_tile() +
        scale_fill_viridis_c() +
        ylab(expression(log(kappa))) +
        xlab(expression(tilde(gamma)))

ggsave(plot = p, filename = "Horseshoe_y.png", dpi=600, units="in", width=5.2/2, height=5.2/2)

p <- ggplot(data=param_grid, aes(x=gamma_tilde, y=kappa, fill=mhs_0)) +
        geom_tile() +
        scale_fill_viridis_c() +
        ylab(expression(log(kappa))) +
        xlab(expression(tilde(gamma)))

ggsave(plot = p, filename = "Modified_horseshoe_zero.png", dpi=600, units="in", width=5.2/2, height=5.2/2)

p <- ggplot(data=param_grid, aes(x=gamma_tilde, y=kappa, fill=mhs_y)) +
        geom_tile() +
        scale_fill_viridis_c() +
        ylab(expression(log(kappa))) +
        xlab(expression(tilde(gamma)))

ggsave(plot = p, filename = "Modified_horseshoe_y.png", dpi=600, units="in", width=5.2/2, height=5.2/2)
