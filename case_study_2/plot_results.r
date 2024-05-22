library(tidyverse)
library(rstan)
library(HDInterval)

calc_hdi <- function(samples, model){
    hdis  <- apply(samples, FUN=function(x) hdi(x, 0.95), MARGIN=2)
    means  <- apply(samples, FUN=mean, MARGIN=2)
    data.frame("mean"    = means,
               "lower"   = hdis[1,],
               "upper"   = hdis[2,],
               "region"  = levels(as.factor(data$regions)),
               "model"   = model)
}

run_t_test <- function(df, g1, g2){
    df_g1 <- filter(df, group==g1)
    df_g2 <- filter(df, group==g2)
    test_res <- t.test(x=df_g1$counts, y=df_g2$counts,var.equal = F, paired=F, conf.level = 0.95, alternative = "two.sided")
    data.frame("mean"      = -diff(test_res$estimate) * log2(exp(1)),
               "p"         = test_res$p.value,
               "lower"     = test_res$conf.int[1] * log2(exp(1)),
               "upper"     = test_res$conf.int[2] * log2(exp(1)),
               "df"        = test_res$parameter,
               "model"     = "t-test CI" )
}

# read data
data <- read_csv("data.csv")
data$regions    <- factor(data$regions)
data$group      <- factor(data$group)
data$animal_idx <- factor(data$animal_idx)

# Bayesian fit
fit_hs      <- readRDS("fits/fit_hs.rds")
fit_inflate <- readRDS("fits/fit_inflate.rds")

mus        <- extract(fit_hs, "theta")$theta
mean_diffs_hs <- (mus[,1,] - mus[,2,] ) * log2(exp(1))

mus        <- extract(fit_inflate, "theta")$theta
mean_diffs_inflate <- (mus[,1,] - mus[,2,] ) * log2(exp(1))

bayesian_res <- bind_rows(calc_hdi(mean_diffs_hs, "hs"), calc_hdi(mean_diffs_inflate, "inflate"))

# t test
data$counts[data$counts == 0] <- 1 # add zeroes, cant log otherwise
data$counts <- log(data$counts)

data$regions <- factor(data$regions)
data$group   <- factor(data$group)

split_data <- split(data, f=data$regions)
freq_res   <- bind_rows(lapply(split_data, FUN=function(x) run_t_test(x,"het", "ko")),.id = "region")

# Plotting
res_df <- bind_rows(freq_res, bayesian_res)
res_df$p[is.na(res_df$p)] <- 1
res_df$region <- reorder(factor(res_df$region), res_df$p, decreasing = T)

# The x-axis is too long to be easily visualised in a few inches. Split it over two facets
res_1 <- res_df %>% filter(region %in% levels(region)[1:25])
res_2 <- res_df %>% filter(region %in% levels(region)[26:50])

theme_set(theme_classic(base_size = 10))
theme_update(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

p <- ggplot(data=res_1, aes(x=region, y=mean, group=model, color=model, fill=model)) +
        geom_hline(yintercept=0, alpha=0.33, linewidth=0.25) +
        geom_crossbar(aes(ymin=lower, ymax=upper), alpha=0.33, linewidth=0.25,position = position_dodge2(), width=0.75)+
        scale_color_manual(values = c("#7570b3", "#e72989", "#d95f01"))+
        scale_fill_manual(values = c("#7570b3", "#e72989", "#d95f01"))+
        ylab(expression(log[2](het/ko)))+
        xlab("brain region") +
        coord_cartesian(ylim = c(-9,6)) +
        scale_y_continuous(breaks = seq(-9,6, by =3)) + theme(legend.position = "none")
ggsave(plot = p, filename = "results/results_row1.png", dpi=600, width=17.8, height=17.8/4, units="cm")

p <- ggplot(data=res_2, aes(x=region, y=mean, group=model, color=model, fill=model)) +
        geom_hline(yintercept=0, alpha=0.33, linewidth=0.25) +
        geom_crossbar(aes(ymin=lower, ymax=upper), alpha=0.33, linewidth=0.25,position = position_dodge2(), width=0.75)+
        scale_color_manual(values = c("#7570b3", "#e72989", "#d95f01"))+
        scale_fill_manual(values = c("#7570b3", "#e72989", "#d95f01"))+
        ylab(expression(log[2](het/ko)))+
        xlab("brain region") +
        coord_cartesian(ylim = c(-9,6)) +
        scale_y_continuous(breaks = seq(-9,6, by =3)) + theme(legend.position = "none")
ggsave(plot = p, filename = "results/results_row2.png", dpi=600, width=17.8, height=17.8/4, units="cm")
