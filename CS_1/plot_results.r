library(tidyverse)
library(rstan)
library(HDInterval)

# Calculate Bayesian HDI
calc_hdi <- function(samples, g1_idx, g2_idx, model){
    # fold changes
    nat_log_ab <- samples[,g1_idx,] - samples[,g2_idx,]
    log_2_fold <- nat_log_ab * log2(exp(1))
    hdis  <- apply(log_2_fold, FUN=function(x) hdi(x, 0.95), MARGIN=2)
    means <- apply(log_2_fold, FUN=mean, MARGIN=2)

    data.frame("mean"    = means,
               "lower"   = hdis[1,],
               "upper"   = hdis[2,],
               "region"  = levels(factor(data$BrainRegion)),
               "model"   = model)
}

# Confidence intervals
run_t_test <- function(df, g1, g2){
    df_g1 <- filter(df, group==g1)
    df_g2 <- filter(df, group==g2)
    test_res <- t.test(x=df_g1$mu_cpr, y=df_g2$mu_cpr,var.equal = F, paired=F, conf.level = 0.95, alternative = "two.sided")

    data.frame("mean"      = -diff(test_res$estimate),
               "p"         = test_res$p.value,
               "lower"     = test_res$conf.int[1],
               "upper"     = test_res$conf.int[2],
               "df"        = test_res$parameter,
               "model"     = "t-test CI" )
}

data <- read_csv("data.csv")
grp_idx <- list("Lesion_Familiar"=1, "Lesion_Novel"=2, "Sham_Familiar"=3, "Sham_Novel"=4)
fit   <- readRDS("fits/fit_poiss.rds")
theta <- extract(fit, "theta")$theta

bayes_res_SFSN <- calc_hdi(extract(fit, "theta")$theta, grp_idx$Sham_Familiar, grp_idx$Sham_Novel, "Normal HDI")
bayes_res_LFLN <- calc_hdi(extract(fit, "theta")$theta, grp_idx$Lesion_Familiar, grp_idx$Lesion_Novel, "Normal HDI")

# Mean for each animal within each group, this will be input to a t-test
sample_data <- data %>%
                group_by(BrainRegion, group, rat_ID) %>%
                summarise(mu_cpr = log2(mean(counts_per_area))) %>%
                ungroup() %>%
                select(-c(rat_ID))

sample_data$BrainRegion <- factor(sample_data$BrainRegion)
sample_data$group       <- factor(sample_data$group)

# list split by brainregions
sample_data <- split(sample_data, f=sample_data$BrainRegion)
freq_res_SFSN <- bind_rows(lapply(sample_data, FUN=function(x) run_t_test(x,"Sham_Familiar", "Sham_Novel")), .id="region")
freq_res_LFLN <- bind_rows(lapply(sample_data, FUN=function(x) run_t_test(x,"Lesion_Familiar", "Lesion_Novel")), .id="region")

##### Plotting
sham_df <- bind_rows(freq_res_SFSN, bayes_res_SFSN)
sham_df$p[is.na(sham_df$p)] <-1
sham_df$region <- reorder(factor(sham_df$region), sham_df$p, decreasing = T)

lesion_df <- bind_rows(freq_res_LFLN, bayes_res_LFLN)
lesion_df$p[is.na(lesion_df$p)] <-1
lesion_df$region <- reorder(factor(lesion_df$region), lesion_df$p, decreasing = T)

theme_set(theme_classic(base_size = 10))
theme_update(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), legend.position="none")

p <- ggplot(data=sham_df, aes(x=region, y=mean, group=model, color=model, fill=model)) +
        geom_hline(yintercept=0, alpha=0.33, linewidth=0.25) +
        geom_crossbar(aes(ymin=lower, ymax=upper), alpha=0.33, linewidth=0.25,position = position_dodge2(), width=0.75)+
        scale_color_brewer(palette = "Dark2") +
        scale_fill_brewer(palette = "Dark2") +
        ylab(expression(log[2](SF/SN)))+
        xlab("brain region") +
        coord_cartesian(ylim = c(-2,1)) +
        scale_y_continuous(breaks = seq(-2,2, by = 1))
ggsave(plot = p, filename = "results/SFSN.png", dpi=600, width=17.8, height=17.8/4, units="cm")

p <- ggplot(data=lesion_df, aes(x=region, y=mean, group=model, color=model, fill=model)) +
        geom_hline(yintercept=0, alpha=0.33, linewidth=0.25) +
        geom_crossbar(aes(ymin=lower, ymax=upper), alpha=0.33, linewidth=0.25,position = position_dodge2(), width=0.75)+
        scale_color_brewer(palette = "Dark2") +
        scale_fill_brewer(palette = "Dark2") +
        ylab(expression(log[2](LF/LN)))+
        xlab("brain region") +
        coord_cartesian(ylim = c(-2,1)) +
        scale_y_continuous(breaks = seq(-2,2, by = 1))
ggsave(plot = p, filename = "results/LFLN.png", dpi=600, width=17.8, height=17.8/4, units="cm")
