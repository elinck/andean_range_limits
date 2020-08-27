# load libraries
library(cowplot)
library(tidyverse)
library(brms)
library(ape)
library(phytools)
library(MCMCglmm)
library(loo)
library(rstan)
library(tidybayes)
library(sjstats)

# run stan faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

### slope data prep

# load tree
slope.tree <- ape::read.tree("~/Dropbox/andean_range_limits/data/blood_slope.tre")

# load data
slope_df  <- read.csv("~/Dropbox/andean_range_limits/data/blood_slopes.csv")

# make variables same order of magnitude; standardardize
slope_df$slope_hb <- slope_df$slope_hb * 1e4
slope_df$error_hb <- slope_df$error_hb * 1e4
slope_df$slope_hct <- slope_df$slope_hct * 1e5
slope_df$error_hct <- slope_df$error_hct * 1e5
slope_df$slope_mchc <- slope_df$slope_mchc * 1e4
slope_df$error_mchc <- slope_df$error_mchc * 1e4
slope_df$elev_range <- ((slope_df$elev_range - mean(slope_df$elev_range)) / sd(slope_df$elev_range))
slope_df$median_elevation <- ((slope_df$median_elevation - mean(slope_df$median_elevation)) / sd(slope_df$median_elevation))
slope_df$sampling_range <- ((slope_df$sampling_range - mean(slope_df$sampling_range)) / sd(slope_df$sampling_range))
slope_df$mass <- ((slope_df$mass - mean(slope_df$mass)) / sd(slope_df$mass))
slope_df$phylo <- slope_df$species
                        
# check normality
blood_tidy <- slope_df %>% pivot_longer(c(slope_hb, slope_hct, slope_mchc), 
                                         names_to = "key", values_to = "value")

s1 <- ggplot(blood_tidy, aes(x=value)) + 
  theme_bw() +
  facet_wrap(~key, scales="free_x") +
  geom_histogram() +
  theme(panel.grid = element_blank())

s2 <- ggplot(blood_tidy, aes(sample=value)) + 
  theme_bw() +
  facet_wrap(~key, scales="free_y") +
  stat_qq() +
  theme(panel.grid = element_blank())

plot_grid(s1,s2,nrow=2) # distribution has heavy tails

# get covariance matrix

A <- ape::vcv.phylo(slope.tree)

### slope models, hemoglobin

# full model
slope_full_hb <- brm(
  formula = bf(slope_hb  ~ 1 + elev_range + sampling_range  + mass + median_elevation +
                 (1 | gr(phylo, cov=A))),
  data = slope_df, 
  family = student(), 
  data2 = list(A = A),
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  prior = c(
    prior(normal(0, 10), "b", coef="elev_range"),
    prior(normal(0, 10), "b", coef="sampling_range"),
    prior(normal(0, 10), "b", coef="mass"),
    prior(normal(0, 10), "b", coef="median_elevation"),
    prior(normal(0, 10), "Intercept"),
    prior(student_t(3, 0, 2), "sd"),
    prior(student_t(3, 0, 2), "sigma")
  )
)

# summarize, visualize
summary(slope_full_hb)
pp_check(slope_full_hb, nsamples = 100)

# simple model
slope_hb <- brm(
  formula = bf(slope_hb  ~ 1 + elev_range + median_elevation + (1 | gr(phylo, cov=A))),
  data = slope_df, 
  family = student(), 
  data2 = list(A = A),
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  prior = c(
    prior(normal(0, 10), "b", coef="elev_range"),
    prior(normal(0, 10), "b", coef="median_elevation"),
    prior(normal(0, 10), "Intercept"),
    prior(student_t(3, 0, 2), "sd"),
    prior(student_t(3, 0, 2), "sigma")
  )
)

summary(slope_hb)
pp_check(slope_hb, nsamples = 100)

# null model, phylogeny only 
slope_null_hb <- brm(
  slope_hb  ~ 0 + (1 | gr(phylo, cov=A)),
  data = slope_df, 
  family = student(), 
  data2 = list(A = A),
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  prior = c(
    prior(student_t(3, 0, 2), "sd"),
    prior(student_t(3, 0, 2), "sigma")
  )
)

# summarize
summary(slope_null_hb)
pp_check(slope_null_hb, nsamples = 100)

# compare
loo_slope_hb <- loo(slope_hb, slope_full_hb, slope_null_hb)  # simple model best, then full, then null

# export full model data for plotting
slope_full_hb %>%
  gather_draws(b_Intercept, b_elev_range, b_median_elevation, b_mass, b_sampling_range) %>%
  write_csv("~/Dropbox/andean_range_limits/data/slope_full_hb.csv")
  
# export looic
write.csv(as.data.frame(loo_slope_hb$diffs), "~/Dropbox/andean_range_limits/data/slope_full_hb_loo_elpd.csv")
write.csv(as.data.frame(loo_slope_hb$ic_diffs__), "~/Dropbox/andean_range_limits/data/slope_full_hb_loo_ic.csv")

### slope models, hct

# full model
slope_full_hct <- brm(
  formula = bf(slope_hct  ~ 1 + elev_range + sampling_range  + mass + median_elevation +
                 (1 | gr(phylo, cov=A))),
  data = slope_df, 
  family = student(), 
  data2 = list(A = A),
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  prior = c(
    prior(normal(0, 10), "b", coef="elev_range"),
    prior(normal(0, 10), "b", coef="sampling_range"),
    prior(normal(0, 10), "b", coef="mass"),
    prior(normal(0, 10), "b", coef="median_elevation"),
    prior(normal(0, 10), "Intercept"),
    prior(student_t(3, 0, 2), "sd"),
    prior(student_t(3, 0, 2), "sigma")
  )
)

# summarize, visualize
summary(slope_full_hct)
pp_check(slope_full_hct, nsamples = 100)

# null model, phylogeny only 
slope_null_hct <- brm(
  slope_hct  ~ 0 + (1 | gr(phylo, cov=A)),
  data = slope_df, 
  family = student(), 
  data2 = list(A = A),
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  prior = c(
    prior(student_t(3, 0, 2), "sd"),
    prior(student_t(3, 0, 2), "sigma")
  )
)

# summarize
summary(slope_null_hct)
pp_check(slope_null_hct, nsamples = 100)

# compare
loo_slope_hct <- loo(slope_full_hct, slope_null_hct) # null model better, but not significantly

# export full model data for plotting
slope_full_hct %>%
  gather_draws(b_Intercept, b_elev_range, b_median_elevation, b_mass, b_sampling_range) %>%
  write_csv("~/Dropbox/andean_range_limits/data/slope_full_hct.csv")

# export looic
write.csv(as.data.frame(loo_slope_hct$diffs), "~/Dropbox/andean_range_limits/data/slope_full_hct_loo_elpd.csv")
write.csv(as.data.frame(loo_slope_hct$ic_diffs__), "~/Dropbox/andean_range_limits/data/slope_full_hct_loo_ic.csv")

### slope models, mchc

# full model
slope_full_mchc <- brm(
  formula = bf(slope_mchc  ~ 1 + elev_range + sampling_range  + mass + median_elevation +
                 (1 | gr(phylo, cov=A))),
  data = slope_df, 
  family = student(), 
  data2 = list(A = A),
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  prior = c(
    prior(normal(0, 10), "b", coef="elev_range"),
    prior(normal(0, 10), "b", coef="median_elevation"),
    prior(normal(0, 10), "Intercept"),
    prior(student_t(3, 0, 2), "sd"),
    prior(student_t(3, 0, 2), "sigma")
  )
)

# summarize, visualize
summary(slope_full_mchc)
pp_check(slope_full_mchc, nsamples = 100)

# null model, phylogeny only 
slope_null_mchc <- brm(
  slope_mchc  ~ 0 + (1 | gr(phylo, cov=A)),
  data = slope_df, 
  family = student(), 
  data2 = list(A = A),
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  sample_prior = TRUE
)

# summarize
summary(slope_null_mchc)
plot(slope_null_mchc)
pp_check(slope_null_mchc, nsamples = 100)
pairs(slope_null_mchc)

# compare
loo_slope_mchc <- loo(slope_full_mchc, slope_null_mchc)

# export full model data for plotting
slope_full_mchc %>%
  gather_draws(b_Intercept, b_elev_range, b_median_elevation, b_mass, b_sampling_range) %>%
  write_csv("~/Dropbox/andean_range_limits/data/slope_full_mchc.csv")

# export looic
write.csv(as.data.frame(loo_slope_mchc$diffs), "~/Dropbox/andean_range_limits/data/slope_full_mchc_loo_elpd.csv")
write.csv(as.data.frame(loo_slope_mchc$ic_diffs__), "~/Dropbox/andean_range_limits/data/slope_full_mchc_loo_ic.csv")

### variance data prep

# load tree
variance.tree <- ape::read.tree("~/Dropbox/andean_range_limits/data/blood_variances.tre")

# load data
variance_df  <- read.csv("~/Dropbox/andean_range_limits/data/blood_variances.csv")

# make variables same order of magnitude; standardardize
variance_df$variance_hb <- variance_df$variance_hb  * 1e1
variance_df$variance_hct <- variance_df$variance_hct  * 1e1
variance_df$variance_mchc <- variance_df$variance_mchc  * 1e1

variance_df$elev_range <- ((variance_df$elev_range - mean(variance_df$elev_range)) / sd(variance_df$elev_range))
variance_df$median_elevation <- ((variance_df$median_elevation - mean(variance_df$median_elevation)) / sd(variance_df$median_elevation))
variance_df$mass <- ((variance_df$mass - mean(variance_df$mass)) / sd(variance_df$mass))
variance_df$range_position <- ((variance_df$range_position - mean(variance_df$range_position)) / sd(variance_df$range_position))
variance_df$edge_distance <- ((variance_df$edge_distance - mean(variance_df$edge_distance)) / sd(variance_df$edge_distance))
variance_df$bin_elevation <- ((variance_df$bin_elevation - mean(variance_df$bin_elevation)) / sd(variance_df$bin_elevation))
variance_df$phylo <- variance_df$species

# check normality
blood_tidy <- variance_df %>% pivot_longer(c(variance_hb, variance_hct, variance_mchc), 
                                        names_to = "key", values_to = "value")

s1 <- ggplot(blood_tidy, aes(x=value)) + 
  theme_bw() +
  facet_wrap(~key, scales="free_x") +
  geom_histogram() +
  theme(panel.grid = element_blank())

s2 <- ggplot(blood_tidy, aes(sample=value)) + 
  theme_bw() +
  facet_wrap(~key, scales="free_y") +
  stat_qq() +
  theme(panel.grid = element_blank())

plot_grid(s1,s2,nrow=2) # distribution has left skew; tails also too fat

# drop outlier8
variance_df <- variance_df[variance_df$variance_hb<1.5,]
variance_df <- variance_df[variance_df$variance_hct<1.5,]
variance_df <- variance_df[variance_df$variance_mchc<1.1,]

# get covariance matrix
B <- ape::vcv.phylo(variance.tree)

### variance models, hemoglobin

# full model
variance_full_hb <- brm(
  formula = bf(variance_hb ~ 1 + range_position + edge_distance + bin_elevation +
                 (1 | gr(phylo, cov=B))),
  data = variance_df, 
  family = student(),
  data2 = list(B = B),
  iter = 10000,
  control = list(adapt_delta = 0.99),
  prior = c(
    prior(normal(0, 2), "b", coef="edge_distance"),
    prior(normal(0, 2), "b", coef="bin_elevation"),
    prior(normal(0, 2), "b", coef="range_position"),
    prior(normal(0, 2), "Intercept"),
    prior(student_t(3, 0, 2), "sd"),
    prior(student_t(3, 0, 2), "sigma")
  )
)

summary(variance_full_hb)
pp_check(variance_full_hb, nsamples = 100)

# simple model
variance_hb <- brm(
  formula = bf(variance_hb ~ edge_distance + (1 | gr(phylo, cov=B))),
  data = variance_df, 
  family = student(),
  data2 = list(B = B),
  iter = 10000,
  control = list(adapt_delta = 0.99),
  save_all_pars = TRUE,
  prior = c(
    prior(normal(0, 2), "b", coef="edge_distance"),
    prior(normal(0, 2), "Intercept"),
    prior(student_t(3, 0, 2), "sd"),
    prior(student_t(3, 0, 2), "sigma")
  )
)

summary(variance_hb)
pp_check(variance_hb, nsamples = 100)

# null model
variance_null_hb <- brm(
  formula = bf(variance_hb ~ 0 + (1 | gr(phylo, cov=B))),
  data = variance_df, 
  family = student(),
  data2 = list(B = B),
  iter = 10000,
  save_all_pars = TRUE,
  prior = c(
    prior(student_t(3, 0, 2), "sd"),
    prior(student_t(3, 0, 2), "sigma")
  )
)

summary(variance_null_hb)
pp_check(variance_null_hb, nsamples = 100)

# compare
loo_variance_hb <- loo(variance_hb, variance_full_hb, variance_null_hb)

# export full model data for plotting
variance_full_hb %>%
  gather_draws(b_Intercept, b_range_position, b_edge_distance, b_bin_elevation) %>%
  write_csv("~/Dropbox/andean_range_limits/data/variance_full_hb.csv")

# export looic
write.csv(as.data.frame(loo_variance_hb$diffs), "~/Dropbox/andean_range_limits/data/variance_full_hb_loo_elpd.csv")
write.csv(as.data.frame(loo_variance_hb$ic_diffs__), "~/Dropbox/andean_range_limits/data/variance_full_hb_loo_ic.csv")

### variance models, hct

# full model
variance_full_hct <- brm(
  formula = bf(variance_hct ~ 1 + range_position + edge_distance + bin_elevation +
                 (1 | gr(phylo, cov=B))),
  data = variance_df, 
  family = student(),
  data2 = list(B = B),
  iter = 10000,
  control = list(adapt_delta = 0.99),
  prior = c(
    prior(normal(0, 2), "b", coef="edge_distance"),
    prior(normal(0, 2), "b", coef="bin_elevation"),
    prior(normal(0, 2), "b", coef="range_position"),
    prior(normal(0, 2), "Intercept"),
    prior(student_t(3, 0, 2), "sd"),
    prior(student_t(3, 0, 2), "sigma")
  )
)

summary(variance_full_hct)
pp_check(variance_full_hct, nsamples = 100)

# null model
variance_null_hct <- brm(
  formula = bf(variance_hct ~ 0 + (1 | gr(phylo, cov=B))),
  data = variance_df, 
  family = student(),
  data2 = list(B = B),
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  prior = c(
    prior(student_t(3, 0, 2), "sd"),
    prior(student_t(3, 0, 2), "sigma")
  )
)


summary(variance_null_hct)
pp_check(variance_null_hct, nsamples = 100)

# compare
loo_variance_hct <- loo(variance_full_hct, variance_null_hct)

# export full model data for plotting
variance_full_hct %>%
  gather_draws(b_Intercept, b_range_position, b_edge_distance, b_bin_elevation) %>%
  write_csv("~/Dropbox/andean_range_limits/data/variance_full_hct.csv")

# export looic
write.csv(as.data.frame(loo_variance_hct$diffs), "~/Dropbox/andean_range_limits/data/variance_full_hct_loo_elpd.csv")
write.csv(as.data.frame(loo_variance_hct$ic_diffs__), "~/Dropbox/andean_range_limits/data/variance_full_hct_loo_ic.csv")

### variance models, mchc

# full model
variance_full_mchc <- brm(
  formula = bf(variance_mchc ~ 1 + range_position + edge_distance + bin_elevation +
                 (1 | gr(phylo, cov=B))),
  data = variance_df, 
  family = student(),
  data2 = list(B = B),
  iter = 10000,
  control = list(adapt_delta = 0.995),
  prior = c(
    prior(normal(0, 1), "b", coef="edge_distance"),
    prior(normal(0, 1), "b", coef="bin_elevation"),
    prior(normal(0, 1), "b", coef="range_position"),
    prior(normal(0, 1), "Intercept"),
    prior(student_t(2, 0, 1), "sd"),
    prior(student_t(2, 0, 1), "sigma")
  )
)

summary(variance_full_mchc)
pp_check(variance_full_mchc, nsamples = 100)

# null model
variance_null_mchc <- brm(
  formula = bf(variance_mchc ~ 0 + (1 | gr(phylo, cov=B))),
  data = variance_df, 
  family = student(),
  data2 = list(B = B),
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  prior = c(
    prior(student_t(2, 0, 1), "sd"),
    prior(student_t(2, 0, 1), "sigma")
  )
)

summary(variance_null_mchc)
pp_check(variance_null_mchc, nsamples = 100)

# compare
loo_variance_mchc <- loo(variance_full_mchc, variance_null_mchc)

# export full model data for plotting
variance_full_mchc %>%
  gather_draws(b_Intercept, b_range_position, b_edge_distance, b_bin_elevation) %>%
  write_csv("~/Dropbox/andean_range_limits/data/variance_full_mchc.csv")

# export looic
write.csv(as.data.frame(loo_variance_mchc$diffs), "~/Dropbox/andean_range_limits/data/variance_full_mchc_loo_elpd.csv")
write.csv(as.data.frame(loo_variance_mchc$ic_diffs__), "~/Dropbox/andean_range_limits/data/variance_full_mchc_loo_ic.csv")
