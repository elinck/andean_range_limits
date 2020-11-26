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
library(modelr)

### set up

# run stan faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# load slope tree
slope.tree <- ape::read.tree("~/Dropbox/andean_range_limits/data/blood_slope.tre")

# load slope data
slope_df  <- read.csv("~/Dropbox/andean_range_limits/data/blood_slopes.csv")

# load variance tree
variance.tree <- ape::read.tree("~/Dropbox/andean_range_limits/data/blood_variances.tre")

# load variance data
variance_df  <- read.csv("~/Dropbox/andean_range_limits/data/blood_variances.csv")

### basic distribution models, slope and variance

slope_hb_dist <- 
  brm(data = slope_df, family = student(),
      slope_hb | se(error_hb, sigma=TRUE) ~ 1,
      inits = inits_list,
      prior = c(
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(slope_hb_dist, nsamples = 100) + xlim(-0.025, 0.025)

# trace and density plots
plot(slope_hb_dist)

# pairwise plots
pairs(slope_hb_dist)

# export predictions
predict(slope_hb_dist, summary=FALSE, nsamples = 100) %>% 
  write.csv("~/Dropbox/andean_range_limits/data/slope_hb_dist_draws.csv")


slope_hct_dist <- 
  brm(data = slope_df, family = student(),
      slope_hct | se(error_hct, sigma=TRUE) ~ 1,
      inits = inits_list,
      prior = c(
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(slope_hct_dist, nsamples = 50) + xlim(-0.001, 0.001)

# trace and density plots
plot(slope_hct_dist)

# pairwise plots
pairs(slope_hct_dist)

# export predictions
predict(slope_hct_dist, summary=FALSE, nsamples = 100) %>% 
  write.csv("~/Dropbox/andean_range_limits/data/slope_hct_dist_draws.csv")

slope_mchc_dist <- 
  brm(data = slope_df, family = student(),
      slope_mchc | se(error_mchc, sigma=TRUE) ~ 1,
      inits = inits_list,
      prior = c(
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(slope_mchc_dist, nsamples = 50) + xlim(-0.025, 0.025)

# trace and density plots
plot(slope_mchc_dist)

# pairwise plots
pairs(slope_mchc_dist)

# export predictions
predict(slope_mchc_dist, summary=FALSE, nsamples = 100) %>% 
  write.csv("~/Dropbox/andean_range_limits/data/slope_mchc_dist_draws.csv")

variance_hb_dist <- 
  brm(data = variance_df, family = lognormal(),
      variance_hb ~ 1,
      inits = inits_list,
      prior = c(
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(variance_hb_dist, nsamples = 100) 

# trace and density plots
plot(variance_hb_dist)

# pairwise plots
pairs(variance_hb_dist)

# export predictions
predict(variance_hb_dist, summary=FALSE, nsamples = 100) %>% 
  write.csv("~/Dropbox/andean_range_limits/data/variance_hb_dist_draws.csv")

variance_hct_dist <- 
  brm(data = variance_df, family = lognormal(),
      variance_hct ~ 1,
      inits = inits_list,
      prior = c(
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(variance_hct_dist, nsamples = 50)

# trace and density plots
plot(variance_hct_dist)

# pairwise plots
pairs(variance_hct_dist)

# export predictions
predict(variance_hct_dist, summary=FALSE, nsamples = 100) %>% 
  write.csv("~/Dropbox/andean_range_limits/data/variance_hct_dist_draws.csv")

variance_mchc_dist <- 
  brm(data = variance_df, family = lognormal(),
      variance_mchc ~ 1,
      inits = inits_list,
      prior = c(
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(variance_mchc_dist, nsamples = 50) 

# trace and density plots
plot(variance_mchc_dist)

# pairwise plots
pairs(variance_mchc_dist)

# export predictions
predict(variance_mchc_dist, summary=FALSE, nsamples = 100) %>% 
  write.csv("~/Dropbox/andean_range_limits/data/variance_mchc_dist_draws.csv")


### slope data prep for predictive models

# check order of magnitude of variables
head(slope_df)

# make variables same order of magnitude, standarize
slope_df <- 
  slope_df %>% mutate(slope_hb = slope_hb*1e3,
                    error_hb = error_hb*1e3,
                    slope_hct = slope_hct*1e5,
                    slope_mchc = slope_hct*1e1,
                    elev_range = elev_range*1e-4,
                    median_elevation = median_elevation*1e-4,
                    mass = mass*1e-2)

# standardize predictors
slope_df <- 
  slope_df %>% mutate(elev_range_s = (elev_range - mean(elev_range)) / sd(elev_range),
                      median_elevation_s = (median_elevation - mean(median_elevation)) / sd(median_elevation),
                      sampling_range_s = (sampling_range - mean(sampling_range)) / sd(median_elevation),
                      mass_s = (mass - mean(mass)) / sd(mass))

# add "phylo" variable
slope_df$phylo <- slope_df$species

# get covariance matrix from phylogeny
A <- ape::vcv.phylo(slope.tree)

# write transformed slope data 
write.csv(slope_df, "~/Dropbox/andean_range_limits/data/blood_slopes_m.csv")

### slope models, hemoglobin

# here we specify the initial (i.e., starting) values for models with a predictor
inits      <- list(Yl = slope_df$slope_hb)
inits_list <- list(inits, inits)

# fit full model with phylogeny
slope_hb_1 <- 
  brm(data = slope_df, family = student(),
      slope_hb | se(error_hb, sigma=TRUE) ~ 1 + elev_range_s + 
        median_elevation_s + mass_s + sampling_range_s + elev_range_s*median_elevation_s + (1 | gr(phylo, cov=A)),
      data2 = list(A = A),
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="elev_range_s"),
        prior(normal(0, 2.5), "b", coef="median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="elev_range_s:median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="sampling_range_s"),
        prior(normal(0, 2.5), "b", coef="mass_s"),
        prior(normal(0, 10), "Intercept"),
        prior(student_t(3, 0, 2), "sd"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(slope_hb_1, nsamples = 50) + xlim(-10, 10)

# trace and density plots
plot(slope_hb_1)

# pairwise plots
pairs(slope_hb_1)

# fit full model without phylogeny
slope_hb_2 <- 
  brm(data = slope_df, family = student(),
      slope_hb | se(error_hb, sigma=TRUE) ~ 1 + elev_range_s + 
        median_elevation_s + mass_s + sampling_range_s + elev_range_s*median_elevation_s,
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="elev_range_s"),
        prior(normal(0, 2.5), "b", coef="median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="elev_range_s:median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="sampling_range_s"),
        prior(normal(0, 2.5), "b", coef="mass_s"),
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(slope_hb_2, nsamples = 50) + xlim(-10, 10)

# trace and density plots
plot(slope_hb_2)

# pairwise plots
pairs(slope_hb_2)

# fit full model without phylogeny or interaction term
slope_hb_3 <- 
  brm(data = slope_df, family = student(),
      slope_hb | se(error_hb, sigma=TRUE) ~ 1 + elev_range_s + 
        median_elevation_s + mass_s + sampling_range_s,
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="elev_range_s"),
        prior(normal(0, 2.5), "b", coef="median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="sampling_range_s"),
        prior(normal(0, 2.5), "b", coef="mass_s"),
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(slope_hb_3, nsamples = 50) + xlim(-10, 10)

# trace and density plots
plot(slope_hb_3)

# pairwise plots
pairs(slope_hb_3)

# fit the null model with phylogeny alone
slope_hb_0 <- 
  brm(data = slope_df, family = student(),
      slope_hb | se(error_hb, sigma=TRUE) ~ 0 + (1 | gr(phylo, cov=A)),
      data2 = list(A = A),
      inits = inits_list,
      prior = c(
        prior(student_t(3, 0, 2), "sd"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 5000, warmup = 2000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(slope_hb_0, nsamples = 50) + xlim(-10, 10)

# trace and density plots
plot(slope_hb_0)

# pairwise plots
pairs(slope_hb_0)

# compare models
loo_slope_hb <- loo(slope_hb_1, slope_hb_2, slope_hb_3, slope_hb_0)  # simple model best, then full, then null

# export counterfactual prediction for effect of elevational range
range_seq <- tibble(elev_range_s = seq(from = -3, to = 3, by = 0.1),
                   error_hb = mean(slope_df$error_hb),
                   median_elevation_s = mean(slope_df$median_elevation_s),
                   mass_s = mean(slope_df$mass_s),
                   sampling_range_s = mean(slope_df$sampling_range_s))

fitted(slope_hb_3, newdata=range_seq) %>%
  as_tibble() %>%
  rename(f_ll=Q2.5,
         f_ul=Q97.5) %>%
  bind_cols(
    predict(slope_hb_3,
          newdata = range_seq) %>%
      as_tibble() %>%
      transmute(p_ll = Q2.5,
                p_ul = Q97.5),
    range_seq) %>%
  write.csv("~/Dropbox/andean_range_limits/data/slope_hb_range_counter.csv")

# export counterfactual prediction for effect of median range elevation
median_seq <- tibble(elev_range_s = mean(slope_df$elev_range_s),
                     error_hb = mean(slope_df$error_hb),
                     median_elevation_s = seq(from = -3, to = 3, by = 0.1),
                     mass_s = mean(slope_df$mass_s),
                     sampling_range_s = mean(slope_df$sampling_range_s))

fitted(slope_hb_3, newdata=median_seq) %>%
  as_tibble() %>%
  rename(f_ll=Q2.5,
         f_ul=Q97.5) %>%
  bind_cols(
    predict(slope_hb_3,
            newdata = median_seq) %>%
      as_tibble() %>%
      transmute(p_ll = Q2.5,
                p_ul = Q97.5),
    median_seq) %>%
  write.csv("~/Dropbox/andean_range_limits/data/slope_hb_median_counter.csv")

# export counterfactual prediction for interaction term
tmp = list()
for(i in -1:1){
  interaction_seq <- tibble(elev_range_s = c(i,i,i),
                            error_hb = mean(slope_df$error_hb),
                            median_elevation_s = seq(from = -1, to = 1, by = 1),
                            mass_s = mean(slope_df$mass_s),
                            sampling_range_s = mean(slope_df$sampling_range_s))
  
  # amake predictions, assign to dataframe
  tmp[[i+2]] <- fitted(slope_hb_2, newdata=interaction_seq) %>%
    as_tibble() %>%
    rename(f_ll=Q2.5,
           f_ul=Q97.5) %>%
    bind_cols(
      predict(slope_hb_2,
              newdata = interaction_seq) %>%
        as_tibble() %>%
        transmute(p_ll = Q2.5,
                  p_ul = Q97.5),
      interaction_seq)
  
  # name panel for faceting
  tmp[[i+2]]$panel <- paste0("panel_",i)
}

# assemble and export
do.call(rbind, tmp) %>% write.csv("~/Dropbox/andean_range_limits/data/slope_hb_interaction.csv")

# export draws from sampled posterior
slope_hb_1 %>%
  gather_draws(b_elev_range_s, b_median_elevation_s, b_mass_s, b_sampling_range_s, `b_elev_range_s:median_elevation_s`) %>%
  write.csv("~/Dropbox/andean_range_limits/data/slope_hb_draws_interaction.csv")

slope_hb_3 %>%
  gather_draws(b_elev_range_s, b_median_elevation_s, b_mass_s, b_sampling_range_s) %>%
  write.csv("~/Dropbox/andean_range_limits/data/slope_hb_draws.csv")


# export looic
write.csv(as.data.frame(loo_slope_hb$diffs), "~/Dropbox/andean_range_limits/data/slope_full_hb_loo_elpd.csv")
write.csv(as.data.frame(loo_slope_hb$ic_diffs__), "~/Dropbox/andean_range_limits/data/slope_full_hb_loo_ic.csv")

### slope models, hematocrit

# here we specify the initial (i.e., starting) values
inits      <- list(Yl = slope_df$slope_hct)
inits_list <- list(inits, inits)

# fit full model with phylogeny
slope_hct_1 <- 
  brm(data = slope_df, family = student(),
      slope_hct | se(error_hct, sigma=TRUE) ~ 1 + elev_range_s + 
        median_elevation_s + mass_s + sampling_range_s + elev_range_s*median_elevation_s + (1 | gr(phylo, cov=A)),
      data2 = list(A = A),
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="elev_range_s"),
        prior(normal(0, 2.5), "b", coef="median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="elev_range_s:median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="sampling_range_s"),
        prior(normal(0, 2.5), "b", coef="mass_s"),
        prior(normal(0, 10), "Intercept"),
        prior(student_t(3, 0, 2), "sd"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(slope_hct_1, nsamples = 50) + xlim(-15, 15)

# trace and density plots
plot(slope_hct_1)

# pairwise plots
pairs(slope_hct_1)

# fit full model without phylogeny
slope_hct_2 <- 
  brm(data = slope_df, family = student(),
      slope_hct | se(error_hct, sigma=TRUE) ~ 1 + elev_range_s + 
        median_elevation_s + mass_s + sampling_range_s + elev_range_s*median_elevation_s,
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="elev_range_s"),
        prior(normal(0, 2.5), "b", coef="median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="elev_range_s:median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="sampling_range_s"),
        prior(normal(0, 2.5), "b", coef="mass_s"),
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(slope_hct_2, nsamples = 50) + xlim(-15, 15)

# trace and density plots
plot(slope_hct_2)

# pairwise plots
pairs(slope_hct_2)

# fit full model without phylogeny or interaction term
slope_hct_3 <- 
  brm(data = slope_df, family = student(),
      slope_hct | se(error_hct, sigma=TRUE) ~ 1 + elev_range_s + 
        median_elevation_s + mass_s + sampling_range_s,
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="elev_range_s"),
        prior(normal(0, 2.5), "b", coef="median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="sampling_range_s"),
        prior(normal(0, 2.5), "b", coef="mass_s"),
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(slope_hct_3, nsamples = 50) + xlim(-15, 15)

# trace and density plots
plot(slope_hct_3)

# pairwise plots
pairs(slope_hct_3)

# fit the null model with phylogeny alone
slope_hct_0 <- 
  brm(data = slope_df, family = student(),
      slope_hct | se(error_hct, sigma=TRUE) ~ 0 + (1 | gr(phylo, cov=A)),
      data2 = list(A = A),
      inits = inits_list,
      prior = c(
        prior(student_t(3, 0, 2), "sd"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 5000, warmup = 2000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(slope_hct_0, nsamples = 50) + xlim(-10, 10)

# trace and density plots
plot(slope_hct_0)

# pairwise plots
pairs(slope_hct_0)

# compare models
loo_slope_hct <- loo(slope_hct_1, slope_hct_2, slope_hct_3, slope_hct_0)  # simple model best, then full, then null

# export draws from sampled posterior
slope_hct_1 %>%
  gather_draws(b_elev_range_s, b_median_elevation_s, b_mass_s, b_sampling_range_s,`b_elev_range_s:median_elevation_s`) %>%
  write.csv("~/Dropbox/andean_range_limits/data/slope_hct_draws.csv")

# export looic
write.csv(as.data.frame(loo_slope_hct$diffs), "~/Dropbox/andean_range_limits/data/slope_full_hct_loo_elpd.csv")
write.csv(as.data.frame(loo_slope_hct$ic_diffs__), "~/Dropbox/andean_range_limits/data/slope_full_hct_loo_ic.csv")

### slope models, mchc

# here we specify the initial (i.e., starting) values
inits      <- list(Yl = slope_df$slope_mchc)
inits_list <- list(inits, inits)

# fit full model with phylogeny
slope_mchc_1 <- 
  brm(data = slope_df, family = student(),
      slope_mchc | se(error_mchc, sigma=TRUE) ~ 1 + elev_range_s + 
        median_elevation_s + mass_s + sampling_range_s + elev_range_s*median_elevation_s + (1 | gr(phylo, cov=A)),
      data2 = list(A = A),
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="elev_range_s"),
        prior(normal(0, 2.5), "b", coef="median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="elev_range_s:median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="sampling_range_s"),
        prior(normal(0, 2.5), "b", coef="mass_s"),
        prior(normal(0, 10), "Intercept"),
        prior(student_t(3, 0, 2), "sd"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(slope_mchc_1, nsamples = 50) + xlim(-10, 10)

# trace and density plots
plot(slope_mchc_1)

# pairwise plots
pairs(slope_mchc_1)

# fit full model without phylogeny
slope_mchc_2 <- 
  brm(data = slope_df, family = student(),
      slope_mchc | se(error_mchc, sigma=TRUE) ~ 1 + elev_range_s + 
        median_elevation_s + mass_s + sampling_range_s + elev_range_s*median_elevation_s,
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="elev_range_s"),
        prior(normal(0, 2.5), "b", coef="median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="elev_range_s:median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="sampling_range_s"),
        prior(normal(0, 2.5), "b", coef="mass_s"),
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(slope_mchc_2, nsamples = 50) + xlim(-10, 10)

# trace and density plots
plot(slope_mchc_2)

# pairwise plots
pairs(slope_mchc_2)

# fit full model without phylogeny or interaction term
slope_mchc_3 <- 
  brm(data = slope_df, family = student(),
      slope_mchc | se(error_mchc, sigma=TRUE) ~ 1 + elev_range_s + 
        median_elevation_s + mass_s + sampling_range_s,
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="elev_range_s"),
        prior(normal(0, 2.5), "b", coef="median_elevation_s"),
        prior(normal(0, 2.5), "b", coef="sampling_range_s"),
        prior(normal(0, 2.5), "b", coef="mass_s"),
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(slope_mchc_3, nsamples = 50) + xlim(-10, 10)

# trace and density plots
plot(slope_mchc_3)

# pairwise plots
pairs(slope_mchc_3)

# fit the null model with phylogeny alone
slope_mchc_0 <- 
  brm(data = slope_df, family = student(),
      slope_mchc | se(error_mchc, sigma=TRUE) ~ 0 + (1 | gr(phylo, cov=A)),
      data2 = list(A = A),
      inits = inits_list,
      prior = c(
        prior(student_t(3, 0, 2), "sd"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 5000, warmup = 2000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# compare models
loo_slope_mchc <- loo(slope_mchc_1, slope_mchc_2, slope_mchc_3, slope_mchc_0)  # simple model best, then full, then null

# export draws from sampled posterior
slope_mchc_1 %>%
  gather_draws(b_elev_range_s, b_median_elevation_s, b_mass_s, b_sampling_range_s,`b_elev_range_s:median_elevation_s`) %>%
  write.csv("~/Dropbox/andean_range_limits/data/slope_mchc_draws.csv")

# export looic
write.csv(as.data.frame(loo_slope_mchc$diffs), "~/Dropbox/andean_range_limits/data/slope_full_mchc_loo_elpd.csv")
write.csv(as.data.frame(loo_slope_mchc$ic_diffs__), "~/Dropbox/andean_range_limits/data/slope_full_mchc_loo_ic.csv")

### variance data prep

# check order of magnitude of variables
head(variance_df)

# make variables same order of magnitude
variance_df <- 
  variance_df %>% mutate(variance_hb = variance_hb*1e1,
                         variance_hct = variance_hct*1e1,
                         variance_mchc = variance_mchc*1e1)

# standardize predictors
variance_df <- 
  variance_df %>% mutate(bin_elevation_s = (bin_elevation - mean(bin_elevation)) / sd(bin_elevation),
                         edge_distance_s = (edge_distance - mean(edge_distance)) / sd(edge_distance))

# add "phylo" variable
variance_df$phylo <- variance_df$species

# get covariance matrix
A <- ape::vcv.phylo(variance.tree)

# write transformed variance
write.csv(variance_df, "~/Dropbox/andean_range_limits/data/blood_variance_m.csv")

### variance models, hemoglobin

# here we specify the initial (i.e., starting) values
inits      <- list(Yl = variance_df$variance_hb)
inits_list <- list(inits, inits)

# fit full model with phylogeny
variance_hb_1 <- 
  brm(data = variance_df, family = lognormal(),
      variance_hb ~ 1 + bin_elevation_s + edge_distance_s + bin_elevation_s*edge_distance_s + (1 | gr(phylo, cov=A)),
      data2 = list(A = A),
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="bin_elevation_s"),
        prior(normal(0, 2.5), "b", coef="edge_distance_s"),
        prior(normal(0, 10), "Intercept"),
        prior(student_t(3, 0, 2), "sd"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(variance_hb_1, nsamples = 50) + xlim(0, 5)

# trace and density plots
plot(variance_hb_1)

# pairwise plots
pairs(variance_hb_1)

# fit full model without phylogeny
variance_hb_2 <- 
  brm(data = variance_df, family = lognormal(),
      variance_hb ~ 1 + bin_elevation_s + edge_distance_s + bin_elevation_s*edge_distance_s,
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="bin_elevation_s"),
        prior(normal(0, 2.5), "b", coef="edge_distance_s"),
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(variance_hb_2, nsamples = 50) + xlim(0, 5)

# trace and density plots
plot(variance_hb_2)

# pairwise plots
pairs(variance_hb_2)

# fit full model without phylogeny or interaction term
variance_hb_3 <- 
  brm(data = variance_df, family = lognormal(),
      variance_hb ~ 1 + bin_elevation_s + edge_distance_s,
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="bin_elevation_s"),
        prior(normal(0, 2.5), "b", coef="edge_distance_s"),
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(variance_hb_3, nsamples = 50) + xlim(0, 5)

# trace and density plots
plot(variance_hb_3)

# pairwise plots
pairs(variance_hb_3)

# fit the null model with phylogeny alone
variance_hb_0 <- 
  brm(data = variance_df, family = lognormal(),
      variance_hb ~ 0 + (1 | gr(phylo, cov=A)),
      data2 = list(A = A),
      inits = inits_list,
      prior = c(
        prior(student_t(3, 0, 2), "sd"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# compare models
loo_variance_hb <- loo(variance_hb_1, variance_hb_2, variance_hb_3, variance_hb_0) 

# export counterfactual prediction for effect of elevational range
bin_seq <- tibble(bin_elevation_s = seq(from = -3, to = 3, by = 0.1),
                    edge_distance_s = mean(variance_df$edge_distance_s))

fitted(variance_hb_2, newdata=bin_seq) %>%
  as_tibble() %>%
  rename(f_ll=Q2.5,
         f_ul=Q97.5) %>%
  bind_cols(
    predict(variance_hb_2,
            newdata = bin_seq) %>%
      as_tibble() %>%
      transmute(p_ll = Q2.5,
                p_ul = Q97.5),
    bin_seq) %>%
  write.csv("~/Dropbox/andean_range_limits/data/variance_hb_bin_counter.csv")

# export counterfactual prediction for effect of edge effects
edge_seq <- tibble(bin_elevation_s = mean(variance_df$bin_elevation_s),
                  edge_distance_s = seq(from = -3, to = 3, by = 0.1))

fitted(variance_hb_2, newdata=edge_seq) %>%
  as_tibble() %>%
  rename(f_ll=Q2.5,
         f_ul=Q97.5) %>%
  bind_cols(
    predict(variance_hb_2,
            newdata = edge_seq) %>%
      as_tibble() %>%
      transmute(p_ll = Q2.5,
                p_ul = Q97.5),
    edge_seq) %>%
  write.csv("~/Dropbox/andean_range_limits/data/variance_hb_edge_counter.csv")

# export counterfactual prediction for interaction term
tmp = list()
for(i in -1:1){
  interaction_seq <- tibble(bin_elevation_s = c(i,i,i),
                            edge_distance_s = seq(from = -1, to = 1, by = 1))
  
  # amake predictions, assign to dataframe
  tmp[[i+2]] <- fitted(variance_hb_2, newdata=interaction_seq) %>%
    as_tibble() %>%
    rename(f_ll=Q2.5,
           f_ul=Q97.5) %>%
    bind_cols(
      predict(variance_hb_2,
              newdata = interaction_seq) %>%
        as_tibble() %>%
        transmute(p_ll = Q2.5,
                  p_ul = Q97.5),
      interaction_seq)
  
  # name panel for faceting
  tmp[[i+2]]$panel <- paste0("panel_",i)
}

# assemble and export
do.call(rbind, tmp) %>% write.csv("~/Dropbox/andean_range_limits/data/variance_hb_interaction.csv")

# export draws from sampled posterior
variance_hb_1 %>%
  gather_draws(b_bin_elevation_s, b_edge_distance_s, `b_bin_elevation_s:edge_distance_s`) %>%
  write.csv("~/Dropbox/andean_range_limits/data/variance_hb_draws.csv")

# export looic
write.csv(as.data.frame(loo_variance_hb$diffs), "~/Dropbox/andean_range_limits/data/variance_full_hb_loo_elpd.csv")
write.csv(as.data.frame(loo_variance_hb$ic_diffs__), "~/Dropbox/andean_range_limits/data/variance_full_hb_loo_ic.csv")

### variance models, hematocrit

# here we specify the initial (i.e., starting) values
inits      <- list(Yl = variance_df$variance_hct)
inits_list <- list(inits, inits)

# fit full model with phylogeny
variance_hct_1 <- 
  brm(data = variance_df, family = lognormal(),
      variance_hct ~ 1 + bin_elevation_s + edge_distance_s + bin_elevation_s*edge_distance_s + (1 | gr(phylo, cov=A)),
      data2 = list(A = A),
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="bin_elevation_s"),
        prior(normal(0, 2.5), "b", coef="edge_distance_s"),
        prior(normal(0, 10), "Intercept"),
        prior(student_t(3, 0, 2), "sd"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(variance_hct_1, nsamples = 50) + xlim(0, 5)

# trace and density plots
plot(variance_hct_1)

# pairwise plots
pairs(variance_hct_1)

# fit full model without phylogeny
variance_hct_2 <- 
  brm(data = variance_df, family = lognormal(),
      variance_hct ~ 1 + bin_elevation_s + edge_distance_s + bin_elevation_s*edge_distance_s,
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="bin_elevation_s"),
        prior(normal(0, 2.5), "b", coef="edge_distance_s"),
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(variance_hct_2, nsamples = 50) + xlim(0, 5)

# trace and density plots
plot(variance_hct_2)

# pairwise plots
pairs(variance_hct_2)

# fit full model without phylogeny or interaction term
variance_hct_3 <- 
  brm(data = variance_df, family = lognormal(),
      variance_hct ~ 1 + bin_elevation_s + edge_distance_s,
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="bin_elevation_s"),
        prior(normal(0, 2.5), "b", coef="edge_distance_s"),
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(variance_hct_3, nsamples = 50) + xlim(0, 5)

# trace and density plots
plot(variance_hct_3)

# pairwise plots
pairs(variance_hct_3)

# fit the null model with phylogeny alone
variance_hct_0 <- 
  brm(data = variance_df, family = lognormal(),
      variance_hct ~ 0 + (1 | gr(phylo, cov=A)),
      data2 = list(A = A),
      inits = inits_list,
      prior = c(
        prior(student_t(3, 0, 2), "sd"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(variance_hct_0, nsamples = 50) + xlim(0, 5)

# trace and density plots
plot(variance_hct_0)

# pairwise plots
pairs(variance_hct_0)

# compare models
loo_variance_hct <- loo(variance_hct_1, variance_hct_2, variance_hct_3, variance_hct_0) 

# export counterfactual prediction for effect of elevational range
bin_seq <- tibble(bin_elevation_s = seq(from = -3, to = 3, by = 0.1),
                  edge_distance_s = mean(variance_df$edge_distance_s))

fitted(variance_hct_2, newdata=bin_seq) %>%
  as_tibble() %>%
  rename(f_ll=Q2.5,
         f_ul=Q97.5) %>%
  bind_cols(
    predict(variance_hct_2,
            newdata = bin_seq) %>%
      as_tibble() %>%
      transmute(p_ll = Q2.5,
                p_ul = Q97.5),
    bin_seq) %>%
  write.csv("~/Dropbox/andean_range_limits/data/variance_hct_bin_counter.csv")

# export counterfactual prediction for effect of edge effects
edge_seq <- tibble(bin_elevation_s = mean(variance_df$bin_elevation_s),
                   edge_distance_s = seq(from = -3, to = 3, by = 0.1))

fitted(variance_hct_2, newdata=edge_seq) %>%
  as_tibble() %>%
  rename(f_ll=Q2.5,
         f_ul=Q97.5) %>%
  bind_cols(
    predict(variance_hct_2,
            newdata = edge_seq) %>%
      as_tibble() %>%
      transmute(p_ll = Q2.5,
                p_ul = Q97.5),
    edge_seq) %>%
  write.csv("~/Dropbox/andean_range_limits/data/variance_hct_edge_counter.csv")

# export counterfactual prediction for interaction term
tmp = list()
for(i in -1:1){
  interaction_seq <- tibble(bin_elevation_s = c(i,i,i),
                            edge_distance_s = seq(from = -1, to = 1, by = 1))
  
  # amake predictions, assign to dataframe
  tmp[[i+2]] <- fitted(variance_hct_2, newdata=interaction_seq) %>%
    as_tibble() %>%
    rename(f_ll=Q2.5,
           f_ul=Q97.5) %>%
    bind_cols(
      predict(variance_hct_2,
              newdata = interaction_seq) %>%
        as_tibble() %>%
        transmute(p_ll = Q2.5,
                  p_ul = Q97.5),
      interaction_seq)
  
  # name panel for faceting
  tmp[[i+2]]$panel <- paste0("panel_",i)
}

# assemble and export
do.call(rbind, tmp) %>% write.csv("~/Dropbox/andean_range_limits/data/variance_hct_interaction.csv")

# export draws from sampled posterior
variance_hct_1 %>%
  gather_draws(b_bin_elevation_s, b_edge_distance_s, `b_bin_elevation_s:edge_distance_s`) %>%
  write.csv("~/Dropbox/andean_range_limits/data/variance_hct_draws.csv")

# export looic
write.csv(as.data.frame(loo_variance_hct$diffs), "~/Dropbox/andean_range_limits/data/variance_full_hct_loo_elpd.csv")
write.csv(as.data.frame(loo_variance_hct$ic_diffs__), "~/Dropbox/andean_range_limits/data/variance_full_hct_loo_ic.csv")

### variance models, mchc

# here we specify the initial (i.e., starting) values
inits      <- list(Yl = variance_df$variance_mchc)
inits_list <- list(inits, inits)

# fit full model with phylogeny
variance_mchc_1 <- 
  brm(data = variance_df, family = lognormal(),
      variance_mchc ~ 1 + bin_elevation_s + edge_distance_s + bin_elevation_s*edge_distance_s + (1 | gr(phylo, cov=A)),
      data2 = list(A = A),
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="bin_elevation_s"),
        prior(normal(0, 2.5), "b", coef="edge_distance_s"),
        prior(normal(0, 10), "Intercept"),
        prior(student_t(3, 0, 2), "sd"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(variance_mchc_1, nsamples = 50) + xlim(0, 5)

# trace and density plots
plot(variance_mchc_1)

# pairwise plots
pairs(variance_mchc_1)

# fit full model without phylogeny
variance_mchc_2 <- 
  brm(data = variance_df, family = lognormal(),
      variance_mchc ~ 1 + bin_elevation_s + edge_distance_s + bin_elevation_s*edge_distance_s,
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="bin_elevation_s"),
        prior(normal(0, 2.5), "b", coef="edge_distance_s"),
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(variance_mchc_2, nsamples = 50) + xlim(0, 5)

# trace and density plots
plot(variance_mchc_2)

# pairwise plots
pairs(variance_mchc_2)

# fit full model without phylogeny or interaction term
variance_mchc_3 <- 
  brm(data = variance_df, family = lognormal(),
      variance_mchc ~ 1 + bin_elevation_s + edge_distance_s,
      inits = inits_list,
      prior = c(
        prior(normal(0, 2.5), "b", coef="bin_elevation_s"),
        prior(normal(0, 2.5), "b", coef="edge_distance_s"),
        prior(normal(0, 10), "Intercept"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(variance_mchc_3, nsamples = 50) + xlim(0, 5)

# trace and density plots
plot(variance_mchc_3)

# pairwise plots
pairs(variance_mchc_3)

# fit the null model with phylogeny alone
variance_mchc_0 <- 
  brm(data = variance_df, family = lognormal(),
      variance_mchc ~ 0 + (1 | gr(phylo, cov=A)),
      data2 = list(A = A),
      inits = inits_list,
      prior = c(
        prior(student_t(3, 0, 2), "sd"),
        prior(cauchy(0, 2.5), "sigma")),
      iter = 10000, warmup = 5000, cores = 2, chains = 2,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 12),
      save_mevars = TRUE
  )

# pp check
pp_check(variance_mchc_0, nsamples = 50) + xlim(0, 5)

# trace and density plots
plot(variance_mchc_0)

# pairwise plots
pairs(variance_mchc_0)

# compare models
loo_variance_mchc <- loo(variance_mchc_1, variance_mchc_2, variance_mchc_3, variance_mchc_0) 

# export draws from sampled posterior
variance_mchc_1 %>%
  gather_draws(b_bin_elevation_s, b_edge_distance_s, `b_bin_elevation_s:edge_distance_s`) %>%
  write.csv("~/Dropbox/andean_range_limits/data/variance_mchc_draws.csv")

# export looic
write.csv(as.data.frame(loo_variance_mchc$diffs), "~/Dropbox/andean_range_limits/data/variance_full_mchc_loo_elpd.csv")
write.csv(as.data.frame(loo_variance_mchc$ic_diffs__), "~/Dropbox/andean_range_limits/data/variance_full_mchc_loo_ic.csv")
