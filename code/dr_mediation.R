#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

outco_list <- 3 #as.numeric(args[1])
media_list <- 2 #as.numeric(args[2])
#delta <- as.numeric(args[3])

## what is the exposure being held at (presume A = 1). Can we set it to A = 0? use 1 - A
## 

packages <- c("data.table","tidyverse","skimr","here","haven","broom")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

remotes::install_github("tlverse/tlverse")

library(tlverse)
library(tmle3)
library(tmle3mediate)
library(sl3)

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

# READ AND CREATE DATA
a <- read_csv(here("data","numom_processed_2022_08_31.csv")) %>% 
  mutate(preeacog_late = pree_acog*(1-preeacog_early))

mean(a$preeacog_late)

dim(a)

names(a)

outcome_list <- c("pree_acog",
                  "preeacog_early",
                  "preeacog_late",
                  "preeacog_sev")

mediator_list <- c("vitc", 
                   "caroten")

exposure_list <- c("fv_totdens_2_5")

covariate_list <-  c("momage",              "bmiprepreg",          "smokecigspre",       
                     "gravidity",           "pa_totmetwk_v1",     
                     "puqe_tot",            "realm_tot",            
                     "momaccult",           "epds_tot_v1",         "stress_tot_v1",
                     "anx_tot",             "walk_nat",            "adi_nat",
                     "povperc",             "momracehisp2",        "momracehisp3",
                     "momracehisp4",        "momracehisp5",        "momracehisp6",
                     "momracehisp7",        "momracehisp8",
                     "smokerpre1",
                     "marital2",            "marital4",            "marital5",
                     "insurance2",          "insurance3",          "momeduc2",
                     "momeduc3",            "momeduc4",            "momeduc5",
                     "momeduc6",            "artcat2",             "artcat3",
                     "alc_bingeprepreg1",   "prediab1",            "prehtn1",
                     "pregplanned1",        "sleepsat1",
                     "hei2015_total_nofv")

covariate_list

# CREATE SUPERLEARNER LIBRARY
# choose base learners

sl3_list_learners("binomial")

lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_gam <- make_learner(Lrnr_gam)

# glmnet learner
grid_params <- seq(0,1,by=.25)
lrnr_glmnet <- vector("list", length = length(grid_params))
for(i in 1:length(grid_params)){
  lrnr_glmnet[[i]] <- make_learner(Lrnr_glmnet, alpha = grid_params[i])
}
lrnr_glmnet <- make_learner(Lrnr_glmnet)

# xgboost learner
grid_params <- list(max_depth = c(2, 5, 8),
                    eta = c(0.01, 0.1, 0.3))
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
lrnr_xgboost <- vector("list", length = nrow(grid))
for(i in 1:nrow(grid)){
  lrnr_xgboost[[i]] <- make_learner(Lrnr_xgboost, max_depth=grid[i,]$max_depth, eta=grid[i,]$eta)
}
lrnr_xgboost <- make_learner(Lrnr_xgboost)

# MARS
grid_params <- c(3,4,5)
lrnr_earth <- vector("list", length = length(grid_params))
for(i in 1:length(grid_params)){
  lrnr_earth[[i]] <- make_learner(Lrnr_earth, degree = grid_params[i])
}
lrnr_earth <- make_learner(Lrnr_earth)

sl_ <- make_learner(Stack, unlist(list(lrnr_xgboost,
                                       lrnr_glmnet,
                                       lrnr_earth,
                                       lrnr_mean,
                                       lrnr_glm), 
                                      recursive = TRUE))
sl_

# we picked mostly regression based learners to 
# avoid problems observed here with trees --> impossibly high variance

# glmnet VarImp learner
lrn_lasso <- Lrnr_glmnet$new(alpha = 1)

lasso_screen <- Lrnr_screener_coefs$new(learner = lrn_lasso,
                                        threshold = 0.2, # this keeps only coefs with \beta \geq exp(.2) in the logit model
                                        max_retain = 15)

LASSOscreen1_top15_stack <- Pipeline$new(lasso_screen, sl_)

sl_2 <- make_learner(Stack, unlist(list(sl_, 
                                        LASSOscreen1_top15_stack),
                                   recursive = T))

sl_2 <- make_learner(Lrnr_glm)

# DEFINE SL_Y AND SL_A 
# We only need one, because they're the same

Q_learner <- Lrnr_sl$new(learners = sl_, 
                         metalearner = Lrnr_nnls$new(convex=T))
g_learner <- Lrnr_sl$new(learners = sl_, 
                         metalearner = Lrnr_nnls$new(convex=T))
learner_list <- list(Y = Q_learner,
                     A = g_learner)

# PREPARE THE THINGS WE WANT TO FEED IN TO TMLE3
tmle_spec_NIE <- tmle_NIE(
  e_learners = Q_learner,
  psi_Z_learners = Q_learner,
  max_iter = 1
)

(covariate_list %in% names(a))
(exposure_list[1] %in% names(a))
(mediator_list[1] %in% names(a))
(outcome_list[1] %in% names(a))

nodes_ <- list(W = covariate_list, 
               A = exposure_list, 
               Z = mediator_list[media_list],
               Y = outcome_list[outco_list])

# RUN TMLE3 
set.seed(1234)
tmle_NIE_fit <- tmle3(
  tmle_spec_NIE, a, nodes_, learner_list
)

influence_function <- as.data.frame(as.matrix(tmle_NIE_fit$estimates[[1]]$IC))
names(influence_function) <- "IC"

write_csv(
  influence_function,
  here("misc", paste0("influence_function_mediation", "_",
                      mediator_list[media_list], "_", 
                      outcome_list[outco_list], ".csv"))
)

# export RD to table
RD <- tmle_NIE_fit$summary$psi_transformed
RD_LCL <- tmle_NIE_fit$summary$lower_transformed
RD_UCL <- tmle_NIE_fit$summary$upper_transformed
round(tibble(RD,RD_LCL,RD_UCL),3)

write_csv(
  round(tibble(RD,RD_LCL,RD_UCL),3),
  here("misc", paste0("tmle_med_results_", "_",
                      mediator_list[media_list], "_", 
                      outcome_list[outco_list], ".csv"))
  )

saveRDS(tmle_NIE_fit, 
        here("misc", paste0("tmle_med_fit_", "_",
                            mediator_list[media_list], "_", 
                            outcome_list[outco_list], ".rds"))
        )

# diagnostics
tmle_task <- tmle_spec_NIE$make_tmle_task(a, nodes_)

initial_likelihood <- tmle_spec_NIE$make_initial_likelihood(
  tmle_task,
  learner_list
)

## save propensity score for diagnosis
propensity_score <- initial_likelihood$get_likelihoods(tmle_task)$A
propensity_score <- propensity_score*a$fv_totdens_2_5 + (1 - propensity_score)*(1 - a$fv_totdens_2_5)

# min and max
print(min(propensity_score))
print(max(propensity_score))

plap_ <- tibble(exposure=a$fv_totdens_2_5,pscore=propensity_score)
plap_ <- plap_ %>% mutate(sw= exposure*(mean(exposure)/propensity_score) + 
                            (1-exposure)*((1-mean(exposure))/(1-propensity_score)))

# distribution of PS and stabilized weights
summary(plap_$sw)
summary(propensity_score)

## initial diagnostic
ggplot(plap_) + 
  geom_jitter(aes(y = pscore, 
                  x = factor(exposure)),
              width=.1, height=0, alpha=.25)

ggsave(here("figures","ps_diagnosis-2022_06_02.pdf"),
       width = 15,
       height = 15, 
       units = "cm")

# ps overlap plot
ggplot(plap_) + geom_histogram(aes(y = ..density.., pscore,fill=factor(exposure)),
                               colour="grey50", 
                               alpha=0.75, bins=50,
                               position="identity") +
  scale_fill_manual(values=c("blue","orange")) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

ggsave(here("figures","ps_overlap_hist-2022_06_02.pdf"),
       width = 15,
       height = 15, 
       units = "cm")

ggplot(plap_) + geom_density(aes(pscore,color=factor(exposure))) +
  scale_fill_manual(values=c("blue","orange")) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

ggsave(here("figures","ps_overlap_dens-2022_06_02.pdf"),
       width = 15,
       height = 15, 
       units = "cm")

# save outcome predictions for diagnosis
outcome_preds <- initial_likelihood$get_likelihoods(tmle_task)$Y

# min and max
print(min(outcome_preds))
print(max(outcome_preds))

opred_ <- tibble(outcome=a$pree_acog,pred=outcome_preds)

## initial diagnostic
ggplot(opred_) + 
  geom_jitter(aes(y = pred, 
                  x = factor(outcome)),
              width=.1, height=0, alpha=.25)

ggsave(here("figures","outcome_diagnosis-2022_06_02.pdf"),
       width = 15,
       height = 15, 
       units = "cm")

# super learner coefficients for PS model
g_fit <- tmle_NIE_fit$likelihood$factor_list[["A"]]$learner
print(g_fit$fit_object$full_fit$learner_fits$Lrnr_nnls_TRUE)

# super learner coefficients for mediator model
E_fit <- tmle_NIE_fit$likelihood$factor_list[["E"]]$learner
print(E_fit$fit_object$full_fit$learner_fits$Lrnr_nnls_TRUE)

# super learner coefficients for outcome model
Q_fit <- tmle_NIE_fit$likelihood$factor_list[["Y"]]$learner
print(Q_fit$fit_object$full_fit$learner_fits$Lrnr_nnls_TRUE)


# PIDE
# set the IPSI multiplicative shift
# 
# Edward notes:
# 1) estimate IPSI effect using npcausal package
# 2)
# 
# 
delta_ipsi <- delta

plap_ <- plap_ %>% 
  mutate(propensity_score_shift = (delta_ipsi*propensity_score) / ((delta_ipsi*propensity_score) + 1 - propensity_score))

write_csv(
  plap_,
  here("misc", paste0("ps_overlap_shift_data", "_",
                      mediator_list[media_list], "_", 
                      outcome_list[outco_list], ".csv"))
)

ggplot(plap_) + 
  geom_density(aes(propensity_score)) +
  geom_density(aes(propensity_score_shift), linetype="dashed") +
  scale_x_continuous(expand=c(0,0), limits = c(0,1)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Propensity Score") + ylab("Density")

ggsave(here("figures","ps_shift_dens-2022_06_02.pdf"),
       width = 15,
       height = 15, 
       units = "cm")

# instantiate tmle3 spec for stochastic mediation
tmle_spec_pie_decomp <- tmle_medshift(
  delta = delta_ipsi,
  e_learners = Q_learner,
  phi_learners = Q_learner
)

# compute the TML estimate
pie_decomp <- tmle3(
  tmle_spec_pie_decomp, a, nodes_, learner_list
)
pie_decomp

# get the PIDE

mean_Y <- a %>% summarize(meanY = mean(get(outcome_list[outco_list])))
mean_Y <- as.numeric(mean_Y)

RD_pie <- pie_decomp$summary$tmle_est - mean_Y
RD_LCL_pie <- pie_decomp$summary$lower - mean_Y # NB: doesn't account for variance in mean_Y!!
RD_UCL_pie <- pie_decomp$summary$upper - mean_Y # NB: doesn't account for variance in mean_Y!!
round(tibble(RD_pie,RD_LCL_pie,RD_UCL_pie),3)

write_csv(
  round(tibble(RD_pie,RD_LCL_pie,RD_UCL_pie),3),
  here("misc", paste0("tmle_pie_results_", exposure_list, "_", outcome_list[outco_list], ".csv"))
)

saveRDS(pie_decomp, 
        here("misc", paste0("pie_decomp", exposure_list, "_", outcome_list[outco_list], ".rds"))
)