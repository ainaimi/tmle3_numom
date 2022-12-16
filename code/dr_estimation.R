#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

for (iii in 1:3){
  for (jjj in 1:6){

# expos_list <- as.numeric(args[1])
# outco_list <- as.numeric(args[2])

expos_list <- iii
outco_list <- jjj

print(expos_list)
print(outco_list)

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
library(sl3)

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

#### TODO:
#### 1) implement screening algorithms
#### 2) generate diagnostic plots for PS and outcome model
#### 3) double check covariates for inclusion in PS and outcome model (w Lisa)
#### 4) data preprocessing (centering and scaling, variable transformations)
#### 5) add seed to command args
#### 6) generalize program for use in multiple projects
#### 7) find the quantile that corresponds to 70 mg for vit c; use this quantile to define cutpoint for total caroten; use half the value of that cutpoint to define caroten dens cutpoit

# READ AND CREATE DATA
a <- read_csv(here("data","numom_processed_2022_08_31.csv"))%>% 
  mutate(preeacog_late = pree_acog*(1-preeacog_early))

mean(a$preeacog_late)

dim(a)

names(a)

outcome_list <- c("pree_acog", 
                  "preeacog_early",
                  "preeacog_late",
                  "preeacog_sev",
                  "vitc",
                  "caroten")

exposure_list <- c("fv_totdens_2_5",
                   "vitc_70",
                   "caroten_7400")

if(expos_list==1){
  covariate_list <-  c("momage",              "bmiprepreg",          "smokecigspre",       
                       "gravidity",           "pa_totmetwk_v1",     
                       "puqe_tot",            "realm_tot",            #"dt_kcal",            
                       "momaccult",           "epds_tot_v1",         "stress_tot_v1",      
                       "anx_tot",             "walk_nat",            "adi_nat",            
                       "povperc",             "momracehisp2",        "momracehisp3",       
                       "momracehisp4",        "momracehisp5",        "momracehisp6",       
                       "momracehisp7",        "momracehisp8",        "smokerpre1",         
                       "marital2",            "marital4",            "marital5",           
                       "insurance2",          "insurance3",          "momeduc2",           
                       "momeduc3",            "momeduc4",            "momeduc5",           
                       "momeduc6",            "artcat2",             "artcat3",            
                       "prediab1",            "prehtn1",             "alc_bingeprepreg1", 
                       "pregplanned1",        "sleepsat1")#, 
                       #"hei2015_total_nofv")
} else{
  covariate_list <-  c("momage",              "bmiprepreg",          "smokecigspre",       
                       "gravidity",           "pa_totmetwk_v1",     
                       "puqe_tot",            "realm_tot",            #"dt_kcal",            
                       "momaccult",           "epds_tot_v1",         "stress_tot_v1",      
                       "anx_tot",             "walk_nat",            "adi_nat",            
                       "povperc",             "momracehisp2",        "momracehisp3",       
                       "momracehisp4",        "momracehisp5",        "momracehisp6",       
                       "momracehisp7",        "momracehisp8",        "smokerpre1",         
                       "marital2",            "marital4",            "marital5",           
                       "insurance2",          "insurance3",          "momeduc2",           
                       "momeduc3",            "momeduc4",            "momeduc5",           
                       "momeduc6",            "artcat2",             "artcat3",            
                       "prediab1",            "prehtn1",             "alc_bingeprepreg1", 
                       "pregplanned1",        "sleepsat1", 
                       "hei2015_total_nofv",  "fv_totdens")
}


# CREATE SUPERLEARNER LIBRARY
# choose base learners

sl3_list_learners("binomial")

lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm)

# ranger learner
grid_params <- list(num.trees = c(250, 500, 1000, 2000),
                    mtry = c(2,4,6),
                    min.node.size = c(50,100))
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
lrnr_ranger <- vector("list", length = nrow(grid))
for(i in 1:nrow(grid)){
  lrnr_ranger[[i]] <- make_learner(Lrnr_ranger, 
                                   num.trees=grid[i,]$num.trees, 
                                   mtry=grid[i,]$mtry,
                                   min.node.size=grid[i,]$min.node.size)
}
lrnr_ranger <- make_learner(Lrnr_ranger)

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

# nnet learner
grid_params <- list(size = c(3, 5, 7, 9),
                    decay = c(0, 0.001, 0.01, 0.1),
                    skip = c(T,F))
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
lrnr_nnet <- vector("list", length = nrow(grid))
for(i in 1:nrow(grid)){
  lrnr_nnet[[i]] <- make_learner(Lrnr_nnet, size=grid[i,]$size, decay=grid[i,]$decay, skip=grid[i,]$skip)
}
lrnr_nnet <- make_learner(Lrnr_nnet)

sl_ <- make_learner(Stack, unlist(list(lrnr_nnet,
                                       lrnr_xgboost,
                                       lrnr_glmnet,
                                       lrnr_ranger,
                                       lrnr_mean,
                                       lrnr_glm), 
                                      recursive = TRUE))
sl_

# ranger VarImp learner
ranger_with_importance <- Lrnr_ranger$new(importance = "impurity_corrected")
RFscreenIC_top10 <- Lrnr_screener_importance$new(
  learner = ranger_with_importance, num_screen = 10
)
RFscreenIC_top10_stack <- Pipeline$new(RFscreenIC_top10, sl_)

ranger_with_importance <- Lrnr_ranger$new(importance = "impurity_corrected")
RFscreenIC_top15 <- Lrnr_screener_importance$new(
  learner = ranger_with_importance, num_screen = 15
)
RFscreenIC_top15_stack <- Pipeline$new(RFscreenIC_top15, sl_)

ranger_with_importance <- Lrnr_ranger$new(importance = "permutation")
RFscreenPerm_top10 <- Lrnr_screener_importance$new(
  learner = ranger_with_importance, num_screen = 10
)
RFscreenPerm_top10_stack <- Pipeline$new(RFscreenPerm_top10, sl_)

ranger_with_importance <- Lrnr_ranger$new(importance = "permutation")
RFscreenPerm_top15 <- Lrnr_screener_importance$new(
  learner = ranger_with_importance, num_screen = 15
)
RFscreenPerm_top15_stack <- Pipeline$new(RFscreenPerm_top15, sl_)

# glmnet VarImp learner
lrn_lasso <- Lrnr_glmnet$new(alpha = 1)
lasso_screen <- Lrnr_screener_coefs$new(learner = lrn_lasso,
                                        threshold = 0.1,
                                        max_retain = 10)
LASSOscreen1_top10_stack <- Pipeline$new(lasso_screen, sl_)

lasso_screen <- Lrnr_screener_coefs$new(learner = lrn_lasso,
                                        threshold = 0.1,
                                        max_retain = 15)
LASSOscreen1_top15_stack <- Pipeline$new(lasso_screen, sl_)

lasso_screen <- Lrnr_screener_coefs$new(learner = lrn_lasso,
                                        threshold = 0.2,
                                        max_retain = 10)
LASSOscreen2_top10_stack <- Pipeline$new(lasso_screen, sl_)

lasso_screen <- Lrnr_screener_coefs$new(learner = lrn_lasso,
                                        threshold = 0.2,
                                        max_retain = 15)
LASSOscreen2_top15_stack <- Pipeline$new(lasso_screen, sl_)

sl_2 <- make_learner(Stack, unlist(list(sl_, 
                                        # RFscreenIC_top10_stack,
                                        # RFscreenIC_top15_stack,
                                        # RFscreenPerm_top10_stack,
                                        # RFscreenPerm_top15_stack,
                                        # LASSOscreen1_top10_stack,
                                        LASSOscreen1_top15_stack),
                                        # LASSOscreen2_top10_stack,
                                        # LASSOscreen2_top15_stack),
                                   recursive = T))

# DEFINE SL_Y AND SL_A 
# We only need one, because they're the same

Q_learner <- Lrnr_sl$new(learners = sl_2, 
                         metalearner = Lrnr_nnls$new(convex=T))
g_learner <- Lrnr_sl$new(learners = sl_2, 
                         metalearner = Lrnr_nnls$new(convex=T))
learner_list <- list(Y = Q_learner,
                     A = g_learner)

# PREPARE THE THINGS WE WANT TO FEED IN TO TMLE3
ate_spec <- tmle_ATE(treatment_level = 1, 
                     control_level = 0)
nodes_ <- list(W = covariate_list, 
               A = exposure_list[expos_list], 
               Y = outcome_list[outco_list])

# RUN TMLE3 
set.seed(123)
tmle_fit <- tmle3(ate_spec, 
                  a, 
                  nodes_, 
                  learner_list)

# export RD to table
RD <- tmle_fit$summary$psi_transformed
RD_LCL <- tmle_fit$summary$lower_transformed
RD_UCL <- tmle_fit$summary$upper_transformed
write_csv(
  round(tibble(RD,RD_LCL,RD_UCL),3),
  here("misc", paste0("tmle_results_", 
                      exposure_list[expos_list], "_", 
                      outcome_list[outco_list], ".csv"))
  )

saveRDS(tmle_fit, 
        here("misc", paste0("tmle_fit_", 
                            exposure_list[expos_list], "_", 
                            outcome_list[outco_list], ".rds"))
        )

  }
}

tmle_task <- ate_spec$make_tmle_task(a, nodes_)

initial_likelihood <- ate_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)

## save propensity score for diagnosis
propensity_score <- initial_likelihood$get_likelihoods(tmle_task)$A
propensity_score <- propensity_score*a$tot_c + (1 - propensity_score)*(1 - a$tot_c)

# min and max
print(min(propensity_score))
print(max(propensity_score))

plap_ <- tibble(exposure=a$tot_c,pscore=propensity_score)
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
g_fit <- tmle_fit$likelihood$factor_list[["A"]]$learner
print(g_fit)
coef(g_fit)

# super learner coefficients for outcome model
Q_fit <- tmle_fit$likelihood$factor_list[["Y"]]$learner
print(Q_fit)
coef(Q_fit)

## AIPW

outcome <- a$pree_acog
exposure <- a$tot_c_weight
covariates.Q <- a %>% select(all_of(covariate_list))
covariates.g <- a %>% 
  select(all_of(covariate_list)) %>% 
  select(-hei2015_total_score, -pa_totmetwk_v1, 
         -dt_kcal)

library(AIPW)

set.seed(123)
AIPW_SL <- AIPW$new(Y = outcome,
                    A = exposure,
                    W.Q = covariates.Q,
                    W.g = covariates.g,
                    Q.SL.library = Q_learner,
                    g.SL.library = g_learner,
                    k_split = 10,
                    save.sl.fit=T,
                    verbose=FALSE)$
  fit()$
  summary(g.bound = 0.025)$ 
  plot.p_score()

AIPW_SL$result[3,]

str(AIPW_SL$obs_est)

plot_dat <- tibble(aipw_eif1 = AIPW_SL$obs_est$aipw_eif1,
                   aipw_eif0 = AIPW_SL$obs_est$aipw_eif0,
                   aipw_eif_diff = aipw_eif1 - aipw_eif0,
                   tmle_eif_diff = tmle_fit_$estimates[[1]]$IC,
                   BMI = a$bmiprepreg)

ggplot(plot_dat) + 
  geom_point(aes(y = tmle_eif_diff, x = BMI)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Scaled Prepregnancy BMI") + ylab("EIF for Y1 - Y0")


p1 <- ggplot(plot_dat) + 
  geom_histogram(aes(aipw_eif1)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Count") + xlab("EIF for Y1")
p2 <- ggplot(plot_dat) + 
  geom_histogram(aes(aipw_eif0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Count") + xlab("EIF for Y0")
p3 <- ggplot(plot_dat) + 
  geom_histogram(aes(aipw_eif_diff)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Count") + xlab("EIF for Y1 - Y0")
p4 <- ggplot(plot_dat) + 
  geom_point(aes(y = aipw_eif_diff, x = BMI)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Scaled Prepregnancy BMI") + ylab("EIF for Y1 - Y0")

plot_eif <- gridExtra::grid.arrange(p1, p2, p3, p4, nrow=2)
ggsave(filename = here("figures","eif_plot.pdf"),
       plot = plot_eif)
  
