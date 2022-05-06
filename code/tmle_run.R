#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

packages <- c("data.table","tidyverse","skimr","here","haven")

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


# READ AND CREATE DATA
a <- read_csv(here("data","Numom_FV_VitC_Mediation_3.17.22.csv"))

dim(a)

skim(a)

# We need to address missing data
a <- a %>% select(vitc_dens, vita_dens, pct_addsug, pct_satfat, 
                  momage, momeduc, married,insurpub,smokerpre,
                  prediab,prehtn,bmiprepreg ,vitcsuppl_yn,pree_acog) %>% 
  na.omit()


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

# glmnet learner
grid_params <- seq(0,1,by=.25)
lrnr_glmnet <- vector("list", length = length(grid_params))
for(i in 1:length(grid_params)){
  lrnr_glmnet[[i]] <- make_learner(Lrnr_glmnet, alpha = grid_params[i])
}

# xgboost learner
grid_params <- list(max_depth = c(2, 5, 8),
                    eta = c(0.01, 0.1, 0.3))
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
lrnr_xgboost <- vector("list", length = nrow(grid))
for(i in 1:nrow(grid)){
  lrnr_xgboost[[i]] <- make_learner(Lrnr_xgboost, max_depth=grid[i,]$max_depth, eta=grid[i,]$eta)
}


learners_ <- make_learner(Stack, unlist(list(lrnr_mean, lrnr_glm,
                                              lrnr_ranger, lrnr_glmnet,
                                              lrnr_xgboost), 
                                         recursive = TRUE))

# default metalearner appropriate to data types
sl_ <- Lrnr_sl$new(
  learners = learners_
)


# DEFINE SL_Y AND SL_A 
# We only need one, because they're the same
sl_lib <- Lrnr_sl$new(learners = sl_, metalearner = make_learner(Lrnr_nnls))


# PREPARE THE THINGS WE WANT TO FEED IN TO TMLE3
ate_spec <- tmle_ATE(treatment_level = 1, control_level = 0)
learner_list <- list(A = sl_, Y = sl_)
nodes_ <- list(W = names(a)[1:12],
                    A = "vitcsuppl_yn", 
                    Y="pree_acog")

# RUN TMLE3 FOR BOTH FRUIT AND VEGETABLES 
tmle_fit_ <- tmle3(ate_spec, a, nodes_fruit, learner_list)
saveRDS(tmle_fit_fruit, here("misc","tmle_fit-preliminary.rds"))