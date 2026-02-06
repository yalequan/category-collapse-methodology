# ==============================================================================
# GRM Category Collapse Simulation Study
# ==============================================================================
# Citation: Quan, Y., & Wang, C. (2025). Collapsing or Not? A Practical Guide 
#           to Handling Sparse Responses for Polytomous Items. Methodology, 
#           21(1), Article e14303. https://doi.org/10.5964/meth.14303
#
# Purpose: Monte Carlo simulation to investigate the impact of collapsing
#          sparse response categories in polytomous items under the Graded
#          Response Model (GRM). Examines parameter recovery, bias, and model
#          fit across different collapse directions and sparseness levels.
# ==============================================================================

remove(list = ls())
library(mirt)
library(tidyverse)
library(tibble)
library(ggpubr)
library(ggplot2)
library(gdata)
library(dplyr)
options(scipen = 999)

# Packages and Options ----

set.seed(1234) # Keep item parameters are the same across replications

# Simulation settings ----

## Assessment setup ----

R = 150
L = 12 #test length
m = 5 # number of categories
n = c(150, 250, 500, 1000, 1500, 2000)

## Generate item parameters ----

# a parameter

a = rlnorm(L,.4,.5)

# Creating modified b parameters

repeat{
  b1=runif(L,-2,-1)
  b2=runif(L,-1,0)
  b3=runif(L,0,1)
  b4=runif(L,1,2)
  b=matrix(c(b1,b2,b3,b4),nrow=L)
  if(all(apply(apply(b, 1, diff),2,min)>0.5)){break} 
  #ensure 0.5 distance between thresholds
}

b = t(apply(b,1,sort)) #make sure b parameters are ordered


# Items 7, 8, 9 will have a 5% response rate for response option 3
# Items 10, 11, 12 will have a 2.5% response rate for response option 3

# 5% endorsement rate
b[7,3] = b[7,4]-0.25
b[8,3] = b[8,4]-0.28
b[9,3] = b[9,4]-0.24

# 2.5% endorsement rate condition
b[10,3] = b[10,4]-0.11
b[11,3] = b[11,4]-0.15
b[12,3] = b[12,4]-0.14

#discrimination and threshold parameter
d = -a*b

# Create matrix of item parameters
grm_ip=cbind(a,d)

# 5% Simulation -----

### Storage -----

Results_item_fit_5 = matrix(,0,8)
Results_item_parameter_5 = matrix(,0,6)
Results_person_parameter_5 = matrix(,0,6)

# Storage for Standard Errors of Parameter Estimates
Baseline_5_Item_Parameters = matrix(,0,5)
Collapse_Up_5_Item_Parameters = matrix(,0,5)
Collapse_Down_5_Item_Parameters = matrix(,0,5)

# Storage for Theta Parameters
Baseline_5_Theta = matrix(,0,2)
Collapse_Down_5_Theta = matrix(,0,2)
Collapse_Up_5_Theta = matrix(,0,2)
Theta_5_true = matrix(,0,1)

# Storage for item fit indicies
Baseline_5_Item_Fit = matrix(,0,6)
Collapse_Down_5_Item_Fit = matrix(,0,6)
Collapse_Up_5_Item_Fit = matrix(,0,6)

for (N in n){
  cat("\n")
  ### Beginning of Simulation -----
  for (r in 1:R){
    cat("\r5% Run ", r, " Sample Size: ", N)
  ### Generate True Data -----
    set.seed(r)
    theta=rnorm(N)
    storage_theta = as.data.frame(cbind(N, theta, c(1:N)))
    names(storage_theta) = c('Sample_Size', 'theta', "ID")
    Theta_5_true = rbind(Theta_5_true, storage_theta)
    resp_full = simdata(matrix(a),d,itemtype = "graded",Theta=matrix(theta))
    
  ### Fit the Baseline model ----
    md = mirt(resp_full,1,itemtype = "graded", SE = TRUE, verbose = F,
              technical=list(NCYCLES=2000, message = F))
    
  #### Storing Baseline item parameter SE ----
    est.se= mirt::coef(md, as.data.frame =  TRUE, printSE = TRUE)
    est.se = as.data.frame(est.se)
    est.se = tibble::rownames_to_column(est.se, "item")
    est.se = cbind(N, est.se)
    Baseline_5_Item_Parameters = rbind(Baseline_5_Item_Parameters, est.se)
    
  #### Storing Baseline Theta parameters -----
    est.theta = fscores(md, full.scores.SE = TRUE)
    est.theta = as.data.frame(est.theta)
    est.theta = cbind(N, est.theta, c(1:N))
    names(est.theta) = c('Sample_Size', 'theta', 'SE', "ID")
    Baseline_5_Theta = rbind(Baseline_5_Theta, est.theta)
    
  #### Baseline Model-item Fit -----
    fit.baseline = cbind(N, itemfit(md, fit_stats = "S_X2")[, c(1,2,5)]) 
    Baseline_5_Item_Fit = rbind(Baseline_5_Item_Fit, fit.baseline)
    
  ### 5% Collapse UP Condition ----
    
    # Collapsing UP items 7-12 (collapsing 3 up into 4)
    resp_collapsed_UP = resp_full
    resp_collapsed_UP[, 7:12][resp_collapsed_UP[, 7:12] == 3] = 4
  #### Fit 5% Collapse UP Model ----
    md2 = mirt(resp_collapsed_UP,1,itemtype = "graded", SE = TRUE,  verbose = F,
               technical=list(NCYCLES=2000, message = F))
    
  #### Storing collapsed item parameter SE ----
    est.se.collapse.5.up = mirt::coef(md2, as.data.frame =  T, printSE = TRUE)
    est.se.collapse.5.up = as.data.frame(est.se.collapse.5.up)
    est.se.collapse.5.up = tibble::rownames_to_column(est.se.collapse.5.up, "item")
    est.se.collapse.5.up = cbind(N, est.se.collapse.5.up)
    Collapse_Up_5_Item_Parameters = rbind(Collapse_Up_5_Item_Parameters, 
                                          est.se.collapse.5.up)
  #### Storing collapsed Theta parameters ----
    est.theta.up = fscores(md2, full.scores.SE = TRUE)
    est.theta.up = as.data.frame(est.theta.up)
    est.theta.up = cbind(N, est.theta.up, c(1:N))
    names(est.theta.up) = c('Sample_Size', 'theta', 'SE', "ID")
    Collapse_Up_5_Theta = rbind(Collapse_Up_5_Theta, est.theta.up)
    
  #### Model-Item Fit ----
    fit.collapse_up = cbind(N, itemfit(md2, fit_stats = "S_X2")[, c(1,2,5)])
    Collapse_Up_5_Item_Fit = rbind(Collapse_Up_5_Item_Fit, fit.collapse_up)
    
  ### 5% Collapse DOWN Condition ----
    
    # Collapsing DOWN items 7-12 (collapsing 3 DOWN into 2)
    resp_collapsed_DOWN = resp_full
    resp_collapsed_DOWN[, 7:12][resp_collapsed_DOWN[, 7:12] == 3] = 2
  #### Fit 5% Collapse DOWN Condition -----
    md3 = mirt(resp_collapsed_DOWN,1,itemtype = "graded", SE = TRUE,  verbose = F,
               technical=list(NCYCLES=2000, message = F))
    
  #### Storing collapsed item parameter SE ----
    est.se.collapse.5.down = mirt::coef(md3, as.data.frame =  T, printSE = TRUE)
    est.se.collapse.5.down = as.data.frame(est.se.collapse.5.down)
    est.se.collapse.5.down = tibble::rownames_to_column(est.se.collapse.5.down,
                                                        "item")
    est.se.collapse.5.down = cbind(N, est.se.collapse.5.down)
    Collapse_Down_5_Item_Parameters = rbind(Collapse_Down_5_Item_Parameters, 
                                            est.se.collapse.5.down)
  #### Storing collapsed Theta parameters ----
    est.theta.down = fscores(md3, full.scores.SE = TRUE)
    est.theta.down = as.data.frame(est.theta.down)
    est.theta.down = cbind(N, est.theta.down, c(1:N))
    names(est.theta.down) = c('Sample_Size', 'theta', 'SE', "ID")
    Collapse_Down_5_Theta = rbind(Collapse_Down_5_Theta, est.theta.down)
    
  #### Model-Item Fit -----
    fit.collapse_down= cbind(N, itemfit(md3, fit_stats = "S_X2")[, c(1,2,5)])
    Collapse_Down_5_Item_Fit = rbind(Collapse_Down_5_Item_Fit, fit.collapse_down)
  }}

### 5% Simulation Cleanup ----

# Collapsed Up item parameters
item_param_5_collapased_up = Collapse_Up_5_Item_Parameters %>%
  drop_na(SE) %>% 
  separate(item, into = c("item", "param"), "[.]")
item_param_5_collapased_up$item = parse_number(item_param_5_collapased_up$item)

# Collapsed Down item parameters
item_param_5_collapased_down = Collapse_Down_5_Item_Parameters %>%
  drop_na(SE) %>% 
  separate(item, into = c("item", "param"), "[.]")
item_param_5_collapased_down$item = parse_number(item_param_5_collapased_down$item)

# Baseline
item_param_5_baseline = Baseline_5_Item_Parameters %>%
  drop_na(SE) %>% 
  separate(item, into = c("item", "param"), "[.]")
item_param_5_baseline$item = parse_number(item_param_5_baseline$item)

# True Item Parameters
true_item_param_5 = as.data.frame(grm_ip)
true_item_param_5 = cbind(rownames(true_item_param_5), data.frame(true_item_param_5, row.names=NULL))
colnames(true_item_param_5) = c("item", "a1", "d1", "d2", "d3", "d4")
true_item_param_5$item = as.numeric(true_item_param_5$item)

### 5% Parameter Recovery -----

#### Item Parameters -----

true_item_param_5_long = gather(true_item_param_5, param, true, 
                                a1:d4, factor_key = F)

## Calculate Relative Bias and Bias Item Parameters
# For the baseline we can compare the true item parameters with the estimated

rel_bias_baseline_5 = item_param_5_baseline %>%
  group_by(N, item, param) %>%
  summarise(avg = mean(par)) 
rel_bias_baseline_5 = rel_bias_baseline_5%>%
  right_join(true_item_param_5_long, by = c("item", "param"))
rel_bias_baseline_5$rel_bias = 
  (rel_bias_baseline_5$avg - rel_bias_baseline_5$true)/rel_bias_baseline_5$true
rel_bias_baseline_5$collapse = "Baseline"

# For the collapsed up we can compare the true item parameters with the
# estimated. Note that d4 is not compared because the d4_hat parameter is
# combined into the d3_hat parameter so we compare d3_hat against the true d3

rel_bias_up_5 = item_param_5_collapased_up %>%
  group_by(N, item, param) %>%
  summarise(avg = mean(par)) 
rel_bias_up_5 = rel_bias_up_5%>%
  right_join(true_item_param_5_long, by = c("item", "param"))
rel_bias_up_5$rel_bias = 
  (rel_bias_up_5$avg - rel_bias_up_5$true)/rel_bias_up_5$true
rel_bias_up_5$collapse = "Collapse Up Condition"


# For collapsed down we compare:
# a1_hat to a1, d1_hat to d1, d2_hat to d2, and d3_hat to d4. So we need
# to set the "true" d3_hat parameter to be the generated d4 for items 7-12

true_item_param_5_long$true_adj_down = true_item_param_5_long$true

true_item_param_5_long[true_item_param_5_long$item == 7 & 
                         true_item_param_5_long$param == "d3", ][4] = 
  true_item_param_5_long[true_item_param_5_long$item == 7 & 
                           true_item_param_5_long$param == "d4", ][3]

true_item_param_5_long[true_item_param_5_long$item == 8 & 
                         true_item_param_5_long$param == "d3", ][4] =
  true_item_param_5_long[true_item_param_5_long$item == 8 & 
                           true_item_param_5_long$param == "d4", ][3]

true_item_param_5_long[true_item_param_5_long$item == 9 & 
                         true_item_param_5_long$param == "d3", ][4] =
  true_item_param_5_long[true_item_param_5_long$item == 9 & 
                           true_item_param_5_long$param == "d4", ][3]

true_item_param_5_long[true_item_param_5_long$item == 10 & 
                         true_item_param_5_long$param == "d3", ][4] =
  true_item_param_5_long[true_item_param_5_long$item == 10 &
                           true_item_param_5_long$param == "d4", ][3]

true_item_param_5_long[true_item_param_5_long$item == 11 & 
                         true_item_param_5_long$param == "d3", ][4] =
  true_item_param_5_long[true_item_param_5_long$item == 11 &
                           true_item_param_5_long$param == "d4", ][3]

true_item_param_5_long[true_item_param_5_long$item == 12 & 
                         true_item_param_5_long$param == "d3", ][4] =
  true_item_param_5_long[true_item_param_5_long$item == 12 & 
                           true_item_param_5_long$param == "d4", ][3]

rel_bias_down_5 = item_param_5_collapased_down %>%
  group_by(N, item, param) %>%
  summarise(avg = mean(par)) 
rel_bias_down_5 = rel_bias_down_5%>%
  right_join(true_item_param_5_long, by = c("item", "param"))
rel_bias_down_5$rel_bias = 
  (rel_bias_down_5$avg - rel_bias_down_5$true_adj_down)/rel_bias_down_5$true_adj_down
rel_bias_down_5$collapse = "Collapse Down Condition"

item_bias_full_df_5 = rbind(dplyr::select(rel_bias_baseline_5, 
                                          c(N, item, param, rel_bias, collapse)),
                            dplyr::select(rel_bias_up_5, 
                                          c(N, item, param, rel_bias, collapse)),
                            dplyr::select(rel_bias_down_5, 
                                          c(N, item, param, rel_bias, collapse)))

#### Theta Parameter Recovery ----

Theta_5_df = as.data.frame(Theta_5_true)
colnames(Theta_5_df)[2] = "True_Theta"
Theta_5_df$Baseline_est = Baseline_5_Theta$theta
Theta_5_df$Collapse_Up_est = Collapse_Up_5_Theta$theta
Theta_5_df$Collapse_Down_est = Collapse_Down_5_Theta$theta

#### Standard errors -----

Theta_5_df = cbind(Theta_5_df,
                   Baseline_5_Theta$SE, 
                   Collapse_Up_5_Theta$SE, 
                   Collapse_Down_5_Theta$SE)
colnames(Theta_5_df)[7:9] = c("Baseline_SE", "Collapse_up_SE", 
                              "Collapse_down_SE")

#### Theta Bias -----

Theta_5_Bias_Baseline = as.data.frame(cbind(Theta_5_df$Sample_Size,
                                            Theta_5_df$Baseline_est-Theta_5_df$True_Theta))%>%
  group_by(V1) %>%
  summarise(avg_bias = mean(V2))
Theta_5_Bias_Baseline = as.data.frame(Theta_5_Bias_Baseline)
colnames(Theta_5_Bias_Baseline) = c("N", "Avg_Bias")
Theta_5_Bias_Baseline$collapse = "Baseline"

Theta_5_Bias_Up = as.data.frame(cbind(Theta_5_df$Sample_Size,
                                      Theta_5_df$Collapse_Up_est-Theta_5_df$True_Theta))%>%
  group_by(V1) %>%
  summarise(avg_bias = mean(V2))
Theta_5_Bias_Up = as.data.frame(Theta_5_Bias_Up)
colnames(Theta_5_Bias_Up) = c("N", "Avg_Bias")
Theta_5_Bias_Up$collapse = "Collapse_Up"


Theta_5_Bias_Down = as.data.frame(cbind(Theta_5_df$Sample_Size,
                                        Theta_5_df$Collapse_Down_est-Theta_5_df$True_Theta))%>%
  group_by(V1) %>%
  summarise(avg_bias = mean(V2))
Theta_5_Bias_Down = as.data.frame(Theta_5_Bias_Down)
colnames(Theta_5_Bias_Down) = c("N", "Avg_Bias")
Theta_5_Bias_Down$collapse = "Collapse_Down"


Theta_5_Bias = as.data.frame(rbind(Theta_5_Bias_Baseline,
                                   Theta_5_Bias_Up,
                                   Theta_5_Bias_Down))

#### Theta SE -----

Theta_5_SE_Baseline = as.data.frame(cbind(Theta_5_df$Sample_Size,
                                          Theta_5_df$Baseline_SE))%>%
  group_by(V1) %>%
  summarise(avg_bias = mean(V2))
Theta_5_SE_Baseline = as.data.frame(Theta_5_SE_Baseline)
colnames(Theta_5_SE_Baseline) = c("N", "Avg_SE")
Theta_5_SE_Baseline$collapse = "Baseline"

Theta_5_SE_Up = as.data.frame(cbind(Theta_5_df$Sample_Size,
                                    Theta_5_df$Collapse_up_SE))%>%
  group_by(V1) %>%
  summarise(avg_bias = mean(V2))
Theta_5_SE_Up = as.data.frame(Theta_5_SE_Up)
colnames(Theta_5_SE_Up) = c("N", "Avg_SE")
Theta_5_SE_Up$collapse = "Collapse_Up"


Theta_5_SE_Down = as.data.frame(cbind(Theta_5_df$Sample_Size,
                                      Theta_5_df$Collapse_down_SE))%>%
  group_by(V1) %>%
  summarise(avg_bias = mean(V2))
Theta_5_SE_Down = as.data.frame(Theta_5_SE_Down)
colnames(Theta_5_SE_Down) = c("N", "Avg_SE")
Theta_5_SE_Down$collapse = "Collapse_Down"


Theta_5_SE = as.data.frame(rbind(Theta_5_SE_Baseline,
                                 Theta_5_SE_Up,
                                 Theta_5_SE_Down))

#### Flagged item fit ----

Baseline_5_Item_Fit = Baseline_5_Item_Fit %>%
  separate(item, into = c("drop", "item"), "[_]") %>%
  dplyr::select(-drop) %>% 
  mutate_if(is.character, as.numeric)
Baseline_5_Item_Fit$collapse = "Baseline"

Collapse_Down_5_Item_Fit = Collapse_Down_5_Item_Fit%>%
  separate(item, into = c("drop", "item"), "[_]") %>%
  dplyr::select(-drop) %>%
  mutate_if(is.character, as.numeric)
Collapse_Down_5_Item_Fit$collapse = "Collapse Down"

Collapse_Up_5_Item_Fit = Collapse_Up_5_Item_Fit%>%
  separate(item, into = c("drop", "item"), "[_]") %>%
  dplyr::select(-drop) %>%
  mutate_if(is.character, as.numeric)
Collapse_Up_5_Item_Fit$collapse = "Collapse Up"

Item_fit_5 = rbind(Baseline_5_Item_Fit, 
                   Collapse_Down_5_Item_Fit, 
                   Collapse_Up_5_Item_Fit)

Item_fit_5$flag_SX2_05 = ifelse(Item_fit_5$p.S_X2 < 0.05, 1, 0)
Item_fit_5$flag_SX2_10 = ifelse(Item_fit_5$p.S_X2 < 0.10, 1, 0)

Item_fit_5 = Item_fit_5 %>%
  group_by(N, item, collapse) %>%
  summarise(percent_flag_SX2_05 = sum(flag_SX2_05, na.rm =T)/(R),
            percent_flag_SX2_10 = sum(flag_SX2_10, na.rm =T)/(R))

#### Keep Results -----

Item_fit_5$Sparse = "5% Collapse Condition"
Theta_5_SE$Sparse = "5% Collapse Condition"
Theta_5_Bias$Sparse = "5% Collapse Condition"
item_bias_full_df_5$Sparse = "5% Collapse Condition"

gdata::keep(Item_fit_5,
            Theta_5_SE,
            Theta_5_Bias,
            item_bias_full_df_5,
            grm_ip,
            a,
            b,
            d,
            sure = T)

# 2.5% Simulation -----

Results_item_fit_2_5 = matrix(,0,8)
Results_item_parameter_2_5 = matrix(,0,6)
Results_person_parameter_2_5 = matrix(,0,6)

# Storage for Standard Errors of Parameter Estimates
Baseline_2_5_Item_Parameters = matrix(,0,5)
Collapse_Up_2_5_Item_Parameters = matrix(,0,5)
Collapse_Down_2_5_Item_Parameters = matrix(,0,5)

# Storage for Theta Parameters
Baseline_2_5_Theta = matrix(,0,2)
Collapse_Down_2_5_Theta = matrix(,0,2)
Collapse_Up_2_5_Theta = matrix(,0,2)
Theta_2_5_true = matrix(,0,1)

# Storage for item fit indicies
Baseline_2_5_Item_Fit = matrix(,0,6)
Collapse_Down_2_5_Item_Fit = matrix(,0,6)
Collapse_Up_2_5_Item_Fit = matrix(,0,6)

R = 100
L = 12 #test length
m = 5 # number of categories
n = c(150, 250, 500, 1000, 1500, 2000)

for (N in n){
  
  # Beginning of Simulation
  for (r in 1:R){
    
    print(paste0("2.5% Run ", r, " Sample Size: ", N))
    
    ## Generate True Data
    set.seed(r)
    theta=rnorm(N)
    storage_theta = as.data.frame(cbind(N, theta))
    Theta_2_5_true = rbind(Theta_2_5_true, storage_theta)
    resp_full = simdata(matrix(a),d,itemtype = "graded",Theta=matrix(theta))
    
    ### Estimate the Baseline model ----
    md = mirt(resp_full,1,itemtype = "graded", SE = TRUE, verbose = F,
              technical=list(NCYCLES=2000, message = F))
    
    ## Storing Baseline item parameter SE
    est.se= mirt::coef(md, as.data.frame =  TRUE, printSE = TRUE)
    est.se = as.data.frame(est.se)
    est.se = tibble::rownames_to_column(est.se, "item")
    est.se = cbind(N, est.se)
    Baseline_2_5_Item_Parameters = rbind(Baseline_2_5_Item_Parameters, est.se)
    
    ## Storing Baseline Theta parameters
    est.theta = fscores(md, full.scores.SE = TRUE)
    est.theta = as.data.frame(est.theta)
    est.theta = cbind(N, est.theta)
    names(est.theta) = c('Sample_Size', 'theta', 'SE')
    Baseline_2_5_Theta = rbind(Baseline_2_5_Theta, est.theta)
    
    # Model-item Fit
    fit.baseline = cbind(N, itemfit(md, fit_stats = "S_X2")[, c(1,2,5)]) 
    Baseline_2_5_Item_Fit = rbind(Baseline_2_5_Item_Fit, fit.baseline)
    
    ## 2.5% Collapse UP Condition ----
    
    # Collapsing UP items 10-12 (collapsing 3 up into 4)
    resp_collapsed_UP = resp_full
    resp_collapsed_UP[, 10:12][resp_collapsed_UP[, 10:12] == 3] = 4
    # Estimate 5% Collapse UP Model
    md2 = mirt(resp_collapsed_UP,1,itemtype = "graded", SE = TRUE,  verbose = F,
               technical=list(NCYCLES=2000, message = F))
    
    # Storing collapsed item parameter SE
    est.se.collapse.2_5.up = mirt::coef(md2, as.data.frame =  T, printSE = TRUE)
    est.se.collapse.2_5.up = as.data.frame(est.se.collapse.2_5.up)
    est.se.collapse.2_5.up = tibble::rownames_to_column(est.se.collapse.2_5.up, "item")
    est.se.collapse.2_5.up = cbind(N, est.se.collapse.2_5.up)
    Collapse_Up_2_5_Item_Parameters = rbind(Collapse_Up_2_5_Item_Parameters, 
                                            est.se.collapse.2_5.up)
    # Storing collapsed Theta parameters
    est.theta.up = fscores(md2, full.scores.SE = TRUE)
    est.theta.up = as.data.frame(est.theta.up)
    est.theta.up = cbind(N, est.theta.up)
    names(est.theta.up) = c('Sample_Size', 'theta', 'SE')
    Collapse_Up_2_5_Theta = rbind(Collapse_Up_2_5_Theta, est.theta.up)
    
    # Model-Item Fit
    fit.collapse_up = cbind(N, itemfit(md2, fit_stats = "S_X2")[, c(1,2,5)])
    Collapse_Up_2_5_Item_Fit = rbind(Collapse_Up_2_5_Item_Fit, fit.collapse_up)
    
    ### 2.5% Collapse DOWN Condition ----
    
    # Collapsing DOWN items 10-12 (collapsing 3 DOWN into 2)
    resp_collapsed_DOWN = resp_full
    resp_collapsed_DOWN[, 10:12][resp_collapsed_DOWN[, 10:12] == 3] = 2
    # Estimating 5% Collapse DOWN Condition 
    md3 = mirt(resp_collapsed_DOWN,1,itemtype = "graded", SE = TRUE,  verbose = F,
               technical=list(NCYCLES=2000, message = F))
    
    # Storing collapsed item parameter SE
    est.se.collapse.2_5.down = mirt::coef(md3, as.data.frame =  T, printSE = TRUE)
    est.se.collapse.2_5.down = as.data.frame(est.se.collapse.2_5.down)
    est.se.collapse.2_5.down = tibble::rownames_to_column(est.se.collapse.2_5.down,
                                                          "item")
    est.se.collapse.2_5.down = cbind(N, est.se.collapse.2_5.down)
    Collapse_Down_2_5_Item_Parameters = rbind(Collapse_Down_2_5_Item_Parameters, 
                                              est.se.collapse.2_5.down)
    # Storing collapsed Theta parameters
    est.theta.down = fscores(md3, full.scores.SE = TRUE)
    est.theta.down = as.data.frame(est.theta.down)
    est.theta.down = cbind(N, est.theta.down)
    names(est.theta.down) = c('Sample_Size', 'theta', 'SE')
    Collapse_Down_2_5_Theta = rbind(Collapse_Down_2_5_Theta, est.theta.down)
    
    # Model-Item Fit
    fit.collapse_down= cbind(N, itemfit(md3, fit_stats = "S_X2")[, c(1,2,5)])
    Collapse_Down_2_5_Item_Fit = rbind(Collapse_Down_2_5_Item_Fit, fit.collapse_down)
  }}

# Cleanup ----

# Collapsed Up item parameters
item_param_2_5_collapased_up = Collapse_Up_2_5_Item_Parameters %>%
  drop_na(SE) %>% 
  separate(item, into = c("item", "param"), "[.]")
item_param_2_5_collapased_up$item = parse_number(item_param_2_5_collapased_up$item)

# Collapsed Down item parameters
item_param_2_5_collapased_down = Collapse_Down_2_5_Item_Parameters %>%
  drop_na(SE) %>% 
  separate(item, into = c("item", "param"), "[.]")
item_param_2_5_collapased_down$item = parse_number(item_param_2_5_collapased_down$item)

# Baseline
item_param_2_5_baseline = Baseline_2_5_Item_Parameters %>%
  drop_na(SE) %>% 
  separate(item, into = c("item", "param"), "[.]")
item_param_2_5_baseline$item = parse_number(item_param_2_5_baseline$item)

# True Item Parameters
true_item_param_2_5 = as.data.frame(grm_ip)
true_item_param_2_5 = cbind(rownames(true_item_param_2_5), data.frame(true_item_param_2_5, row.names=NULL))
colnames(true_item_param_2_5) = c("item", "a1", "d1", "d2", "d3", "d4")
true_item_param_2_5$item = as.numeric(true_item_param_2_5$item)

### 2.5% Collapse Item Parameter Recovery ----

true_item_param_2_5_long = gather(true_item_param_2_5, param, true, 
                                  a1:d4, factor_key = F)

## Calculate Relative Bias and Bias Item Parameters
# For the baseline we can compare the true item parameters with the estimated

rel_bias_baseline_2_5 = item_param_2_5_baseline %>%
  group_by(N, item, param) %>%
  summarise(avg = mean(par)) 
rel_bias_baseline_2_5 = rel_bias_baseline_2_5%>%
  right_join(true_item_param_2_5_long, by = c("item", "param"))
rel_bias_baseline_2_5$rel_bias = 
  (rel_bias_baseline_2_5$avg - rel_bias_baseline_2_5$true)/rel_bias_baseline_2_5$true
rel_bias_baseline_2_5$collapse = "Baseline"

# For the collapsed up we can compare the true item parameters with the
# estimated. Note that d4 is not compared because the d4_hat parameter is
# combined into the d3_hat parameter so we compare d3_hat against the true d3

rel_bias_up_2_5 = item_param_2_5_collapased_up %>%
  group_by(N, item, param) %>%
  summarise(avg = mean(par)) 
rel_bias_up_2_5 = rel_bias_up_2_5%>%
  right_join(true_item_param_2_5_long, by = c("item", "param"))
rel_bias_up_2_5$rel_bias = 
  (rel_bias_up_2_5$avg - rel_bias_up_2_5$true)/rel_bias_up_2_5$true
rel_bias_up_2_5$collapse = "Collapse Up Condition"


# For collapsed down we compare:
# a1_hat to a1, d1_hat to d1, d2_hat to d2, and d3_hat to d4. So we need
# to set the "true" d3_hat parameter to be the generated d4 for items 7-12

true_item_param_2_5_long$true_adj_down = true_item_param_2_5_long$true

true_item_param_2_5_long[true_item_param_2_5_long$item == 7 & 
                           true_item_param_2_5_long$param == "d3", ][4] = 
  true_item_param_2_5_long[true_item_param_2_5_long$item == 7 & 
                             true_item_param_2_5_long$param == "d4", ][3]

true_item_param_2_5_long[true_item_param_2_5_long$item == 8 & 
                           true_item_param_2_5_long$param == "d3", ][4] =
  true_item_param_2_5_long[true_item_param_2_5_long$item == 8 & 
                             true_item_param_2_5_long$param == "d4", ][3]

true_item_param_2_5_long[true_item_param_2_5_long$item == 9 & 
                           true_item_param_2_5_long$param == "d3", ][4] =
  true_item_param_2_5_long[true_item_param_2_5_long$item == 9 & 
                             true_item_param_2_5_long$param == "d4", ][3]

true_item_param_2_5_long[true_item_param_2_5_long$item == 10 & 
                           true_item_param_2_5_long$param == "d3", ][4] =
  true_item_param_2_5_long[true_item_param_2_5_long$item == 10 &
                             true_item_param_2_5_long$param == "d4", ][3]

true_item_param_2_5_long[true_item_param_2_5_long$item == 11 & 
                           true_item_param_2_5_long$param == "d3", ][4] =
  true_item_param_2_5_long[true_item_param_2_5_long$item == 11 &
                             true_item_param_2_5_long$param == "d4", ][3]

true_item_param_2_5_long[true_item_param_2_5_long$item == 12 & 
                           true_item_param_2_5_long$param == "d3", ][4] =
  true_item_param_2_5_long[true_item_param_2_5_long$item == 12 & 
                             true_item_param_2_5_long$param == "d4", ][3]

rel_bias_down_2_5 = item_param_2_5_collapased_down %>%
  group_by(N, item, param) %>%
  summarise(avg = mean(par)) 
rel_bias_down_2_5 = rel_bias_down_2_5%>%
  right_join(true_item_param_2_5_long, by = c("item", "param"))
rel_bias_down_2_5$rel_bias = 
  (rel_bias_down_2_5$avg - rel_bias_down_2_5$true_adj_down)/rel_bias_down_2_5$true_adj_down
rel_bias_down_2_5$collapse = "Collapse Down Condition"

item_bias_full_df_2_5 = rbind(dplyr::select(rel_bias_baseline_2_5, 
                                            c(N, item, param, rel_bias, collapse)),
                              dplyr::select(rel_bias_up_2_5, 
                                            c(N, item, param, rel_bias, collapse)),
                              dplyr::select(rel_bias_down_2_5, 
                                            c(N, item, param, rel_bias, collapse)))

# Theta Parameter Recovery ----

Theta_2_5_df = as.data.frame(Theta_2_5_true)
colnames(Theta_2_5_df)[2] = "True_Theta"

Theta_2_5_df$Baseline_est = Baseline_2_5_Theta$theta
Theta_2_5_df$Collapse_Up_est = Collapse_Up_2_5_Theta$theta
Theta_2_5_df$Collapse_Down_est = Collapse_Down_2_5_Theta$theta

# Standard errors

Theta_2_5_df = cbind(Theta_2_5_df,
                     Baseline_2_5_Theta$SE, 
                     Collapse_Up_2_5_Theta$SE, 
                     Collapse_Down_2_5_Theta$SE)
colnames(Theta_2_5_df)[6:8] = c("Baseline_SE", "Collapse_up_SE", 
                                "Collapse_down_SE")

# Theta Bias

Theta_2_5_Bias_Baseline = as.data.frame(cbind(Theta_2_5_df$N,
                                              Theta_2_5_df$Baseline_est-Theta_2_5_df$True_Theta))%>%
  group_by(V1) %>%
  summarise(avg_bias = mean(V2))
Theta_2_5_Bias_Baseline = as.data.frame(Theta_2_5_Bias_Baseline)
colnames(Theta_2_5_Bias_Baseline) = c("N", "Avg_Bias")
Theta_2_5_Bias_Baseline$collapse = "Baseline"

Theta_2_5_Bias_Up = as.data.frame(cbind(Theta_2_5_df$N,
                                        Theta_2_5_df$Collapse_Up_est-Theta_2_5_df$True_Theta))%>%
  group_by(V1) %>%
  summarise(avg_bias = mean(V2))
Theta_2_5_Bias_Up = as.data.frame(Theta_2_5_Bias_Up)
colnames(Theta_2_5_Bias_Up) = c("N", "Avg_Bias")
Theta_2_5_Bias_Up$collapse = "Collapse_Up"


Theta_2_5_Bias_Down = as.data.frame(cbind(Theta_2_5_df$N,
                                          Theta_2_5_df$Collapse_Down_est-Theta_2_5_df$True_Theta))%>%
  group_by(V1) %>%
  summarise(avg_bias = mean(V2))
Theta_2_5_Bias_Down = as.data.frame(Theta_2_5_Bias_Down)
colnames(Theta_2_5_Bias_Down) = c("N", "Avg_Bias")
Theta_2_5_Bias_Down$collapse = "Collapse_Down"


Theta_2_5_Bias = as.data.frame(rbind(Theta_2_5_Bias_Baseline,
                                     Theta_2_5_Bias_Up,
                                     Theta_2_5_Bias_Down))

# Theta SE

Theta_2_5_SE_Baseline = as.data.frame(cbind(Theta_2_5_df$N,
                                            Theta_2_5_df$Baseline_SE))%>%
  group_by(V1) %>%
  summarise(avg_bias = mean(V2))
Theta_2_5_SE_Baseline = as.data.frame(Theta_2_5_SE_Baseline)
colnames(Theta_2_5_SE_Baseline) = c("N", "Avg_SE")
Theta_2_5_SE_Baseline$collapse = "Baseline"

Theta_2_5_SE_Up = as.data.frame(cbind(Theta_2_5_df$N,
                                      Theta_2_5_df$Collapse_up_SE))%>%
  group_by(V1) %>%
  summarise(avg_bias = mean(V2))
Theta_2_5_SE_Up = as.data.frame(Theta_2_5_SE_Up)
colnames(Theta_2_5_SE_Up) = c("N", "Avg_SE")
Theta_2_5_SE_Up$collapse = "Collapse_Up"


Theta_2_5_SE_Down = as.data.frame(cbind(Theta_2_5_df$N,
                                        Theta_2_5_df$Collapse_down_SE))%>%
  group_by(V1) %>%
  summarise(avg_bias = mean(V2))
Theta_2_5_SE_Down = as.data.frame(Theta_2_5_SE_Down)
colnames(Theta_2_5_SE_Down) = c("N", "Avg_SE")
Theta_2_5_SE_Down$collapse = "Collapse_Down"


Theta_2_5_SE = as.data.frame(rbind(Theta_2_5_SE_Baseline,
                                   Theta_2_5_SE_Up,
                                   Theta_2_5_SE_Down))

# Flagged item fit

Baseline_2_5_Item_Fit = Baseline_2_5_Item_Fit %>%
  separate(item, into = c("drop", "item"), "[_]") %>%
  dplyr::select(-drop) %>% 
  mutate_if(is.character, as.numeric)
Baseline_2_5_Item_Fit$collapse = "Baseline"

Collapse_Down_2_5_Item_Fit = Collapse_Down_2_5_Item_Fit%>%
  separate(item, into = c("drop", "item"), "[_]") %>%
  dplyr::select(-drop) %>%
  mutate_if(is.character, as.numeric)
Collapse_Down_2_5_Item_Fit$collapse = "Collapse Down"

Collapse_Up_2_5_Item_Fit = Collapse_Up_2_5_Item_Fit%>%
  separate(item, into = c("drop", "item"), "[_]") %>%
  dplyr::select(-drop) %>%
  mutate_if(is.character, as.numeric)
Collapse_Up_2_5_Item_Fit$collapse = "Collapse Up"

Item_fit_2_5 = rbind(Baseline_2_5_Item_Fit, 
                     Collapse_Down_2_5_Item_Fit, 
                     Collapse_Up_2_5_Item_Fit)

Item_fit_2_5$flag_SX2_05 = ifelse(Item_fit_2_5$p.S_X2 < 0.05, 1, 0)
Item_fit_2_5$flag_SX2_10 = ifelse(Item_fit_2_5$p.S_X2 < 0.10, 1, 0)

Item_fit_2_5 = Item_fit_2_5 %>%
  group_by(N, item, collapse) %>%
  summarise(percent_flag_SX2_05 = sum(flag_SX2_05, na.rm =T)/(R),
            percent_flag_SX2_10 = sum(flag_SX2_10, na.rm =T)/(R))

# Keep Results

Item_fit_2_5$Sparse = "2.5% Collapse Condition"
Theta_2_5_SE$Sparse = "2.5% Collapse Condition"
Theta_2_5_Bias$Sparse = "2.5% Collapse Condition"
item_bias_full_df_2_5$Sparse = "2.5% Collapse Condition"

# Save Data ----
gdata::keep(Item_fit_5,
            Theta_5_SE,
            Theta_5_Bias,
            item_bias_full_df_5,
            Item_fit_2_5,
            Theta_2_5_SE,
            Theta_2_5_Bias,
            item_bias_full_df_2_5,
            grm_ip,
            sure = T)

#save.image(file = "Data_GRM.Rdata")

# Plotting -----
library(mirt)
library(tidyverse)
library(tibble)
library(ggpubr)
library(ggplot2)
library(gdata)
library(dplyr)
options(scipen = 999)
load(file = "Data_GRM.Rdata")

plotting_item_bias= rbind(item_bias_full_df_5,
                          item_bias_full_df_2_5)

plotting_theta_bias = rbind(Theta_5_Bias,
                            Theta_2_5_Bias)

plotting_theta_SE = rbind(Theta_5_SE,
                          Theta_2_5_SE)

plotting_item_fit = rbind(Item_fit_5,
                          Item_fit_2_5)

# Person Parameter SE
ggplot(data = plotting_theta_SE, aes(x = factor(N), y = Avg_SE, fill = collapse)) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  xlab("Sample Size") + 
  ylab("Average Standard Error of Recovered Person Parameters") + 
  scale_fill_manual(name = "Collapse Direction", 
                    labels = c("Baseline",
                               "Collapse Down",
                               "Collapse Up"),
                    values = c("black", "darkgray", "grey")) +
  facet_wrap(~Sparse, nrow = 2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  theme(axis.title = element_text(size = 20))

ggplot(data = plotting_theta_SE, aes(x = factor(N), y = Avg_SE, fill = collapse)) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  xlab("Sample Size") + 
  ylab("Average Standard Error of Recovered Person Parameters") + 
  ggtitle("Standard Errors of Recovered GRM Person Parameters") +
  scale_fill_manual(name = "Collapse Direction", 
                    labels = c("Baseline",
                               "Collapse Down",
                               "Collapse Up"),
                    values = c("purple", "chocolate", "lightblue")) +
  facet_wrap(~Sparse, nrow = 2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 30, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  theme(axis.title = element_text(size = 20))


# Person Parameter Bias
ggplot(data = plotting_theta_bias, aes(x = factor(N), y = Avg_Bias, fill = collapse)) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  xlab("Sample Size") + 
  ylab("Average Recovered Person Parameter Bias") + 
  scale_fill_manual(name = "Collapse Direction", 
                    labels = c("Baseline",
                               "Collapse Down",
                               "Collapse Up"),
                    values = c("black", "darkgray", "grey")) +
  facet_wrap(~Sparse, nrow = 2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  ylim(-0.005, 0.005) +
  theme(axis.title = element_text(size = 20))

ggplot(data = plotting_theta_bias, aes(x = factor(N), y = Avg_Bias, fill = collapse)) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  xlab("Sample Size") + 
  ylab("Average Recovered Person Parameter Bias") + 
  ggtitle("Bias of Recovered GRM Person Parameters") +
  scale_fill_manual(name = "Collapse Direction", 
                    labels = c("Baseline",
                               "Collapse Down",
                               "Collapse Up"),
                    values = c("purple", "chocolate", "lightblue")) +
  facet_wrap(~Sparse, nrow = 2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 30, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  ylim(-0.005, 0.005) +
  theme(axis.title = element_text(size = 20))

# Item Parameter Bias
plotting_item_bias = dplyr::filter(plotting_item_bias, item >=7)
plotting_item_bias$item <- sub("^", "Item ", plotting_item_bias$item )
plotting_item_bias$item = as.factor(plotting_item_bias$item)
plotting_item_bias$item = factor(plotting_item_bias$item,
                                 levels = c("Item 7","Item 8",
                                            "Item 9","Item 10","Item 11",
                                            "Item 12"))

plotting_item_bias = plotting_item_bias %>% na.omit()

ggplot(data = dplyr::filter(plotting_item_bias, Sparse == "5% Collapse Condition"), 
       aes(x = factor(N), y = rel_bias, fill = collapse)) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  xlab("Sample Size") + 
  ylab("Average Parameter Bias") + 
  scale_fill_manual(name = "Collapse Direction", 
                      labels = c("Baseline",
                                 "Collapse Down",
                                 "Collapse Up"),
                      values = c("black", "darkgray", "grey")) +
  facet_grid(param ~ item) + 
  ylim(-0.10, 0.10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(strip.text.y = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  theme(axis.title = element_text(size = 20)) + 
  geom_hline(yintercept = 0.05, linetype = "dashed") + 
  geom_hline(yintercept = -0.05, linetype = "dashed")

ggplot(data = dplyr::filter(plotting_item_bias, Sparse == "5% Collapse Condition"), 
       aes(x = factor(N), y = rel_bias, fill = collapse)) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  xlab("Sample Size") + 
  ylab("Average Parameter Bias") + 
  ggtitle("Relative Bias of Recovered GRM Item Parameters",
          subtitle = "5% Sparseness Condition") +
  scale_fill_manual(name = "Collapse Direction", 
                    labels = c("Baseline",
                               "Collapse Down",
                               "Collapse Up"),
                    values = c("purple", "chocolate", "lightblue")) +
  facet_grid(param ~ item) + 
  ylim(-0.10, 0.10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 30, face = "bold"),
        plot.subtitle = element_text(size = 25, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(strip.text.y = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  theme(axis.title = element_text(size = 20)) + 
  geom_hline(yintercept = 0.05, linetype = "dashed") + 
  geom_hline(yintercept = -0.05, linetype = "dashed")

ggplot(data = dplyr::filter(plotting_item_bias, Sparse == "2.5% Collapse Condition"), 
       aes(x = factor(N), y = rel_bias, fill = collapse)) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  xlab("Sample Size") + 
  ylab("Average Parameter Bias") + 
  scale_fill_manual(name = "Collapse Direction", 
                    labels = c("Baseline",
                               "Collapse Down",
                               "Collapse Up"),
                    values = c("black", "darkgray", "grey")) +
  facet_grid(param ~ item) + 
  #ylim(-0.10, 0.10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 25, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(strip.text.y = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  theme(axis.title = element_text(size = 20)) + 
  geom_hline(yintercept = 0.05, linetype = "dashed") + 
  geom_hline(yintercept = -0.05, linetype = "dashed")

ggplot(data = dplyr::filter(plotting_item_bias, Sparse == "2.5% Collapse Condition"), 
       aes(x = factor(N), y = rel_bias, fill = collapse)) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  xlab("Sample Size") + 
  ylab("Average Parameter Bias") + 
  ggtitle("Relative Bias of Recovered GRM Item Parameters",
          subtitle = "2.5% Sparseness Condition") +
  scale_fill_manual(name = "Collapse Direction", 
                    labels = c("Baseline",
                               "Collapse Down",
                               "Collapse Up"),
                    values = c("purple", "chocolate", "lightblue")) +
  facet_grid(param ~ item) + 
  #ylim(-0.10, 0.10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 30, face = "bold"),
        plot.subtitle = element_text(size = 25, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(strip.text.y = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  theme(axis.title = element_text(size = 20)) + 
  geom_hline(yintercept = 0.05, linetype = "dashed") + 
  geom_hline(yintercept = -0.05, linetype = "dashed")

## Item Fit ----

plotting_item_fit1 = dplyr::filter(plotting_item_fit, item >=7)
plotting_item_fit1$item <- sub("^", "Item ", plotting_item_fit1$item )
plotting_item_fit1$item = as.factor(plotting_item_fit1$item)
plotting_item_fit1$item = factor(plotting_item_fit1$item,
                                 levels = c("Item 7","Item 8",
                                            "Item 9","Item 10","Item 11",
                                            "Item 12"))

plotting_item_fit_2 = dplyr::filter(plotting_item_fit, item < 7)
plotting_item_fit_2$item <- sub("^", "Item ", plotting_item_fit_2$item )
plotting_item_fit_2$item = as.factor(plotting_item_fit_2$item)
plotting_item_fit_2$item = factor(plotting_item_fit_2$item,
                                levels = c("Item 1","Item 2",
                                           "Item 3","Item 4","Item 5",
                                           "Item 6"))

ggplot(data = plotting_item_fit1, 
       aes(x = factor(N), y = percent_flag_SX2_05, fill = collapse)) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  xlab("Sample Size") + 
  ylab("Percentage of Replications with Data Model Misift") + 
  scale_fill_manual(name = "Collapse Direction", 
                      labels = c("Baseline",
                                 "Collapse Down",
                                 "Collapse Up"),
                     values = c("black", "darkgray", "grey")) +
  facet_grid(Sparse ~ item) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(strip.text.y = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  theme(axis.title = element_text(size = 20)) + 
  scale_y_continuous(breaks = scales::breaks_width(0.05),
                     labels = scales::percent)
### Lineplot

ggplot(data = plotting_item_fit1, 
       aes(x = factor(N), y = percent_flag_SX2_05, group = collapse)) + 
  geom_point()+
  geom_line(aes(linetype = collapse)) +
  xlab("Sample Size") + 
  ylab("Percentage of Replications with Data Model Misift") + 
  scale_linetype_manual("Collapse Direction",
                        values = c('solid', 'dashed', 'dotted')) +
  facet_grid(Sparse ~ item) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(strip.text.y = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  theme(axis.title = element_text(size = 20)) + 
  scale_y_continuous(breaks = scales::breaks_width(0.05),
                     limits = c(0, 0.15),
                     labels = scales::percent)

ggplot(data = plotting_item_fit_2, 
       aes(x = factor(N), y = percent_flag_SX2_05, group = collapse)) + 
  geom_point()+
  geom_line(aes(linetype = collapse)) +
  xlab("Sample Size") + 
  ylab("Percentage of Replications with Data Model Misift") + 
  scale_linetype_manual("Collapse Direction",
                        values = c('solid', 'dashed', 'dotted')) +
  facet_grid(Sparse ~ item) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(strip.text.y = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  theme(axis.title = element_text(size = 20)) + 
  scale_y_continuous(breaks = scales::breaks_width(0.05),
                     limits = c(0, 0.15),
                     labels = scales::percent)


ggplot(data = filter(plotting_item_fit1, collapse == "Baseline"), 
       aes(x = factor(N), y = percent_flag_SX2_05, group=1)) + 
  geom_point() + 
  geom_line() +
  xlab("Sample Size") + 
  ylab("Percentage of Replications with Data Model Misift") + 
  facet_grid(Sparse ~ item) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(strip.text.y = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  theme(axis.title = element_text(size = 20)) + 
  scale_y_continuous(breaks = scales::breaks_width(0.05),
                     limits = c(0, 0.15),
                     labels = scales::percent)

ggplot(data = filter(plotting_item_fit1, collapse == "Collapse Down"), 
       aes(x = factor(N), y = percent_flag_SX2_05, group=1)) + 
  geom_point() + 
  geom_line() +
  xlab("Sample Size") + 
  ylab("Percentage of Replications with Data Model Misift") + 
  facet_grid(Sparse ~ item) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(strip.text.y = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  theme(axis.title = element_text(size = 20)) + 
  scale_y_continuous(breaks = scales::breaks_width(0.05),
                     labels = scales::percent)

ggplot(data = filter(plotting_item_fit1, collapse == "Collapse Up"), 
       aes(x = factor(N), y = percent_flag_SX2_05, group=1)) + 
  geom_point() + 
  geom_line() +
  xlab("Sample Size") + 
  ylab("Percentage of Replications with Data Model Misift") + 
  facet_grid(Sparse ~ item) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(strip.text.y = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  theme(axis.title = element_text(size = 20)) + 
  scale_y_continuous(breaks = scales::breaks_width(0.05),
                     limits = c(0, 0.15),
                     labels = scales::percent)

###
ggplot(data = plotting_item_fit1, 
       aes(x = factor(N), y = percent_flag_SX2_05, fill = collapse)) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  xlab("Sample Size") + 
  ylab("Percentage of Replications with Data Model Misift") + 
  ggtitle("GRM Cateogry Collapse Data Model Misfit") +
  scale_fill_manual(name = "Collapse Direction", 
                    labels = c("Baseline",
                               "Collapse Down",
                               "Collapse Up"),
                    values = c("purple", "chocolate", "lightblue")) +
  facet_grid(Sparse ~ item) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20)) + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme(plot.title = element_text(size = 30, face = "bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(strip.text.y = element_text(size = 15)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor=element_line(colour="grey", linetype = "dashed"),
        panel.grid.major=element_line(colour="grey", linetype = "dashed")) +
  theme(axis.title = element_text(size = 20)) + 
  scale_y_continuous(breaks = scales::breaks_width(0.05),
                     labels = scales::percent)
