###############################################################
# ControlThompsaon_EpsilonGreedy.R
# Yuki Atsusaka
# Created: 9/12/2021
# Last update: 9/28/2021
# Aim: to implement the Control-augmented Thompson Sampling
# and Epsilon Greedy algorithms for adaptive experiments
###############################################################

library(tidyverse)
rm(list=ls())

###############################################
# FIRST BATCH
###############################################

# Using U.S. Video Treatments as an example
dt <- read.csv("analysis_ready.csv", sep=";")
dt_batch1 <- dt[1:800,]     # Simulating Batch 1
dt_batch2 <- dt[801:1609,]  # Simulating Batch 2 (This is not observed before computing epsilon in reality)

table(dt_batch1$clicked)     # See the outcome
table(dt_batch1$Treatment2)  # See the treatment arms

K <- length(unique(dt_batch1$Treatment2)) # Number of Treatment Arms

# Just to rename and clean data
dt_batch1 <- dt_batch1 %>% 
             mutate(outcome = clicked, #
                    treatment = Treatment2) %>% 
              select(outcome, treatment)   

dt_cum <- rbind(dt_batch1)  # Simulating Cumulative Sample

# Estimate the mean outcome values
# (This needs to be modified depending on the outcome type)

table(dt_cum$treatment) # Use this order (Cash, CDC, Lottery)

mean_cash <- mean(dt_cum$outcome[dt_cum$treatment=="Cash Voucher Incentive"])
mean_cdc <- mean(dt_cum$outcome[dt_cum$treatment=="CDC Health Information"])
mean_lottery <- mean(dt_cum$outcome[dt_cum$treatment=="Lottery Incentive"])

effect_vec <- c(mean_cash, mean_cdc, mean_lottery)
arm_vec <- c("Cash Voucher Incentive", "CDC Health Information", "Lottery Incentive")

#++++++++++++++++++++++++++++++++++++++++++++++++++#
# THE CONTROL-AUGMENTED THOMPSON SAMPLING ALGORITHM
#++++++++++++++++++++++++++++++++++++++++++++++++++#

# THOMPSON SAMPLING #############################################
S <- 10000 # Number of Monte Carlo simulations (Not sample size!)
# K # Number of treatment arms defined above

# Number of 1s
successes <- table(dt_cum$outcome, dt_cum$treatment)[2,] 

# Number of samples in each group
trials <- c(sum(table(dt_cum$outcome, dt_cum$treatment)[,1]),
            sum(table(dt_cum$outcome, dt_cum$treatment)[,2]),
            sum(table(dt_cum$outcome, dt_cum$treatment)[,3])) 
  
# Sampling from a Beta-Binomial posterior density S times
set.seed(142)   
draws <- replicate(n=S, rbeta(K, successes+1, trials-successes+1))
    
# Computing how many times k-th arm had the highest posterior
argmax <- apply(X=draws, MARGIN=2, FUN=which.max)

# Tallying up the maximal probabilities = Thompson sampling probabilities
assign_probs_ts <- table(argmax) / S # Normalizing with S

# We will update the above probabilities below following 
# the control-augmented design proposed in Offert-Westort et al. (2021)
################################################################

# Name of the current best arm
current_best_arm <- arm_vec[which.max(assign_probs_ts)]

# Cumulative sample size for each group
cum_sample_size <- table(dt_cum$treatment)

# N(Current Best Arm) - N(Control) [in Cumulative Sample]
control_best_diff <- table(dt_cum$treatment)[current_best_arm] - table(dt_cum$treatment)["CDC Health Information"] 

# Whether N(Control) is smaller than N(Current Best Arm)
is_control_smaller <- control_best_diff > 0

# Part I: Reserving samples for the control, if necessary
if(is_control_smaller==TRUE){

expected_batch_size <- 200 # Expected N of the next batch [Tuning Parameter]
prob_ceiling <- 0.9        # Arbitrary ceiling to reserve non-zero prob for non-controls [Tuning Parameter]

control_rsv <- control_best_diff/expected_batch_size # Reserved prob for the control
control_rsv <- min(control_rsv, prob_ceiling)        # Choose whichever is small

remaining_prob <- 1 - control_rsv # The remaining assignment probability

}else{
  
remaining_prob <- 1 # Initialized to one  
}

# Part II: Compute Thompson Sampling posterior probabilities
control_prob_K <- remaining_prob/K           # Assignment prob for the control group (no exploration)
control_prob <- control_rsv + control_prob_K # Add these probs up

remaining_prob_ts <- 1 - control_prob        # The remaining prob that can be use for the "final" Thompson Sampling

assign_probs_ts_no_control <- assign_probs_ts[-2] # 2 = Control [This part must be generalized]

ts_prob <- remaining_prob_ts*assign_probs_ts_no_control # Weighting the remaining probs with the initial posteriors

assign_probs_batch2_ts <- c(ts_prob[1], control_prob, ts_prob[2]) # [This part must be generalized]

print(assign_probs_batch2_ts)
print(sum(assign_probs_batch2_ts))           # This must sum up to one

#++++++++++++++++++++++++++++++++++++++++++++++++++#
# THE EPSILON GREEDY ALGORITHM
#++++++++++++++++++++++++++++++++++++++++++++++++++#

greedy_arm <- arm_vec[which.max(effect_vec)] # Arm with the highest mean outcome

t = 2 # Second batch

epsilon <- 1/K * 1/(t^0.5) # Time-varying epsilon (t=time count, 0.5=decay rate)

# CHECK HOW EPSILON BEHAVES IN A TIME-VARYING SETTING
# t_check <- seq(2, 40, by = 1) # Simulating 20 batches
# epsilon_check1 <- 1/(t_check^0.5)
# epsilon_check2 <- 1/K * 1/(t_check^0.5)
# plot(t_check, epsilon_check1, xlim=c(2,40), ylim=c(0,1), xlab="Number of Batches (t)",
#      yaxt="n", xaxt="n", ylab="Epsilon = Assignment Probability for Non-Best Arms", pch=16, col="firebrick4")
# points(t_check, epsilon_check2, col="navy", pch=16)
# text(x=11, y=0.6, labels="epsilon = 1/(t^0.5)", col="firebrick4", cex=1)
# text(x=9, y=0.04, labels="epsilon = 1/K * 1/(t^0.5)", col="navy", cex=1)
# abline(h=0.1, lty=2, col="gray50")
# axis(1, at=seq(2, 40, by = 1), las=1)
# axis(2, at=seq(0, 1, by = 0.1), las=1)

#epsilon <- 0.1           # Time-fixed Epsilon 
c_epsilon <- 1 - epsilon  # Compliment of epsilon

# Assigning (1-epsilon) to the greedy arm (with the highest mean outcome)
# while assigning epsilon/(K-1) to the remaining arms
assign_probs_batch2_eg <- ifelse(greedy_arm==arm_vec, c_epsilon, epsilon/(K-1))


probs <- rbind(assign_probs_batch2_ts, # Thompson Sampling
               assign_probs_batch2_eg) # Epsilon Greedy
colnames(probs) <- c("Cash", "CDC", "Lottery")
print(probs)

############################
# HYPOTHETICAL NEW SUBJECTS
############################
next_obs <- seq(400) # Next 400 observations (next batch)
L <- length(next_obs)

# Next Treatment Assignment
set.seed(1234)
next_treat_ts <- sample(arm_vec, size=L, replace=T, prob=assign_probs_batch2_ts)
next_treat_eg <- sample(arm_vec, size=L, replace=T, prob=assign_probs_batch2_eg)

# Transform these into a data.frame
new_ts <- data.frame(next_obs, next_treat_ts)
head(new_ts, n=20)

new_eg <- data.frame(next_obs, next_treat_eg)
head(new_eg, n=20)

# Compare the mean outcome values, and repeat.

###############################################
# SECOND BATCH
###############################################
# Assume that we obtain Batch 2 with the above probabilities
head(dt_batch2)
dt_batch2 <- dt_batch2 %>% 
             mutate(outcome = clicked, #
                    treatment = Treatment2) %>% 
              select(outcome, treatment)   

# Just to rename and clean data
dt_cum <- rbind(dt_batch1, dt_batch2)  # Cumulative Sample

# Estimate the mean outcome values

mean_cash <- mean(dt_cum$outcome[dt_cum$treatment=="Cash Voucher Incentive"])
mean_cdc <- mean(dt_cum$outcome[dt_cum$treatment=="CDC Health Information"])
mean_lottery <- mean(dt_cum$outcome[dt_cum$treatment=="Lottery Incentive"])
rbind(mean_cash, mean_cdc, mean_lottery)

#########################################################
# FIGURE OUT HOW TO INCLUDE WEIGHTS
#########################################################
# # Given that treatment assignment is not "random" as in uniform, our mean estimates are biased
# # To fix this issue, we will rely on the standard inverse probability weighting approach
# # with stabilized weights proposed in Hajek (1971)
# mean_cash_ts_ipw <- mean(dt_cum$outcome[dt_cum$treatment=="Cash Voucher Incentive"]/assign_probs_batch2_ts[1])
# mean_cdc_ts_ipw <- mean(dt_cum$outcome[dt_cum$treatment=="CDC Health Information"]/assign_probs_batch2_ts[2])
# mean_lottery_ts_ipw <- mean(dt_cum$outcome[dt_cum$treatment=="Lottery Incentive"]/assign_probs_batch2_ts[3])
# rbind(mean_cash_ts_ipw, mean_cdc_ts_ipw, mean_lottery_ts_ipw)
# # Extremely small assignment probs may cause some issues
# 
# # Let's use the Hajekc stabilized weights instead to avoid extreme weights 
# mean_cash_ts_haj <- mean_cash_ts_ipw/

#++++++++++++++++++++++++++++++++++++++++++++++++++#
# THE CONTROL-AUGMENTED THOMPSON SAMPLING ALGORITHM
#++++++++++++++++++++++++++++++++++++++++++++++++++#

# THOMPSON SAMPLING #############################################
S <- 10000 # Number of Monte Carlo simulations (Not sample size!)
# K # Number of treatment arms defined above

# Number of 1s
successes <- table(dt_cum$outcome, dt_cum$treatment)[2,] 

# Number of samples in each group
trials <- c(sum(table(dt_cum$outcome, dt_cum$treatment)[,1]),
            sum(table(dt_cum$outcome, dt_cum$treatment)[,2]),
            sum(table(dt_cum$outcome, dt_cum$treatment)[,3])) 
  
# Sampling from a Beta-Binomial posterior density S times
set.seed(142)   
draws <- replicate(n=S, rbeta(K, successes+1, trials-successes+1))
    
# Computing how many times k-th arm had the highest posterior
argmax <- apply(X=draws, MARGIN=2, FUN=which.max)

# Tallying up the maximal probabilities = Thompson sampling probabilities
assign_probs_ts <- table(argmax) / S # Normalizing with S

# We will update the above probabilities below following 
# the control-augmented design proposed in Offert-Westort et al. (2021)
################################################################

# Name of the current best arm
current_best_arm <- arm_vec[which.max(assign_probs_ts)]

# Cumulative sample size for each group
cum_sample_size <- table(dt_cum$treatment)

# N(Current Best Arm) - N(Control) [in Cumulative Sample]
control_best_diff <- table(dt_cum$treatment)[current_best_arm] - table(dt_cum$treatment)["CDC Health Information"] 

# Whether N(Control) is smaller than N(Current Best Arm)
is_control_smaller <- control_best_diff > 0

# Part I: Reserving samples for the control, if necessary
if(is_control_smaller==TRUE){

expected_batch_size <- 200 # Expected N of the next batch [Tuning Parameter]
prob_ceiling <- 0.9        # Arbitrary ceiling to reserve non-zero prob for non-controls [Tuning Parameter]

control_rsv <- control_best_diff/expected_batch_size # Reserved prob for the control
control_rsv <- min(control_rsv, prob_ceiling)        # Choose whichever is small

remaining_prob <- 1 - control_rsv # The remaining assignment probability

}else{
  
remaining_prob <- 1 # Initialized to one  
}

# Part II: Compute Thompson Sampling posterior probabilities
control_prob_K <- remaining_prob/K           # Assignment prob for the control group (no exploration)
control_prob <- control_rsv + control_prob_K # Add these probs up

remaining_prob_ts <- 1 - control_prob        # The remaining prob that can be use for the "final" Thompson Sampling

assign_probs_ts_no_control <- assign_probs_ts[-2] # 2 = Control [This part must be generalized]

ts_prob <- remaining_prob_ts*assign_probs_ts_no_control # Weighting the remaining probs with the initial posteriors

assign_probs_batch3_ts <- c(ts_prob[1], control_prob, ts_prob[2]) # [This part must be generalized]

print(assign_probs_batch3_ts)
print(sum(assign_probs_batch3_ts))           # This must sum up to one

#++++++++++++++++++++++++++++++++++++++++++++++++++#
# THE EPSILON GREEDY ALGORITHM
#++++++++++++++++++++++++++++++++++++++++++++++++++#

effect_vec_ts <- c(mean_cdc, mean_cash, mean_lottery)
arm_vec <- c("CDC Health Information", "Cash Voucher Incentive", "Lottery Incentive")
            
greedy_arm <- arm_vec[which.max(effect_vec)] # Arm with the highest mean outcome

t = 3 # Third batch next time

epsilon <- 1/K * 1/(t^0.5) # Time-varying epsilon (t=time count, 0.5=decay rate)
c_epsilon <- 1 - epsilon # Compliment of epsilon

# Assigning (1-epsilon) to the greedy arm (with the highest mean outcome)
# while assigning epsilon/(K-1) to the remaining arms
assign_probs_batch3_eg <- ifelse(greedy_arm==arm_vec, c_epsilon, epsilon/(K-1))
assign_probs_batch3_eg
print(sum(assign_probs_batch3_eg)) # This must be one

# COMPARE
rbind(assign_probs_batch2_ts, 
      assign_probs_batch3_ts,
      assign_probs_batch2_eg, 
      assign_probs_batch3_eg)

# NEXT STEP
# SIMULATE COVARIATE DISTRIBTION  IN EACH BATCH IN THE PRESENCE OF HETEROGENEOUS TREATMETNE EFFECTS

###############################################################
# END OF THIS R SOURCE FILE
###############################################################

