##This script will reconstruct the trial data and 2 models, one for all-cause mortality (ACM) and another for cardiovascular mortality (CVM)##

#Load the necessary packages (ensure that they are installed on your machine beforehand)
library(brms)
library(dplyr)
library(tidyr)
library(bayestestR)
library(marginaleffects)
library(rstan)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(ggdist)

#Save the stanmodel to hard disk
rstan_options(auto_write = TRUE)
#Allow Stan to use more cores in order to increase computational speed
options(mc.cores = parallel::detectCores())

###Record the number of patients and of events in the primary and secondary prevention subgroups
#Of patients (Table 1 in the NEJM publication)
ba_total_pri <- 2100
placebo_total_pri <- 2106
ba_total_sec <- 4892
placebo_total_sec <- 4872
#of all-cause mortality (ACM) events
ba_acm_pri <- 75 #Abstract of JAMA paper 
placebo_acm_pri <- 109 #Abstract of JAMA paper
ba_acm_sec <- 434 - 75 #Table 2 of NEJM - Abstract of JAMA paper
placebo_acm_sec <- 420 - 109 #Table 2 of NEJM - Abstract of JAMA paper
#of CVM events (if one would like to calculate effect on CVM instead of ACM)
ba_cvm_pri <- 37 #Abstract of JAMA paper 
placebo_cvm_pri <- 65 #Abstract of JAMA paper
ba_cvm_sec <- 269 - 37 #Table 2 of NEJM - Abstract of JAMA paper
placebo_cvm_sec <- 257 - 65 #Table 2 of NEJM - Abstract of JAMA paper



###Construct two tables reflecting the above information for ACM
#Primary Prevention Table
primary_acm_data <- data.frame(
  prevention_status = "Primary", #Specify that this is for primary prevention
  death = c(1,  1, 0, 0), #Create 2 rows for those who died and 2 for those who lived (stratified by treatment arm)
  ba = c(1, 0, 1, 0), #Create treatment arms (stratified by those who died and those who lived)
  n = c(ba_acm_pri, placebo_acm_pri, #Number of all-cause deaths in the BA and palcebo arms
        ba_total_pri - ba_acm_pri, placebo_total_pri - placebo_acm_pri) #Number of people alive in the BA and placebo arms (total - death)
)
#Secondary Prevention Table
secondary_acm_data <- data.frame(
  prevention_status = "Secondary", #Specify that this is for secondary prevention
  death = c(1,  1, 0, 0),#Create 2 rows for those who died and 2 for those who lived (stratified by treatment arm)
  ba = c(1, 0, 1, 0), #Create treatment arms (stratified by those who died and those who lived)
  n = c(ba_acm_sec, placebo_acm_sec, #Number of all-cause deaths in the BA and palcebo arms
        ba_total_sec - ba_acm_sec, placebo_total_sec - placebo_acm_sec) #Number of people alive in the BA and placebo arms (total - death)
)



#Merge the two tables to recreate the entire trial dataset stratified into primary and secondary prevention
acm_data <- bind_rows(primary_acm_data, secondary_acm_data)
#Expand our aggregated data (so that each row represents a single patient)
acm_data <- acm_data %>% uncount(n)


###Construct two tables reflecting the above information for cardiovascular mortality (CVM)
#Primary Prevention Table
primary_cvm_data <- data.frame(
  prevention_status = "Primary", #Specify that this is for primary prevention
  death = c(1,  1, 0, 0), #Create 2 rows for those who died and 2 for those who lived (stratified by treatment arm)
  ba = c(1, 0, 1, 0), #Create treatment arms (stratified by those who died and those who lived)
  n = c(ba_cvm_pri, placebo_cvm_pri, #Number of all-cause deaths in the BA and palcebo arms
        ba_total_pri - ba_cvm_pri, placebo_total_pri - placebo_cvm_pri) #Number of people alive in the BA and placebo arms (total - death)
)
#Secondary Prevention Table
secondary_cvm_data <- data.frame(
  prevention_status = "Secondary", #Specify that this is for secondary prevention
  death = c(1,  1, 0, 0),#Create 2 rows for those who died and 2 for those who lived (stratified by treatment arm)
  ba = c(1, 0, 1, 0), #Create treatment arms (stratified by those who died and those who lived)
  n = c(ba_cvm_sec, placebo_cvm_sec, #Number of all-cause deaths in the BA and palcebo arms
        ba_total_sec - ba_cvm_sec, placebo_total_sec - placebo_cvm_sec) #Number of people alive in the BA and placebo arms (total - death)
)

#Merge the two tables
cvm_data <- bind_rows(primary_cvm_data, secondary_cvm_data)
#Expand our data (so that each row represents a single patient)
cvm_data <- cvm_data %>% uncount(n)



#Create formula
formula <- bf( #Set a formula such that:
  death ~ #Death is modeled as a function of:
    ba + #Assignment to bempedoic acid
    prevention_status #Whether the patient falls under the primary or secondary prevention categories
  + ba*prevention_status #An interaction between the two
)


#Set weak priors (can be changed as necessary; in this example, weak priors were used so as to allow the data from
#this trial only to shape the posterior estimates)
weak_priors <- c(
  prior(normal(0, 10), class = b, coef = "prevention_statusSecondary"), #Prior about BA's Treatment effect. This assumes we have no prior information on the effects of previous CVD on mortality rates and that we estimate it completely from the data.
  prior(normal(0, 10), class = b, coef = "ba"), #Prior about BA's Treatment effect (on primary prevention). This assumes we have no clear sense of the likely range of reductions or increases in mortality likely to be observed.
  prior(normal(0, 10), class = b, coef = "ba:prevention_statusSecondary") #Prior about how secondary prevention status changes BA's effect. This assumes we have no clear sense of the likely extent of interactions likely to be observed.
)



#Run model with the priors of your choice to observe the range of values deemed plausible
#Note that this simulation only takes the priors into account and does not yet incorporate the data
priors_model <- brm(data = acm_data, #Use the dataset for ACM (switch to CVM if desired)
               family = bernoulli, #Using a bernoulli distribution for the binary outcome (either dead or alive)
               formula = formula, #Use the aforementioned formula
               seed = 100, #Set seed for reproducibility
               prior = weak_priors, #Use a weak prior (can be changed as desired)
               sample_prior = "only") #Ignore the data and sample from the prior only

#Simulate predictions to see what the resultant risk ratios are before the data (based on our priors)
avg_comparisons(priors_model, 
                by = "prevention_status", 
                variables = "ba",
                comparison = "ratio")

#Run the model incorporating both the prior you chose and the data at hand (for all-cause mortality).
b_model_acm <- brm(data = acm_data, #Use the dataset for ACM
               family = bernoulli, #Using a bernoulli distribution for the outcome (either dead or alive)
               formula = formula, #Use our formula
               seed = 100, #Set seed for reproducibility
               prior = weak_priors, #Use a weak prior (can be changed as desired)
               iter = 5000) #Use 5000 iterations per chain

#Run the model for cardiovascular mortality as well
b_model_cvm <- brm(data = cvm_data, #Use the dataset for CVM
                   family = bernoulli, #Using a bernoulli distribution for the outcome (either dead or alive)
                   formula = formula, #Use our formula
                   seed = 100, #Set seed for reproducibility
                   prior = weak_priors, #Use a weak prior (can be changed as desired)
                   iter = 5000) #Use 5000 iterations per chain



