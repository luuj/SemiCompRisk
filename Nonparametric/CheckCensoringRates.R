####################################################################################
####################################################################################
# CausalSemiComp
# Censoring rates for the nonparameteric simulations
####################################################################################
rm(list = ls())
library(dplyr)
library(survival)
library(CausalSemiComp)
library(boot)

set.seed(17)
#####################################################################
no.protected <- T
no.large <- T
n.sample <- 100000
#####################################################################

# For Scenario 10752
params <- GetScenarioParams(scenario.num = 10752)
cens.rate <- 0.0053 #

sim.df <- SimDataWeibFrail(n.sample = n.sample, params, no.protected = no.protected, cens.exp.rate = cens.rate)
mean(sim.df$delta2)
mean(sim.df$delta1)

cens.rate <- 0.0185 #
sim.df <- SimDataWeibFrail(n.sample = n.sample, params, no.protected = no.protected, cens.exp.rate = cens.rate)
mean(sim.df$delta2)
mean(sim.df$delta1)

# For Scenario 2
params <- GetScenarioParams(scenario.num = 20001)
cens.rate <- 0.007 #

sim.df <- SimDataWeibFrail(n.sample = n.sample, params, no.protected = no.protected, cens.exp.rate = cens.rate)
mean(sim.df$delta2)
mean(sim.df$delta1)

cens.rate <- 0.025 #
sim.df <- SimDataWeibFrail(n.sample = n.sample, params, no.protected = no.protected, cens.exp.rate = cens.rate)
mean(sim.df$delta2)
mean(sim.df$delta1)

# For Scenario 3
params <- GetScenarioParams(scenario.num = 30001)
cens.rate <- 0.0055 #

sim.df <- SimDataWeibFrail(n.sample = n.sample, params, no.protected = no.protected, cens.exp.rate = cens.rate)
mean(sim.df$delta2)
mean(sim.df$delta1)

cens.rate <- 0.020 #
sim.df <- SimDataWeibFrail(n.sample = n.sample, params, no.protected = no.protected, cens.exp.rate = cens.rate)
mean(sim.df$delta2)
mean(sim.df$delta1)



