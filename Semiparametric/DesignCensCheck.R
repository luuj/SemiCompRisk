####################################################################################
####################################################################################
# CausalSemiComp
# Calculate parameters of the censoring distribution to achieve desired censoring rates.
####################################################################################
rm(list = ls())
library(dplyr)
library(survival)
library(CausalSemiComp)
library(boot)

set.seed(17)

####################################################################################

cens.rate.L <- c(0.0045, 0.0037, 0.00001, 0.0045, 0.0037, 0.00001, 0.0045, 0.0037, 0.00001)
cens.rate.M <- c(0.0295, 0.026, 0.0165, 0.0295, 0.026, 0.0165, 0.0295, 0.026, 0.0165)
scenarios <- c(50001, 50002, 50003, 50501, 50502, 50503, 51001, 51002, 51003)
Res <- as.data.frame(cbind(scenarios, cens.rate.L, cens.rate.M))
colnames(Res) <- c("Scenario", "exp.rate.L", "exp.rate.M", "rho", "theta")
for (j in 1:9)
{
  Daniel::CatIndex(j)
  params <- GetScenarioParams(scenario.num = scenarios[j])
  Res$rho[j] <- params$rho
  Res$theta[j] <- params$theta
  sim.df.L <- SimDataWeibFrail(n.sample = 2000000, params = params, no.protected = F, no.large = F,
                               cens.exp.rate = cens.rate.L[j], cens.admin = 100)
  sim.df.M <- SimDataWeibFrail(n.sample = 2000000, params = params, no.protected = F, no.large = F,
                               cens.exp.rate = cens.rate.M[j], cens.admin = 100)
  Res$`Cens-L-T1`[j] <- mean(1 - sim.df.L$delta1) %>% round(3)
  Res$`Cens-L-T2`[j] <- mean(1 - sim.df.L$delta2) %>% round(3)
  Res$`Cens-M-T1`[j] <- mean(1 - sim.df.M$delta1) %>% round(3)
  Res$`Cens-M-T2`[j] <- mean(1 - sim.df.M$delta2) %>% round(3)
}

Res
Res %>% round(3)
