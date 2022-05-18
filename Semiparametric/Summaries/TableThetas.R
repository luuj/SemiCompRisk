########################################################################################################
# Table 7 for "Causal Inference for Semi-competing risks data" by Nevo and Gorfine
########################################################################################################
library(dplyr)
library(xtable)
library(Daniel)
rm(list = ls())

setwd("/Users/danielnevo/Dropbox/CausalSemiComp/Sims/SemiParam/Summaries/AllResults/")

#### Load a single table with simulation results and true values
# See AllCombine.R
df <- read.csv("AllSimResults.csv")

OneEst <- function(theta0 = NULL, theta1 = NULL, se0  = NULL, se1  = NULL, ret = c("est", "se"))
{
  if(ret =="est"){
    est <- 0.5 * theta0 + 0.5 * theta1
  return(est)}
  if(ret =="se"){
  se <- sqrt(0.25 * se0^2 + 0.25 * se1^2)
  return(se)}
}

#### 50001
T50001.r1 <- df %>% filter(scenario==50001) %>% group_by(cens) %>%
  summarise(m = mean(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T50001.r2 <- df %>% filter(scenario==50001) %>% group_by(cens) %>%
  summarise(s = sd(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T50001.r3 <- df %>% filter(scenario==50001) %>% group_by(cens) %>%
  summarise(m.s = mean(OneEst(se0 =  SE.theta0, se1 =  SE.theta1, ret = "se"))) %>%
  select(-cens) %>% unlist

T50001.r4 <- df %>% filter(scenario==50001) %>% group_by(cens) %>%
  summarise(CP = mean(CoverCInorm(est = OneEst(theta0 = theta0, theta1 = theta1, ret = "est"),
                                  se = OneEst(se0 = SE.theta0, se1 = SE.theta1, ret = "se"),
                                  truth = 2/3))) %>%
  select(-cens) %>% unlist


#### 50002
T50002.r1 <- df %>% filter(scenario==50002) %>% group_by(cens) %>%
  summarise(m = mean(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T50002.r2 <- df %>% filter(scenario==50002) %>% group_by(cens) %>%
  summarise(s = sd(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T50002.r3 <- df %>% filter(scenario==50002) %>% group_by(cens) %>%
  summarise(m.s = mean(OneEst(se0 =  SE.theta0, se1 =  SE.theta1, ret = "se"))) %>%
  select(-cens) %>% unlist

T50002.r4 <- df %>% filter(scenario==50002) %>% group_by(cens) %>%
  summarise(CP = mean(CoverCInorm(est = OneEst(theta0 = theta0, theta1 = theta1, ret = "est"),
                                  se = OneEst(se0 = SE.theta0, se1 = SE.theta1, ret = "se"),
                                  truth = 1))) %>%
  select(-cens) %>% unlist

#### 50003
T50003.r1 <- df %>% filter(scenario==50003) %>% group_by(cens) %>%
  summarise(m = mean(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T50003.r2 <- df %>% filter(scenario==50003) %>% group_by(cens) %>%
  summarise(s = sd(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T50003.r3 <- df %>% filter(scenario==50003) %>% group_by(cens) %>%
  summarise(m.s = mean(OneEst(se0 =  SE.theta0, se1 =  SE.theta1, ret = "se"))) %>%
  select(-cens) %>% unlist

T50003.r4 <- df %>% filter(scenario==50003) %>% group_by(cens) %>%
  summarise(CP = mean(CoverCInorm(est = OneEst(theta0 = theta0, theta1 = theta1, ret = "est"),
                                  se = OneEst(se0 = SE.theta0, se1 = SE.theta1, ret = "se"),
                                  truth = 2))) %>%
  select(-cens) %>% unlist


#### 50501
T50501.r1 <- df %>% filter(scenario==50501) %>% group_by(cens) %>%
  summarise(m = mean(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T50501.r2 <- df %>% filter(scenario==50501) %>% group_by(cens) %>%
  summarise(s = sd(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T50501.r3 <- df %>% filter(scenario==50501) %>% group_by(cens) %>%
  summarise(m.s = mean(OneEst(se0 =  SE.theta0, se1 =  SE.theta1, ret = "se"))) %>%
  select(-cens) %>% unlist

T50501.r4 <- df %>% filter(scenario==50501) %>% group_by(cens) %>%
  summarise(CP = mean(CoverCInorm(est = OneEst(theta0 = theta0, theta1 = theta1, ret = "est"),
                                  se = OneEst(se0 = SE.theta0, se1 = SE.theta1, ret = "se"),
                                  truth = 2/3))) %>%
  select(-cens) %>% unlist


#### 50502
T50502.r1 <- df %>% filter(scenario==50502) %>% group_by(cens) %>%
  summarise(m = mean(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T50502.r2 <- df %>% filter(scenario==50502) %>% group_by(cens) %>%
  summarise(s = sd(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T50502.r3 <- df %>% filter(scenario==50502) %>% group_by(cens) %>%
  summarise(m.s = mean(OneEst(se0 =  SE.theta0, se1 =  SE.theta1, ret = "se"))) %>%
  select(-cens) %>% unlist

T50502.r4 <- df %>% filter(scenario==50502) %>% group_by(cens) %>%
  summarise(CP = mean(CoverCInorm(est = OneEst(theta0 = theta0, theta1 = theta1, ret = "est"),
                                  se = OneEst(se0 = SE.theta0, se1 = SE.theta1, ret = "se"),
                                  truth = 1))) %>%
  select(-cens) %>% unlist

#### 50503
T50503.r1 <- df %>% filter(scenario==50503) %>% group_by(cens) %>%
  summarise(m = mean(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T50503.r2 <- df %>% filter(scenario==50503) %>% group_by(cens) %>%
  summarise(s = sd(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T50503.r3 <- df %>% filter(scenario==50503) %>% group_by(cens) %>%
  summarise(m.s = mean(OneEst(se0 =  SE.theta0, se1 =  SE.theta1, ret = "se"))) %>%
  select(-cens) %>% unlist

T50503.r4 <- df %>% filter(scenario==50503) %>% group_by(cens) %>%
  summarise(CP = mean(CoverCInorm(est = OneEst(theta0 = theta0, theta1 = theta1, ret = "est"),
                                  se = OneEst(se0 = SE.theta0, se1 = SE.theta1, ret = "se"),
                                  truth = 2))) %>%
  select(-cens) %>% unlist


#### 51001
T51001.r1 <- df %>% filter(scenario==51001) %>% group_by(cens) %>%
  summarise(m = mean(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T51001.r2 <- df %>% filter(scenario==51001) %>% group_by(cens) %>%
  summarise(s = sd(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T51001.r3 <- df %>% filter(scenario==51001) %>% group_by(cens) %>%
  summarise(m.s = mean(OneEst(se0 =  SE.theta0, se1 =  SE.theta1, ret = "se"))) %>%
  select(-cens) %>% unlist

T51001.r4 <- df %>% filter(scenario==51001) %>% group_by(cens) %>%
  summarise(CP = mean(CoverCInorm(est = OneEst(theta0 = theta0, theta1 = theta1, ret = "est"),
                                  se = OneEst(se0 = SE.theta0, se1 = SE.theta1, ret = "se"),
                                  truth = 2/3))) %>%
  select(-cens) %>% unlist


#### 51002
T51002.r1 <- df %>% filter(scenario==51002) %>% group_by(cens) %>%
  summarise(m = mean(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T51002.r2 <- df %>% filter(scenario==51002) %>% group_by(cens) %>%
  summarise(s = sd(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T51002.r3 <- df %>% filter(scenario==51002) %>% group_by(cens) %>%
  summarise(m.s = mean(OneEst(se0 =  SE.theta0, se1 =  SE.theta1, ret = "se"))) %>%
  select(-cens) %>% unlist

T51002.r4 <- df %>% filter(scenario==51002) %>% group_by(cens) %>%
  summarise(CP = mean(CoverCInorm(est = OneEst(theta0 = theta0, theta1 = theta1, ret = "est"),
                                  se = OneEst(se0 = SE.theta0, se1 = SE.theta1, ret = "se"),
                                  truth = 1))) %>%
  select(-cens) %>% unlist

#### 51003

T51003.r1 <- df %>% filter(scenario==51003) %>% group_by(cens) %>%
  summarise(m = mean(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T51003.r2 <- df %>% filter(scenario==51003) %>% group_by(cens) %>%
  summarise(s = sd(OneEst(theta0 = theta0, theta1 = theta1, ret = "est"))) %>%
  select(-cens) %>% unlist

T51003.r3 <- df %>% filter(scenario==51003) %>% group_by(cens) %>%
  summarise(m.s = mean(OneEst(se0 =  SE.theta0, se1 =  SE.theta1, ret = "se"))) %>%
  select(-cens) %>% unlist

T51003.r4 <- df %>% filter(scenario==51003) %>% group_by(cens) %>%
  summarise(CP = mean(CoverCInorm(est = OneEst(theta0 = theta0, theta1 = theta1, ret = "est"),
                                  se = OneEst(se0 = SE.theta0, se1 = SE.theta1, ret = "se"), truth = 2))) %>%
  select(-cens) %>% unlist




c11 <- rbind(T50001.r1, T50001.r2, T50001.r3, T50001.r4)
c12 <- rbind(T50501.r1, T50501.r2, T50501.r3, T50501.r4)
c13 <- rbind(T51001.r1, T51001.r2, T51001.r3, T51001.r4)

c21 <- rbind(T50002.r1, T50002.r2, T50002.r3, T50002.r4)
c22 <- rbind(T50502.r1, T50502.r2, T50502.r3, T50502.r4)
c23 <- rbind(T51002.r1, T51002.r2, T51002.r3, T51002.r4)


c31 <- rbind(T50003.r1, T50003.r2, T50003.r3, T50003.r4)
c32 <- rbind(T50503.r1, T50503.r2, T50503.r3, T50503.r4)
c33 <- rbind(T51003.r1, T51003.r2, T51003.r3, T51003.r4)


rbind(cbind(c11, c12, c13),
      cbind(c21, c22, c23),
      cbind(c31, c32, c33)) %>% round(2) %>% xtable


#%>% round(2) %>% xtable
rbind(T50502b.r1, T50502b.r2, T50502b.r3, T50502b.r4) %>% round(2) %>% xtable
