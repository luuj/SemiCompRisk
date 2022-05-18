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

params <- GetScenarioParams(50001)
true.betas <- as.data.frame(t(c(params$beta.a0.01, params$beta.a0.02, params$beta.a0.12,
                                params$beta.a1.01, params$beta.a1.02, params$beta.a1.12)))
beta.names <- c("beta.a0.01.1", "beta.a0.01.2", "beta.a0.02.1", "beta.a0.02.2", "beta.a0.12.1", "beta.a0.12.2",
                "beta.a1.01.1", "beta.a1.01.2", "beta.a1.02.1", "beta.a1.02.2", "beta.a1.12.1", "beta.a1.12.2")

colnames(true.betas) <- paste0("true.", beta.names)


#### 50501
T50501.r1 <- df %>% filter(scenario==50501) %>% group_by(cens) %>%
  summarise(m0.01.1 = mean(beta.a1.01.1),
            m0.01.2 = mean(beta.a1.01.2),
            m0.02.1 = mean(beta.a1.02.1),
            m0.02.2 = mean(beta.a1.02.2),
            m0.12.1 = mean(beta.a1.12.1),
            m0.12.2 = mean(beta.a1.12.2)) %>%
  select(-cens) %>% unlist

T50501.r2 <- df %>% filter(scenario==50501) %>% group_by(cens) %>%
  summarise(s0.01.1 = sd(beta.a1.01.1),
            s0.01.2 = sd(beta.a1.01.2),
            s0.02.1 = sd(beta.a1.02.1),
            s0.02.2 = sd(beta.a1.02.2),
            s0.12.1 = sd(beta.a1.12.1),
            s0.12.2 = sd(beta.a1.12.2)) %>%
  select(-cens) %>% unlist

T50501.r3 <- df %>% filter(scenario==50501) %>% group_by(cens) %>%
  summarise(s0.01.1 = mean(SE.beta.a1.01.1),
            s0.01.2 = mean(SE.beta.a1.01.2),
            s0.02.1 = mean(SE.beta.a1.02.1),
            s0.02.2 = mean(SE.beta.a1.02.2),
            s0.12.1 = mean(SE.beta.a1.12.1),
            s0.12.2 = mean(SE.beta.a1.12.2)) %>%
  select(-cens) %>% unlist

T50501.r4 <- df %>% filter(scenario==50501) %>% group_by(cens) %>%
  summarise(cp0.01.1 = mean(CoverCInorm(est = beta.a1.01.1, se = SE.beta.a1.01.1,
                                        truth = params$beta.a1.01[1])),
            cp0.01.2 = mean(CoverCInorm(est = beta.a1.01.2, se = SE.beta.a1.01.2,
                                        truth = params$beta.a1.01[2])),
            cp0.02.1 = mean(CoverCInorm(est = beta.a1.02.1, se = SE.beta.a1.02.1,
                                        truth = params$beta.a1.02[1])),
            cp0.02.2 = mean(CoverCInorm(est = beta.a1.02.2, se = SE.beta.a1.02.2,
                                        truth = params$beta.a1.02[2])),
            cp0.12.1 = mean(CoverCInorm(est = beta.a1.12.1, se = SE.beta.a1.12.1,
                                        truth = params$beta.a1.12[1])),
            cp0.12.2 = mean(CoverCInorm(est = beta.a1.12.2, se = SE.beta.a1.12.2,
                                        truth = params$beta.a1.12[2]))) %>%
  select(-cens) %>% unlist

rbind(T50501.r1, T50501.r2, T50501.r3, T50501.r4) %>% round(2) %>% xtable


#### 50502
T50502.r1 <- df %>% filter(scenario==50502) %>% group_by(cens) %>%
  summarise(m0.01.1 = mean(beta.a1.01.1),
            m0.01.2 = mean(beta.a1.01.2),
            m0.02.1 = mean(beta.a1.02.1),
            m0.02.2 = mean(beta.a1.02.2),
            m0.12.1 = mean(beta.a1.12.1),
            m0.12.2 = mean(beta.a1.12.2)) %>%
  select(-cens) %>% unlist

T50502.r2 <- df %>% filter(scenario==50502) %>% group_by(cens) %>%
  summarise(s0.01.1 = sd(beta.a1.01.1),
            s0.01.2 = sd(beta.a1.01.2),
            s0.02.1 = sd(beta.a1.02.1),
            s0.02.2 = sd(beta.a1.02.2),
            s0.12.1 = sd(beta.a1.12.1),
            s0.12.2 = sd(beta.a1.12.2)) %>%
  select(-cens) %>% unlist

T50502.r3 <- df %>% filter(scenario==50502) %>% group_by(cens) %>%
  summarise(s0.01.1 = mean(SE.beta.a1.01.1),
            s0.01.2 = mean(SE.beta.a1.01.2),
            s0.02.1 = mean(SE.beta.a1.02.1),
            s0.02.2 = mean(SE.beta.a1.02.2),
            s0.12.1 = mean(SE.beta.a1.12.1),
            s0.12.2 = mean(SE.beta.a1.12.2)) %>%
  select(-cens) %>% unlist

T50502.r4 <- df %>% filter(scenario==50502) %>% group_by(cens) %>%
  summarise(cp0.01.1 = mean(CoverCInorm(est = beta.a1.01.1, se = SE.beta.a1.01.1,
                                        truth = params$beta.a1.01[1])),
            cp0.01.2 = mean(CoverCInorm(est = beta.a1.01.2, se = SE.beta.a1.01.2,
                                        truth = params$beta.a1.01[2])),
            cp0.02.1 = mean(CoverCInorm(est = beta.a1.02.1, se = SE.beta.a1.02.1,
                                        truth = params$beta.a1.02[1])),
            cp0.02.2 = mean(CoverCInorm(est = beta.a1.02.2, se = SE.beta.a1.02.2,
                                        truth = params$beta.a1.02[2])),
            cp0.12.1 = mean(CoverCInorm(est = beta.a1.12.1, se = SE.beta.a1.12.1,
                                        truth = params$beta.a1.12[1])),
            cp0.12.2 = mean(CoverCInorm(est = beta.a1.12.2, se = SE.beta.a1.12.2,
                                        truth = params$beta.a1.12[2]))) %>%
  select(-cens) %>% unlist

rbind(T50502.r1, T50502.r2, T50502.r3, T50502.r4) %>% round(2) %>% xtable

#### 50503
T50503.r1 <- df %>% filter(scenario==50503) %>% group_by(cens) %>%
  summarise(m0.01.1 = mean(beta.a1.01.1),
            m0.01.2 = mean(beta.a1.01.2),
            m0.02.1 = mean(beta.a1.02.1),
            m0.02.2 = mean(beta.a1.02.2),
            m0.12.1 = mean(beta.a1.12.1),
            m0.12.2 = mean(beta.a1.12.2)) %>%
  select(-cens) %>% unlist

T50503.r2 <- df %>% filter(scenario==50503) %>% group_by(cens) %>%
  summarise(s0.01.1 = sd(beta.a1.01.1),
            s0.01.2 = sd(beta.a1.01.2),
            s0.02.1 = sd(beta.a1.02.1),
            s0.02.2 = sd(beta.a1.02.2),
            s0.12.1 = sd(beta.a1.12.1),
            s0.12.2 = sd(beta.a1.12.2)) %>%
  select(-cens) %>% unlist

T50503.r3 <- df %>% filter(scenario==50503) %>% group_by(cens) %>%
  summarise(s0.01.1 = mean(SE.beta.a1.01.1),
            s0.01.2 = mean(SE.beta.a1.01.2),
            s0.02.1 = mean(SE.beta.a1.02.1),
            s0.02.2 = mean(SE.beta.a1.02.2),
            s0.12.1 = mean(SE.beta.a1.12.1),
            s0.12.2 = mean(SE.beta.a1.12.2)) %>%
  select(-cens) %>% unlist

T50503.r4 <- df %>% filter(scenario==50503) %>% group_by(cens) %>%
  summarise(cp0.01.1 = mean(CoverCInorm(est = beta.a1.01.1, se = SE.beta.a1.01.1,
                                        truth = params$beta.a1.01[1])),
            cp0.01.2 = mean(CoverCInorm(est = beta.a1.01.2, se = SE.beta.a1.01.2,
                                        truth = params$beta.a1.01[2])),
            cp0.02.1 = mean(CoverCInorm(est = beta.a1.02.1, se = SE.beta.a1.02.1,
                                        truth = params$beta.a1.02[1])),
            cp0.02.2 = mean(CoverCInorm(est = beta.a1.02.2, se = SE.beta.a1.02.2,
                                        truth = params$beta.a1.02[2])),
            cp0.12.1 = mean(CoverCInorm(est = beta.a1.12.1, se = SE.beta.a1.12.1,
                                        truth = params$beta.a1.12[1])),
            cp0.12.2 = mean(CoverCInorm(est = beta.a1.12.2, se = SE.beta.a1.12.2,
                                        truth = params$beta.a1.12[2]))) %>%
  select(-cens) %>% unlist

rbind(T50501.r1, T50501.r2, T50501.r3, T50501.r4) %>% round(2) %>% xtable
rbind(T50502.r1, T50502.r2, T50502.r3, T50502.r4) %>% round(2) %>% xtable
rbind(T50503.r1, T50503.r2, T50503.r3, T50503.r4) %>% round(2) %>% xtable

