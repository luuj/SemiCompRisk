########################################################################################################
# Table 3 for "Causal Inference for Semi-competing risks data" by Nevo and Gorfine
########################################################################################################
library(dplyr)
library(xtable)
library(Daniel)
rm(list = ls())

setwd("/Users/danielnevo/Dropbox/CausalSemiComp/Sims/SemiParam/Summaries/AllResults/")

#### Load a single table with simulation results and true values
# See AllCombine.R
df <- read.csv("AllSimResults.csv")

#### Scenario 50501
Table3.50501.r1 <- df %>% filter(scenario==50501) %>% group_by(cens) %>% summarise(m.ATE.T1.ad = mean(ATE.T1.ad),
                                                                med.ATE.T1.ad = mean(med.ATE.T1.ad),
                                                                m.ATE.T2.ad = mean(ATE.T2.ad),
                                                                med.ATE.T2.ad = mean(med.ATE.T2.ad),
                                                                m.ATE.T2.nd = mean(ATE.T2.nd),
                                                                med.ATE.T2.nd = mean(med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist
Table3.50501.r2 <- df %>% filter(scenario==50501) %>% group_by(cens) %>% summarise(emp.s.ATE.T1.ad = sd(ATE.T1.ad),
                                                                emp.s.med.ATE.T1.ad = sd(med.ATE.T1.ad),
                                                                emp.s.ATE.T2.ad = sd(ATE.T2.ad),
                                                                emp.s.med.ATE.T2.ad = sd(med.ATE.T2.ad),
                                                                emp.s.ATE.T2.nd = sd(ATE.T2.nd),
                                                                emp.s.med.ATE.T2.nd = sd(med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist

Table3.50501.r3 <- df %>% filter(scenario==50501) %>% group_by(cens) %>% summarise(est.s.ATE.T1.ad = mean(SE.ATE.T1.ad),
                                                                                   est.s.med.ATE.T1.ad = mean(SE.med.ATE.T1.ad),
                                                                                   est.s.ATE.T2.ad = mean(SE.ATE.T2.ad),
                                                                                   est.s.med.ATE.T2.ad = mean(SE.med.ATE.T2.ad),
                                                                                   est.s.ATE.T2.nd = mean(SE.ATE.T2.nd),
                                                                                   est.s.med.ATE.T2.nd = mean(SE.med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist


Table3.50501.r4 <- df %>% filter(scenario==50501) %>% group_by(cens) %>%
  summarise(cp.ATE.T1.ad = mean(CoverCInorm(est = ATE.T1.ad, se = SE.ATE.T1.ad, truth = true.ATE.T1.ad)),
            cp.med.ATE.T1.ad = mean(CoverCInorm(est = med.ATE.T1.ad, se = SE.med.ATE.T1.ad, truth = true.med.ATE.T1.ad)),
            cp.ATE.T2.ad = mean(CoverCInorm(est = ATE.T2.ad, se = SE.ATE.T2.ad, truth = true.ATE.T2.ad)),
            cp.med.ATE.T2.ad = mean(CoverCInorm(est = med.ATE.T2.ad, se = SE.med.ATE.T2.ad, truth = true.med.ATE.T2.ad)),
            cp.ATE.T2.nd = mean(CoverCInorm(est = ATE.T2.nd, se = SE.ATE.T2.nd, truth = true.ATE.T2.nd)),
            cp.med.ATE.T2.nd = mean(CoverCInorm(est = med.ATE.T2.nd, se = SE.med.ATE.T2.nd, truth = true.med.ATE.T2.nd))) %>%
  select(-cens) %>% unlist

rbind(Table3.50501.r1, Table3.50501.r2, Table3.50501.r3, Table3.50501.r4)  %>% round(2) %>% xtable

####### 50502
Table3.50502.r1 <- df %>% filter(scenario==50502) %>% group_by(cens) %>% summarise(m.ATE.T1.ad = mean(ATE.T1.ad),
                                                                                   med.ATE.T1.ad = mean(med.ATE.T1.ad),
                                                                                   m.ATE.T2.ad = mean(ATE.T2.ad),
                                                                                   med.ATE.T2.ad = mean(med.ATE.T2.ad),
                                                                                   m.ATE.T2.nd = mean(ATE.T2.nd),
                                                                                   med.ATE.T2.nd = mean(med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist


Table3.50502.r2 <- df %>% filter(scenario==50502) %>% group_by(cens) %>% summarise(emp.s.ATE.T1.ad = sd(ATE.T1.ad),
                                                                                   emp.s.med.ATE.T1.ad = sd(med.ATE.T1.ad),
                                                                                   emp.s.ATE.T2.ad = sd(ATE.T2.ad),
                                                                                   emp.s.med.ATE.T2.ad = sd(med.ATE.T2.ad),
                                                                                   emp.s.ATE.T2.nd = sd(ATE.T2.nd),
                                                                                   emp.s.med.ATE.T2.nd = sd(med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist

Table3.50502.r3 <- df %>% filter(scenario==50502) %>% group_by(cens) %>% summarise(est.s.ATE.T1.ad = mean(SE.ATE.T1.ad),
                                                                                   est.s.med.ATE.T1.ad = mean(SE.med.ATE.T1.ad),
                                                                                   est.s.ATE.T2.ad = mean(SE.ATE.T2.ad),
                                                                                   est.s.med.ATE.T2.ad = mean(SE.med.ATE.T2.ad),
                                                                                   est.s.ATE.T2.nd = mean(SE.ATE.T2.nd),
                                                                                   est.s.med.ATE.T2.nd = mean(SE.med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist


Table3.50502.r4 <- df %>% filter(scenario==50502) %>% group_by(cens) %>%
  summarise(cp.ATE.T1.ad = mean(CoverCInorm(est = ATE.T1.ad, se = SE.ATE.T1.ad, truth = true.ATE.T1.ad)),
            cp.med.ATE.T1.ad = mean(CoverCInorm(est = med.ATE.T1.ad, se = SE.med.ATE.T1.ad, truth = true.med.ATE.T1.ad)),
            cp.ATE.T2.ad = mean(CoverCInorm(est = ATE.T2.ad, se = SE.ATE.T2.ad, truth = true.ATE.T2.ad)),
            cp.med.ATE.T2.ad = mean(CoverCInorm(est = med.ATE.T2.ad, se = SE.med.ATE.T2.ad, truth = true.med.ATE.T2.ad)),
            cp.ATE.T2.nd = mean(CoverCInorm(est = ATE.T2.nd, se = SE.ATE.T2.nd, truth = true.ATE.T2.nd)),
            cp.med.ATE.T2.nd = mean(CoverCInorm(est = med.ATE.T2.nd, se = SE.med.ATE.T2.nd, truth = true.med.ATE.T2.nd))) %>%
  select(-cens) %>% unlist


rbind(Table3.50502.r1, Table3.50502.r2, Table3.50502.r3, Table3.50502.r4) %>% round(2) %>% xtable

####### 50503
Table3.50503.r1 <- df %>% filter(scenario==50503) %>% group_by(cens) %>% summarise(m.ATE.T1.ad = mean(ATE.T1.ad),
                                                                                   med.ATE.T1.ad = mean(med.ATE.T1.ad),
                                                                                   m.ATE.T2.ad = mean(ATE.T2.ad),
                                                                                   med.ATE.T2.ad = mean(med.ATE.T2.ad),
                                                                                   m.ATE.T2.nd = mean(ATE.T2.nd),
                                                                                   med.ATE.T2.nd = mean(med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist


Table3.50503.r2 <- df %>% filter(scenario==50503) %>% group_by(cens) %>% summarise(emp.s.ATE.T1.ad = sd(ATE.T1.ad),
                                                                                   emp.s.med.ATE.T1.ad = sd(med.ATE.T1.ad),
                                                                                   emp.s.ATE.T2.ad = sd(ATE.T2.ad),
                                                                                   emp.s.med.ATE.T2.ad = sd(med.ATE.T2.ad),
                                                                                   emp.s.ATE.T2.nd = sd(ATE.T2.nd),
                                                                                   emp.s.med.ATE.T2.nd = sd(med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist

Table3.50503.r3 <- df %>% filter(scenario==50503) %>% group_by(cens) %>% summarise(est.s.ATE.T1.ad = mean(SE.ATE.T1.ad),
                                                                                   est.s.med.ATE.T1.ad = mean(SE.med.ATE.T1.ad),
                                                                                   est.s.ATE.T2.ad = mean(SE.ATE.T2.ad),
                                                                                   est.s.med.ATE.T2.ad = mean(SE.med.ATE.T2.ad),
                                                                                   est.s.ATE.T2.nd = mean(SE.ATE.T2.nd),
                                                                                   est.s.med.ATE.T2.nd = mean(SE.med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist


Table3.50503.r4 <- df %>% filter(scenario==50503) %>% group_by(cens) %>%
  summarise(cp.ATE.T1.ad = mean(CoverCInorm(est = ATE.T1.ad, se = SE.ATE.T1.ad, truth = true.ATE.T1.ad)),
            cp.med.ATE.T1.ad = mean(CoverCInorm(est = med.ATE.T1.ad, se = SE.med.ATE.T1.ad, truth = true.med.ATE.T1.ad)),
            cp.ATE.T2.ad = mean(CoverCInorm(est = ATE.T2.ad, se = SE.ATE.T2.ad, truth = true.ATE.T2.ad)),
            cp.med.ATE.T2.ad = mean(CoverCInorm(est = med.ATE.T2.ad, se = SE.med.ATE.T2.ad, truth = true.med.ATE.T2.ad)),
            cp.ATE.T2.nd = mean(CoverCInorm(est = ATE.T2.nd, se = SE.ATE.T2.nd, truth = true.ATE.T2.nd)),
            cp.med.ATE.T2.nd = mean(CoverCInorm(est = med.ATE.T2.nd, se = SE.med.ATE.T2.nd, truth = true.med.ATE.T2.nd))) %>%
  select(-cens) %>% unlist


rbind(Table3.50503.r1, Table3.50503.r2, Table3.50503.r3, Table3.50503.r4) %>% round(2) %>% xtable



