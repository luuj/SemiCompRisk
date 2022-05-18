########################################################################################################
# Table summarizing beta estimation for "Causal Inference for Semi-competing risks data" by Nevo and Gorfine
########################################################################################################
library(dplyr)
library(xtable)
library(Daniel)
rm(list = ls())

setwd("/Users/danielnevo/Dropbox/CausalSemiComp/Sims/SemiParam/Summaries/AllResults/")

#### Load a single table with simulation results and true values
# See AllCombine.R
df <- read.csv("AllSimResults.csv")

#### Scenario 51001

Table3.51001.r1 <- df %>% filter(scenario==51001) %>% group_by(cens) %>% summarise(m.ATE.T1.ad = mean(ATE.T1.ad),
                                                                med.ATE.T1.ad = mean(med.ATE.T1.ad),
                                                                m.ATE.T2.ad = mean(ATE.T2.ad),
                                                                med.ATE.T2.ad = mean(med.ATE.T2.ad),
                                                                m.ATE.T2.nd = mean(ATE.T2.nd),
                                                                med.ATE.T2.nd = mean(med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist


Table3.51001.r2 <- df %>% filter(scenario==51001) %>% group_by(cens) %>% summarise(emp.s.ATE.T1.ad = sd(ATE.T1.ad),
                                                                emp.s.med.ATE.T1.ad = sd(med.ATE.T1.ad),
                                                                emp.s.ATE.T2.ad = sd(ATE.T2.ad),
                                                                emp.s.med.ATE.T2.ad = sd(med.ATE.T2.ad),
                                                                emp.s.ATE.T2.nd = sd(ATE.T2.nd),
                                                                emp.s.med.ATE.T2.nd = sd(med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist

Table3.51001.r3 <- df %>% filter(scenario==51001) %>% group_by(cens) %>% summarise(est.s.ATE.T1.ad = mean(SE.ATE.T1.ad),
                                                                                   est.s.med.ATE.T1.ad = mean(SE.med.ATE.T1.ad),
                                                                                   est.s.ATE.T2.ad = mean(SE.ATE.T2.ad),
                                                                                   est.s.med.ATE.T2.ad = mean(SE.med.ATE.T2.ad),
                                                                                   est.s.ATE.T2.nd = mean(SE.ATE.T2.nd),
                                                                                   est.s.med.ATE.T2.nd = mean(SE.med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist


Table3.51001.r4 <- df %>% filter(scenario==51001) %>% group_by(cens) %>%
  summarise(cp.ATE.T1.ad = mean(CoverCInorm(est = ATE.T1.ad, se = SE.ATE.T1.ad, truth = true.ATE.T1.ad)),
            cp.med.ATE.T1.ad = mean(CoverCInorm(est = med.ATE.T1.ad, se = SE.med.ATE.T1.ad, truth = true.med.ATE.T1.ad)),
            cp.ATE.T2.ad = mean(CoverCInorm(est = ATE.T2.ad, se = SE.ATE.T2.ad, truth = true.ATE.T2.ad)),
            cp.med.ATE.T2.ad = mean(CoverCInorm(est = med.ATE.T2.ad, se = SE.med.ATE.T2.ad, truth = true.med.ATE.T2.ad)),
            cp.ATE.T2.nd = mean(CoverCInorm(est = ATE.T2.nd, se = SE.ATE.T2.nd, truth = true.ATE.T2.nd)),
            cp.med.ATE.T2.nd = mean(CoverCInorm(est = med.ATE.T2.nd, se = SE.med.ATE.T2.nd, truth = true.med.ATE.T2.nd))) %>%
  select(-cens) %>% unlist


rbind(Table3.51001.r1, Table3.51001.r2, Table3.51001.r3, Table3.51001.r4)  %>% round(2) %>% xtable



####### 51002
Table3.51002.r1 <- df %>% filter(scenario==51002) %>% group_by(cens) %>% summarise(m.ATE.T1.ad = mean(ATE.T1.ad),
                                                                                   med.ATE.T1.ad = mean(med.ATE.T1.ad),
                                                                                   m.ATE.T2.ad = mean(ATE.T2.ad),
                                                                                   med.ATE.T2.ad = mean(med.ATE.T2.ad),
                                                                                   m.ATE.T2.nd = mean(ATE.T2.nd),
                                                                                   med.ATE.T2.nd = mean(med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist


Table3.51002.r2 <- df %>% filter(scenario==51002) %>% group_by(cens) %>% summarise(emp.s.ATE.T1.ad = sd(ATE.T1.ad),
                                                                                   emp.s.med.ATE.T1.ad = sd(med.ATE.T1.ad),
                                                                                   emp.s.ATE.T2.ad = sd(ATE.T2.ad),
                                                                                   emp.s.med.ATE.T2.ad = sd(med.ATE.T2.ad),
                                                                                   emp.s.ATE.T2.nd = sd(ATE.T2.nd),
                                                                                   emp.s.med.ATE.T2.nd = sd(med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist

Table3.51002.r3 <- df %>% filter(scenario==51002) %>% group_by(cens) %>% summarise(est.s.ATE.T1.ad = mean(SE.ATE.T1.ad),
                                                                                   est.s.med.ATE.T1.ad = mean(SE.med.ATE.T1.ad),
                                                                                   est.s.ATE.T2.ad = mean(SE.ATE.T2.ad),
                                                                                   est.s.med.ATE.T2.ad = mean(SE.med.ATE.T2.ad),
                                                                                   est.s.ATE.T2.nd = mean(SE.ATE.T2.nd),
                                                                                   est.s.med.ATE.T2.nd = mean(SE.med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist


Table3.51002.r4 <- df %>% filter(scenario==51002) %>% group_by(cens) %>%
  summarise(cp.ATE.T1.ad = mean(CoverCInorm(est = ATE.T1.ad, se = SE.ATE.T1.ad, truth = true.ATE.T1.ad)),
            cp.med.ATE.T1.ad = mean(CoverCInorm(est = med.ATE.T1.ad, se = SE.med.ATE.T1.ad, truth = true.med.ATE.T1.ad)),
            cp.ATE.T2.ad = mean(CoverCInorm(est = ATE.T2.ad, se = SE.ATE.T2.ad, truth = true.ATE.T2.ad)),
            cp.med.ATE.T2.ad = mean(CoverCInorm(est = med.ATE.T2.ad, se = SE.med.ATE.T2.ad, truth = true.med.ATE.T2.ad)),
            cp.ATE.T2.nd = mean(CoverCInorm(est = ATE.T2.nd, se = SE.ATE.T2.nd, truth = true.ATE.T2.nd)),
            cp.med.ATE.T2.nd = mean(CoverCInorm(est = med.ATE.T2.nd, se = SE.med.ATE.T2.nd, truth = true.med.ATE.T2.nd))) %>%
  select(-cens) %>% unlist


rbind(Table3.51002.r1, Table3.51002.r2, Table3.51002.r3, Table3.51002.r4) %>% round(2) %>% xtable

####### 51003
Table3.51003.r1 <- df %>% filter(scenario==51003) %>% group_by(cens) %>% summarise(m.ATE.T1.ad = mean(ATE.T1.ad),
                                                                                   med.ATE.T1.ad = mean(med.ATE.T1.ad),
                                                                                   m.ATE.T2.ad = mean(ATE.T2.ad),
                                                                                   med.ATE.T2.ad = mean(med.ATE.T2.ad),
                                                                                   m.ATE.T2.nd = mean(ATE.T2.nd),
                                                                                   med.ATE.T2.nd = mean(med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist


Table3.51003.r2 <- df %>% filter(scenario==51003) %>% group_by(cens) %>% summarise(emp.s.ATE.T1.ad = sd(ATE.T1.ad),
                                                                                   emp.s.med.ATE.T1.ad = sd(med.ATE.T1.ad),
                                                                                   emp.s.ATE.T2.ad = sd(ATE.T2.ad),
                                                                                   emp.s.med.ATE.T2.ad = sd(med.ATE.T2.ad),
                                                                                   emp.s.ATE.T2.nd = sd(ATE.T2.nd),
                                                                                   emp.s.med.ATE.T2.nd = sd(med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist

Table3.51003.r3 <- df %>% filter(scenario==51003) %>% group_by(cens) %>% summarise(est.s.ATE.T1.ad = mean(SE.ATE.T1.ad),
                                                                                   est.s.med.ATE.T1.ad = mean(SE.med.ATE.T1.ad),
                                                                                   est.s.ATE.T2.ad = mean(SE.ATE.T2.ad),
                                                                                   est.s.med.ATE.T2.ad = mean(SE.med.ATE.T2.ad),
                                                                                   est.s.ATE.T2.nd = mean(SE.ATE.T2.nd),
                                                                                   est.s.med.ATE.T2.nd = mean(SE.med.ATE.T2.nd)) %>%
  select(-cens) %>% unlist


Table3.51003.r4 <- df %>% filter(scenario==51003) %>% group_by(cens) %>%
  summarise(cp.ATE.T1.ad = mean(CoverCInorm(est = ATE.T1.ad, se = SE.ATE.T1.ad, truth = true.ATE.T1.ad)),
            cp.med.ATE.T1.ad = mean(CoverCInorm(est = med.ATE.T1.ad, se = SE.med.ATE.T1.ad, truth = true.med.ATE.T1.ad)),
            cp.ATE.T2.ad = mean(CoverCInorm(est = ATE.T2.ad, se = SE.ATE.T2.ad, truth = true.ATE.T2.ad)),
            cp.med.ATE.T2.ad = mean(CoverCInorm(est = med.ATE.T2.ad, se = SE.med.ATE.T2.ad, truth = true.med.ATE.T2.ad)),
            cp.ATE.T2.nd = mean(CoverCInorm(est = ATE.T2.nd, se = SE.ATE.T2.nd, truth = true.ATE.T2.nd)),
            cp.med.ATE.T2.nd = mean(CoverCInorm(est = med.ATE.T2.nd, se = SE.med.ATE.T2.nd, truth = true.med.ATE.T2.nd))) %>%
  select(-cens) %>% unlist


rbind(Table3.51003.r1, Table3.51003.r2, Table3.51003.r3, Table3.51003.r4) %>% round(2) %>% xtable



