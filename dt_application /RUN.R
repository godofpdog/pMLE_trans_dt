####################################################
# pMLE of transformation                           #
# model with doubly truncated data by EM algorithm #
# Application : AIDS data                          #
####################################################


rm(list=ls())
library(nleqslv)

## working direction setting

setwd("/Users/liuyi/survival/dt_application")

## model selection    ("ph","prop odds","log trans","Box Cox")

md = "ph"

rho=0.7

ini_beta=0

##  method of est Gnt  ##  ##IPW or cMLE

Gnt_est = "cMLE"


## execution

options(digits=6)

source("programs/main_pgm.r")

   