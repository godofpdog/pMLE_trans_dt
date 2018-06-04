
####################################################
# pMLE of transformation                           #
# model with doubly truncated data by EM algorithm #
####################################################


rm(list = ls())
library(nleqslv)


############################################ 
                                           #
setwd("/Users/liuyi/survival/simulation/") # 
                                           #
############################################



########  transformation functions  ######## 
                                           #
  md = "log trans"                         #  ("ph","prop odds","log trans","Box Cox")
                                           #
############################################



#############   setup   ####################  
                                           #
  REP=100                                  #
  nn=n=100                                 #
  stheta=0.25                              #
  d0=12                                    #
  jack=n                                   #                                   
  sbeta=c(-2,-3)                           #                               
  rho=0.5                                  #
                                           #
############################################


options(digits=6)


##  inverse function of baseline harzard  ##  

R.in=function(x) log(1+x)

##  method of est Gnt  ##  ##IPW or cMLE

Gnt_est = "cMLE"

## main program  ##

source("programs/main_pgm.r")
