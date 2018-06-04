
####################################################  ## parameter test
# MLE of transformation                            #
# model with doubly truncated data by EM algorithm #
####################################################


rm(list = ls())
library(nleqslv)


############################################  ##設定路徑
                                           #
setwd("F:\\semi_dt_em")   # 
                                           #
############################################



########  transformation functions  ########  ##選擇模型 
                                           #
  md = "log trans"                                #  ("ph","prop odds","log trans","Box Cox")
                                           #
############################################



#############   setup   ####################  ##參數設定
                                           #
  REP=10                                   #
  nn=n=100                                 #
  theta=sp=0.25                            #
  d0=12                                    #
  jack=n                                   #
  sbeta=c(-2,-3)                           #
  coef=1.96                                #
  rho=0.5                                  #
                                           #
############################################


options(digits=6)


##  inverse function of baseline harzard  ##  

R.in=function(x) log(1+x)

##  method of est Gnt  ##  ##IPW or cMLE

Gnt_est = "IPW"


## main program  ##

source(file.path("programs/parameter_test.r"))