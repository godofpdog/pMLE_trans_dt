B=matrix(NA,REP,2)
I=rep(NA,REP)
ST0.2=ST0.5=ST0.8=rep(NA,REP)
surv_ft=matrix(c(ST0.2,ST0.5,ST0.8),REP,3)
colnames(surv_ft)=c("S0.2","S0.5","S0.8")
DL=NX=matrix(NA,REP,n)
B=JB=matrix(NA,REP,2)
JB=matrix(NA,jack,2)
I=rep(NA,REP)
DL=NX=matrix(NA,REP,n)
nZ=matrix(NA,REP,2)
Q=rep(NA,REP)
ST0.2J=ST0.5J=ST0.8J=rep(NA,jack)
JsdB1=JsdB2=Jsd0.2=Jsd0.5=Jsd0.8=rep(NA,REP)
surv_ft_jack=matrix(c(ST0.2J,ST0.5J,ST0.8J),REP,3)
colnames(surv_ft_jack)=c("S0.2J","S0.5J","S0.8J")
beta_range1=matrix(NA,REP,2)
beta_range2=matrix(NA,REP,2)
st0.2_range=st0.5_range=st0.8_range=matrix(NA,REP,2)

length_beta=length(sbeta)
nbeta_range1=nbeta_range2=matrix(NA,REP,2)
nst0.2_range=nst0.5_range=nst0.8_range=matrix(NA,REP,2)


#initial value

beta_cover_cnt1=beta_cover_cnt2=st0.2_cover_cnt=st0.5_cover_cnt=st0.8_cover_cnt=0