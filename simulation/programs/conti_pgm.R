## conti pgm ##

R=read.table(file.path("output/replications.txt"))[1,1]

result_tmp=read.table(file.path("output/result_tmp.txt"))

B=matrix(c(result_tmp[,1],result_tmp[,2]),REP,length_beta)

ST0.2=result_tmp[,3]
ST0.5=result_tmp[,4]
ST0.8=result_tmp[,5]

while(R<=REP){

theta=sp=stheta

cat("replication : ",R,"\n")

iter_rep=data.frame(iter=R)

write.table(iter_rep,file.path("output/replications.txt"),quote=F)

source(file.path("programs/data_gen.r"))

if(Gnt_est=="IPW") source(file.path("programs/ipw_em.r")) 

if(Gnt_est=="cMLE") source(file.path("programs/cmle_em.r"))

#source(file.path("programs/em.r"))

result_tmp=data.frame(beta1=B[,1],beta2=B[,2],st0.2=ST0.2,st0.5=ST0.5,st0.8=ST0.8)

write.table(result_tmp,file.path("output/result_tmp.txt"),quote=F)

R=R+1

}

if(R==(REP+1)) source(file.path("programs/result.r"))
