





##  result  ##

mean_beta1=mean(B[,1])
mean_beta2=mean(B[,2])
bias_beta1=mean_beta1-sbeta[1]
bias_beta2=mean_beta2-sbeta[2]
sd_beta1=sd(B[,1])
sd_beta2=sd(B[,2])
r1=(bias_beta1^2+sd_beta1^2)^0.5
r2=(bias_beta2^2+sd_beta2^2)^0.5
mean_ST0.2=mean(ST0.2)
mean_ST0.5=mean(ST0.5)
mean_ST0.8=mean(ST0.8)
bias_ST0.2=mean_ST0.2-0.2
bias_ST0.5=mean_ST0.5-0.5
bias_ST0.8=mean_ST0.8-0.8
sd_ST0.2=sd(ST0.2)
sd_ST0.5=sd(ST0.5)
sd_ST0.8=sd(ST0.8)
rs2=(bias_ST0.2^2+sd_ST0.2^2)^0.5
rs5=(bias_ST0.5^2+sd_ST0.5^2)^0.5
rs8=(bias_ST0.8^2+sd_ST0.8^2)^0.5
b1=c(mean_beta1,bias_beta1,sd_beta1,r1)
b2=c(mean_beta2,bias_beta2,sd_beta2,r2)
s2=c(mean_ST0.2,bias_ST0.2,sd_ST0.2,rs2)
s5=c(mean_ST0.5,bias_ST0.5,sd_ST0.5,rs5)
s8=c(mean_ST0.8,bias_ST0.8,sd_ST0.8,rs8)

cat("## summary ##","\n\n")
cat(model,"\n")
cat("the average truncated rate is : ",mean(Q),"\n")
cat("the average iterations is :",mean(I),"\n")
cat("the sample size is : ",n,"\n")
cat("use",Gnt_est,"to est Gnt","\n\n")

summary=matrix(c(b1,b2,s2,s5,s8),4,5)
rownames(summary)=c("mean","bias","sd","rmse")
colnames(summary)=c("beta1","beta2","st0.2","st0.5","st0.8")
print(summary)
write.table(summary,file.path("output/summary.txt"),quote=F)
cat("\n\n\n")
cat("t0.2 is",t0.2,"\n")
cat("t0.5 is",t0.5,"\n")
cat("t0.8 is",t0.8,"\n")