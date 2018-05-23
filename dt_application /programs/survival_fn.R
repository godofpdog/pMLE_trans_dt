n=length(nX)
time_points=NULL
NN=N1=N0=NULL
nX_z1=nX[which(z1==1)]
nX_z0=nX[which(z1==0)]
nU_z1=nU[which(z1==1)]
nU_z0=nU[which(z1==0)]
n1=length(nX_z1)
n0=length(nX_z0)
i=0
j=1
k=0
repeat{
   i=i+1
   if(i>n) break
   if(i!=n){
      if(nX[i]!=nX[i+1]){
         time_points[j]=nX[i]
         NN[j]=length(which(nX==nX[i]))
         j=j+1
      }
   else{
      time_points[j]=nX[i] 
      NN[j]=length(which(nX==nX[i]))  
   }
   }
}
i=0
j=1
k=0
repeat{
   i=i+1
   if(i>n0) break
   if(i!=n0){
      if(nX_z0[i]!=nX_z0[i+1]){
         #time_points[j]=nX[i]
         N0[j]=length(which(nX_z0==nX_z0[i]))
         j=j+1
      }
   else{
      #time_points[j]=nX[i] 
      N0[j]=length(which(nX_z0==nX_z0[i]))  
   }
   }
}
i=0
j=1
k=0
repeat{
   i=i+1
   if(i>n1) break
   if(i!=n1){
      if(nX_z1[i]!=nX_z1[i+1]){
         #time_points[j]=nX[i]
         N1[j]=length(which(nX_z1==nX_z1[i]))
         j=j+1
      }
   }
   else{
      #time_points[j]=nX[i] 
      N1[j]=length(which(nX_z1==nX_z1[i]))  
   }
}
Ri0=Ri1=NULL
time_set.fn=function(x){
   result=NULL
   j=1
   result[1]=x[1]
   for(i in 2:length(x)){
      if(x[i]!=x[i-1]){
         j=j+1
         result[j]=x[i]
      }
   }
   return(result)
}
Ri0=Ri1=NULL
nX_set_z1=time_set.fn(nX_z1)
nX_set_z0=time_set.fn(nX_z0)
nU_set_z1=time_set.fn(nU_z1)
nU_set_z0=time_set.fn(nU_z0)
for(i in 1:length(nX_set_z1)){
   Ri1[i]=sum(1*(nX_set_z1[i]<=nX_z1)*1*(nU_z1<=nX_set_z1[i])*1*(nX_set_z1[i]<=nX_z1+d0))
}
for(i in 1:length(nX_set_z0)){
   Ri0[i]=sum(1*(nX_set_z0[i]<=nX_z0)*1*(nU_z0<=nX_set_z0[i])*1*(nX_set_z0[i]<=nX_z0+d0))
}
#for(i in 1:length(nX_z1)){
#   Ri1[i]=sum(1*(nX_z1[i]<=nX_z1)*1)
#}
#for(i in 1:length(nX_z0)){
#   Ri0[i]=sum(1*(nX_z0[i]<=nX_z0)*1)
#}

Snt_z0=Snt_z1=Snt_model_z0=Snt_model_z1=rep(1,length(time_points))
ind_z0=ind_z1=ind_all=0
prod_z0=prod_z1=prod_all=NULL
cumsum_ind=0
dL_cumsum=0
for(i in 2:length(time_points)){
   if(time_points[i]%in%nX_z1){
      ind_z1=ind_z1+1
      prod_z1[i]=1-N1[ind_z1]/Ri1[ind_z1]
      Snt_z1[i]=Snt_z1[i-1]*prod_z1[i]
   }
   else{
      Snt_z1[i]=Snt_z1[i-1]
   }
   if(time_points[i]%in%nX_z0){
      ind_z0=ind_z0+1
      prod_z0[i]=1-N0[ind_z0]/Ri0[ind_z0]
      Snt_z0[i]=Snt_z0[i-1]*prod_z0[i]
   }
   else{
      Snt_z0[i]=Snt_z0[i-1]
   }
   #ind_all=ind_all+1
   if(i==1){
      cumsum_ind=NN[i]+NN[i-1]
   }
   else{
      cumsum_ind=NN[i]+cumsum_ind
   }
   dL_cumsum=dL[min(which(time_points[i]==nX))]+dL_cumsum
   Snt_model_z1[i]=exp(-G.fn(dL_cumsum*exp(beta)))
   Snt_model_z0[i]=exp(-G.fn(dL_cumsum))
   Snt_model_z1[i]=exp(-G.fn(sum(dL[1:cumsum_ind])*exp(beta)))
   Snt_model_z0[i]=exp(-G.fn(sum(dL[1:cumsum_ind])*exp(0)))
}

plot(time_points,Snt_model_z0,type="l",ylab="Estimated survival functions",xlab="Age (in months)")
lines(time_points,Snt_model_z1,col=2)
lines(time_points,Snt_z0,col=1,lty=2)
lines(time_points,Snt_z1,col=2,lty=2)
legend(30,0.85,c("PMLE Z=0","PMLE Z=1","NPMLE Z=0","NPMLE Z=1"),col=c(1,2,1,2),lty=c(1,1,2,2))

sup_z1=max(Snt_model_z1-Snt_z1)
sup_z0=max(Snt_model_z0-Snt_z0)
cat("sup(Snt_model_z1-Snt_z1) is ",sup_z1,"\n")
cat("sup(Snt_model_z0-Snt_z0) is ",sup_z0,"\n")
#print(KSZ0)
#print(KSZ1)


