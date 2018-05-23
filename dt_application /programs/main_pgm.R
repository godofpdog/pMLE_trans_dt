source("programs/data.r")
source("programs/models.r")
library(nleqslv)
iterem=0
qt=1/295
n=294
d0=54
p=0.8
theta=4
beta=ini_beta

obs=data[,1]
z1=data[,2]
z2=data[,3]
u=data[,4]
v=data[,5]
tp=data[,6]

bind=cbind(obs,z1,z2,u,v,tp)
T.sort=bind[order(bind[,6]),]

  obs=T.sort[,1]
  z1=T.sort[,2] 
  z2=T.sort[,3]
  u=T.sort[,4]
  v=T.sort[,5]
  tp=T.sort[,6]

  ti=tp

  Gn.Ti=1-exp(-theta*ti)
  Gn.Ti.d0=1-exp(-theta*(ti-d0))

  GG=Gn.Ti-Gn.Ti.d0
  GGP=GG[GG>0]
  n=length(GGP)
 
  obs=obs[1:n]
  z1=z1[1:n] 
  z2=z2[1:n]
  u=u[1:n]
  v=v[1:n]
  tp=tp[1:n]

  ui=u
  ti=tp

  t=-1/theta*log(1-p)

  IT=ifelse(ti<=t,1,0)
  Fn.t=(1/sum(1/(GGP)))*sum(IT/(GGP))
  Fn.t

  Fn.Ui.d0=Fn.Ui=NULL
  IT1=IT2=matrix(0,n,n)

for(i in 1:n)
{

IT1[i,]=ifelse(ti<=v[i],1,0)
Fn.Ui.d0[i]=(1/sum(1/(GGP)))*sum(IT1[i,]/(GGP))
IT2[i,]=ifelse(ti<=u[i],1,0)
Fn.Ui[i]=(1/sum(1/(GGP)))*sum(IT2[i,]/(GGP))

}

IU=ifelse(ui<=t,1,0)
Gn.t=(1/sum(1/(Fn.Ui.d0-Fn.Ui)))*sum(IU/(Fn.Ui.d0-Fn.Ui))
Gn.t

   Gn.Ti=Gn.Ti.d0=NULL
   IU1=IU2=matrix(0,n,n) 
 
for(i in 1:n)
{

IU1[i,]=ifelse(ui<=tp[i],1,0)
Gn.Ti[i]=(1/sum(1/(Fn.Ui.d0-Fn.Ui)))*sum(IU1[i,]/(Fn.Ui.d0-Fn.Ui))
IU2[i,]=ifelse(ui<=tp[i]-d0,1,0)
Gn.Ti.d0[i]=(1/sum(1/(Fn.Ui.d0-Fn.Ui)))*sum(IU2[i,]/(Fn.Ui.d0-Fn.Ui))

}

  
  G=1-exp(-theta*t)

  iter=0

repeat
{
  iter=iter+1

  GG=Gn.Ti-Gn.Ti.d0
  GGP=GG[GG>0]
  n=length(GGP)

 obs=obs[1:n]
 z1=z1[1:n] 
 z2=z2[1:n]
 u=u[1:n]
 v=v[1:n]
 tp=tp[1:n]

 ui=u
 ti=tp


  IT=ifelse(ti<=t,1,0)
  Fn.t=(1/sum(1/(GGP)))*sum(IT/(GGP))

  Fn.Ui.d0=Fn.Ui=NULL
  IT1=IT2=matrix(0,n,n)

  for(i in 1:n)
    {
 
     IT1[i,]=ifelse(ti<=v[i],1,0)
     Fn.Ui.d0[i]=(1/sum(1/(GGP)))*sum(IT1[i,]/(GGP))
     IT2[i,]=ifelse(ti<=u[i],1,0)
     Fn.Ui[i]=(1/sum(1/(GGP)))*sum(IT2[i,]/(GGP))

    }

  IU=ifelse(ui<=t,1,0)
  Gn.t=(1/sum(1/(Fn.Ui.d0-Fn.Ui)))*sum(IU/(Fn.Ui.d0-Fn.Ui))

   Gn.Ti=Gn.Ti.d0=NULL
   IU1=IU2=matrix(0,n,n)


  for(i in 1:n)
    {

     IU1[i,]=ifelse(ui<=tp[i],1,0)
     Gn.Ti[i]=(1/sum(1/(Fn.Ui.d0-Fn.Ui)))*sum(IU1[i,]/(Fn.Ui.d0-Fn.Ui))
     IU2[i,]=ifelse(ui<=tp[i]-d0,1,0)
     Gn.Ti.d0[i]=(1/sum(1/(Fn.Ui.d0-Fn.Ui)))*sum(IU2[i,]/(Fn.Ui.d0-Fn.Ui))

    }

  dev=abs(Gn.t-G)
  if (dev<1/10000|iter>100) break

}

z1=z1[1:n]
nT=nX=nX[1:n]
Z=z1
Gnt=Gn.Ti
cat("estmating beta ...","\n")

iter=0

dL=dLd=rep(NA,n)
dL[1]=nT[1]/10
for(i in 2:n){
   dL[i]=(nT[i]-nT[i-1])/10
}
dLd[1]=(nT[1]+d0-nT[n])/10
for(i in 2:n){
   dLd[i]=(nT[i]-nT[i-1])/10
}

for(i in 2:n){
   dL[i]=exp(nT[i])-exp(nT[i-1])
   dLd[i]=exp(nT[i]+d0)-exp(nT[i-1]+d0)
}
dL[1]=exp(nT[1])
dLd[1]=exp(nT[1]+d0)-exp(d0)
Gnd=Gn.Ti.d0
nbeta=0
repeat{

      #dL=dL/sum(dL*exp(-cumsum(dL)))
      iter=iter+1
      bz=Z*beta
      cat("iter=",iter,"\n")

      G=prod=alpha=rep(NA,n)
      Ltz=Ldz=Szd=matrix(NA,n,n)

      G=Gnt
      Lt=cumsum(dL)
      Ld=cumsum(dLd)
      
      for(i in 1:n){         
         Ltz[i,]=Lt*exp(bz[i])
         Ldz[i,]=Ld*exp(bz[i])  
      }

      Stz=exp(-G.fn(Ltz))
      Sdz=exp(-G.fn(Ldz))
      Sz=Stz-Sdz

      for(i in 1:n){
         Szd[i,1]=1-Sz[i,1]
         prod[1]=G[1]*Szd[i,1]
         for(j in 2:n){
            Szd[i,j]=Sz[i,j-1]-Sz[i,j]
            prod[j]=G[j]*Szd[i,j]
         }
         alpha[i]=sum(prod)
      }

      w=matrix(NA,n,n)
      pij=p=matrix(NA,n,n)

      for(i in 1:n){
         p[i,]=dL*exp(bz[i])*g.fn(cumsum(dL)*exp(bz[i]))*exp(-G.fn(cumsum(dL)*exp(bz[i])))
         w[i,]=1*(nT[i]==nT)+p[i,]*(1-Gnt+Gnd)/alpha[i]
      }

      wpi=rowSums(w)
      wpj=colSums(w)
    
      repeat{

      #dL=dL/sum(dL*exp(-cumsum(dL)))

      for(l in 1:n){
         tmpl=matrix(NA,n,n)
         for(j in 1:n){
            LT=(exp(bz))*sum(dL[1:j])
            Ai=g.dn(LT)/g.fn(LT)
            tmpl[,j]=w[,j]*exp(bz)*(g.fn(LT)-Ai)*ifelse(nT[j]>=nT[l],1,0)
         }
         dL[l]=wpj[l]/sum(tmpl)
      }

      SB=function(x){
         sb=matrix(NA,n,1)   
         for(i in 1:n){
            LT=cumsum(dL)*exp(Z[i]*x)
            sb[i,]=Z[i]*sum((1+g.dn(LT)/g.fn(LT)*cumsum(dL)*exp(Z[i]*x)-g.fn(LT)*cumsum(dL)*exp(Z[i]*x))*w[i,])
         }
         return(sum(sb))
      }

      tmpp=nleqslv(x=beta,SB)$x
      db=(abs(beta-tmpp))
      beta=tmpp
      #cat("beta=",beta,"db=",db,"\n")
   
      if(iter==200|db<=1e-2) break
      }
      db2=(abs(beta-nbeta))
      nbeta=beta
      #bz=Z%*%beta
      #dL=dL/sum(dL*exp(-cumsum(dL)))
      #cat("nbeta=",nbeta,"\n")
      if(iterem==200|db2<=1e-4) break
   }
   cat("est beta is",beta,"\n")
   #cat("est average un-truncated probability is",mean(alpha),"\n")
   cat(iter,"iterations","\n\n\n")
   source("programs/survival_fn.r")


