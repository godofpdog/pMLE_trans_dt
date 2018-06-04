   #n=nn
   #Q=rep(NA,REP)
   #theta=stheta
   #thetad=sthetad

   #cat("replication:",R,"\n")
   beta=sbeta
   N=0
   z2=sample(1:4,nn,TRUE,c(0.25,0.25,0.25,0.25))
   z1=rbinom(nn,1,0.5)
   Z=matrix(c(z1,z2),nn,2)
   bz=Z%*%beta
   C=T=V=X=D=U=rep(NA,nn)
   TT=rep(NA,N)

   for(i in 1:nn){

      repeat{

         N=N+1
         u=runif(1)
         eps=log(u/(1-u))
         hTT =-bz[i]+eps
         TT=10*exp(hTT)
         TT=R.in(G.in(-log(1-u))/exp(bz[i]))
         #t=-(log(1-u)/exp(-1*z1[i]-1*z2[i]))
         #DD=rexp(1,1/thetad)
         #DD=runif(1,0,theta)
         UU=rexp(1,theta)
         VV=UU+d0
         if(UU<=TT&&TT<=VV) break
      }
   
      T[i]=TT
      V[i]=VV
      #D[i]=DD
      U[i]=UU
   }


   q=(N-n)/N  
   Q[R]=q

   cat("the truncated rate is ",q,"\n")

   k=nn

   Order=order(T)
   nT=rnT=T[Order]
   nU=rnU=U[Order]
   nV=rnV=V[Order]   
   Z1=rZ1=z1[Order]
   Z2=rZ2=z2[Order]
   Z=matrix(c(Z1,Z2),nn,2)