   k=nn
   n=nn
   cat("IPW to est Gnt","\n")
   iter=0
   iterem=0
   nbeta=sbeta
   G.en=function(x) 1-exp(-theta*x)
   nTd=ifelse(nT-d0<0,0,nT-d0)
   #nT-d0
   Gnt=pexp(sp,nT)
   Gnd=pexp(sp,nTd)

   repeat{

      iter=iter+1

      cat("iter= ",iter,"\n")

      Fh=function(Gnt,Gnd,dd){
         Fh=rep(NA,n)
         for(i in 1:n){
            id=ifelse(nT<=nU[i]+dd,1,0)
            Fh[i]=1/sum(1/(Gnt-Gnd))*sum(id/(Gnt-Gnd)) 
         }
         return(Fh)
      }

      Gh=function(Fnt,Fnd,dd){
         Gh=rep(NA,n)
         for(i in 1:n){
            id=ifelse(nU<=nT[i]-dd,1,0)
            Gh[i]=1/sum(1/(Fnd-Fnt))*sum(id/(Fnd-Fnt))
         }
         return(Gh)
      }

      Fnt=Fh(Gnt,Gnd,0)
      Fnd=Fh(Gnt,Gnd,d0)
      if(max(abs(Gnt-Gh(Fnt,Fnd,0)))<=0.0001) break
      Gnt=Gh(Fnt,Fnd,0)
      Gnd=Gh(Fnt,Fnd,d0)

   }

   plot(Gnt,type="l",xlab="time points",ylab="CDFs")
   lines(pexp(nT,sp),col=2)
   ln=n*0.65
   legend(ln,0.55,c("est_Gnt","true_Gnt"),col=c(1:2),lty=1)

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

   repeat{

      #dL=dL/sum(dL*exp(-cumsum(dL)))
      iter=iter+1
      bz=Z%*%beta
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

      #for(i in 1:n){
       #  for(j in 1:n){
        #    p[i,j]=dL[j]*exp(bz[i])*g.fn(sum(dL[1:j])*exp(bz[i]))*exp(-G.fn(sum(dL[1:j]*exp(bz[i]))))
         #   w[i,j]=ifelse(  (nT[j]-d0)>=0 ,ifelse(nT[i]==nT[j],1,0)+p[i,j]* (1-G.en(nT[j])+G.en(nT[j]-d0) )/alpha[i], ifelse(nT[i]==nT[j],1,0)+p[i,j]* (1-G.en(nT[j])) /alpha[i] )      #   }
     # }}

      wpi=rowSums(w)
      wpj=colSums(w)
    
      repeat{

      #dL=dL/sum(dL*exp(-cumsum(dL)))

      for(l in 1:n){
         tmpl=matrix(NA,n,n)
         for(j in 1:n){
            LT=(exp(bz))%*%sum(dL[1:j])
            Ai=g.dn(LT)/g.fn(LT)
            tmpl[,j]=w[,j]*exp(bz)*(g.fn(LT)-Ai)*ifelse(nT[j]>=nT[l],1,0)
         }
         dL[l]=wpj[l]/sum(tmpl)
      }

      SB=function(x){
         sb=matrix(NA,n,2)   
         for(i in 1:n){
            LT=cumsum(dL)*exp(Z[i,]%*%x)
            sb[i,]=Z[i,]*sum((1+g.dn(LT)/g.fn(LT)*cumsum(dL)*exp(Z[i,]%*%x)-g.fn(LT)*cumsum(dL)*exp(Z[i,]%*%x))*w[i,])
         }
         return(colSums(sb))
      }

      tmpp=nleqslv(x=beta,SB)$x
      db=max(abs(beta-tmpp))
      beta=tmpp
      #cat("beta=",beta,"db=",db,"\n")
   
      if(iter==200|db<=1e-2) break
      }
      db2=max(abs(beta-nbeta))
      nbeta=beta
      #bz=Z%*%beta
      #dL=dL/sum(dL*exp(-cumsum(dL)))
      #cat("nbeta=",nbeta,"\n")
      if(iterem==200|db2<=1e-2) break



   }
   cat("est beta is",beta,"\n")
   cat("est average un-truncated probability is",mean(alpha),"\n")
   cat(iter,"iterations","\n\n\n")

   nX=nT
   B[R,]=beta
   I[R]=iter
   BZ=c(1,1)%*%beta
   #DL[R,]=dL
   #NX[R,]=nT
   BZ=(c(1,1)%*%beta)[1,1]
   t0.8=R.in(G.in(-log(1-0.2))/exp(BZ));t0.5=R.in(G.in(-log(1-0.5))/exp(BZ));t0.2=R.in(G.in(-log(1-0.8))/exp(BZ)); 
   id0.8=sum(ifelse(nX<=t0.8,1,0));id0.5=sum(ifelse(nX<=t0.5,1,0));id0.2=sum(ifelse(nX<=t0.2,1,0))
   ST0.2[R]=exp(-G.fn(sum(dL[1:id0.2])*exp(BZ)))
   ST0.5[R]=exp(-G.fn(sum(dL[1:id0.5])*exp(BZ)))
   ST0.8[R]=exp(-G.fn(sum(dL[1:id0.8])*exp(BZ)))
   