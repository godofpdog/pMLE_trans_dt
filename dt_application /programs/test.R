nX_z1=nX[which(z1==1)]
nX_z0=nX[which(z1==0)]
nU_z1=nU[which(z1==1)]
nU_z0=nU[which(z1==0)]

Ri0=Ri1=NULL

for(i in 1:n){
   Ri1[i]=sum(1*(nX_z1[i]<=nX_z1)*1*(nU_z1<=nX_z1[i])*1*(nX_z1[i]<=nX_z1+d0))
}
for(i in 1:n){
   Ri0[i]=sum(1*(nX_z0[i]<=nX_z0)*1*(nU_z0<=nX_z0[i])*1*(nX_z0[i]<=nX_z0+d0))
}

Fnt1=cumprod(1-1/Ri1)
Fnt0=cumprod(1-1/Ri0)

RR=NULL
for(i in 1:n){
   RR[i]=sum(1*(nX[i]<=nX))
}
Fnt=cumprod(1-1/RR)


plot(Fnt)
dL=lambda
lines(exp(-cumsum(G.fn(dL*exp(1.67375)))),col=2)
lines(exp(-cumsum(G.fn(dL*exp(0)))),col=3)
lines(exp(-cumsum(dL)))