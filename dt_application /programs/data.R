



data=read.table(file.path("AIDS.txt"),header=T)
d0=54
n=294
sp=theta=10
sbeta=beta=0.5
obs=data[,1]
z1=data[,2]
z2=data[,3]
U=data[,4]
V=data[,5]
T=data[,6]
bind=cbind(obs,z1,z2,U,V,T)
data_sort=bind[order(bind[,6]),]
obs=data_sort[,1]
z1=data_sort[,2]
z2=data_sort[,3]
nU=data_sort[,4]
nV=data_sort[,5]
nT=nX=data_sort[,6]



