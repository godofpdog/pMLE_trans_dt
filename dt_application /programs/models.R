link.fn=function(lfn){
   if(lfn=="ph"){
      model="Cox ph model"
      G.fn=function(x) x
      G.in=function(x) x
      g.fn=function(x) 1
      g.dn=function(x) 0
   }
   if(lfn=="prop odds"){
      model="proportional odds model"
      G.fn=function(x) log(1+x)
      G.in=function(x) exp(x)-1    
      g.fn=function(x) 1/(1+x)
      g.dn=function(x) 1/(1+x)^2*(-1)
   }
   if(lfn=="log trans"){
      model="logarithm transformation model"
      G.fn=function(x) log(1+rho*x)/rho     
      G.in=function(x) (exp(rho*x)-1)/rho 
      g.fn=function(x) 1/(1+rho*x)
      g.dn=function(x) rho/(1+rho*x)^2*(-1)
   }
   if(lfn=="Box Cox"){
      model="Box-Cox transformation model"
      G.fn=function(x) ((1+x)^rho-1)/rho  
      G.in=function(x) (1+rho*x)^(1/rho)-1 
      g.fn=function(x) (1+x)^(rho-1)
      g.dn=function(x) (rho-1)*(1+x)^(rho-2)
   }
   list(model=model,G.fn=G.fn,G.in=G.in,g.fn=g.fn,g.dn=g.dn)
}

G.fn=link.fn(paste(md))$G.fn
g.fn=link.fn(paste(md))$g.fn
g.dn=link.fn(paste(md))$g.dn
G.in=link.fn(paste(md))$G.in
model=link.fn(paste(md))$model

cover_fn=function(x,lower_bound,upper_bound){
   result=rep(0,length(x))
   for(i in 1:length(x)){
      if(x[i]>=lower_bound && x<=upper_bound) result[i]=1
   }
   return(result)
}
lol=rho