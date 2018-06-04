

source(file.path("programs/space.r"))

for(R in 1:1000){

   source(file.path("programs/data_gen.r"))

}


cat("the average truncated rate is : ",mean(Q))
   