
##Function to Compute convolution density of Exponentials with varying rate -------------
## Input: Vector of Rate parameters
## Output: Value of density function

dexpconv=function(lambdas,t){
  library(combinat)
  lambdas = sort(lambdas)
  ulambdas = unique(lambdas)
  r = length(lambdas)
  n = length(ulambdas)
  ki = as.numeric(table(lambdas))
  
  
  sum(ulambdas^ki*exp(-t*ulambdas)* 
    sapply(
      sapply(1:n,function(i){
        sapply(1:ki[i],function(y){
          (((-1)^(ki[i]-y))/factorial(y-1))*t^(y-1)*
            sum(
              xsimplex(n-1,ki[i]-y,function(x){
                prod(choose(ki[-i]+x-1,x)*(ulambdas[-i]^(ki[-i])/(ulambdas[-i]-ulambdas[i])^(ki[-i]+x)))
              })
            )
        })
      })
      ,sum
    )
  )
}

qexpconv=function(p,lambdas){
  uniroot(function(x){
    integrate(Vectorize(function(y){
      dexpconv(lambdas,y)
      }),0,x)$value-p
    },interval=c(0,100000))$root
}

qmexpconv=function(p,lambdas,prob){
  ##Supply lambda and probs in lists
  uniroot(function(x){
    integrate(Vectorize(function(y){
      sum(sapply(lambdas,function(z){dexpconv(z,y)})*(sapply(prob,prod)/sum(sapply(prob,prod))))
    }),0,x)$value-p
  },interval=c(0,100000),
  extendInt = "yes")$root
}

