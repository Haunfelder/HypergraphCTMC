library(dplyr)
library(ggplot2)
library(igraph)
library(lubridate)
library(plyr)
library(dplyr)
library(Matrix)
library(knitr)
library(xtable)
library(doParallel)
library(parallel)
library(foreach)
library(combinat)
library(data.table)

setwd("C:/Users/haunf/Documents/GitHub/HypergraphCTMC/Code")

##Sim 1 Parameters
lambda_true = c(1/250, 1/200, 1/200, 1/300, 1/300, 1/150 , 1/150)
names(lambda_true) = c("lambda00","lambda10","lambda20","lambda30","lambda01","lambda11","lambda21")

trans_probs_true = c(0.45, 0.45, 0.2, 0.45, 0.25, 0.3, 0.35, 0.35, 0.2, 0.3)
names(trans_probs_true) = c("p01|00","p10|00","p11|10","p20|10","p11|01",
                          "p21|11","p31|21","p31|30",
                          "p21|20","p30|20")
pars1 = c(lambda_true,trans_probs_true)



##Sim 2 Parameters
lambda_true = c(1/500, 1/100, 1/100, 1/150, 1/150, 1/200 , 1/200)
names(lambda_true) = c("lambda00","lambda10","lambda20","lambda30","lambda01","lambda11","lambda21")

trans_probs_true = c(0.10, 0.30, 0.35, 0.6, 0.9, 0.9, 0.8, 0.8, 0.25, 0.65)
names(trans_probs_true) = c("p01|00","p10|00","p11|10","p20|10","p11|01",
                          "p21|11","p31|21","p31|30",
                          "p21|20","p30|20")
pars2 = c(lambda_true,trans_probs_true)


#Set Parameters (Simulation 3)
lambda_true = c(1/200, 1/100, 1/100, 1/150, 1/400, 1/450 , 1/450)
names(lambda_true) = c("lambda00","lambda10","lambda20","lambda30","lambda01","lambda11","lambda21")

trans_probs_true = c(0.10, 0.4, 0.05, 0.7, 0.5, 0.5, 0.5, 0.8, 0.05, 0.75)
names(trans_probs_true) = c("p01|00","p10|00","p11|10","p20|10","p11|01",
                          "p21|11","p31|21","p31|30",
                          "p21|20","p30|20")
pars3 = c(lambda_true,trans_probs_true)

#Set Parameters (Simulation 4)
lambda_true = c(1/200, 1/400, 1/400, 1/450, 1/100, 1/150 , 1/150)
names(lambda_true) = c("lambda00","lambda10","lambda20","lambda30","lambda01","lambda11","lambda21")

trans_probs_true = c(0.60, 0.2, 0.6, 0.1, 0.7, 0.7, 0.7, 0.6, 0.5, 0.4)
names(trans_probs_true) = c("p01|00","p10|00","p11|10","p20|10","p11|01",
                          "p21|11","p31|21","p31|30",
                          "p21|20","p30|20")
pars4 = c(lambda_true,trans_probs_true)






nsim = 200


em_lambda = matrix(0,nrow=nsim,ncol=length(lambda_true))
em_trans_probs = matrix(0,nrow=nsim,ncol=length(trans_probs_true))




source("expconv.R")

quantiles=c(0.50,0.75,0.90,0.95)
bigt=c(sapply(quantiles,function(x){qexpconv(x,pars[c("lambda00","lambda10")])}),
       sapply(quantiles,function(x){qexpconv(x,pars[c("lambda00","lambda01")])}),
       sapply(quantiles,function(x){qexpconv(x,pars[c("lambda00","lambda10","lambda20")])}),
       sapply(quantiles,function(x){qexpconv(x,pars[c("lambda00","lambda10","lambda20","lambda30")])}),
       sapply(quantiles,function(x){qmexpconv(x, list(pars[c("lambda00","lambda01","lambda11")],pars[c("lambda00","lambda10","lambda11")]), list(pars[c("p01|00","p11|01")],pars[c("p10|00","p11|10")]))}),
       sapply(quantiles,function(x){qmexpconv(x, list(pars[c("lambda00","lambda01","lambda11","lambda21")],
                                                      pars[c("lambda00","lambda10","lambda11","lambda21")],
                                                      pars[c("lambda00","lambda10","lambda20","lambda21")]),
                                              list(pars[c("p01|00","p11|01","p21|11")],
                                                   pars[c("p10|00","p11|10","p21|11")],
                                                   pars[c("p10|00","p20|10","p21|20")]))
       })
)
q = length(quantiles)
bigtdf = data.frame(dyads=rep(c(1,0,2,3,1,2),each=q),
                    triads=rep(c(0,1,0,0,1,1),each=q), 
                    quantiles=rep(quantiles,6),
                    bigt=bigt)

n = seq(500,2000,300)
tandn = expand.grid(bigt, n)



#library(doSNOW)
cl <- makeCluster(3)
#registerDoSNOW(cl)
clusterExport(cl,varlist=c("triad.sim","pars","bigt","tandn"))
clusterEvalQ(cl,library(data.table))
clusterEvalQ(cl,library(dplyr))
registerDoParallel(cl)

#pb <- txtProgressBar(max = dim(tandn)[1], style = 3)
#progress <- function(n) setTxtProgressBar(pb, n)
#opts <- list(progress = progress)

sim_results=foreach (i = 1:dim(tandn)[1],.packages=c("dplyr","data.table")) %:%
  foreach( j = 1:nsim,.combine='rbind',.packages=c("dplyr","data.table")) %dopar%{
    ctmc_sim = rbindlist(lapply(1:tandn[i,2],function(x){data.frame(triad_sim(pars,tandn[i,1]),TrioID=x)}))
    #ctmc.sim=rbindlist(parLapplyLB(cl,1:1000,function(x){data.frame(triad.sim(pars,bigt[i]),TrioID=x)}))
    colnames(ctmc_sim) = c("dyads","triads","TimeBetween","DateSubmitted","TrioID")
    ctmc_sim = data.frame(ctmc_sim)
    
    both=EM_alg(ctmc_sim)
    c(both[[1]], both[[2]])
  }


sim_results=list()
for(i in 1:dim(tandn)[1]){
  sim_results[[i]] = matrix(0,nrow=nsim,ncol=length(pars))
  for(j in 1:nsim){
    clusterExport(cl,varlist=c("i","j","tandn","triad.sim"))
    #ctmc.sim=rbindlist(lapply(1:tandn[i,2],function(x){data.frame(triad.sim(pars,tandn[i,1]),TrioID=x)}))
    ctmc_sim = rbindlist(parLapplyLB(cl,1:tandn[i,2],function(x){data.frame(triad.sim(pars,tandn[i,1]),TrioID=x)}))
    colnames(ctmc_sim) = c("dyads","triads","TimeBetween","DateSubmitted","TrioID")
    ctmc_sim = data.frame(ctmc_sim)
    
    both = EM_alg(ctmc_sim)
    sim_results[[i]][j,] = c(both[[1]],  both[[2]])
  }
  print((i/dim(tandn)[1])*100)
}



#close(pb)
#stopCluster(cl) 

##Compute MSE
pars = c(lambda_true,trans_probs_true)
mses = lapply(sim_results,function(x){
  apply(
    sweep(x,2,pars[match(colnames(x),names(pars))],"-"),2,function(y){
      mean(y[y<Inf & y>-Inf]^2,na.rm=T)
    })
})

##Construct dataframe of mse results
df = data.frame(T=tandn,Reduce(rbind,mses))


##Make Plot
df_lambda = melt(df[,c(1,grep("lambda",colnames(df)))],id="T")
df_trans_probs = melt(df[,c(1,grep("\\.",colnames(df)))],id="T")
df_trans_probs$variable=gsub("\\.","\\|",df_trans_probs$variable)

ggplot(df_lambda, aes(x=T,y=log(value),group=variable,colour=variable))+
  geom_line() +
  #coord_cartesian(ylim=c(0,0.0025))+
  ylab("log MSE")+
  theme(legend.title=element_blank())

ggplot(df_trans_probs, aes(x=T,y=log(value),group=variable,colour=variable))+
  geom_line()+
  #coord_cartesian(ylim=c(0,0.0025))+
  ylab("log MSE")+
  theme(legend.title=element_blank())

##Surface as a function of t and n
df_lambda=melt(df[,c(1,grep("lambda",colnames(df)))],id="T")
df_trans_probs=melt(df[,c(1,grep("\\.",colnames(df)))],id="T")
df_trans_probs$variable=gsub("\\.","\\|",df_trans_probs$variable)

ggplot(df_lambda,aes(x=bigt,y=n,fill=log(value)))+
  geom_raster()+
  scale_fill_gradient(low="red",high="yellow")



stat_density2d(aes(fill=..level..),geom="polygon")
