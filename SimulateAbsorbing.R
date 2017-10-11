

## Simulate Absorbing Model --------------------------------------------------------------

#lambda_true=rep(1/500,7)
lambda_true=rep(500,7)
names(lambda_true)=c("lambda00","lambda10","lambda20","lambda30","lambda01","lambda11","lambda21")

trans_probs_true=c(rep(1/3,4),rep(0.5,4),rep(1/3,2))
names(trans_probs_true)=c("p01|00","p10|00","p11|10","p20|10","p11|01",
                     "p21|11","p31|21","p31|30",
                     "p21|20","p30|20")
pars=c(lambda_true,trans_probs_true)
bigT=2000
n=10000


triad.sim=function(pars,bigT){
  pars["p00|00"]=1-pars[["p10|00"]]-pars[["p01|00"]]
  pars["p10|10"]=1-pars[["p20|10"]]-pars[["p11|10"]]
  pars["p01|01"]=1-pars[["p11|01"]]
  pars["p11|11"]=1-pars[["p21|11"]]
  pars["p21|21"]=1-pars[["p31|21"]]
  pars["p30|30"]=1-pars[["p31|30"]]
  pars["p20|20"]=1-pars[["p30|20"]]-pars[["p21|20"]]
  cumtime=0
  timebetween=NULL
  dyads=0
  triads=0
  i=1
  while(cumtime[i]<bigT & !(dyads[i]==3 & triads[i]==1)){
    timebetween=c(timebetween,rexp(1,1/pars[paste("lambda",dyads[i],triads[i],sep="")]))
    cumtime=c(cumtime,sum(timebetween))
    ind.p=which(gsub("p(\\w+)\\|(\\w+)","\\2",names(pars)) %in% paste(dyads[i],triads[i],sep=""))
    new.motif=gsub("p(\\w+)\\|(\\w+)","\\1",names(pars))[ind.p][sample(1:length(ind.p),1,prob=pars[ind.p])]
    dyads=c(dyads,strsplit(new.motif,split="")[[1]][1])
    triads=c(triads,strsplit(new.motif,split="")[[1]][2])
    i=i+1
    if(dyads[i]==dyads[i-1] & triads[i]==triads[i-1]){cumtime[i]=Inf}
  }
  
  if(dyads[i]==3 & triads[i]==1 & (cumtime[i]+timebetween[i-1])<bigT){
    timebetween[i]=bigT-sum(timebetween[1:(i-1)])
    cumtime[i+1]=bigT
    return(data.frame(dyads=dyads,triads=triads,TimeBetween=timebetween,cumtime=cumtime[-1]))
  }else{
    if(i>2){
      timebetween[i-1]=bigT-sum(timebetween[1:(i-2)])
    }else{
      timebetween[i-1]=bigT
    }
    return(data.frame(dyads=dyads[-i],triads=triads[-i],TimeBetween=timebetween[-i],cumtime=cumtime[-1]))
  }
}

triad.sim(pars,bigT)



library(data.table)
ctmc.sim=rbindlist(lapply(1:1000,function(x){data.frame(triad.sim(pars,40000),TrioID=x)}))
colnames(ctmc.sim)=c("dyads","triads","TimeBetween","DateSubmitted","TrioID")
ctmc.sim=data.frame(ctmc.sim)

ctmc.sim %>% group_by(dyads,triads) %>% tally()


## Estimating Parameters from Simulated Data ------------------------------------------

EM_alg=function(ctmc.sim){
  
  niter=100
  lambda=matrix(0,nrow=niter,ncol=7)
  lambda[1,]=rep(10,7)
  colnames(lambda)=c("lambda00","lambda01","lambda10","lambda20","lambda30","lambda11","lambda21")
  
  trans_probs=matrix(0,nrow=niter,ncol=10)
  trans_probs[1,]=rep(0.5,10)
  colnames(trans_probs)=c("p01|00","p10|00","p11|10",
                          "p20|10","p11|01","p21|11",
                          "p31|21","p31|30","p21|20",
                          "p30|20")
  
  lamb_temp=ctmc.sim %>%
    dplyr::group_by(TrioID) %>%
    dplyr::filter(DateSubmitted==max(DateSubmitted)) %>%
    #dplyr::group_by(dyads,triads) %>%
    dplyr::summarize(ind.l=match(paste("lambda",dyads,triads,sep=""),colnames(lambda)),
                     TimeBetween,DateSubmitted,
                     dyads,triads) %>%
    dplyr::filter(!is.na(ind.l))
  
  tal_temp=ctmc.sim %>% 
    dplyr::group_by(TrioID) %>%
    dplyr::filter(DateSubmitted<max(DateSubmitted)) %>%
    dplyr::group_by(dyads,triads) %>% 
    dplyr::tally() %>%
    dplyr::rename(n.trans=n)
  
  max_tal_temp=ctmc.sim %>%
    dplyr::group_by(TrioID) %>%
    dplyr::filter(DateSubmitted==max(DateSubmitted)) %>%
    dplyr::group_by(dyads,triads) %>%
    tally() %>%
    dplyr::rename(n.cens=n)
  
  sum_temp=ctmc.sim %>% 
    dplyr::group_by(TrioID) %>%
    dplyr::filter(DateSubmitted<max(DateSubmitted)) %>%
    dplyr::group_by(dyads,triads) %>%
    dplyr::summarize(sumtime=sum(as.numeric(TimeBetween))) %>%
    filter(!(dyads==3 & triads==1))
  
  trans_temp=ctmc.sim %>%
    dplyr::group_by(TrioID) %>%
    dplyr::mutate(transition=paste(lead(dyads),lead(triads),dyads,triads,sep="")) %>%
    dplyr::filter(DateSubmitted<max(DateSubmitted)) %>%
    dplyr::group_by(dyads,triads,transition) %>%
    tally()
  
  lamb_trans_temp=ctmc.sim %>%
    dplyr::group_by(TrioID) %>%
    dplyr::filter(DateSubmitted==max(DateSubmitted)) %>%
    #dplyr::group_by(dyads,triads) %>%
    dplyr::summarize(ind.l=match(paste("lambda",dyads,triads,sep=""),colnames(lambda)),
                     TimeBetween,dyads,triads) %>%
    dplyr::filter(!is.na(ind.l)) 
  
  
  for(i in 2:niter){
    # lambda_update=sum_temp %>%
    #   left_join(
    #     #Subract 1 for each censored observation, which is a TrioID ending in a motif other than m31  
    #     #Count by TrioID
    #     tal_temp %>%
    #       dplyr::left_join(lamb_temp %>%
    #                          dplyr::group_by(dyads,triads) %>%
    #                          dplyr::summarize(num=sum(1-exp(-lambda[i-1,ind.l]*as.numeric(TimeBetween))),
    #                                           den=sum(-2*as.numeric(TimeBetween)*exp(-lambda[i-1,ind.l]*as.numeric(TimeBetween))-
    #                                                     (1/lambda[i-1,ind.l])*exp(-lambda[i-1,ind.l]*as.numeric(TimeBetween))+
    #                                                     1/lambda[i-1,ind.l]))
    #                        ,
    #                        by=c("dyads"="dyads","triads"="triads")
    #       ) %>% dplyr::left_join(max_tal_temp,
    #                              by=c("dyads"="dyads","triads"="triads")
    #       ) %>%
    #       dplyr::group_by(dyads,triads) %>% 
    #       dplyr::mutate(num.mod=ifelse(!is.na(num),num,0),den.mod=ifelse(!is.na(den),den,0),
    #                     n.mod=ifelse(!is.na(n.y),n.x-n.y,n.x)),
    #     by=c("dyads"="dyads","triads"="triads")
    #   ) %>%
    #   dplyr::mutate(lambdahat=(num.mod+n.mod)/(den.mod+sumtime))
    
    lambda_update=sum_temp %>%
      left_join(
        #Subract 1 for each censored observation, which is a TrioID ending in a motif other than m31  
        #Count by TrioID
        tal_temp %>%
          dplyr::left_join(lamb_temp %>%
                             dplyr::group_by(dyads,triads) %>%
                             dplyr::summarize(num=sum((lambda[i-1,ind.l]-(as.numeric(TimeBetween)+lambda[i-1,ind.l])*exp(-as.numeric(TimeBetween)/lambda[i-1,ind.l]))/(1-exp(-as.numeric(TimeBetween)/lambda[i-1,ind.l]))+as.numeric(TimeBetween)*exp(-as.numeric(TimeBetween)/lambda[i-1,ind.l])),
                                            #num=sum(lambda[i-1,ind.l]*(1-exp(-as.numeric(TimeBetween)/lambda[i-1,ind.l]))+as.numeric(TimeBetween)*exp(-as.numeric(TimeBetween)/lambda[i-1,ind.l])),
                                            den=sum((1-exp(-as.numeric(TimeBetween)/lambda[i-1,ind.l]))))
                           ,
                           by=c("dyads"="dyads","triads"="triads")
          ) %>% dplyr::left_join(max_tal_temp,
                                 by=c("dyads"="dyads","triads"="triads")
          ) %>%
          dplyr::group_by(dyads,triads) %>% 
          dplyr::mutate(num.mod=ifelse(!is.na(num),num,0),den.mod=ifelse(!is.na(den),den,0),
                        n.mod=ifelse(!is.na(n.cens),n.trans-n.cens,n.trans)),
        by=c("dyads"="dyads","triads"="triads")
      ) %>%
      dplyr::mutate(lambdahat=(num.mod+sumtime)/(den.mod+n.trans))
    lambda[i,match(paste("lambda",lambda_update %>% .$dyads,lambda_update %>% .$triads,sep=""),colnames(lambda))]=(lambda_update %>% .$lambdahat)
    
    trans_prob_update=trans_temp %>% 
      left_join(lamb_trans_temp%>%
                  dplyr::group_by(dyads,triads) %>%
                  dplyr::summarize(sum.e.part=sum(1-exp(-as.numeric(TimeBetween)/lambda[i-1,ind.l]))),
                by=c("dyads"="dyads","triads"="triads")
      ) 
    ind.p=match(trans_prob_update %>% .$transition,gsub("p","",gsub("\\|","",colnames(trans_probs))))
    n=trans_prob_update %>% .$n
    names(n)=trans_prob_update %>% .$transition
    e.part=trans_prob_update %>% .$sum.e.part
    e.part[is.na(e.part)]=1
    names(e.part)=trans_prob_update %>% .$transition
    
    # ## 3 possible transitions
    # trans_probs[i,"p01|00"]=(n["0100"]*e.part["0100"])/((n["0100"]+e.part["0100"])*(n["1000"]+e.part["1000"])-n["0100"]*n["1000"])
    # trans_probs[i,"p10|00"]=(n["1000"]*e.part["1000"])/((n["1000"]+e.part["1000"])*(n["0100"]+e.part["0100"])-n["1000"]*n["0100"])
    # trans_probs[i,"p11|10"]=(n["1110"]*e.part["1110"])/((n["1110"]+e.part["1110"])*(n["2010"]+e.part["2010"])-n["1110"]*n["2010"])
    # trans_probs[i,"p20|10"]=(n["2010"]*e.part["2010"])/((n["2010"]+e.part["2010"])*(n["1110"]+e.part["1110"])-n["2010"]*n["1110"])
    # trans_probs[i,"p21|20"]=(n["2120"]*e.part["2120"])/((n["2120"]+e.part["2120"])*(n["3020"]+e.part["3020"])-n["2120"]*n["3020"])
    # trans_probs[i,"p30|20"]=(n["3020"]*e.part["3020"])/((n["3020"]+e.part["3020"])*(n["2120"]+e.part["2120"])-n["3020"]*n["2120"])
    # 
    # 
    # ## 2 possible transitions
    # trans_probs[i,"p11|01"]=n["1101"]/(e.part["1101"]+n["1101"])
    # trans_probs[i,"p21|11"]=n["2111"]/(e.part["2111"]+n["2111"])
    # trans_probs[i,"p31|21"]=n["3121"]/(e.part["3121"]+n["3121"])
    # trans_probs[i,"p31|30"]=n["3130"]/(e.part["3130"]+n["3130"])
    # 
    ## 3 possible transitions
    trans_probs[i,"p01|00"]=(n["0100"])/(e.part["0100"]+n["1000"]+n["0100"])
    trans_probs[i,"p10|00"]=(n["1000"])/(e.part["1000"]+n["1000"]+n["0100"])
    trans_probs[i,"p11|10"]=(n["1110"])/(e.part["1110"]+n["1110"]+n["2010"])
    trans_probs[i,"p20|10"]=(n["2010"])/(e.part["2010"]+n["2010"]+n["1110"])
    trans_probs[i,"p21|20"]=(n["2120"])/(e.part["2120"]+n["2120"]+n["3020"])
    trans_probs[i,"p30|20"]=(n["3020"])/(e.part["3020"]+n["3020"]+n["2120"])
    
    
    ## 2 possible transitions
    trans_probs[i,"p11|01"]=n["1101"]/(e.part["1101"]+n["1101"])
    trans_probs[i,"p21|11"]=n["2111"]/(e.part["2111"]+n["2111"])
    trans_probs[i,"p31|21"]=n["3121"]/(e.part["3121"]+n["3121"])
    trans_probs[i,"p31|30"]=n["3130"]/(e.part["3130"]+n["3130"])
    
    
    
  }
  return(list(lambda[niter,],trans_probs[niter,]))
  
}


##Careful here, timebetween can be infinite if it was actually censored.  Need to change triad.sim function to only report Inf if it is truly censored if you want to use this
ctmc.sim %>% 
  dplyr::group_by(TrioID) %>%
  dplyr::filter(DateSubmitted==max(DateSubmitted)) %>%
  dplyr::group_by(dyads,triads) %>%
  dplyr::summarize(n.abs=sum(DateSubmitted==Inf),n.truecens=sum(DateSubmitted<Inf))


##  Simulation using functions ------------------------------------------------------------------------------------------------


library(parallel)
library(foreach)
cl <- makeCluster(4)
clusterExport(cl,varlist=c("triad.sim","pars"))

nsim=100
em_lambda=matrix(0,nrow=nsim,ncol=length(lambda_true))
em_trans_probs=matrix(0,nrow=nsim,ncol=length(trans_probs_true))
for(j in 1:nsim){
  ctmc.sim=rbindlist(lapply(1:500,function(x){data.frame(triad.sim(pars[[i]],1187.8956),TrioID=x)}))
  #rbindlist(parLapplyLB(cl,1:1000,function(x){data.frame(triad.sim(pars,5000),TrioID=x)}))
  colnames(ctmc.sim)=c("dyads","triads","TimeBetween","DateSubmitted","TrioID")
  ctmc.sim=data.frame(ctmc.sim)
  
  both=EM_alg(ctmc.sim)
  em_lambda[j,]=both[[1]]
  em_trans_probs[j,]=both[[2]]
  print(j)
}
colnames(em_lambda)=names(lambda_true)
colnames(em_trans_probs)=names(trans_probs_true)


ggplot(data.frame(melt(em_lambda)),aes(x=value))+
  geom_histogram()+
  facet_grid(.~Var2)

ggplot(data.frame(melt(em_trans_probs)),aes(x=value))+
  geom_histogram()+
  facet_grid(.~Var2)


