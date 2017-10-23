
## Markov chain simulation function

triad_sim=function(pars,bigT){
  pars["p00|00"]=1-pars[["p10|00"]]-pars[["p01|00"]]
  pars["p10|10"]=1-pars[["p20|10"]]-pars[["p11|10"]]
  pars["p01|01"]=1-pars[["p11|01"]]
  pars["p11|11"]=1-pars[["p21|11"]]
  pars["p21|21"]=1-pars[["p31|21"]]
  pars["p30|30"]=1-pars[["p31|30"]]
  pars["p20|20"]=1-pars[["p30|20"]]-pars[["p21|20"]]
  cumtime = 0
  timebetween = NULL
  dyads = 0
  triads = 0
  i = 1
  while(cumtime[i] < bigT & !(dyads[i]==3 & triads[i]==1)){
    timebetween = c(timebetween,rexp(1, pars[paste("lambda",dyads[i],triads[i],sep="")]))
    cumtime = c(cumtime,sum(timebetween))
    ind.p = which(gsub("p(\\w+)\\|(\\w+)","\\2",names(pars)) %in% paste(dyads[i],triads[i],sep=""))
    new.motif = gsub("p(\\w+)\\|(\\w+)","\\1",names(pars))[ind.p][sample(1:length(ind.p),1,prob=pars[ind.p])]
    dyads = c(dyads,strsplit(new.motif,split="")[[1]][1])
    triads = c(triads,strsplit(new.motif,split="")[[1]][2])
    i = i+1
    if(dyads[i]==dyads[i-1] & triads[i]==triads[i-1]){cumtime[i]=Inf}
  }
  
  if(dyads[i]==3 & triads[i]==1 & (cumtime[i]+timebetween[i-1])<bigT){
    timebetween[i] = bigT-sum(timebetween[1:(i-1)])
    cumtime[i+1] = bigT
    return(data.frame(dyads=dyads,triads=triads,TimeBetween=timebetween,cumtime=cumtime[-1]))
  }else{
    if(i>2){
      timebetween[i-1] = bigT-sum(timebetween[1:(i-2)])
    }else{
      timebetween[i-1] = bigT
    }
    return(data.frame(dyads = dyads[-i],
                      triads = triads[-i],
                      TimeBetween = timebetween[-i],
                      cumtime = cumtime[-1]))
  }
}