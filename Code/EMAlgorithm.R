
## This function returns Maximum Likelihood Estimates for holding parameters and 
## transition probabilities using an EM algorithm 
# Input: A hypergraph data frame, number of iterations
# Output:  A named vector of maximum likelihood estimates

EM_alg=function(ctmc.sim,niter=30){
  
  lambda=matrix(0,nrow=niter,ncol=6)
  lambda[1,]=rep(500,6)
  colnames(lambda)=c("lambda01","lambda10","lambda20","lambda30","lambda11","lambda21")
  
  trans_probs=matrix(0,nrow=niter,ncol=10)
  trans_probs[1,]=rep(0.333333,10)
  colnames(trans_probs)=c("p01|00","p10|00","p11|10",
                          "p20|10","p11|01","p21|11",
                          "p31|21","p31|30","p21|20",
                          "p30|20")
  
  ##Initialize holding time parameter vector
  lamb_temp=ctmc.sim %>%
    dplyr::group_by(TrioID) %>%
    dplyr::slice(which(DateSubmitted==max(DateSubmitted))) %>%
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
  
  
  M_step=function(l,tp){
    test_lamb=lamb_temp %>%
      dplyr::group_by(TrioID) %>%
       dplyr::mutate(pabs=(1-sum(tp[which(gsub("(\\w+)\\|(\\w+)","\\2",names(tp)) %in% paste(dyads,triads,sep=""))])),
         den=((1-exp(-TimeBetween/l[ind.l]))*pabs+exp(-TimeBetween/l[ind.l])),
         delta0=(1-exp(-TimeBetween/l[ind.l]))*pabs/den,
         delta1=exp(-TimeBetween/l[ind.l])/den)
         #delta0=(1-exp(-TimeBetween/l[ind.l])),
         #delta1=exp(-TimeBetween/l[ind.l]))
        # 
  
    
    # q.iter=function(lamb.uni,curr.par){
    #   lamb=rep(10,7)
    #   names(lamb)=c("lambda00","lambda01","lambda10","lambda20","lambda30","lambda11","lambda21")
    #   lamb[curr.par]=lamb.uni
    #   q.final=sum_temp %>%
    #     left_join(
    #       tal_temp %>%
    #         dplyr::left_join(test_lamb %>%
    #                            dplyr::group_by(dyads,triads) %>%
    #                            dplyr::summarize(q=sum(-delta0*TimeBetween/(exp(TimeBetween*lamb[ind.l])-1)+delta1*TimeBetween))
    #                          ,
    #                          by=c("dyads"="dyads","triads"="triads")
    #         ) %>% dplyr::left_join(max_tal_temp,
    #                                by=c("dyads"="dyads","triads"="triads")
    #         ) %>%
    #         dplyr::group_by(dyads,triads) %>%
    #         dplyr::mutate(q.mod=ifelse(!is.na(q),q,0)),
    #       by=c("dyads"="dyads","triads"="triads")
    #     )  %>%
    #     dplyr::mutate(ind.l=match(paste("lambda",dyads,triads,sep=""),names(lamb))) %>%
    #     dplyr::ungroup() %>%
    #     dplyr::mutate(q.final=n.trans/(q.mod+sumtime)) %>% select(ind.l,q.final)
    #   return(q.final$q.final[q.final$ind.l==match(curr.par,names(lamb))])
    # }
    # 
    # fp=function(curr.par){
    #   x=numeric()
    #   x[1]=1/500
    #   x[2]=q.iter(x,curr.par)
    #   i=2
    #   while(abs((x[i]-x[i-1])/x[i])>0.01 & i<30){
    #     x[i+1]=q.iter(x[i],curr.par)
    #     i=i+1
    #   }
    #   if(i==30){return(NA)}else{
    #     return(1/x[i])}
    # }
    # l=sapply(names(l),function(x){fp(curr.par=x)})
    # 
    q.func=function(lamb.uni,curr.par){
      lamb=l
      names(lamb)=c("lambda01","lambda10","lambda20","lambda30","lambda11","lambda21")
      lamb[curr.par]=lamb.uni
      q.final=sum_temp %>%
        filter(!(dyads==0 & triads==0)) %>% 
        left_join(
          tal_temp %>%
            dplyr::left_join(test_lamb %>%
                               dplyr::group_by(dyads,triads) %>%
                               dplyr::summarize(q=sum(delta0*log((1-exp(-TimeBetween/lamb[ind.l]))*pabs)-delta1*TimeBetween/lamb[ind.l]))
                             ,
                             by=c("dyads"="dyads","triads"="triads")
            ) %>% dplyr::left_join(max_tal_temp,
                                   by=c("dyads"="dyads","triads"="triads")
            ) %>%
            dplyr::group_by(dyads,triads) %>%
            dplyr::mutate(q.mod=ifelse(!is.na(q),q,0)),
          by=c("dyads"="dyads","triads"="triads")
        )  %>%
        dplyr::mutate(ind.l=match(paste("lambda",dyads,triads,sep=""),names(lamb))) %>%
        dplyr::ungroup() %>%
        #dplyr::summarize(q.final=sum(q.mod-n.trans*log(lamb[ind.l])-sumtime/lamb[ind.l])) %>% .$q.final
        dplyr::mutate(q.final=q.mod-n.trans*log(lamb[ind.l])-sumtime/lamb[ind.l]) %>% select(ind.l,q.final)
      return(-q.final$q.final[q.final$ind.l==match(curr.par,names(lamb))])
    }
    
    
    #lambda[i,]=sapply(colnames(lambda),function(x){optimize(q.func,c(1,5000),curr.par=x)$minimum})
    #l[is.na(l)]=sapply(names(l)[is.na(l)],function(x){optimize(q.func,c(1,5000),curr.par=x)$minimum})
    l=sapply(names(l),function(x){optimize(q.func,c(1,5000),curr.par=x)$minimum})
    
    trans_prob_update=trans_temp %>% 
      left_join(test_lamb%>%
                  dplyr::group_by(dyads,triads) %>%
                  dplyr::summarize(d0=sum(delta0)),
                by=c("dyads"="dyads","triads"="triads")
      ) 
    ind.p=match(trans_prob_update %>% .$transition,gsub("p","",gsub("\\|","",colnames(trans_probs))))
    n=trans_prob_update %>% .$n
    names(n)=trans_prob_update %>% .$transition
    e.part=trans_prob_update %>% .$d0
    e.part[is.na(e.part)]=0
    names(e.part)=trans_prob_update %>% .$transition
    
    ## 3 possible transitions
    tp["p01|00"]=(n["0100"])/(e.part["0100"]+n["1000"]+n["0100"])
    tp["p10|00"]=(n["1000"])/(e.part["1000"]+n["1000"]+n["0100"])
    tp["p11|10"]=(n["1110"])/(e.part["1110"]+n["1110"]+n["2010"])
    tp["p20|10"]=(n["2010"])/(e.part["2010"]+n["2010"]+n["1110"])
    tp["p21|20"]=(n["2120"])/(e.part["2120"]+n["2120"]+n["3020"])
    tp["p30|20"]=(n["3020"])/(e.part["3020"]+n["3020"]+n["2120"])
    
    
    ## 2 possible transitions
    tp["p11|01"]=n["1101"]/(e.part["1101"]+n["1101"])
    tp["p21|11"]=n["2111"]/(e.part["2111"]+n["2111"])
    tp["p31|21"]=n["3121"]/(e.part["3121"]+n["3121"])
    tp["p31|30"]=n["3130"]/(e.part["3130"]+n["3130"]) 
    
    return(list(l,tp))
  }
  
  update=M_step(lambda[1,],trans_probs[1,])
  lambda[2,]=update[[1]]
  trans_probs[2,]=update[[2]]
  i=2
  while(max(abs((lambda[i,]-lambda[i-1,])/lambda[i,]))>0.01 & i<niter){
    update=M_step(lambda[i,],trans_probs[i,])
    lambda[i+1,]=update[[1]]
    trans_probs[i+1,]=update[[2]]
    i=i+1
  }
  
  return(list(lambda[i,],trans_probs[i,]))
  
}

EM_alg(ctmc.sim)
both.obs=EM_alg(ctmc.em.obs,niter=100)




dq.func=function(lamb.uni,curr.par){
  lamb=rep(10,7)
  names(lamb)=c("lambda00","lambda01","lambda10","lambda20","lambda30","lambda11","lambda21")
  lamb[curr.par]=lamb.uni
  q.final=sum_temp %>%
    left_join(
      #Subract 1 for each censored observation, which is a TrioID ending in a motif other than m31  
      #Count by TrioID
      tal_temp %>%
        dplyr::left_join(test_lamb %>%
                           dplyr::group_by(dyads,triads) %>% 
                           dplyr::summarize(q=sum((-delta0*TimeBetween*exp(-TimeBetween/lamb[ind.l]))/((1-exp(-TimeBetween/lamb[ind.l]))*lamb[ind.l]^2)+delta1*TimeBetween/lamb[ind.l]^2))
                         ,
                         by=c("dyads"="dyads","triads"="triads")
        ) %>% dplyr::left_join(max_tal_temp,
                               by=c("dyads"="dyads","triads"="triads")
        ) %>%
        dplyr::group_by(dyads,triads) %>% 
        dplyr::mutate(q.mod=ifelse(!is.na(q),q,0)),
      by=c("dyads"="dyads","triads"="triads")
    )  %>%
    dplyr::mutate(ind.l=match(paste("lambda",dyads,triads,sep=""),names(lamb))) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(q.final=q.mod-n.trans/(lamb[ind.l])+sumtime/lamb[ind.l]^2) %>% .$q.final           
  return(q.final[curr.par])
}                     
dq.func(20,"lambda00")
uniroot(dq.func,c(1,10000),curr.par="lambda30")

sapply(colnames(lambda),function(x){uniroot(dq.func,c(1,100000),curr.par=x)$root})
lambda[i,]=optim(rep(10,7),q.func)$par


# 
# dq.func=function(lamb.uni,curr.par){
#   lamb=rep(10,7)
#   names(lamb)=c("lambda00","lambda01","lambda10","lambda20","lambda30","lambda11","lambda21")
#   lamb[curr.par]=lamb.uni
#   q.final=sum_temp %>%
#     left_join(
#       #Subract 1 for each censored observation, which is a TrioID ending in a motif other than m31  
#       #Count by TrioID
#       tal_temp %>%
#         dplyr::left_join(test_lamb %>%
#                            dplyr::group_by(dyads,triads) %>% 
#                            dplyr::summarize(q=sum(delta0*exp(-TimeBetween/lamb[ind.l])/((1-exp(-TimeBetween/lamb[ind.l]))*lamb[ind.l])+delta1*TimeBetween/lamb[ind.l]^2))
#                              #q=sum((-delta0*TimeBetween*exp(-TimeBetween/lamb[ind.l]))/((1-exp(-TimeBetween/lamb[ind.l]))*lamb[ind.l]^2)+delta1*TimeBetween/lamb[ind.l]^2))
#                          ,
#                          by=c("dyads"="dyads","triads"="triads")
#         ) %>% dplyr::left_join(max_tal_temp,
#                                by=c("dyads"="dyads","triads"="triads")
#         ) %>%
#         dplyr::group_by(dyads,triads) %>% 
#         dplyr::mutate(q.mod=ifelse(!is.na(q),q,0)),
#       by=c("dyads"="dyads","triads"="triads")
#     )  %>%
#     dplyr::mutate(ind.l=match(paste("lambda",dyads,triads,sep=""),names(lamb))) %>% 
#     dplyr::ungroup() %>% 
#     dplyr::mutate(q.final=q.mod-n.trans/(lamb[ind.l])+sumtime/lamb[ind.l]^2) %>% .$q.final           
#   return(q.final[curr.par])
# }    


#lambda[i,]=sapply(colnames(lambda),function(x){uniroot(dq.func,c(1,1000),extendInt="yes",curr.par=x)$root})




##Averages
#Transitions
ctmc.em.obs %>% 
  dplyr::group_by(TrioID) %>%
  dplyr::filter(DateSubmitted<max(DateSubmitted)) %>%
  dplyr::group_by(dyads,triads) %>%
  dplyr::summarize(mean(TimeBetween))

ctmc.em.obs %>% 
  dplyr::group_by(TrioID) %>%
  dplyr::filter(DateSubmitted==max(DateSubmitted)) %>%
  dplyr::group_by(dyads,triads) %>%
  dplyr::summarize(mean(TimeBetween))




## Standard Error Estimates --------------------------------------------------------------------------------------------

# Computes the hessian and returns standard error estimates
lambda_mle=c(251.2032, 283.7163, 202.1596, 197.9379, 135.5240, 150.1860, 104.8401)
transprobs_mle=c(0.4499191, 0.4518757, 0.1983434, 0.4532005, 0.2482894, 0.3039186, 0.2837092, 0.1996541, 0.1959190, 0.2999592)
lambda_mle=both.obs[[1]]
transprobs_mle=both.obs[[2]]
se_est=function(lambda_mle,transprobs_mle){
  
  lamb_temp=ctmc.sim %>%
    dplyr::group_by(TrioID) %>%
    dplyr::filter(DateSubmitted==max(DateSubmitted)) %>%
    #dplyr::group_by(dyads,triads) %>%
    dplyr::summarize(ind.l=match(paste("lambda",dyads,triads,sep=""),names(lambda_mle)),
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
    dplyr::summarize(ind.l=match(paste("lambda",dyads,triads,sep=""),names(lambda_mle)),
                     TimeBetween,dyads,triads) %>%
    dplyr::filter(!is.na(ind.l)) 

    
   test_lamb=lamb_temp %>%
      dplyr::group_by(TrioID) %>%
      dplyr::mutate(pabs=(1-sum(transprobs_mle[which(gsub("(\\w+)\\|(\\w+)","\\2",names(transprobs_mle)) %in% paste(dyads,triads,sep=""))])),
        #den=((1-exp(-TimeBetween/lambda[i-1,ind.l]))*pabs+exp(-TimeBetween/lambda[i-1,ind.l])),
        #delta0=(1-exp(-TimeBetween/lambda[i-1,ind.l]))*pabs/den,
        #delta1=exp(-TimeBetween/lambda[i-1,ind.l])/den)
        delta0=(1-exp(-TimeBetween/lambda_mle[ind.l])),
        delta1=exp(-TimeBetween/lambda_mle[ind.l]))
    
   lamb_2d=sum_temp %>%
     left_join(
       tal_temp %>%
         dplyr::left_join(test_lamb %>% 
                            dplyr::group_by(dyads,triads) %>% 
                            # dplyr::summarize(h=sum(-(exp(-TimeBetween/lambda_mle[ind.l])*TimeBetween/lambda_mle[ind.l]^2-(exp(-TimeBetween/lambda_mle[ind.l])*pabs*TimeBetween)/lambda_mle[ind.l]^2)^2/
                            #                    (exp(-TimeBetween/lambda_mle[ind.l])+(1-exp(-TimeBetween/lambda_mle[ind.l]))*pabs)^2+
                            #                    ((-2*exp(-TimeBetween/lambda_mle[ind.l])*TimeBetween)/lambda_mle[ind.l]^3+(2*exp(-TimeBetween/lambda_mle[ind.l])*pabs*TimeBetween)/lambda_mle[ind.l]^3+
                            #                    (exp(-TimeBetween/lambda_mle[ind.l])*TimeBetween^2)/lambda_mle[ind.l]^4+(exp(-TimeBetween/lambda_mle[ind.l])*pabs*TimeBetween^2)/lambda_mle[ind.l]^4)/
                            #                    (exp(-TimeBetween/lambda_mle[ind.l])+(1-exp(-TimeBetween/lambda_mle[ind.l]))*pabs))),
                            dplyr::summarize(h=sum(((-1+pabs)*TimeBetween*(exp(TimeBetween/lambda_mle[ind.l])*pabs*(2*lambda_mle[ind.l]-TimeBetween)+2*(-1+pabs)*(-lambda_mle[ind.l]+TimeBetween)))/
                                               (lambda_mle[ind.l]^4*(1+(-1+exp(TimeBetween/lambda_mle[ind.l]))*pabs)^2))
                                               
                                               
                            ),
                          by=c("dyads"="dyads","triads"="triads")
         ) %>% dplyr::left_join(max_tal_temp,
                                by=c("dyads"="dyads","triads"="triads")
         ) %>%
         dplyr::group_by(dyads,triads) %>%
         dplyr::mutate(h.mod=ifelse(!is.na(h),h,0)),
       by=c("dyads"="dyads","triads"="triads")
     )  %>%
     dplyr::mutate(ind.l=match(paste("lambda",dyads,triads,sep=""),names(lambda_mle)),
                   h.final=n.trans/lambda_mle[ind.l]^2-2*sumtime/lambda_mle[ind.l]^3+h.mod)
    
   tp_2d= trans_temp %>%
     left_join(
       tal_temp %>%
         dplyr::left_join(test_lamb %>% 
                            dplyr::group_by(dyads,triads) %>% 
                            dplyr::summarize(h=sum(-(1-exp(-TimeBetween/lambda_mle[ind.l]))^2/((1-exp(-TimeBetween/lambda_mle[ind.l]))*pabs+exp(-TimeBetween/lambda_mle[ind.l]))^2),
                                             dldp=sum(exp(TimeBetween/lambda_mle[ind.l])/(lambda_mle[ind.l]^2*((1-exp(-TimeBetween/lambda_mle[ind.l]))*pabs+exp(-TimeBetween/lambda_mle[ind.l]))^2)),
                                             dp1dp2=sum(-(1-exp(-TimeBetween/lambda_mle[ind.l]))^2/(((1-exp(-TimeBetween/lambda_mle[ind.l]))*pabs+exp(-TimeBetween/lambda_mle[ind.l]))^2)))
                          ,
                          by=c("dyads"="dyads","triads"="triads")
         ) %>% dplyr::left_join(max_tal_temp,
                                by=c("dyads"="dyads","triads"="triads")
         ) %>%
         dplyr::group_by(dyads,triads) %>%
         dplyr::mutate(h.mod=ifelse(!is.na(h),h,0),
                       dldp=ifelse(!is.na(dldp),dldp,0),
                       dp1dp2=ifelse(!is.na(dp1dp2),dp1dp2,0)),
       by=c("dyads"="dyads","triads"="triads")
     )  %>%
     dplyr::mutate(ind.l=match(paste("lambda",dyads,triads,sep=""),names(lambda_mle)),
                   ind.p=match(transition,gsub("p","",gsub("\\|","",names(transprobs_mle)))),
                   h.final=(h.mod-n/transprobs_mle[ind.p]^2))
    lamb_2d=lamb_2d %>% filter(!(dyads==0 & triads==0))
    hessian=diag(c(lamb_2d %>% .$h.final, tp_2d %>% .$h.final))
    rownames(hessian)=colnames(hessian)=c(names(lambda_mle),names(transprobs_mle))
    tp_2d= tp_2d %>% filter(!(dyads==0 & triads==0))
    lamb.names=paste("lambda",tp_2d %>%.$dyads,tp_2d %>%.$triads,sep="")
    tp.names=gsub("(\\d\\d)(\\d\\d)","p\\1\\|\\2",tp_2d %>%.$transition)
    hessian[cbind(lamb.names,tp.names)]=tp_2d %>% .$dldp
    hessian[cbind(tp.names,lamb.names)]=tp_2d %>% .$dldp
    hessian["p21|20","p30|20"]=tp_2d %>% dplyr::filter(transition=="2120") %>%  .$dp1dp2
    hessian["p30|20","p21|20"]=tp_2d %>% dplyr::filter(transition=="2120") %>%  .$dp1dp2
    
  return(sqrt(diag(solve(-hessian))))
  
}


##Error Bar Plots -----------------------------------------------------------

#Holding Parameters
df.error=data.frame(par=gsub("(\\w+)(\\d\\d)","\\1\\[\\'\\2\\'\\]",names(both.obs[[1]])),lamb=both.obs[[1]],se=sqrt(diag(solve(-hessian)))[grep("lambda",names(sqrt(diag(solve(-hessian)))))])

#tikzDevice::tikz("HoldingEstimates.tikz",width=6,height=3)
ggplot(df.error, aes(x=par,y=lamb/365))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=(lamb-2*se)/365,ymax=(lamb+2*se)/365),width=0.3,lwd=1.2)+
  scale_x_discrete("",labels=parse(text=levels(df.error$par)))+
  ggtitle("Holding Parameters")+
  #theme(text=element_text(size=20,family="serif"))+
  ylab("Years")+
  theme_hc(base_size=25,base_family="serif")
#dev.off()
ggsave(file="C:\\Users\\haunf\\Google Drive\\Academic\\Research\\HyperGraph\\CSUInternalNetworkSummary\\MotifModelPaper\\Images\\Figure11a.png",
       width=12.44,
       height=4.22)
       
#Transition Probabilities
df.error=data.frame(par=gsub("(\\w)(\\d+)\\|(\\d+)","m\\[\\'\\2\\'\\]\\*\\'\\|\\'\\*m\\[\\'\\3\\'\\]",names(both.obs[[2]])),
                    tp=both.obs[[2]],se=sqrt(diag(solve(-hessian)))[grep("p",names(sqrt(diag(solve(-hessian)))))])
df.error=df.error %>%
  dplyr::mutate(To=gsub("(\\w)(\\d+)\\|(\\d+)","m\\[\\'\\2\\'\\]",names(both.obs[[2]])),
                From=gsub("(\\w)(\\d+)\\|(\\d+)","m\\[\\'\\3\\'\\]",names(both.obs[[2]])))

df.error=df.error %>% bind_rows(df.error %>%
  dplyr::group_by(From) %>%
  dplyr::summarize(tp=1-sum(tp)) %>%
  dplyr::mutate(To="Absorbing",
                se=0,
                par=paste("m['",From,"'*'|'*m['",From,"']",sep="")))

ggplot(df.error, aes(x=From,y=tp,fill=To))+
  geom_bar(position="stack",stat="identity")+
  scale_fill_discrete(name="Transition To",breaks=levels(as.factor(df.error$To)),labels=parse(text=levels(as.factor(df.error$To))))+
  #geom_errorbar(aes(ymin=tp-2*se,ymax=tp+2*se),width=0.3,lwd=1.2)+
  scale_x_discrete("",labels=parse(text=levels(as.factor(df.error$From))))+
  ggtitle("Transition Probabilities")+
  #theme_hc()+
  theme_minimal()+
  theme(text=element_text(size=20,family="serif"),
        legend.text.align=0)+
  ylab("")
ggsave(file="C:\\Users\\haunf\\Google Drive\\Academic\\Research\\HyperGraph\\CSUInternalNetworkSummary\\MotifModelPaper\\Images\\TPEstimates.png",
       width=12.44,
       height=4.22)
