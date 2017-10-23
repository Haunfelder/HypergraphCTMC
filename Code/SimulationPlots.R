

library(lattice)

##Compute MSE
pars=c(lambda_true,trans_probs_true)
mses=lapply(sim_results,function(x){
  apply(
    sweep(x,2,pars[match(colnames(x),names(pars))],"-"),2,function(y){
      mean(y[y<Inf & y>-Inf]^2,na.rm=T)
    })
})

names(tandn)=c("T","n")

##Construct dataframe of mse results
df=data.frame(tandn,Reduce(rbind,mses))
df_ests=data.frame(rep(tandn[,1],each=nsim),rep(tandn[,2],each=nsim),Reduce(rbind,sim_results))
names(df_ests)[1:2]=c("Var1","Var2")
names(df_ests)[-2:-1]=names(pars)

##Make Plot
df_lambda=melt(df[,c(1,grep("lambda",colnames(df)))],id="T")
df_trans_probs=melt(df[,c(1,grep("\\.",colnames(df)))],id="T")
df_trans_probs$variable=gsub("\\.","\\|",df_trans_probs$variable)

ggplot(df_lambda, aes(x=T,y=log(value),group=variable,colour=variable))+
  geom_line()+
  #coord_cartesian(ylim=c(0,0.0025))+
  ylab("log MSE")+
  theme(legend.title=element_blank())

ggplot(df_trans_probs, aes(x=T,y=log(value),group=variable,colour=variable))+
  geom_line()+
  #coord_cartesian(ylim=c(0,0.0025))+
  ylab("log MSE")+
  theme(legend.title=element_blank())

##Surface as a function of t and n
df_lambda=melt(df[,c(1,2,grep("lambda",colnames(df)))],id=c("Var1","Var2"))
df_lambda_ests=melt(df_ests[,c(1,2,grep("lambda",colnames(df)))],id=c("Var1","Var2"))
df_trans_probs=melt(df[,c(1,2,grep("\\.",colnames(df)))],id=c("Var1","Var2"))
df_trans_probs$variable=gsub("\\.","\\|",df_trans_probs$variable)
df_trans_probs_ests=melt(df_ests[,c(1,2,grep("\\.",colnames(df)))],id=c("Var1","Var2"))
df_trans_probs_ests$variable=gsub("\\.","\\|",df_trans_probs_ests$variable)
df_real=data.frame(variable=names(pars),value=pars)

##Add dyads and triads
df_lambda_ests=df_lambda_ests %>%
  mutate(dyads=sapply(str_split(gsub("lambda","",variable),""),function(x){x[1]}),
         triads=sapply(str_split(gsub("lambda","",variable),""),function(x){x[2]}))

df_real_lambdas= df_real %>%
  filter(grepl("lambda",variable)) %>%
  mutate(dyads=sapply(str_split(gsub("lambda","",variable),""),function(x){x[1]}),
         triads=sapply(str_split(gsub("lambda","",variable),""),function(x){x[2]}))

df_lambda_ests=df_lambda_ests %>% mutate(bigt.factor=factor(as.character(Var1)))


ggplot(df_lambda,aes(x=Var1,y=Var2,fill=log(value)))+
  geom_raster()+
  scale_fill_gradient(low="red",high="yellow")

ggplot(df_lambda_ests,aes(x=factor(Var1),y=value))+
  geom_boxplot()+
  facet_grid(Var2~dyads+triads)+
  geom_hline(data=df_real[grep("lambda",df_real$variable),],aes(yintercept=value),colour="red",lty=2)+
  coord_cartesian(ylim=c(0,0.04))


ggplot(df_lambda_ests,aes(x=Var1,group=cut_interval(x=Var1,length=100),y=value))+
  geom_boxplot(colour="black",notch=FALSE)+
  facet_grid(Var2~dyads+triads)+
  geom_hline(data=df_real_lambdas,aes(yintercept=value),colour="red",lty=2)+
  geom_vline(data=bigtdf,aes(xintercept=bigt,colour=factor(quantiles)))+
  coord_cartesian(ylim=c(0,1000))+
  xlab("T")+
  ylab("Estimate")+
  guides(colour=guide_legend(title="Quantiles"))+
  theme(text = element_text(size=20))








ggplot(df_trans_probs_ests,aes(x=Var2,y=value,group=Var2))+
  geom_boxplot()+
  facet_grid(.~variable)+
  geom_hline(data=df_real[grep("\\|",df_real$variable),],aes(yintercept=value),colour="red",lty=2)
ggplot(df_lambda,aes(x=Var1,))


stat_density2d(aes(fill=..level..),geom="polygon")


##Plots of bias and variance
i=1
for(i in 1:4){
  #sim_results[[i]]=sapply(sim_results[[i]],function(x){rownames(x)=c(); return(x)})
  
  df_ests=data.frame(rep(tandn[[i]][,1],times=sapply(sim_results[[i]],function(x){dim(x)[1]})),
                     rep(tandn[[i]][,2],times=sapply(sim_results[[i]],function(x){dim(x)[1]})),
                     do.call(rbind,sim_results[[i]]))
  names(df_ests)[1:2]=c("Var1","Var2")
  names(df_ests)[-2:-1]=colnames(sim_results[[i]][[1]])
  
  ##Surface as a function of t and n
  df_lambda=melt(df[,c(1,2,grep("lambda",colnames(df)))],id=c("Var1","Var2"))
  df_lambda_ests=melt(df_ests[,c(1,2,grep("lambda",colnames(df)))],id=c("Var1","Var2"))
  df_trans_probs=melt(df[,c(1,2,grep("\\.",colnames(df)))],id=c("Var1","Var2"))
  df_trans_probs$variable=gsub("\\.","\\|",df_trans_probs$variable)
  df_trans_probs_ests=melt(df_ests[,c(1,2,grep("\\.",colnames(df)))],id=c("Var1","Var2"))
  df_trans_probs_ests$variable=gsub("\\.","\\|",df_trans_probs_ests$variable)
  df_real=data.frame(variable=names(pars[[i]]),value=pars[[i]])
  
  ##Add dyads and triads
  df_lambda_ests=df_lambda_ests %>%
    mutate(dyads=sapply(str_split(gsub("lambda","",variable),""),function(x){x[1]}),
           triads=sapply(str_split(gsub("lambda","",variable),""),function(x){x[2]}))
  
  df_real_lambdas= df_real %>%
    filter(grepl("lambda",variable)) %>%
    mutate(dyads=sapply(str_split(gsub("lambda","",variable),""),function(x){x[1]}),
           triads=sapply(str_split(gsub("lambda","",variable),""),function(x){x[2]}))
  
  df_lambda_ests=df_lambda_ests %>% mutate(bigt.factor=factor(as.character(Var1)))
  
  
  
  df_biasvar=df_lambda_ests %>%
    dplyr::filter(value<4900) %>%
    dplyr::group_by(Var1,Var2,dyads,triads) %>%
    dplyr::summarize(Expected=mean(value[value<Inf],na.rm=T),Variance=sd(value[value<Inf],na.rm=T)^2) %>%
    dplyr::left_join(df_real_lambdas,by=c("dyads"="dyads","triads"="triads")) %>%
    dplyr::mutate(bias=Expected-value)
  
  # ggplot(df_biasvar %>% filter(dyads==0,triads==0), aes(x=as.numeric(Var1),y=as.numeric(Var2),colour=as.numeric(bias))) +
  #   geom_point()
  # 
  # 
  # z=spread(df_biasvar%>% ungroup() %>% filter(dyads==0,triads==0) %>% select(Var1,Var2,bias),Var1,bias)
  # persp(unique(df_biasvar %>% filter(dyads==0,triads==0) %>% .$Var1),
  #       unique(df_biasvar %>% filter(dyads==0,triads==0) %>% .$Var2),
  #       t(as.matrix(z[,-1])),xlab="T",ylab="N",zlab="Bias")
  # 
  # 
  # 
  # df_wf=df_biasvar %>% dplyr::filter(dyads==0,triads==0)
  # wireframe(bias~Var1*Var2,data=df_wf,xlab="T",ylab="n")
  
  df_wf=df_biasvar
  df_wf=df_wf %>% dplyr::mutate(Exp=paste(gsub("lambda","m[",variable),"]",sep=""))
  
  png(paste("C:\\Users\\haunf\\Google Drive\\Academic\\Research\\HyperGraph\\CSUInternalNetworkSummary\\MotifModelPaper\\Images\\SimResults",i,".png",sep=""),width=1667,height=300)
  print(wireframe(bias~Var1*Var2|Exp,data=df_wf,xlab="T",ylab="n",zlab="",
                  strip=strip.custom(factor.levels=expression(m["00"],m["01"],m["10"],m[11],m[20],m[21],m[30])),
                  scales=list(arrows=F),
                  par.settings = list(strip.background=list(col="lightgrey")),
                  par.strip.text=list(cex=3)))
  dev.off()
  
}
wireframe(bias~Var1*Var2|dyads*triads,data=df_wf,xlab="T",ylab="n",zlab="",colorkey=T,scales=list(arrows=F))
wireframe(Variance~Var1*Var2|dyads*triads,data=df_wf,xlab="T",ylab="n",colorkey=T)
wireframe(Variance~Var1*Var2|dyads*triads,data=df_wf,xlab="T",ylab="n",colorkey=T)


ggplot(df_wf, aes(x=Var1, y=Var2, fill=Exp)) +
  geom_tile()



