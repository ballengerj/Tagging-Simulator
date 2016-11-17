setwd(wd)

source('Simulation_control.r')
source('Brownie_control.r')

setwd(paste0(wd,"/Figures/Report Files"))

out =readList("Brownie1.rep")
par_names=c('Obj_Func','SIM_M','SIM_T_region','SIM_F','SIM_report_rate','M','T_region','F','report_rate_CNST','Tag_Loss','Hand_Mort')
result=out[par_names]
age_max_recap<-out$max_age_recap
M_true=result$SIM_M
F_true=result$SIM_F
T_true=result$SIM_T_region
B_true=result$SIM_report_rate
for(i in 1:nstocks)
{
  for(y in 1:nyrs)
  {
    for(a in 1:nages)
    {
      if(a>age_max_recap[y])
      {
        if(out$F_switch==1)
        {
        F_true[(y+nyrs*(i-1)),a]=0
        }
        if(out$M_switch>1)
        {
        M_true[(y+nyrs*(i-1)),a]=0
        }
      }
    }
  }
}

tags=function(iteration1)
{
  if(file.exists(paste0(wd,"/Figures/Report Files/Brownie",iteration1,".cor",sep=""))==TRUE)
  {
  Brownie_out1 =readList(paste(wd,"/Figures/Report Files/Brownie",iteration1,".rep",sep=""))
    par_names1=c('OBS_total_rec',"total_rec")
    results1=Brownie_out1[par_names1]
      c(Tag_recap_true=as.vector(results1$OBS_total_rec),
     Tag_recap=as.vector(results1$total_rec))
  }
  else
  {
    Brownie_out1 =readList(paste(wd,"/Figures/Report Files/Brownie",iteration1,".rep",sep=""))
    par_names1=c('OBS_total_rec',"total_rec")
    results1=Brownie_out1[par_names1]
      c(Tag_recap_true=as.vector(results1$OBS_total_rec),
      Tag_recap=as.vector(rep(NA, times=length(results1$total_rec))))
  }
}

recaps=cbind(sapply(1:ntrials,function(i)tags(i)))
recaps_true=recaps[c(1:(nstocks*ncohorts*nages)),]
recaps_pred=recaps[c((nstocks*ncohorts*nages+1):(2*nstocks*ncohorts*nages)),]

tag_bias<-matrix(NA,nrow=length(recaps_true[,1]),ncol=ntrials)

for(i in 1:length(recaps_true[,1]))  
{
  for(j in 1:ntrials)
  {
    if(is.na(recaps_pred[i,j]))
    {tag_bias[i,j]=NA}
    else{tag_bias[i,j]=recaps_pred[i,j]-recaps_true[i,j]}
  }
}

tag_percent_bias<-matrix(NA,nrow=length(recaps_true[,1]),ncol=ntrials)

for(i in 1:length(recaps_true[,1]))  
{
  for(j in 1:ntrials)
  {
    if(recaps_true[i,j]==0 || is.na(tag_bias[i,j]))
    {tag_percent_bias[i,j]=NA}
    else{tag_percent_bias[i,j]=tag_bias[i,j]/recaps_true[i,j]*100}
  }
}

results=function(iteration)
{
  if(file.exists(paste0(wd,"/Figures/Report Files/Brownie",iteration,".cor",sep=""))==TRUE)
  {
    Brownie_out =readList(paste(wd,"/Figures/Report Files/Brownie",iteration,".rep",sep=""))
    par_names=c('Obj_Func','SIM_M','SIM_T_region','SIM_F','SIM_report_rate','M','T_region','F','report_rate_CNST','Tag_Loss','Hand_Mort')
    results=Brownie_out[par_names]
    c(rep=iteration,
      objective=results$Obj_Func,
      M_est=as.vector(results$M),
      F_est=as.vector(results$F),
      T_est=as.vector(results$T_region),
      B_est=as.vector(results$report_rate_CNST)
    )
  }
  #if(file.exists(paste0(wd,"/Program_files/Brownie",iteration,".cor",sep=""))==FALSE)
  else
  {
    Brownie_out =readList(paste(wd,"/Figures/Report Files/Brownie",iteration,".rep",sep=""))
    par_names=c('Obj_Func','SIM_M','SIM_T_region','SIM_F','SIM_report_rate','M','T_region','F','report_rate_CNST','Tag_Loss','Hand_Mort')
    results=Brownie_out[par_names]
    c(rep=iteration,
      objective=NA,
      M_est=as.vector(rep(NA,times=length(results$M))),
      F_est=as.vector(rep(NA,times=length(results$F))),
      T_est=as.vector(rep(NA, times=length(results$T_region))),
      B_est=as.vector(rep(NA, times=length(results$report_rate_CNST))))
  }
}

true=c(0,0,as.vector(M_true),as.vector(F_true),as.vector(T_true),as.vector(B_true))
est=cbind(sapply(1:ntrials,function(i)results(i)))


Converged<-(ntrials-length(which(is.na(est[2,]))))
percent.converged<-(Converged/ntrials)*100

bias<-matrix(NA,nrow=length(true),ncol=ntrials)

for(i in 1:length(true))  
{
  for(j in 1:ntrials)
  {
    if(is.na(est[i,j]))
    {bias[i,j]=NA}
    else{bias[i,j]=est[i,j]-true[i]}
  }
}

percent_bias<-matrix(NA,nrow=length(true),ncol=ntrials)

for(i in 1:length(true))  
{
  for(j in 1:ntrials)
  {
    if(true[i]==0 || is.na(bias[i,j]))
    {percent_bias[i,j]=NA}
    else{percent_bias[i,j]=bias[i,j]/true[i]*100}
  }
}

##############################################################################################################################################################

plot.beans.combined<-function(parameter.name,data.frame.name,...)
{
  par(mfrow=c(1,1),  mgp=c(2.5,.6,.5), mar=c(5, 4, 3, .5) + 0.1)
  
  beanplot(data.frame.name,bw="nrd0",boxwex=1.2,what=c(0,1,0,0),col=c('grey','black','black','white'),
           beanlines="median",main=paste0(parameter.name," Percent Bias"), ylab="Percent Bias", mgp=c(2.75,.75,.5),
           sub=paste0(Converged, " runs converged out of ", ntrials," (", percent.converged,"% Convergence)"),las=1,
           #cutmin=max(cut,min(as.numeric(unlist(data.frame.name)),na.rm=TRUE)),cutmax=min(cut,max(as.numeric(unlist(data.frame.name)),na.rm=TRUE)), #xlim=c(1,16),
           cex.axis=.5,na.rm=TRUE,...)
  abline(h=0,lty=1,lwd=2)
  abline(h=median(as.numeric(unlist(data.frame.name)),na.rm=TRUE),lty=2,col="red",lwd=2)
  legend("top", c(paste0("Median Bias ", signif(median(as.numeric(unlist(data.frame.name)),na.rm=TRUE),3),"%"),
                  paste0("Max Bias ", signif(max(abs(as.numeric(unlist(data.frame.name))),na.rm=TRUE),3),"%")),lty=c(4,-1),cex=.5,ncol=2,bg='white')
  
}


plot.beans.multi<-function(parameter.name,data.frame.name,plot.number,multiplier1,multiplier2,text.label,...)
{
  if(plot.number<2)
  {
    par(mfrow=c(1,1),  mgp=c(1.85,.6,.5), mar=c(3, 2.5, 1, .5) + 0.1)
  }
  if(plot.number>1 && plot.number<5)
  {
    par(mfrow=c(2,ceiling(plot.number/2)),  mgp=c(1.85,.6,.5), mar=c(3, 2.5, 1, .5) + 0.1)
  }
  if(plot.number>4 && plot.number<7)
  {
    par(mfrow=c(2,3),  mgp=c(1.85,.6,.5), mar=c(3, 2.5, 1, .5) + 0.1)
  }
  if(plot.number>6 && plot.number<10)
  {
    par(mfrow=c(3,3),  mgp=c(1.85,.6,.5), mar=c(3, 2.5, 1, .5) + 0.1)
  }
  if(plot.number>9 && plot.number<13)
  {
    par(mfrow=c(4,3),  mgp=c(1.85,.6,.5), mar=c(3, 2.5, 1, .5) + 0.1)
  }
  if(plot.number>12)
  {
    par(mfrow=c(2,3),  mgp=c(1.85,.6,.5), mar=c(3, 2.5, 1, .5) + 0.1)
  }
  for(n in 1:1)
  {
    for(i in 1:plot.number)
    {
      if(multiplier2==0)
      {
        beanplot(data.frame.name[,seq((1+(i-1)*multiplier1),(i*multiplier1),1)],bw="nrd0",boxwex=1.2,what=c(0,1,0,0),col=c('grey','black','black','white'),
                 beanlines="median",main=paste0(parameter.name,i), ylab="Percent Bias", mgp=c(2.75,.6,.5),las=1,
                 names=c(1:multiplier1),cex.main=.75,cex.lab=.75,
                 #cutmin=max(cut,min(as.numeric(unlist(data.frame.name[,seq((1+(i-1)*multiplier1),(i*multiplier1),1)])),na.rm=TRUE)),
                 #cutmax=min(cut,max(as.numeric(unlist(data.frame.name[,seq((1+(i-1)*multiplier1),(i*multiplier1),1)])),na.rm=TRUE)), #xlim=c(1,16),
                 cex.axis=.5,na.rm=TRUE,...)
        abline(h=0,lty=1,lwd=2)
        abline(h=median(as.numeric(unlist(data.frame.name[,seq((1+(i-1)*multiplier1),(i*multiplier1),1)])),na.rm=TRUE),lty=2,col="red",lwd=2)
        legend("top", c(paste0("Median Bias ", signif(median(as.numeric(unlist(data.frame.name[,seq((1+(i-1)*multiplier1),(i*multiplier1),1)])),na.rm=TRUE),3),"%"),
                        paste0("Max Bias ", signif(max(abs(as.numeric(unlist(data.frame.name[,seq((1+(i-1)*multiplier1),(i*multiplier1),1)]))),na.rm=TRUE),3),"%")),lty=c(4,-1),cex=.5,ncol=2,bg='white')
      }
      if(multiplier2>0)
      {
        beanplot(data.frame.name[,seq((1+(i-1)*multiplier1*multiplier2),(i*multiplier1*multiplier2),1)],bw="nrd0",boxwex=1.2,what=c(0,1,0,0),col=c('grey','black','black','white'),
                 beanlines="median",main=paste0(parameter.name,i), ylab="Percent Bias", mgp=c(2.75,.6,.5),las=1,
                 names=rep(c(1:multiplier2),times=multiplier1),cex.main=.75,cex.lab=.75,
                 #cutmin=max(cut,max(as.numeric(unlist(data.frame.name[,seq((1+(i-1)*multiplier1*multiplier2),(i*multiplier1*multiplier2),1)])),na.rm=TRUE)),
                 #cutmax=min(cut,max(as.numeric(unlist(data.frame.name[,seq((1+(i-1)*multiplier1*multiplier2),(i*multiplier1*multiplier2),1)])),na.rm=TRUE)), #xlim=c(1,16),
                 cex.axis=.5,na.rm=TRUE,...)
        abline(h=0,lty=1,lwd=2)
        abline(h=median(as.numeric(unlist(data.frame.name[,seq((1+(i-1)*multiplier1*multiplier2),(i*multiplier1*multiplier2),1)])),na.rm=TRUE),lty=2,col="red",lwd=2)
        abline(v=seq(multiplier2+0.5,(multiplier1-1)*multiplier2+0.5,multiplier2),lty=1)
        shadowtext(c(seq(multiplier2/2+0.5,multiplier1*multiplier2-multiplier2/2+0.5,multiplier2)),c(rep(0,times=multiplier1)),sapply(seq(1:multiplier1),function(s)c(paste0(text.label,s))),col="black",bg="grey",r=bg.text,cex=cex.text.sub)
        legend("top", c(paste0("Median Bias ", signif(median(as.numeric(unlist(data.frame.name[,seq((1+(i-1)*multiplier1*multiplier2),(i*multiplier1*multiplier2),1)])),na.rm=TRUE),3),"%"),
                        paste0("Max Bias ", signif(max(abs(as.numeric(unlist(data.frame.name[,seq((1+(i-1)*multiplier1*multiplier2),(i*multiplier1*multiplier2),1)]))),na.rm=TRUE),3),"%")),lty=c(4,-1),cex=.5,ncol=2,bg='white')
      }
    }
  }
}


################################################################################################################################
################################################################################################################################
################################################################################################################################
###############################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
time_D<-format(Sys.time(),"%e") 
time_M<-format(Sys.time(),"%m") 
time_Y<-format(Sys.time(),"%Y") 
time_H<-format(Sys.time(),"%H") 
time_Mi<-format(Sys.time(),"%M") 

pdf(file=paste0(wd,"/Figures/parameter bias all_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".pdf",sep=""))

  if(out$report_rate_switch==1)
   {  
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    B_bias=as.data.frame(t(percent_bias[seq(1,nstocks,1)+nages*nstocks^2+2+2*nages*nstocks*nyrs,]))
    colnames(B_bias)<-sapply(seq(1:nstocks),function(i)c(paste0("Region-",i)))
    cut_min_B=(-4)*max(rowSds(as.matrix(B_bias), na.rm=TRUE))
    cut_max_B=(4)*max(rowSds(as.matrix(B_bias), na.rm=TRUE))
    plot.beans.combined('Reporting Rate',B_bias,xlab="Recapture Region")
   }
  
  if(out$movement_switch==1 || out$movement_switch==2 || out$movement_switch==5)
   {
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    T_bias_region=as.data.frame(t(percent_bias[seq(1,nages*nstocks^2,nages)+2+2*nages*nstocks*nyrs,]))
    colnames(T_bias_region)<-sapply(rep(seq(1:nstocks),times=nstocks),function(i)c(paste0("Reg-",i)))
    cut_min_T_region=(-4)*max(rowSds(as.matrix(T_bias_region), na.rm=TRUE))
    cut_max_T_region=(4)*max(rowSds(as.matrix(T_bias_region), na.rm=TRUE))
    plot.beans.combined('Movement Rate',T_bias_region,xlab="Source Region") #,xlim=c(1,ncol(T_bias_region)))
    abline(v=c(seq(nstocks+0.5,nstocks*(nstocks-1)+0.5,nstocks)),lty=1)
    plot.beans.combined('Movement Rate',T_bias_region,xlab="Source Region",add=TRUE) #,xlim=c(1,ncol(T_bias_region)))
    shadowtext(c(seq(nstocks-nstocks/2+.5,nstocks*(nstocks)-nstocks/2+.5,nstocks)),c(rep(0,times=nstocks)),sapply(seq(1:nstocks),function(i)c(paste0("To Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.text.main)
   }
  
  if(out$movement_switch==3 || out$movement_switch==4 || out$movement_switch==6 || out$movement_switch==7)
   {
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    T_bias_age=as.data.frame(t(percent_bias[seq(1,nages*nstocks^2,1)+2+2*nages*nstocks*nyrs,]))
    colnames(T_bias_age)<-sapply(rep(seq(1:nages),times=nstocks*nstocks),function(i)c(paste0(i)))
    min.loc<-(min(as.numeric(unlist(T_bias_age)),na.rm=TRUE)/3)
    cut_min_T_age=(-4)*max(rowSds(as.matrix(T_bias_age), na.rm=TRUE))
    cut_max_T_age=(4)*max(rowSds(as.matrix(T_bias_age), na.rm=TRUE))
    plot.beans.combined('Movement Rate',T_bias_age,xlab="Age") #,xlim=c(4,ncol(T_bias_age)-4))
    abline(v=c(seq(nages*nstocks+0.5,nages*nstocks*(nstocks-1)+0.5,nstocks*nages)),lty=1)
    abline(v=c(seq(nages+0.5,length(T_bias_age[1,])-nages+0.5,nages)),lty=2)
    plot.beans.combined('Movement Rate',T_bias_age,xlab="Age") #,xlim=c(4,ncol(T_bias_age)-4),add=TRUE)
    shadowtext(c(seq(nstocks*nages/2+0.5,nstocks*nages*(nstocks)-(nstocks*nages)/2+0.5,nstocks*nages)),c(rep(0,times=nstocks)),sapply(seq(1:nstocks),function(i)c(paste0("To Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.text.main)
    shadowtext(c(seq(nages/2+0.5,nstocks*nages*(nstocks)-(nages)/2+0.5,nages)),c(rep(min.loc,times=nstocks*nstocks)),sapply(rep(seq(1:nstocks),times=nstocks),function(i)c(paste0("Reg",i))),col="black",bg="grey",r=bg.text,cex=cex.text.sub)

    plot.beans.multi("Movement to Reg-",T_bias_age,nstocks,nstocks,nages,"From Reg-",xlab="Age")
    
    for(i in 1:nstocks)
    {
      par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
      plot.beans.combined(paste0('Movement Rate To Region ',i),T_bias_age[,c(seq((1+(i-1)*nstocks*nages),i*nages*nstocks,1))],xlab="Age" #,xlim=c(1,ncol(T_bias_age[,seq((1+(i-1)*nstocks*nages),(i*nages*nstocks),1)]))
                          ,names=rep(c(1:nages),nstocks))
      abline(v=c(seq(nages+0.5,(nstocks-1)*nages+0.5,nages)),lty=2)
      plot.beans.combined(paste0('Movement Rate To Region ',i),T_bias_age[,seq((1+(i-1)*nstocks*nages),(i*nages*nstocks),1)],xlab="Age", #xlim=c(1,ncol(T_bias_age[,seq((1+(i-1)*nstocks*nages),(i*nages*nstocks),1)]))
                          add=TRUE,names=rep(c(1:nages),nstocks))
      shadowtext(c(seq(nages/2+0.5,nstocks*nages-(nages)/2+0.5,nages)),c(rep(0,times=nstocks)),sapply(seq(1:nstocks),function(i)c(paste0("From Region",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
    }

   }
    
  if(out$M_switch==1)
   {
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    
    M_bias_CNST=as.data.frame(percent_bias[3,])
    colnames(M_bias_CNST)<-"M"
    cut_min_M_CNST=(-4)*max(rowSds(as.matrix(M_bias_CNST), na.rm=TRUE))
    cut_max_M_CNST=(4)*max(rowSds(as.matrix(M_bias_CNST), na.rm=TRUE))
    plot.beans.combined("Natural Mortality",M_bias_CNST,xlab="Constant Natural Mortality")
   }

  if(out$M_switch==2)
   {
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    
    M_bias_age=as.data.frame(t(percent_bias[seq(1,nages*nstocks*nyrs,nstocks*nyrs)+2,]))
    colnames(M_bias_age)<-sapply(seq(1:nages),function(i)c(paste0("Age-",i)))
    cut_min_M_age=(-4)*max(rowSds(as.matrix(M_bias_age), na.rm=TRUE))
    cut_max_M_age=(4)*max(rowSds(as.matrix(M_bias_age), na.rm=TRUE))
    plot.beans.combined("Natural Mortality",M_bias_age,xlab="Age")
   }

  if(out$M_switch==3)
   {
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    
    M_bias_region_age=as.data.frame(t(percent_bias[seq(1,nages*nstocks*nyrs,nyrs)+2,]))
    colnames(M_bias_region_age)<-sapply(rep(seq(1:nstocks),times=nages),function(i)c(paste0("Reg-",i)))
    cut_min_M_region_age=(-4)*max(rowSds(as.matrix(M_bias_region_age), na.rm=TRUE))
    cut_max_M_region_age=(4)*max(rowSds(as.matrix(M_bias_region_age), na.rm=TRUE))
    plot.beans.combined("Natural Mortality",M_bias_region_age #,xlim=c(1,ncol(M_bias_region_age))
                        ,xlab="Region")
    abline(v=c(seq(nstocks+0.5,length(M_bias_region_age[1,])-nstocks+0.5,nstocks)),lty=1)
    plot.beans.combined("Natural Mortality",M_bias_region_age #,xlim=c(1,ncol(M_bias_region_age))
                        ,add=TRUE,xlab="Region")
    shadowtext(c(seq(nstocks/2+0.5,length(M_bias_region_age[1,])-nstocks/2+0.5,nstocks)),c(rep(0,times=nages)),sapply(seq(1:nages),function(i)c(paste0("Age-",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)

    plot.beans.multi("M Age-",M_bias_region_age,nages,nstocks,0,"Region",xlab="Region")
    
        for(i in 1:nages)
    {
      par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
      
      plot.beans.combined(paste0('Natural Mortality Rate Age-',i),M_bias_region_age[,seq((1+(i-1)*nstocks),(i*nstocks),1)],xlab="Region",names=c(1:nstocks))
    }
    
    }

  if(out$M_switch==4)
   {
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    
    M_bias_region_year_age=as.data.frame(t(percent_bias[seq(1,nages*nstocks*nyrs,1)+2,]))
    colnames(M_bias_region_year_age)<-sapply(rep(seq(1:nyrs),times=nages*nstocks),function(i)c(paste0(i)))
    min.loc<-(min(as.numeric(unlist(M_bias_region_year_age)),na.rm=TRUE)/3)
    cut_min_M_region_year_age=(-4)*max(rowSds(as.matrix(M_bias_region_year_age), na.rm=TRUE))
    cut_max_M_region_year_age=(4)*max(rowSds(as.matrix(M_bias_region_year_age), na.rm=TRUE))
    plot.beans.combined("Natural Mortality",M_bias_region_year_age #,xlim=c(4,ncol(M_bias_region_year_age)-4)
                        ,xlab="Year")
    abline(v=c(seq(nyrs*nstocks+0.5,length(M_bias_region_year_age[1,])-nyrs*nstocks+0.5,nyrs*nstocks)),lty=1)
    abline(v=c(seq(nyrs+0.5,length(M_bias_region_year_age[1,])-nyrs+0.5,nyrs)),lty=2)
    plot.beans.combined("Natural Mortality",M_bias_region_year_age #,xlim=c(4,ncol(M_bias_region_year_age)-4)
                        ,xlab="Year",add=TRUE)
    shadowtext(c(seq(nyrs*nstocks/2+0.5,length(M_bias_region_year_age[1,])-nyrs*nstocks/2+0.5,nyrs*nstocks)),c(rep(0,times=nages)),sapply(seq(1:nages),function(i)c(paste0("Age-",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
    shadowtext(c(seq(nyrs/2+0.5,length(M_bias_region_year_age[1,])-nyrs/2+0.5,nyrs)),c(rep(min.loc,times=nstocks*nages)),sapply(rep(seq(1:nstocks),times=nages),function(i)c(paste0("Reg",i))),col="black",bg="grey",r=bg.text,cex=cex.text.sub)

    plot.beans.multi("M Age-",M_bias_region_year_age,nages,nstocks,nyrs,"Region",xlab="Year")
    
    for(i in 1:nages)
    {
      par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
      
      plot.beans.combined(paste0('Natural Mortality Rate Age-',i),M_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)],xlab="Year",names=rep(c(1:nyrs),nstocks))
      abline(v=c(seq(nyrs+0.5,length(M_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)])-nyrs+0.5,nyrs)),lty=2)
      plot.beans.combined(paste0('Natural Mortality Rate Age-',i),M_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)],xlab="Year",names=rep(c(1:nyrs),nstocks),add=TRUE)
      shadowtext(c(seq(nyrs/2+0.5,length(M_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)])-nyrs/2+0.5,nyrs)),c(rep(0,times=nstocks)),sapply(1:nstocks,function(i)c(paste0("Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
    }
  }

  if(out$F_switch>0)
  {
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    
    F_bias_region_year_age=as.data.frame(t(percent_bias[seq(1,nages*nstocks*nyrs,1)+2+nages*nstocks*nyrs,]))
    colnames(F_bias_region_year_age)<-sapply(rep(seq(1:nyrs),times=nages*nstocks),function(i)c(paste0(i)))
    min.loc<-(min(as.numeric(unlist(F_bias_region_year_age)),na.rm=TRUE)/3)
    cut_min_F_region_year_age=(-4)*max(rowSds(as.matrix(F_bias_region_year_age), na.rm=TRUE))
    cut_max_F_region_year_age=(4)*max(rowSds(as.matrix(F_bias_region_year_age), na.rm=TRUE))
    plot.beans.combined("Fishing Mortality",F_bias_region_year_age #,xlim=c(4,ncol(F_bias_region_year_age)-4)
                        ,xlab="Year")
    abline(v=c(seq(nyrs*nstocks+0.5,length(F_bias_region_year_age[1,])-nyrs*nstocks+0.5,nyrs*nstocks)),lty=1)
    abline(v=c(seq(nyrs+0.5,length(F_bias_region_year_age[1,])-nyrs+0.5,nyrs)),lty=2)
    plot.beans.combined("Fishing Mortality",F_bias_region_year_age #,xlim=c(4,ncol(F_bias_region_year_age)-4)
                        ,xlab="Year",add=TRUE)
    shadowtext(c(seq(nyrs*nstocks/2+0.5,length(F_bias_region_year_age[1,])-nyrs*nstocks/2+0.5,nyrs*nstocks)),c(rep(0,times=nages)),sapply(seq(1:nages),function(i)c(paste0("Age-",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
    shadowtext(c(seq(nyrs/2+0.5,length(F_bias_region_year_age[1,])-nyrs/2+0.5,nyrs)),c(rep(min.loc,times=nstocks*nages)),sapply(rep(seq(1:nstocks),times=nages),function(i)c(paste0("Reg",i))),col="black",bg="grey",r=bg.text,cex=cex.text.sub)
   
     plot.beans.multi("F Age-",F_bias_region_year_age,nages,nstocks,nyrs,"Region",xlab="Year")
    
    for(i in 1:nages)
    {
      par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
      
      plot.beans.combined(paste0('Fishing Mortality Rate Age-',i),F_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)],xlab="Year",names=rep(c(1:nyrs),nstocks))
      abline(v=c(seq(nyrs+0.5,length(F_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)])-nyrs+0.5,nyrs)),lty=2)
      plot.beans.combined(paste0('Fishing Mortality Rate Age-',i),F_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)],xlab="Year",names=rep(c(1:nyrs),nstocks),add=TRUE)
      shadowtext(c(seq(nyrs/2+0.5,length(F_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)])-nyrs/2+0.5,nyrs)),c(rep(0,times=nstocks)),sapply(1:nstocks,function(i)c(paste0("Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
    }
    }

#############################################################################################################################################################################
############## Percent Bias for Total Tag Recaptures by Release Events (stock/cohort/age of release combinations), good summary statistic of bias associated with tag recaptures
########################################################################################################################################################################
  
#############################################################################################################################################################################################
### Ages with no releases are not included because total recaptures are summarized by release event not recapture event, therefore if no releases for a given age then total recaps==0
### Using min(age_max_recap) excludes ages greater than the oldest age of release (since these ages have no releases)
#############################################################################################################################
  recap_bias=as.data.frame(t(tag_percent_bias[c(1:(min(age_max_recap)*ncohorts*nstocks)),]))
##########################################################################################################################
  colnames(recap_bias)<-sapply(rep(seq(1:ncohorts),times=(min(age_max_recap)*nstocks)),function(i)c(paste0(i)))
  min.loc<-(min(as.numeric(unlist(recap_bias)),na.rm=TRUE)/3)
  cut_min_recap_bias=(-4)*max(rowSds(as.matrix(recap_bias), na.rm=TRUE))
  cut_max_recap_bias=(4)*max(rowSds(as.matrix(recap_bias), na.rm=TRUE))
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined("Total Tags Recaptured by Release Cohort",recap_bias #,xlim=c(4,ncol(recap_bias)-4)
                      ,xlab="Release Year")
  abline(v=c(seq(ncohorts*nstocks+0.5,(min(age_max_recap)*ncohorts*nstocks)-ncohorts*nstocks+0.5,ncohorts*nstocks)),lty=1)
  abline(v=c(seq(ncohorts+0.5,(min(age_max_recap)*ncohorts*nstocks)-ncohorts+0.5,ncohorts)),lty=2)
  plot.beans.combined("Total Tags Recaptured by Release Cohort",recap_bias #,xlim=c(4,ncol(recap_bias)-4)
                      ,xlab="Release Year",add=TRUE)
  shadowtext(c(seq(ncohorts*nstocks/2+0.5,(min(age_max_recap)*ncohorts*nstocks)-ncohorts*nstocks/2+0.5,ncohorts*nstocks)),c(rep(0,times=min(age_max_recap))),sapply(seq(1:min(age_max_recap)),function(i)c(paste0("Age-",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
  shadowtext(c(seq(ncohorts/2+0.5,(min(age_max_recap)*ncohorts*nstocks)-ncohorts/2+0.5,ncohorts)),c(rep(min.loc,times=nstocks*min(age_max_recap))),sapply(rep(seq(1:nstocks),times=min(age_max_recap)),function(i)c(paste0("Reg",i))),col="black",bg="grey",r=bg.text,cex=cex.text.sub)

  plot.beans.multi("Recaptures Age-",recap_bias,min(age_max_recap),nstocks,ncohorts,"Region",xlab="Year")
  
  for(i in 1:min(age_max_recap))
  {
    plot.beans.combined(paste0('Total Tags Recaptured for Age-',i),recap_bias[,seq((1+(i-1)*nstocks*ncohorts),(i*nstocks*ncohorts),1)],xlab="Year",names=rep(c(1:ncohorts),nstocks))
    abline(v=c(seq(ncohorts+0.5,length(recap_bias[,seq((1+(i-1)*nstocks*ncohorts),(i*nstocks*ncohorts),1)])-ncohorts+0.5,ncohorts)),lty=2)
    plot.beans.combined(paste0('Total Tags Recaptured for Age-',i),recap_bias[,seq((1+(i-1)*nstocks*ncohorts),(i*nstocks*ncohorts),1)],xlab="Year",names=rep(c(1:ncohorts),nstocks),add=TRUE)
    shadowtext(c(seq(ncohorts/2+0.5,length(recap_bias[,seq((1+(i-1)*nstocks*ncohorts),(i*nstocks*ncohorts),1)])-ncohorts/2+0.5,nyrs)),c(rep(0,times=nstocks)),sapply(1:nstocks,function(i)c(paste0("Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
  }

dev.off()



################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
pdf(file=paste0(wd,"/Figures/Summary parameter bias_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".pdf",sep=""))

if(out$report_rate_switch==1)
{  
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined('Reporting Rate',B_bias,xlab="Recapture Region")
}

if(out$movement_switch==1 || out$movement_switch==2 || out$movement_switch==5)
{
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined('Movement Rate',T_bias_region,xlab="Source Region") #,xlim=c(1,ncol(T_bias_region)))
  abline(v=c(seq(nstocks+0.5,nstocks*(nstocks-1)+0.5,nstocks)),lty=1)
  plot.beans.combined('Movement Rate',T_bias_region,xlab="Source Region",add=TRUE) #,xlim=c(1,ncol(T_bias_region)))
  shadowtext(c(seq(nstocks-nstocks/2+.5,nstocks*(nstocks)-nstocks/2+.5,nstocks)),c(rep(0,times=nstocks)),sapply(seq(1:nstocks),function(i)c(paste0("To Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.text.main)
}

if(out$movement_switch==3 || out$movement_switch==4 || out$movement_switch==6 || out$movement_switch==7)
{

  plot.beans.multi("Movement to Reg-",T_bias_age,nstocks,nstocks,nages,"From Reg-",xlab="Age")

}

if(out$M_switch==1)
{
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined("Natural Mortality",M_bias_CNST,xlab="Constant Natural Mortality")
}

if(out$M_switch==2)
{
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined("Natural Mortality",M_bias_age,xlab="Age")
}

if(out$M_switch==3)
{
  plot.beans.multi("M Age-",M_bias_region_age,nages,nstocks,0,"Region",xlab="Region")
}

if(out$M_switch==4)
{
  plot.beans.multi("M Age-",M_bias_region_year_age,nages,nstocks,nyrs,"Region",xlab="Year")
}

if(out$F_switch>0)
{
  plot.beans.multi("F Age-",F_bias_region_year_age,nages,nstocks,nyrs,"Region",xlab="Year")
}

plot.beans.multi("Recaptures Age-",recap_bias,min(age_max_recap),nstocks,ncohorts,"Region",xlab="Year")

dev.off()




################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
pdf(file=paste0(wd,"/Figures/Summary + Ind Plot parameter bias_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".pdf",sep=""))

if(out$report_rate_switch==1)
{  
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined('Reporting Rate',B_bias,xlab="Recapture Region")
}

if(out$movement_switch==1 || out$movement_switch==2 || out$movement_switch==5)
{
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined('Movement Rate',T_bias_region,xlab="Source Region") #,xlim=c(1,ncol(T_bias_region)))
  abline(v=c(seq(nstocks+0.5,nstocks*(nstocks-1)+0.5,nstocks)),lty=1)
  plot.beans.combined('Movement Rate',T_bias_region,xlab="Source Region",add=TRUE) #,xlim=c(1,ncol(T_bias_region)))
  shadowtext(c(seq(nstocks-nstocks/2+.5,nstocks*(nstocks)-nstocks/2+.5,nstocks)),c(rep(0,times=nstocks)),sapply(seq(1:nstocks),function(i)c(paste0("To Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.text.main)
}

if(out$movement_switch==3 || out$movement_switch==4 || out$movement_switch==6 || out$movement_switch==7)
{
  plot.beans.multi("Movement to Reg-",T_bias_age,nstocks,nstocks,nages,"From Reg-",xlab="Age")
  
  for(i in 1:nstocks)
  {
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    plot.beans.combined(paste0('Movement Rate To Region ',i),T_bias_age[,c(seq((1+(i-1)*nstocks*nages),i*nages*nstocks,1))],xlab="Age" #,xlim=c(1,ncol(T_bias_age[,seq((1+(i-1)*nstocks*nages),(i*nages*nstocks),1)]))
                        ,names=rep(c(1:nages),nstocks))
    abline(v=c(seq(nages+0.5,(nstocks-1)*nages+0.5,nages)),lty=2)
    plot.beans.combined(paste0('Movement Rate To Region ',i),T_bias_age[,seq((1+(i-1)*nstocks*nages),(i*nages*nstocks),1)],xlab="Age" #,xlim=c(1,ncol(T_bias_age[,seq((1+(i-1)*nstocks*nages),(i*nages*nstocks),1)]))
                        ,add=TRUE,names=rep(c(1:nages),nstocks))
    shadowtext(c(seq(nages/2+0.5,nstocks*nages-(nages)/2+0.5,nages)),c(rep(0,times=nstocks)),sapply(seq(1:nstocks),function(i)c(paste0("From Region",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
  }
}

if(out$M_switch==1)
{
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined("Natural Mortality",M_bias_CNST,xlab="Constant Natural Mortality")
}

if(out$M_switch==2)
{
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined("Natural Mortality",M_bias_age,xlab="Age")
}

if(out$M_switch==3)
{
  plot.beans.multi("M Age-",M_bias_region_age,nages,nstocks,0,"Region",xlab="Region")
  
  for(i in 1:nages)
  {
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    plot.beans.combined(paste0('Natural Mortality Rate Age-',i),M_bias_region_age[,seq((1+(i-1)*nstocks),(i*nstocks),1)],xlab="Region",names=c(1:nstocks))
  }
  
}

if(out$M_switch==4)
{
  plot.beans.multi("M Age-",M_bias_region_year_age,nages,nstocks,nyrs,"Region",xlab="Year")
  
  for(i in 1:nages)
  {
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    plot.beans.combined(paste0('Natural Mortality Rate Age-',i),M_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)],xlab="Year",names=rep(c(1:nyrs),nstocks))
    abline(v=c(seq(nyrs+0.5,length(M_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)])-nyrs+0.5,nyrs)),lty=2)
    plot.beans.combined(paste0('Natural Mortality Rate Age-',i),M_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)],xlab="Year",names=rep(c(1:nyrs),nstocks),add=TRUE)
    shadowtext(c(seq(nyrs/2+0.5,length(M_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)])-nyrs/2+0.5,nyrs)),c(rep(0,times=nstocks)),sapply(1:nstocks,function(i)c(paste0("Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
  }
}

if(out$F_switch>0)
{
  plot.beans.multi("F Age-",F_bias_region_year_age,nages,nstocks,nyrs,"Region",xlab="Year")
  
  for(i in 1:nages)
  {
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    plot.beans.combined(paste0('Fishing Mortality Rate Age-',i),F_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)],xlab="Year",names=rep(c(1:nyrs),nstocks))
    abline(v=c(seq(nyrs+0.5,length(F_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)])-nyrs+0.5,nyrs)),lty=2)
    plot.beans.combined(paste0('Fishing Mortality Rate Age-',i),F_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)],xlab="Year",names=rep(c(1:nyrs),nstocks),add=TRUE)
    shadowtext(c(seq(nyrs/2+0.5,length(F_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)])-nyrs/2+0.5,nyrs)),c(rep(0,times=nstocks)),sapply(1:nstocks,function(i)c(paste0("Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
  }
}

plot.beans.multi("Recaptures Age-",recap_bias,min(age_max_recap),nstocks,ncohorts,"Region",xlab="Year")

for(i in 1:min(age_max_recap))
{
  plot.beans.combined(paste0('Total Tags Recaptured for Age-',i),recap_bias[,seq((1+(i-1)*nstocks*ncohorts),(i*nstocks*ncohorts),1)],xlab="Year",names=rep(c(1:ncohorts),nstocks))
  abline(v=c(seq(ncohorts+0.5,length(recap_bias[,seq((1+(i-1)*nstocks*ncohorts),(i*nstocks*ncohorts),1)])-ncohorts+0.5,ncohorts)),lty=2)
  plot.beans.combined(paste0('Total Tags Recaptured for Age-',i),recap_bias[,seq((1+(i-1)*nstocks*ncohorts),(i*nstocks*ncohorts),1)],xlab="Year",names=rep(c(1:ncohorts),nstocks),add=TRUE)
  shadowtext(c(seq(ncohorts/2+0.5,length(recap_bias[,seq((1+(i-1)*nstocks*ncohorts),(i*nstocks*ncohorts),1)])-ncohorts/2+0.5,nyrs)),c(rep(0,times=nstocks)),sapply(1:nstocks,function(i)c(paste0("Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
}

dev.off()



######################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################

if(out$report_rate_switch==1)
{  
  png(file=paste0(wd,"/Figures/Reporting Rate Bias_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined('Reporting Rate',B_bias,xlab="Recapture Region")
  dev.off()
}

if(out$movement_switch==1 || out$movement_switch==2 || out$movement_switch==5)
{
  png(file=paste0(wd,"/Figures/Movement Rate Bias_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined('Movement Rate',T_bias_region,xlab="Source Region") #,xlim=c(1,ncol(T_bias_region)))
  abline(v=c(seq(nstocks+0.5,nstocks*(nstocks-1)+0.5,nstocks)),lty=1)
  plot.beans.combined('Movement Rate',T_bias_region,xlab="Source Region",add=TRUE) #,xlim=c(1,ncol(T_bias_region)))
  shadowtext(c(seq(nstocks-nstocks/2+.5,nstocks*(nstocks)-nstocks/2+.5,nstocks)),c(rep(0,times=nstocks)),sapply(seq(1:nstocks),function(i)c(paste0("To Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.text.main)
  dev.off()
  }

if(out$movement_switch==3 || out$movement_switch==4 || out$movement_switch==6 || out$movement_switch==7)
{
  png(file=paste0(wd,"/Figures/Movement Rate Bias By Age Combined_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
  
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined('Movement Rate',T_bias_age,xlab="Age") #,xlim=c(4,ncol(T_bias_age)-4))
  abline(v=c(seq(nages*nstocks+0.5,nages*nstocks*(nstocks-1)+0.5,nstocks*nages)),lty=1)
  abline(v=c(seq(nages+0.5,length(T_bias_age[1,])-nages+0.5,nages)),lty=2)
  plot.beans.combined('Movement Rate',T_bias_age,xlab="Age") #,xlim=c(4,ncol(T_bias_age)-4),add=TRUE)
  shadowtext(c(seq(nstocks*nages/2+0.5,nstocks*nages*(nstocks)-(nstocks*nages)/2+0.5,nstocks*nages)),c(rep(0,times=nstocks)),sapply(seq(1:nstocks),function(i)c(paste0("To Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.text.main)
  shadowtext(c(seq(nages/2+0.5,nstocks*nages*(nstocks)-(nages)/2+0.5,nages)),c(rep(min.loc,times=nstocks*nstocks)),sapply(rep(seq(1:nstocks),times=nstocks),function(i)c(paste0("Reg",i))),col="black",bg="grey",r=bg.text,cex=cex.text.sub)
  dev.off()
  
  png(file=paste0(wd,"/Figures/Movement Rate Bias by Age Compare_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
  
  plot.beans.multi("Movement to Reg-",T_bias_age,nstocks,nstocks,nages,"From Reg-",xlab="Age")
  dev.off()
  for(i in 1:nstocks)
  {
    png(file=paste0(wd,"/Figures/Movement Rate Bias by Age Region",i,"_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
    
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    plot.beans.combined(paste0('Movement Rate To Region ',i),T_bias_age[,c(seq((1+(i-1)*nstocks*nages),i*nages*nstocks,1))],xlab="Age" #,xlim=c(1,ncol(T_bias_age[,seq((1+(i-1)*nstocks*nages),(i*nages*nstocks),1)]))
                        ,names=rep(c(1:nages),nstocks))
    abline(v=c(seq(nages+0.5,(nstocks-1)*nages+0.5,nages)),lty=2)
    plot.beans.combined(paste0('Movement Rate To Region ',i),T_bias_age[,seq((1+(i-1)*nstocks*nages),(i*nages*nstocks),1)],xlab="Age" #,xlim=c(1,ncol(T_bias_age[,seq((1+(i-1)*nstocks*nages),(i*nages*nstocks),1)]))
                        ,add=TRUE,names=rep(c(1:nages),nstocks))
    shadowtext(c(seq(nages/2+0.5,nstocks*nages-(nages)/2+0.5,nages)),c(rep(0,times=nstocks)),sapply(seq(1:nstocks),function(i)c(paste0("From Region",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
  dev.off()
     }
  
}

if(out$M_switch==1)
{
  png(file=paste0(wd,"/Figures/Natural Mortality Bias Constant M_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
  
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined("Natural Mortality",M_bias_CNST,xlab="Constant Natural Mortality")
  dev.off()
}

if(out$M_switch==2)
{
  png(file=paste0(wd,"/Figures/Natural Mortality Bias Age-Specific M_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
  
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined("Natural Mortality",M_bias_age,xlab="Age")
  dev.off()
}

if(out$M_switch==3)
{
  png(file=paste0(wd,"/Figures/Natural Mortality Bias by Age and Region Combined_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
  
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined("Natural Mortality",M_bias_region_age #,xlim=c(1,ncol(M_bias_region_age))
                      ,xlab="Region")
  abline(v=c(seq(nstocks+0.5,length(M_bias_region_age[1,])-nstocks+0.5,nstocks)),lty=1)
  plot.beans.combined("Natural Mortality",M_bias_region_age #,xlim=c(1,ncol(M_bias_region_age))
                      ,add=TRUE,xlab="Region")
  shadowtext(c(seq(nstocks/2+0.5,length(M_bias_region_age[1,])-nstocks/2+0.5,nstocks)),c(rep(0,times=nages)),sapply(seq(1:nages),function(i)c(paste0("Age-",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
  dev.off()
  png(file=paste0(wd,"/Figures/Natural Mortality Bias by Age and Region Comparison_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
  
  plot.beans.multi("M Age-",M_bias_region_age,nages,nstocks,0,"Region",xlab="Region")
  dev.off()
  for(i in 1:nages)
  {
    png(file=paste0(wd,"/Figures/Natural Mortality Bias by Region,  Age",i,"_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
    
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    
    plot.beans.combined(paste0('Natural Mortality Rate Age-',i),M_bias_region_age[,seq((1+(i-1)*nstocks),(i*nstocks),1)],xlab="Region",names=c(1:nstocks))
    dev.off()
  }
  
}

if(out$M_switch==4)
{
  png(file=paste0(wd,"/Figures/Natural Mortality Bias by Age, Region, Year Combined_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
  
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined("Natural Mortality",M_bias_region_year_age #,xlim=c(4,ncol(M_bias_region_year_age)-4)
                      ,xlab="Year")
  abline(v=c(seq(nyrs*nstocks+0.5,length(M_bias_region_year_age[1,])-nyrs*nstocks+0.5,nyrs*nstocks)),lty=1)
  abline(v=c(seq(nyrs+0.5,length(M_bias_region_year_age[1,])-nyrs+0.5,nyrs)),lty=2)
  plot.beans.combined("Natural Mortality",M_bias_region_year_age #,xlim=c(4,ncol(M_bias_region_year_age)-4)
                      ,xlab="Year",add=TRUE)
  shadowtext(c(seq(nyrs*nstocks/2+0.5,length(M_bias_region_year_age[1,])-nyrs*nstocks/2+0.5,nyrs*nstocks)),c(rep(0,times=nages)),sapply(seq(1:nages),function(i)c(paste0("Age-",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
  shadowtext(c(seq(nyrs/2+0.5,length(M_bias_region_year_age[1,])-nyrs/2+0.5,nyrs)),c(rep(min.loc,times=nstocks*nages)),sapply(rep(seq(1:nstocks),times=nages),function(i)c(paste0("Reg",i))),col="black",bg="grey",r=bg.text,cex=cex.text.sub)
  dev.off()
  
  png(file=paste0(wd,"/Figures/Natural Mortality Bias by Age, Region, Year Comparison_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
  
  plot.beans.multi("M Age-",M_bias_region_year_age,nages,nstocks,nyrs,"Region",xlab="Year")
  dev.off()
  for(i in 1:nages)
  {
    png(file=paste0(wd,"/Figures/Natural Mortality Bias by Region, Year, Age-",i,"_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
    
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    
    plot.beans.combined(paste0('Natural Mortality Rate Age-',i),M_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)],xlab="Year",names=rep(c(1:nyrs),nstocks))
    abline(v=c(seq(nyrs+0.5,length(M_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)])-nyrs+0.5,nyrs)),lty=2)
    plot.beans.combined(paste0('Natural Mortality Rate Age-',i),M_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)],xlab="Year",names=rep(c(1:nyrs),nstocks),add=TRUE)
    shadowtext(c(seq(nyrs/2+0.5,length(M_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)])-nyrs/2+0.5,nyrs)),c(rep(0,times=nstocks)),sapply(1:nstocks,function(i)c(paste0("Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
  dev.off()
  }
}

if(out$F_switch>0)
{
  png(file=paste0(wd,"/Figures/Fishing Mortality Bias Combined_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
  
  par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
  plot.beans.combined("Fishing Mortality",F_bias_region_year_age #,xlim=c(4,ncol(F_bias_region_year_age)-4)
                      ,xlab="Year")
  abline(v=c(seq(nyrs*nstocks+0.5,length(F_bias_region_year_age[1,])-nyrs*nstocks+0.5,nyrs*nstocks)),lty=1)
  abline(v=c(seq(nyrs+0.5,length(F_bias_region_year_age[1,])-nyrs+0.5,nyrs)),lty=2)
  plot.beans.combined("Fishing Mortality",F_bias_region_year_age #,xlim=c(4,ncol(F_bias_region_year_age)-4)
                      ,xlab="Year",add=TRUE)
  shadowtext(c(seq(nyrs*nstocks/2+0.5,length(F_bias_region_year_age[1,])-nyrs*nstocks/2+0.5,nyrs*nstocks)),c(rep(0,times=nages)),sapply(seq(1:nages),function(i)c(paste0("Age-",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
  shadowtext(c(seq(nyrs/2+0.5,length(F_bias_region_year_age[1,])-nyrs/2+0.5,nyrs)),c(rep(min.loc,times=nstocks*nages)),sapply(rep(seq(1:nstocks),times=nages),function(i)c(paste0("Reg",i))),col="black",bg="grey",r=bg.text,cex=cex.text.sub)
  dev.off()
  png(file=paste0(wd,"/Figures/Fishing Mortality Bias Comparison_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
  
  plot.beans.multi("F Age-",F_bias_region_year_age,nages,nstocks,nyrs,"Region",xlab="Year")
  dev.off()
  for(i in 1:nages)
  {
    png(file=paste0(wd,"/Figures/Fishing Mortality Bias Age-",i,"_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
    
    par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
    
    plot.beans.combined(paste0('Fishing Mortality Rate Age-',i),F_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)],xlab="Year",names=rep(c(1:nyrs),nstocks))
    abline(v=c(seq(nyrs+0.5,length(F_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)])-nyrs+0.5,nyrs)),lty=2)
    plot.beans.combined(paste0('Fishing Mortality Rate Age-',i),F_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)],xlab="Year",names=rep(c(1:nyrs),nstocks),add=TRUE)
    shadowtext(c(seq(nyrs/2+0.5,length(F_bias_region_year_age[,seq((1+(i-1)*nstocks*nyrs),(i*nstocks*nyrs),1)])-nyrs/2+0.5,nyrs)),c(rep(0,times=nstocks)),sapply(1:nstocks,function(i)c(paste0("Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
  dev.off()
  }
}

png(file=paste0(wd,"/Figures/Total Tag Recaptures Bias Combined_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))

par(mfrow=c(1,1),  mgp=c(2.75,.6,.5))
plot.beans.combined("Total Tags Recaptured by Release Cohort",recap_bias #,xlim=c(4,ncol(recap_bias)-4)
                    ,xlab="Release Year")
abline(v=c(seq(ncohorts*nstocks+0.5,(min(age_max_recap)*ncohorts*nstocks)-ncohorts*nstocks+0.5,ncohorts*nstocks)),lty=1)
abline(v=c(seq(ncohorts+0.5,(min(age_max_recap)*ncohorts*nstocks)-ncohorts+0.5,ncohorts)),lty=2)
plot.beans.combined("Total Tags Recaptured by Release Cohort",recap_bias #,xlim=c(4,ncol(recap_bias)-4)
                    ,xlab="Release Year",add=TRUE)
shadowtext(c(seq(ncohorts*nstocks/2+0.5,(min(age_max_recap)*ncohorts*nstocks)-ncohorts*nstocks/2+0.5,ncohorts*nstocks)),c(rep(0,times=min(age_max_recap))),sapply(seq(1:min(age_max_recap)),function(i)c(paste0("Age-",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
shadowtext(c(seq(ncohorts/2+0.5,(min(age_max_recap)*ncohorts*nstocks)-ncohorts/2+0.5,ncohorts)),c(rep(min.loc,times=nstocks*min(age_max_recap))),sapply(rep(seq(1:nstocks),times=min(age_max_recap)),function(i)c(paste0("Reg",i))),col="black",bg="grey",r=bg.text,cex=cex.text.sub)
dev.off()
png(file=paste0(wd,"/Figures/Total Tag Recaptures Bias Comparison_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))

plot.beans.multi("Recaptures Age-",recap_bias,min(age_max_recap),nstocks,ncohorts,"Region",xlab="Year")
dev.off()
for(i in 1:min(age_max_recap))
{
  png(file=paste0(wd,"/Figures/Total Tag Recaptures Bias Age-",i,"_",time_M,"_",time_D,"_",time_Y,"_",time_H,"_",time_Mi,".png",sep=""))
  
  plot.beans.combined(paste0('Total Tags Recaptured for Age-',i),recap_bias[,seq((1+(i-1)*nstocks*ncohorts),(i*nstocks*ncohorts),1)],xlab="Year",names=rep(c(1:ncohorts),nstocks))
  abline(v=c(seq(ncohorts+0.5,length(recap_bias[,seq((1+(i-1)*nstocks*ncohorts),(i*nstocks*ncohorts),1)])-ncohorts+0.5,ncohorts)),lty=2)
  plot.beans.combined(paste0('Total Tags Recaptured for Age-',i),recap_bias[,seq((1+(i-1)*nstocks*ncohorts),(i*nstocks*ncohorts),1)],xlab="Year",names=rep(c(1:ncohorts),nstocks),add=TRUE)
  shadowtext(c(seq(ncohorts/2+0.5,length(recap_bias[,seq((1+(i-1)*nstocks*ncohorts),(i*nstocks*ncohorts),1)])-ncohorts/2+0.5,nyrs)),c(rep(0,times=nstocks)),sapply(1:nstocks,function(i)c(paste0("Region ",i))),col="black",bg="grey",r=bg.text,cex=cex.medium)
 dev.off()
  }

