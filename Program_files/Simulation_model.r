#################################################################################
###		INDIVIDUAL-BASED TAG SIMULATION MODEL 
#################################################################################

### ESTIMATION OF REGIONAL REPORTING RATES
reporting=sapply(1:nstocks,function(s)sum(report_rates*(catches[,s]/sum(catches[,s]))))

### MIXING RATES
if(mixing == 0)
{
  T=array(0,dim=c(nstocks,nstocks,nages))
  for(s in 1:nstocks)
  {
    for(a in 1:nages)
    {
      T[s,s,a]=1
    }
  }
}
if(mixing == 1)
{
  T=array(1/nstocks,dim=c(nstocks,nstocks,nages))
}
if(mixing == 2)
{
  T=array(0,dim=c(nstocks,nstocks,nages))
  for(s in 1:nstocks)
  {
    for(a in 1:nages)
    {
      for(j in 1:nstocks)
      {
        if(s==j)
        {
          T[s,j,a]=residency[s]
        }
        if(s!=j)
        {
          T[s,j,a]=((1-residency[s])/(nstocks-1))
        }
      }
    }
  }
}
if(mixing == 3)
{
  mix=read.csv(paste(wd,'/Move_Fract.csv',sep=''),row.names=1)
  T=array(0,dim=c(nstocks,nstocks,nages))
  for(a in 1:nages)
  {	
    for(i in 1:nstocks)
    {
      for(j in 1:nstocks)
      {
        T[i,j,a]=mix[i,j]*residency[s]
      }
    }
  }
}

if(mixing == 4)
{
  mix=read.csv(paste(wd,'/mixing_matrix.csv',sep=''),row.names=1)
  T=array(0,dim=c(nstocks,nstocks,nages))
  for(a in 1:nages)
  {	
    for(i in 1:nstocks)
    {
      for(j in 1:nstocks)
      {
        T[i,j,a]=mix[(a*nstocks)-nstocks+i,j]
      }
    }
  }
}
### INDIVIDUAL SPAWN DATE AND AGE AT RELEASE
age_at_release=function(mark_region)
	{
	if(spawn_dist == 0)
		spawned<<-round(runif(1,min=spawn_par[1],max=spawn_par[2]))
	if(spawn_dist == 1)
		spawned<<-round(rnorm(1,spawn_par),0)	
	if(spawn_dist == 2)
		spawned<<-rpois(1,spawn_par[1])

	if(tag_dist[mark_region] == 0)
		tagged<<-round(runif(1,min=tag_par[mark_region,1],max=tag_par[mark_region,2]))
	if(tag_dist[mark_region] == 1)
		tagged<<-round(rnorm(1,tag_par[mark_region,]),0)	
	if(tag_dist[mark_region] == 2)
		tagged<<-rpois(1,tag_par[mark_region,1])
	return(which(rmultinom(1,1,reg_age[mark_region,])==1)+(tagged-spawned)/365)
	}

### REGIONAL FISHERY SAMPLING DISTRIBUTIONS
effort=function(region)
	{
	if(fishery_dist[region] == 0)
		return(rep(dunif(1:365,fishery_par[region,1],fishery_par[region,2])/sum(dunif(1:365,fishery_par[region,1],fishery_par[region,2])),nyrs))
	if(fishery_dist[region] == 1)
		return(rep(dnorm(1:365,fishery_par[region,])/sum(dnorm(1:365,fishery_par[region,])),nyrs))
	if(fishery_dist[region] == 2)
		return(rep(dpois(1:365,fishery_par[region,1])/sum(dpois(1:365,fishery_par[region,1])),nyrs))
	}

### INDIVIDUAL TRAJECTORIES, CURRENTLY PROGRAMMED FOR ANNUAL TIMESTEP
fate=function(mark_region,cohort)
	{
	age1=age_at_release(mark_region)
	time1=(tagged-1)/365
	time2=ifelse(tagged<fishery_par[mark_region,1],(fishery_par[mark_region,1]-tagged)/365,0)
	time3=ifelse(tagged<fishery_par[mark_region,1],(fishery_par[mark_region,2]-fishery_par[mark_region,1]+1)/365,
		ifelse(tagged>fishery_par[mark_region,2],0,(fishery_par[mark_region,2]-tagged+1)/365))
	time4=ifelse(tagged>fishery_par[mark_region,2],(365-tagged+1)/365,(365-fishery_par[mark_region,2])/365)
#	if(tagged<tstep)
#		{
		timesteps=trunc((1+nyrs-cohort)*365/tstep)
#		}
#	if(tagged>tstep)
#		{
#		timesteps=trunc((1+nyrs-cohort-time1)*365/tstep)+1
#		}
	cap_hist=matrix(NA,timesteps,17)
	colnames(cap_hist)=c('age_tagged','cohort','time','age','mark_region','region','hand_mort','doubletag','tag1','tag2','lambda','survive','capture','report','tag1_post','tag2_post','lambda_post')
	region1=ifelse(mix_event1 == 1, which(rmultinom(1,1,T[mark_region,,max(1,age1)])==1), mark_region)
	handle1=rbinom(1,1,handling_mort[trunc(max(1,age1))])
	doubletag1=rbinom(1,1,prop_double)
	loss1=rbinom(1,1,tag_loss1[trunc(max(1,age1))])
	loss2=ifelse(doubletag1==0,1,rbinom(1,1,tag_loss1[trunc(max(1,age1))]))
	chronic1=rbinom(1,1,1-exp(-tag_loss2*(time2+time3/2)))
	chronic2=ifelse(doubletag1==0,1,rbinom(1,1,1-exp(-tag_loss2*(time2+time3/2))))
	chronic3=ifelse(chronic1==1,1,rbinom(1,1,1-exp(-tag_loss2*(time3/2+time4))))
	chronic4=ifelse(chronic2==1,1,
		ifelse(doubletag1==0,1,rbinom(1,1,1-exp(-tag_loss2*(time3/2+time4)))))
	lambda1=ifelse(doubletag1==0,(1-handle1)*(1-loss1)*(1-chronic1),(1-handle1)*(1-loss1*loss2)*(1-chronic1*chronic2))
	lambda2=ifelse(lambda1==0,0,
		ifelse(doubletag1==0,1-chronic3,1-chronic3*chronic4))
	M_before=M[max(1,age1)]*time2
	M_during=M[max(1,age1)]*time3
	M_after=M[max(1,age1)]*time4
	F_1=F[(nyrs*(region1-1)+cohort),trunc(max(1,age1))]*sum(effort(region1)[tagged:tstep])
	if(model_type == 1)
		{
		survive1=ifelse(handle1==1,0,rbinom(1,1,exp(-(M_before))))
		survive2=ifelse(survive1==0,0,rbinom(1,1,exp(-(F_1+M_during))))
		capture1=ifelse(survive1==1&survive2==0,rbinom(1,1,F_1/(F_1+M_during)),0)
		survive3=ifelse(survive2==0,0,rbinom(1,1,exp(-(M_after))))
		}
	if(model_type == 2)
		{
		survive1=ifelse(handle1==1,0,rbinom(1,1,exp(-(M_before+M_during/2))))
		capture1=ifelse(survive1==0,0,rbinom(1,1,1-exp(-F_1)))
		survive2=ifelse(survive1==0,0,rbinom(1,1,exp(-(M_during/2+M_after))))
		}
	report1=ifelse(lambda1==0,0,ifelse(capture1==0,0,rbinom(1,1,reporting[mark_region])))
	cap_hist[1,]=c(age1,cohort,time1,trunc(age1),mark_region,region1,handle1,doubletag1,(1-loss1)*(1-chronic1),(1-loss2)*(1-chronic2),lambda1,survive1*survive2*survive3,capture1,report1,(1-loss1)*(1-chronic1)*(1-chronic3),(1-loss2)*(1-chronic2)*(1-chronic4),lambda1*lambda2)	
	cap_hist[,1]=max(1,trunc(age1,2))
	cap_hist[,2]=cohort
	cap_hist[,3]=round(time1+(0:(timesteps-1))*tstep/365,2)
	cap_hist[,5]=mark_region
	cap_hist[,7]=handle1
	cap_hist[,8]=doubletag1
	if(timesteps>1)
		{
		for(t in 2:timesteps)
			{
			transient=rbinom(1,1,pr_transient)
			cap_hist[t,4]=trunc(cap_hist[t-1,4]+tstep/365)
			cap_hist[t,6]=ifelse(natal_homing==1&transient==0,
				which(rmultinom(1,1,T[mark_region,,min(max_age,trunc(cap_hist[t,4]))])==1),
				which(rmultinom(1,1,T[cap_hist[t-1,6],,min(max_age,trunc(cap_hist[t,4]))])==1))
			t_before=(fishery_par[cap_hist[t,6],1]-1)/tstep
			t_during=(fishery_par[cap_hist[t,6],2]-fishery_par[cap_hist[t,6],1]+1)/tstep
			t_after=(tstep-fishery_par[cap_hist[t,6],2])/tstep
			cap_hist[t,9]=rbinom(1,1,exp(-tag_loss2*(t_before+t_during/2)))*cap_hist[t-1,15]
			cap_hist[t,10]=ifelse(cap_hist[t,8]==1,rbinom(1,1,exp(-tag_loss2*(t_before+t_during/2))),0)*cap_hist[t-1,16]
			cap_hist[t,11]=min(1,cap_hist[t,9]+cap_hist[t,10])*cap_hist[t-1,17]
			M_prior=M[min(max_age,trunc(cap_hist[t,4]))]*t_before
			M_t=M[min(max_age,trunc(cap_hist[t,4]))]*t_during
			M_post=M[min(max_age,trunc(cap_hist[t,4]))]*t_after
			F_t=F[(nyrs*(cap_hist[t,6]-1)+(cohort+trunc(cap_hist[t,3]))),min(max_age,trunc(cap_hist[t,4]))] #*sum(effort(cap_hist[t,6]))[((time1+1+(t-1)*tstep):(time1+(t*tstep)))])
			if(model_type == 1)
				{
				survive_prior=rbinom(1,1,exp(-M_prior))*cap_hist[t-1,12]
				survive_during=ifelse(survive_prior==0,0,rbinom(1,1,exp(-(M_t+F_t))))
				survive_post=ifelse(survive_during==0,0,rbinom(1,1,exp(-M_post)))
				cap_hist[t,12]=survive_prior*survive_during*survive_post
				cap_hist[t,13]=ifelse(survive_prior==1&survive_during==0,rbinom(1,1,F_t/(F_t+M_t)),0)
				}
			if(model_type == 2)
				{
				survive_prior=rbinom(1,1,exp(-(M_prior+M_during/2)))*cap_hist[t-1,12]
				captured_during=ifelse(survive_prior==0,0,rbinom(1,1,1-exp(-F_t)))
				survive_post=ifelse(survive_prior==0,0,rbinom(1,1,exp(-(M_during/2+M_post))))
				cap_hist[t,12]=survive_prior*survive_post
				cap_hist[t,13]=captured_during
				}
			cap_hist[t,14]=ifelse(cap_hist[t,11]==0,0,ifelse(cap_hist[t,13]==0,0,rbinom(1,1,reporting[cap_hist[t,6]])))
			cap_hist[t,15]=rbinom(1,1,exp(-tag_loss2*(t_during/2+t_after)))*cap_hist[t,9]
			cap_hist[t,16]=rbinom(1,1,exp(-tag_loss2*(t_during/2+t_after)))*cap_hist[t,10]
			cap_hist[t,17]=min(1,cap_hist[t,15]+cap_hist[t,16])*cap_hist[t,11]
			}
		}
	if(model_type == 1)
		{
		if(sum(cap_hist[,12])==timesteps)
			{
			return(cap_hist[timesteps,])
			}else
				{
				return(cap_hist[which(sapply(1:timesteps,function(i)cap_hist[i,12])==0)[1],])
				}
		}	
	if(model_type == 2)
		{
		return(c(cap_hist[1,c(1,2,5)],rep(0,cohort-1),sapply(1:timesteps,function(e)cap_hist[e,6]*cap_hist[e,14])))
		}
	}

#pb <- winProgressBar(title=paste("Simulation Progress Bar - Iteration",V,"of",ntrials), label="0% done", min=0, max=100, initial=0)

if(model_type == 1)
	{
	data=matrix(NA,sum(marks)*ncohorts,17)
	colnames(data)=c('age_tagged','cohort','time','age','mark_region','region','hand_mort','doubletag','tag1','tag2','lambda','survive','capture','report','tag1_post','tag2_post','lambda2')
	for(i in 1:ncohorts)
		{
		for(j in 1:nstocks)
			{
			if(i==1&j==1)
				data[1:marks[j],]=t(replicate(marks[j],fate(j,i)))
			if(i==1&j>1)
				data[(sum(marks[1:(j-1)])+1):sum(marks[1:j]),]=t(replicate(marks[j],fate(j,i)))
			if(i>1&j==1)
				data[((i-1)*sum(marks)+1):((i-1)*sum(marks)+sum(marks[1])),]=t(replicate(marks[j],fate(j,i)))
			if(i>1&j>1)
				data[((i-1)*sum(marks)+sum(marks[1:(j-1)])+1):((i-1)*sum(marks)+sum(marks[1:j])),]=t(replicate(marks[j],fate(j,i)))
	       #info <- sprintf("%d%% done", round(((nstocks*(i-1)+j)/(ncohorts*nstocks))*100))
	       #setWinProgressBar(pb, ((nstocks*(i-1)+j)/(ncohorts*nstocks))*100, label=info)
			}
		}
	}
if(model_type == 2)
	{
	data=matrix(NA,sum(marks)*ncohorts,3+nyrs*tstep/365)
	colnames(data)=c('age_tagged','cohort','mark_region',paste('event',1:(nyrs*tstep/365),sep=''))
	for(i in 1:ncohorts)
		{
		for(j in 1:nstocks)
			{
			if(i==1&j==1)
				data[1:marks[j],]=t(replicate(marks[j],fate(j,i)))
			if(i==1&j>1)
				data[(sum(marks[1:(j-1)])+1):sum(marks[1:j]),]=t(replicate(marks[j],fate(j,i)))
			if(i>1&j==1)
				data[((i-1)*sum(marks)+1):((i-1)*sum(marks)+sum(marks[1])),]=t(replicate(marks[j],fate(j,i)))
			if(i>1&j>1)
				data[((i-1)*sum(marks)+sum(marks[1:(j-1)])+1):((i-1)*sum(marks)+sum(marks[1:j])),]=t(replicate(marks[j],fate(j,i)))
	       #info <- sprintf("%d%% done", round(((nstocks*(i-1)+j)/(ncohorts*nstocks))*100))
	       #setWinProgressBar(pb, ((nstocks*(i-1)+j)/(ncohorts*nstocks))*100, label=info)
			}
		}
	}
#close(pb)
#write.csv(data,paste(wd,'/simulated_data.csv',sep=''),row.names=FALSE)

### ESTIMATION OF TAG LOSS FROM DOUBLE-TAGGING
if(prop_double == 0)
	{
	tag_loss=matrix(0,nages,2)
	colnames(tag_loss)=c('double_tagged','single_return')
	tag_loss_vector=rep(0,nages)
	}else
	{
		doubletags=as.data.frame(subset(data,data[,8]==1&data[,14]==1))
		doubletags$tagloss=ifelse(doubletags[,9]==0|doubletags[,10]==0,1,0)
		if(sum(tag_loss1,tag_loss2) == 0)
			{
			doubles=table(trunc(doubletags$age_tagged),doubletags$tagloss)
			tag_loss=cbind(doubles,rep(0,length(doubles[,1])))
			colnames(tag_loss)=c('double_tags','single_return')
			tag_loss_vector=sapply(1:length(tag_loss[,1]),function(a)tag_loss[a,2]/sum(tag_loss[a,1:2]))
			}else
				{
				doubletags$tagloss=ifelse(doubletags[,9]==0|doubletags[,10]==0,1,0)
				tag_loss=table(trunc(doubletags$age_tagged),doubletags$tagloss)
				colnames(tag_loss)=c('double_tagged','single_return')
				tag_loss_vector=sapply(1:length(tag_loss[,1]),function(a)tag_loss[a,2]/sum(tag_loss[a,1:2]))
				}
	}

#################################################################################
###		Summary Matrices for Estimation Model
#################################################################################

releases=aggregate(report~trunc(age_tagged)+cohort+mark_region,data=data,length)[,c(3,2,1,4)]
colnames(releases)=c('mark_region','year_release','age_tagged','number_of_marks')
#write.csv(releases,paste(wd,'tag_releases.csv',sep=''),row.names=FALSE)
#write.csv(releases,'release_matrix.csv',row.names=FALSE)

returns=aggregate(report~region+trunc(time+cohort)+trunc(age_tagged)+cohort+mark_region,data=data,sum)[,c(5,4,3,2,1,6)]
colnames(returns)=c('mark_region','year_release','age_tagged','year_recapture','recapture_region','number_of_returns')
#write.csv(returns,paste(wd,'tag_returns.csv',sep=''),row.names=FALSE)
	
full_matrix=expand.grid(mark_region=1:nstocks,year_release=1:ncohorts,year_recapture=1:nyrs,age_tagged=1:nages,recapture_region=1:nstocks)
return_mat=merge(returns,full_matrix,all=TRUE)
return_mat[is.na(return_mat)]=0
return_matrix=sapply(1:nstocks,function(s)subset(return_mat,recapture_region==s)$number_of_returns)
#write.csv(return_matrix,'return_matrix.csv',row.names=FALSE)
	

