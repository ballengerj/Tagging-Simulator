run.SIM<-function(ntrial,WD,estimation_model_switch1)
{

  setwd(WD)
  source("Simulation_control.r")
  source('Brownie_control.r')
  setwd("Program_files")
  V<<-ntrial
  source("Simulation_model.r")
  
  write.csv(data,paste("individual_fates",ntrial,".csv",sep=""),row.names=FALSE)
  write.csv(releases,paste("SIM_Releases",ntrial,".csv",sep=""),row.names=FALSE)
  write.csv(returns,paste("SIM_Returns",ntrial,".csv",sep=""),row.names=FALSE)
  write.csv(cbind(tag_loss,tag_loss_est=tag_loss_vector),paste("OBS_tag_loss",ntrial,".csv",sep=""),row.names=FALSE)
  write.csv(return_matrix,paste("SIM_return_matrix",ntrial,".csv",sep=""),row.names=FALSE)
  
  write.csv(F,paste("SIM_F",ntrial,".csv",sep=""),row.names=FALSE)		
  write.csv(M,paste("SIM_M",ntrial,".csv",sep=""),row.names=FALSE)		
  write.csv(T,paste("SIM_T",ntrial,".csv",sep=""),row.names=FALSE)		
  write.csv(reporting,paste("SIM_B",ntrial,".csv",sep=""),row.names=FALSE)
  
  write.csv(tag_loss_vector, paste("SIM_Tag_Loss",ntrial,".csv",sep=""), row.names=FALSE)
  write.csv(handling_mort, paste("SIM_Handling_Mortality",ntrial,".csv",sep=""),row.names=FALSE)
  write.csv(neff_value, paste("SIM_neff",ntrial,".csv",sep=""),row.names=FALSE) 
  
  if(estimation_model_switch==0)
  {
    print(paste0("Simulation Model Run ", ntrial," Completed......No Estimation Model Run"))
  }


setwd(WD)

rm(list=(ls()[!ls() %in% c('wd','estimation_model_switch')]))
close(pb)
}

