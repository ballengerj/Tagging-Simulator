write.data<-function(ntrial,WD)
{

  setwd(WD)
  source("Simulation_control.r")
  setwd("Program_files")
  invisible(shell("del Brownie.dat"))
  source(paste(WD,"/Brownie_control.r",sep=""))
  
  ###############################################################################################################################
  ## Read Sim Values 
  ###############################################################################################################################
  
  SIM_F<<- read.csv(paste("SIM_F",ntrial,".csv",sep=""))		
  M<<-read.csv(paste("SIM_M",ntrial,".csv",sep=""))		
  T1<<- read.csv(paste("SIM_T",ntrial,".csv",sep=""))		
  SIM_report_rate<<-read.csv(paste("SIM_B",ntrial,".csv",sep=""))	
  OBS_Tag_Loss<<- read.csv(paste("SIM_Tag_Loss",ntrial,".csv",sep=""))	
  OBS_Hand_Mort<<- read.csv(paste("SIM_Handling_Mortality",ntrial,".csv",sep=""))	
  releases<<- read.csv(paste("SIM_Releases",ntrial,".csv",sep=""))
  OBS_tag_rec<<- read.csv(paste("SIM_return_matrix",ntrial,".csv",sep=""))
  neff_value<<-read.csv(paste("SIM_neff",ntrial,".csv",sep=""))	
  
  source("Brownie_input.r")
  invisible(lapply(Brownie,function(x)write.table(x,paste0("Brownie",ntrial,".dat",sep=""),append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)))
  
  invisible(shell(paste("del SIM_F",ntrial,".csv",sep="")))		
  invisible(shell(paste("del SIM_M",ntrial,".csv",sep="")))		
  invisible(shell(paste("del SIM_T",ntrial,".csv",sep="")))		
  invisible(shell(paste("del SIM_B",ntrial,".csv",sep="")))	
  invisible(shell(paste("del SIM_Tag_Loss",ntrial,".csv",sep="")))	
  invisible(shell(paste("del SIM_Handling_Mortality",ntrial,".csv",sep="")))	
  invisible(shell(paste("del SIM_neff",ntrial,".csv",sep="")))	
  #shell(paste("copy Brownie.dat Brownie",ntrial,".dat", sep=""))
  setwd(WD)
rm(list=(ls()[!ls() %in% c('wd','estimation_model_switch')]))


}

