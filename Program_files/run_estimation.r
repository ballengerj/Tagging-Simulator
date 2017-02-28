run.estimation<-function(ntrial,WD,estimation_model_switch1,...)
{

  setwd(WD)
  source("Simulation_control.r")
  setwd("Program_files")
  source(paste(WD,"/Brownie_control.r",sep=""))
  dir.create(paste0(WD,"/Simulation Results",sep=""))
  dir.create(paste0(WD,"/Simulation Results/Run",ntrial,sep=""))
  dir.create(paste0(WD,"/Figures"))
  dir.create(paste0(WD,"/Figures/Report Files"))
  
  invisible(file.copy(from=paste0(WD,"/Program_files/Brownie.exe",sep=""),to=paste0(WD,"/Simulation Results/Run",ntrial,"/Brownie.exe",sep="")))
  invisible(file.rename(from=paste0(WD,"/Program_files/Brownie",ntrial,".dat",sep=""),to=paste0(WD,"/Simulation Results/Run",ntrial,"/Brownie",ntrial,".dat",sep="")))
  invisible(file.rename(from=paste(WD,"/Program_files/individual_fates",ntrial,".csv",sep=""),to=paste0(WD,"/Simulation Results/Run",ntrial,"/individual_fates",ntrial,".csv",sep="")))
  invisible(file.rename(from=paste(WD,"/Program_files/SIM_Releases",ntrial,".csv",sep=""),to=paste0(WD,"/Simulation Results/Run",ntrial,"/SIM_Releases",ntrial,".csv",sep="")))
  invisible(file.rename(from=paste(WD,"/Program_files/SIM_Returns",ntrial,".csv",sep=""),to=paste0(WD,"/Simulation Results/Run",ntrial,"/SIM_Returns",ntrial,".csv",sep="")))
  invisible(file.rename(from=paste(WD,"/Program_files/OBS_tag_loss",ntrial,".csv",sep=""),to=paste0(WD,"/Simulation Results/Run",ntrial,"/OBS_tag_loss",ntrial,".csv",sep="")))
  invisible(file.rename(from=paste(WD,"/Program_files/SIM_return_matrix",ntrial,".csv",sep=""),to=paste0(WD,"/Simulation Results/Run",ntrial,"/SIM_return_matrix",ntrial,".csv",sep="")))
  
  setwd(paste0(WD,"/Simulation Results/Run",ntrial,sep=""))
  system(paste0("Brownie -nox -ind Brownie",ntrial,".dat",sep=""),wait=TRUE,show.output.on.console=FALSE)  
  
  invisible(shell(paste("copy Brownie.par Brownie",ntrial,".par", sep="")))
  invisible(shell(paste("copy Brownie.rep Brownie",ntrial,".rep", sep="")))
  invisible(shell(paste("copy Brownie.std Brownie",ntrial,".std", sep="")))
  invisible(shell(paste("copy Brownie.cor Brownie",ntrial,".cor", sep="")))
  
  invisible(shell(paste("del Brownie.par", sep="")))
  invisible(shell(paste("del Brownie.rep", sep="")))
  invisible(shell(paste("del Brownie.std", sep="")))
  invisible(shell(paste("del Brownie.cor", sep="")))
  
  invisible(file.copy(from=paste0(WD,"/Simulation Results/Run",ntrial,"/Brownie",ntrial,".rep",sep=""),to=paste0(WD,"/Figures/Report Files/Brownie",ntrial,".rep",sep="")))
  invisible(file.copy(from=paste0(WD,"/Simulation Results/Run",ntrial,"/Brownie",ntrial,".cor",sep=""),to=paste0(WD,"/Figures/Report Files/Brownie",ntrial,".cor",sep="")))
  
   invisible(shell("del admodel.cov"))
   invisible(shell("del admodel.dep"))
   invisible(shell("del admodel.hes")) 
   invisible(shell(paste("del Brownie.b01",sep="")))
   invisible(shell(paste("del Brownie.P01",sep="")))
   invisible(shell(paste("del Brownie.R01",sep="")))
   invisible(shell(paste("del Brownie.b02",sep="")))
   invisible(shell(paste("del Brownie.P02",sep="")))
   invisible(shell(paste("del Brownie.R02",sep="")))
   invisible(shell(paste("del Brownie.b03",sep="")))
   invisible(shell(paste("del Brownie.P03",sep="")))
   invisible(shell(paste("del Brownie.R03",sep="")))
   invisible(shell(paste("del Brownie.bar",sep="")))
   invisible(shell(paste("del Brownie.eva",sep="")))
   invisible(shell(paste("del Brownie.log",sep="")))
   invisible(shell("del fmin.log"))
   invisible(shell("del variance"))

    if(file.exists("Brownie.cor")==TRUE)
    {
     print(paste0("Model Run ",ntrial," Converged"))
    }
    if(file.exists("Brownie.cor")==FALSE)
    {
      print(paste0("Model Run ",ntrial," Did Not Converge"))
    }

setwd(WD)

rm(list=(ls()[!ls() %in% c('wd','estimation_model_switch')]))
}

