#rm(list=(ls()))

#################################################################################
###		DEFINE WORKING DIRECTORY 
#################################################################################

#wd='C:/users/matthew.lauretta/desktop/Yellowfin_Base/yellowfin_continuous_model/'
#wd<<-getwd()

setwd(wd)

#################################################################################
###		LOAD REQUIRED PACKAGES 
#################################################################################

suppressWarnings(suppressMessages(require(PBSmodelling)))
suppressWarnings(suppressMessages(require(beanplot)))
suppressWarnings(suppressMessages(require(matrixStats)))
suppressWarnings(suppressMessages(require(TeachingDemos)))
suppressWarnings(suppressMessages(require(snowfall)))
suppressWarnings(suppressMessages(library(parallel)))
suppressWarnings(suppressMessages(library('snow')))
suppressWarnings(suppressMessages(library(foreach)))
suppressWarnings(suppressMessages(library(doSNOW)))

#################################################################################
###		SET UP TAGGING SIMULATION AND TEST INITIAL RUN
#################################################################################
source('Simulation_control.r')
if(run_trial_SIM==TRUE)
{
source("Program_files/Simulation_model_test.r")
}
#################################################################################
###		RUN MULTIPLE ITERATIONS WITH BROWNIE ESTIMATION MODEL
#################################################################################

source('Brownie_control.r')
source("Program_files/Run_SIM.r")
source("Program_files/write_data.r")
source("Program_files/run_estimation.r")

if(do_parallel==TRUE)
{
  
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd","estimation_model_switch","run.SIM","write.data","run.estimation"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste("Simulation Model Progress Bar"), label=paste("Simulation Run 0 of ",ntrials,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/ntrials*100),label=paste("Simulation Run", n,"of", ntrials,"Completed"))
  opts<-list(progress=progress)
  estimation_model_switch <<-estimation_model_switch
  
  t<- foreach(i=1:ntrials,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
          #,.export=c("wd","estimation_model_switch","run.SIM","write.data","run.estimation")
          ) %dopar% {
            run.SIM(i,wd,estimation_model_switch) 
            write.data(i,wd)
            run.estimation(i,wd,estimation_model_switch) 
          }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()

}else
{
 source("Program_files/Run_SIM.r")
 quiet<-lapply(1:ntrials,function(i)run.SIM(i,wd,estimation_model_switch))
 print("Simulation Runs Finished")
 source("Program_files/write_data.r")
 quiet<-lapply(1:ntrials,function(i)write.data(i,wd))
 print("Data Files Written")
 if(estimation_model_switch==1)
 {      
   source("Program_files/run_estimation.r")
   quiet<-lapply(1:ntrials,function(i)run.estimation(i,wd,estimation_model_switch))
   print("Estimation Model Runs Finished")
   
 }
}

#################################################################################
###		COMPILE RESULTS AND CREATE SUMMARY PLOTS 
#################################################################################

source("Program_files/Simulation_summary.r")
print("Graphics Finished....Simulation Complete...Its Time for a Margarita")

#################################################################################
###		END OF RUN FILE 
#################################################################################

