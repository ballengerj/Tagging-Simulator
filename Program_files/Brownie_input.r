

###########################################################################################
###		Set Up Input Data for ADMB Estimation Model
###########################################################################################
setwd(wd)
source("Simulation_control.r")
setwd("Program_files")
source(paste(wd,"/Brownie_control.r",sep=""))

###############################################################################################################################
## Data Manipulation
###############################################################################################################################

dummy<-(-1345)
max_age_rel=max_age
ages=seq(1,nages,1)	

max_age_recap=rep(nages,times=nyrs)
  if(any(reg_age[1,]==0))
   {
      for(i in 1:nyrs)
       { 
        max_age_recap[i]=min((ages[(which.max(reg_age[1,]==0))-1]+(i-1)),nages)
       }
    }

yrs=seq(1,nyrs,1)

SIM_T1<-matrix(t(T1),ncol=nstocks,byrow=TRUE)
SIM_T1<-matrix(sapply(1:(nstocks*nages),function(i)rep(SIM_T1[i,],times=nyrs)),byrow=TRUE,ncol=nstocks)
SIM_T2=matrix(NA,nstocks*ncohorts*nages*nyrs,nstocks)
for(i in 1:nstocks)
{
  for(a in 1:ncohorts)
  {
   SIM_T2[c(((i-1)*nyrs*nages*ncohorts+(a-1)*nyrs*nages+1):((i-1)*nyrs*nages*ncohorts+(a)*nyrs*nages)),]=SIM_T1[c(((i-1)*nyrs*nages+1):((i)*nyrs*nages)),]
  }
}

SIM_M<-matrix(rbind(do.call(rbind,replicate(nyrs*nstocks,list(M)))),
nrow=nstocks*nyrs)

ntags<-matrix(releases[,4],ncol=1)

ntags1<-matrix(NA,ncohorts*nstocks,nages)
ntags2<-matrix(NA,ncohorts*nstocks*nages)
Reg_Age<-matrix(NA,ncohorts*nstocks*nages)

	for(i in 1:nstocks)
	{
        Reg_Age[(i*nages*ncohorts-(nages*ncohorts-1)):(i*nages*ncohorts)]<-c(rep(reg_age[i,],times=ncohorts))
       } 
count<-0
for(i in 1:(nages*ncohorts*nstocks))
   {
    if(Reg_Age[i]==0) {
        ntags2[i]=0
        count<-count+1
      } else {
        ntags2[i]=ntags[i-(count)]
      }
   }

ntags1<-matrix(ntags2,byrow=TRUE,ncol=nages)
total_recap1<-matrix(NA,ncohorts*nstocks*nages)
for(i in 1:(nages*ncohorts*nstocks)){
total_recap1[i]<-sum(OBS_tag_rec[(i+(i-1)*(nyrs-1)):((i+(i-1)*(nyrs-1))+(nyrs-1)),1:nstocks])}
total_recap<-matrix(total_recap1,byrow=TRUE,ncol=nages)
OBS_tag_not_rec<-ntags1-total_recap


Move_Fract<-matrix(0,nrow=nstocks,ncol=nstocks)
if(movement_switch>3)
  {
   Move_Fract = read.csv(paste(wd,'/Move_Fract.csv',sep=''),row.names=1)
  }


Fixed_T<-matrix(0,nrow=nstocks*ncohorts*nages*nyrs,ncol=nstocks)
if(movement_switch==(-2))
  {
   Fixed_T = read.csv(paste(wd,'/Fixed_T.csv',sep=''),row.names=1)
  }

Fixed_F<-matrix(0,nrow=nstocks*nyrs,ncol=nages)
if(F_switch==(-2))
{
  Fixed_F = read.csv(paste(wd,'/Fixed_F.csv',sep=''),row.names=1)
}

Fixed_M_region<-matrix(0,nrow=nstocks*nyrs,ncol=nages)
if(M_switch==(-3))
{
  Fixed_M_region = read.csv(paste(wd,'/Fixed_M_region.csv',sep=''),row.names=1)
}

neff = matrix(neff_value,ncohorts*nstocks,nages)

###########################################################################################
############### create list to use as ADMB .dat file including basic inputs ###############
###########################################################################################

Brownie<-list(
########## Basic Inputs
"#ncohorts",
ncohorts,
"#nages",
nages,
"#max_age_release",
max_age_rel,
"#nyrs",
nyrs,
"#nstocks",
nstocks,
"#yrs",
yrs,
"#ages",
ages,
"#tagging_month #==(0) assume tagging occurs on Jan 1st (i.e., full year of mortality in release year), ==1-11 mortality is discoundted by (1-tagging_month/12), ==(12) no mortality in year of release",
tagging_month,
"#fish_start_month #==(0) fishing starts on Jan 1st, ==(12) pulse fishery on Dec 31st  //must be between 0 and 12, and less than fish_end_month",
fish_start_month,
"#fish_end_month #==(0) pulse fishery on Jan 1st, ==(12) fishing ends on Dec 31st  //must be between 0 and 12, and greater than fish_start_month",
fish_end_month,
"#max_age_recap # vector defining max age for which recaptures are possible in each year ==nages if plus_group==TRUE",
max_age_recap,
############# Parameter Estimation Phases
"#phase_T_est_CNST //use if movement_switch==1",
phase_T_est_CNST,
"#phase_T_est //use if movement_switch==2",
phase_T_est,
"#phase_T_est_AGE //use if movement_switch==3",
phase_T_est_AGE,
"#phase_T_est_AGE_no_age1  //use if no_move_first_yr==1 and movement_switch==4 (year by region and age)",
phase_T_est_AGE_no_age1,
"#phase_T_res //use if movement_switch==5",
phase_T_res,
"#phase_T_res_age //use if movement_switch==6",
phase_T_res_age,
"#phase_T_res_no_age1 //use if movement_switch==7",
phase_T_res_no_age1,
"#phase_report_rate",
phase_report_rate,
"#phase_F",
phase_F,
"#phase_M_region",
phase_M_region,
"#phase_M_age",
phase_M_age,
"#phase_M_CNST",
phase_M_CNST,
"#phase_M_year",
phase_M_year,
############## Movement Switches and Pens
"#report_rate_switch //#==(-2) fix at Fixed_report_rate (i.e., fix at assumed reporting rate), ==(-1) fix at SIM_reporting rate (i.e., fix at true reporting rate by region, which is the catch weighted fleet specific reporting by region), ==(1)  estimate time and age-invariant region-specific reporting rate",
report_rate_switch,
"#initial_report_rate",
initial_report_rate,
"#report_rate_fixed",
report_rate_fixed,
"#report_rate_pen_switch    //#==0 no prior, ==1 normal prior (use Report_rate_sigma, report_rate_ave), ==2 penalty (use report_rate_pen_low, report_rate_pen_hi)",
report_rate_pen_switch,
"#report_rate_sigma",
report_rate_sigma,
"#report_Rate_ave",
report_rate_ave,
"#report_rate_pen_hi     //# value above which reporting rate pen activated ",
report_rate_pen_hi,
"#report_rate_pen_low  //# value below which reporting rate pen activated ",
report_rate_pen_low,
"#report_rate_pen_mult",
report_rate_pen_mult,
"#F_switch  //#==(-2) fix at Fixed_F (i.e., fix at assumed F),requires that a Fixed_F.csv file is in WD (same dimensions as SIM_F); ==(-1) fix at SIM_F (i.e., fix at true F), ==(0) F estimated for all ages/years/regions (if no observed recaptures for a given age/year/region combo, then F estimates may become highly biased; e.g., if total nages > number of ages of released fish or if only a few fish released at a given age), ==(1) (DO NOT USE; NOT YET IMPLEMNTED) holes in F estimates (need ragged array) ",
F_switch,
"#initial F",
initial_F,
"#Fixed_F",
Fixed_F,
"#F_pen_switch    // # ==0 no prior, ==1 normal prior (use F_sigma, F_ave), ==2 penalty (F_pen_hi)",
F_pen_switch,
"#F_sigma",
F_sigma,
"#F_ave",
F_ave,
"#F_pen_hi     //# value above which F pen activated",
F_pen_hi,
"#F_pen_low     //# value below which F pen activated",
F_pen_low,
"#F_pen_mult",
F_pen_mult,
"#M_switch  //#==(-3) fix at region-specific Fixed_M_region (i.e., fix at assumed M that varies by region), requires that a Fixed_M_region.csv file is in WD (dimensions of nstocks X nages), ==(-2) fix at region-invariant Fixed_M (i.e., fix at assumed M), ==(-1) fix at SIM_M (i.e., fix at true M), ==(1)  estimate one M for all years, ages, and regions (use M_est), ==(2)  estimate time-invariant M by age that is constant by region (use M_age), ==(3)  estimate time-invariant M by age and region (use M_region)",
M_switch,
"#initial_M",
initial_M,
"#Fixed_M",
Fixed_M,
"#Fixed_M_region",
Fixed_M_region,
"#M_pen_switch // # ==0 no prior, ==1 normal prior (use M_sigma, M_ave), ==2 penalty (M_pen_high, M_pen_Low)",
M_pen_switch,
"#M_sigma",
M_sigma,
"#M_ave",
M_ave,
"#M_pen_low",
M_pen_low,
"#M_pen_high",
M_pen_high,
"#M_pen_mult",
M_pen_mult,
"#movement_switch  // #==(-2) fix at assumed (i.e., fix at Fixed_T), requires that a Fixed_T.csv file is in WD (same dimensions as SIM_T), ==(-1) fix at SIM_T value (i.e., fix at true value), ==(0) no movement (i.e., set residency to 100%), ==(1) est time and age-invariant stymmetric movement (use T_est_CNST), ==(2) est time and age-invariant region-specific T (use T_EST), ==(3) est region-specific and age-specific T (use T_est_AGE), ==(4) Estimate time and age-invariant Residency (use T_res), must have a Move_Fract.csv in working directory (nstock X nstock matrix) to define off-diagnol movement, ==(5) Estimate age-varying Residency (use T_res_age), must have a Move_Fract.csv in working directory (nstock X nstock matrix) to define off-diagnol movement",
movement_switch,
"#natal_homing    //# ==0 no natal homing (markovian movement==random, memoryless), ==1 natal homing (follow movement pattern of natal stock at all times)",
natal_homing,
"#no_move_first_yr  //==0 allow movement in event of release, ==1 do not allow movement in first event",
no_move_first_yr,
"#Move_Fract  //User Defined Movement Proportions for Residency Movement switches (determines how non-resident movement is split among regions; rows should sum to 2 with diagnols==1)",
Move_Fract,
"#initial T",
initial_T,
"#Fixed_T",
Fixed_T,
"#move_pen_switch  // ==0 no prior, ==1 normal prior (use T_ave and T_sigma), ==2 penalty (use T_low, T_hi)",
move_pen_switch,
"#T_sigma",
T_sigma,
"#T_ave",
T_ave,
"#T_low",
T_low,
"#T_hi",
T_hi,
"#Tpen_mult",
Tpen_mult,
"#tag_neff_switch  // #==0 effective sample size is constant across cohorts and regions (use neff_value; typically value~100, if too small penalty/priors may be given too much influence), ==1 set effective sample size==ntags (true multinomial, but may overemphasize larger release events), Probably needs further investigation for straight tagging models (setting multinomial sample size equal to # tags vs. constant had little impact in tag-integrated models; see Goethel et al., 2015)",
tag_neff_switch,
"#neff",
neff,
"#TagCst",
TagCst,
"#tag_loss_switch  //==0 tag loss==0, ==1, tag loss based on input parameter",
tag_loss_switch,
"#handling_mortality_switch  //==0 handling mort==0, ==1, based on handling mort parameter",
handling_mortality_switch,
################# Data Inputs
"#OBS_Tag_Loss",
OBS_Tag_Loss,
"#OBS_Hand_Mort",
OBS_Hand_Mort,
############### Constants
"#ntags1",
ntags1,
######################### Observed Data
"#OBS_tag_rec",
OBS_tag_rec,
"#OBS_total_tag_recap",
total_recap,
"#OBS_tag_not_rec",
OBS_tag_not_rec,
########################## Simulation 'True' Values
"#SIM_Mixing",
mixing,
"#SIM_T",
SIM_T2,
"#SIM_M",
SIM_M,
"#SIM_F",
SIM_F,
"#SIM_report_rate",
SIM_report_rate,
"#dummy",
dummy,
"#phase_dummy1",
phase_dummy1
)

###############################################################################################
##### NEED TO DELETE .DAT FILE EACH TIME RUN R CODE OTHERWISE JUST APPENDS EXISTING FILE ######
###############################################################################################

