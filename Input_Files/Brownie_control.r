##########################################################################################
###		ADMB Inputs and Controls
##########################################################################################

################### NOTES ##################################################################################################################################################################################
#Estimation model is run in ADMB, this file is used to parametrize and initialize the estimation model
#The estimation model is run automatically unless estimation_model_switch==0
#The simulator code automatically creates the ADMB .dat file and runs the estimation model
#The ADMB .tpl is provided, but should not be altered except by those familiar with ADMB
#ADMB estimation is achieved in phases where the various parameters become active based on their phase (AT LEAST ONE PARAMETER MUST BE ACTIVE IN PHASE 1)
#The current phase setup is optimal given exploratory tests, but may differ by model (CHANGING PHASES SHOULD BE UNDERTAKEN CAUTIOUSLY BY THOSE UNFAMILIAR WITH ADMB AS DIFFERENT SOLUTIONS MAY BE OBTAINED WITH DIFFERENT PHASE SELECTIONS ESPECIALLY FOR UNSTABLE MODEL PARAMETRIZATIONS)
#Negative phases mean that the paramters will not be estimated (IF A GIVEN PARAMETER SWITCH IS NEGATIVE THE ASSOCIATED PARAMETER PHASE SHOULD BE NEGATIVE)
#A parameter should only be active (POSITIVE PHASE) if the associated parmater switch is activated for that parameter
#Parameters are estimated in log space for estimation stability and initial values provided should also be in log space [i.e., initial parameter value=exp(initial_input); initial_input=ln(desired initial value)]
#Provided initial inputs are uninformed, providing best guesses for each parameter will help stabilize model estimation and improve performance
#The option for a normal prior or penalty term is provided for each paramter, but each should be used carefully
#Normal priors may not be appropriate where little information is available a priori regarding parameter value (i.e., it is only suggested for M and reporting rate where external information (e.g., high reward tagging) might be available to narrow down the potential range of values)
#Penalty terms are better used to stabilize model performance by avoiding unreasonable parameter values (e.g., fishing mortality above some upper limit such 3.0) or flat response surfaces (e.g., movement rates below some minimum threshold that is essentially zero such as 0.00000001)
#Note that both movement rate and reporting rate are estimated as logit tranforms in order to impose a natural bound of (0,1) (exclusive)
#By doing so neither parmater can be exactly 0 or 1
#The penalty term can be useful to avoid flat response surfaces near the bounds (particularly for near 0 movement rates; it is suggested to impose a penalty at an extremely small value that is essentially 0 if it is believed either parameter may venture into this solution space)
#The ability to estimate all parameters simultaneously is limited due to confounding and is highly dependent on data availability and model parametrization (attempting to estimate age-specific M, fishing mortality, reporting rate, and movement is unlikely to lead to stable model estimation; see Lauretta and Goethel, 2016 for more details)
#######################################################################
##F is estimated for each year and age combination of the model (except year/age combinations with no releases, ie age composition==0)
## There must be releases/recaptures in each year and age combination except for ages greater than the max age of release (as defined by the first age with a 0 in the age composition vector)
## Consequently ncohorts MUST EQUAL nyrs, otherwise there will be no releases at the youngest age(s) when year>ncohorts and F will be estimated with no recaptures (the model will be overparametrized)
## Same if there are any years with 0's in the age composition of tags less than the max age of release or if no tags are released in a stock AND fish do not move or are not allowed to move (e.g., if movement_first_yr=0)
#######################################################################################################################################################################################################################################################################################################################################################################

############## Define Assumed Tag Timing #####################################################################################################
tagging_month = c(0,0,0,0,0)
#assumed month that tagging takes place on average for each region
#defines fraction of year that natural mortality occurs in year of release, must be 0-12
#==(0) assume tagging occurs on Jan 1st (i.e., full year of mortality in release year)
#==1-11 mortality is discoundted by (1-tagging_month/12)
#==(12) no mortality in year of release

fish_start_month = c(0,0,0,0,0)
#assumed start of fishing season for each region
########################################################
#CANNOT BE GREATER THAN fish_end_month #################
########################################################
#==(0) assume fishing starts Jan 1
#==(1-11) assume fishing season starts at some point within the year 
#==(12) assume pulse fishery on Dec 31st

fish_end_month = c(12,12,12,12,12)
#assumed start of fishing season for each region
########################################################
#CANNOT BE LESS THAN fish_end_month #################
########################################################
#==(0) assume pulse fishery Jan 1st
#==(1-11) assume fishing season ends at some point within the year 
#==(12) assume fishing ends on Dec 31st
#######################################################################################################################################################################################################################################################################################################################################################################

######### Define How Movement WIll be Estimated ############################################################################################################################################
movement_switch = (2)		
#==(-2) fix at assumed (i.e., fix at Fixed_T), requires that a Fixed_T.csv file is in WD (dimensions of (nstocks*ncohorts*nages*nyrs) X nstocks)
#==(-1) fix at SIM_T value (i.e., fix at true value)
#==(0) no movement (i.e., set residency to 100%)
#==(1) est time and age-invariant SYMMETRIC movement (use T_est_CNST)
#==(2) est time and age-invariant region-specific T (use T_EST)
#==(3) est region-specific and age-specific T (use T_est_AGE)
#==(4) est region-specific and age-specific T, but fish aren't allowed to move in first year (no_move_first_yr==1) of release (use T_est_no_age1); if don't move in the first year then there is no estimate of movement for the first age of release/recapture because fish never move at that age (first movement is at youngest age of release+1)
#==(5) Estimate time and age-invariant Residency (use T_res), must have a Move_Fract.csv in working directory (nstock X nstock matrix) to define off-diagnol movement
#==(6) Estimate age-varying Residency (use T_res_age), must have a Move_Fract.csv in working directory (nstock X nstock matrix) to define off-diagnol movement
#==(7) Estimate age-varying Residency, but fish aren't allowed to move in first year (no_move_first_yr==1) of release (use T_res_no_age1), must have a Move_Fract.csv in working directory (nstock X nstock matrix) to define off-diagnol movement
##############################################################################################################################################

############# Define Active Movement Paramters (see definitions of movement_switch to determine which parameter should be active, all others should be negative)######################################                                                                           
#If movement_switch is negative all following phases should be negative as well
#Reminder: movement_switch=(-2) requires a Fixed_T.csv in the WD
phase_T_est_CNST = (-1)
phase_T_est = (1)
phase_T_est_AGE = (-1)
phase_T_est_AGE_no_age1 = (-1)
phase_T_res = (-1)
phase_T_res_age = (-1)
phase_T_res_no_age1 = (-1)
##############################################################################################################################################################

######## Movement Specifications ##################################################################################################################################
initial_T = (-2)        		                # starting T estimate in log space
natal_homing = (1)		                          # ==0 no natal homing (markovian movement==random, memoryless), ==1 natal homing (follow movement pattern of natal stock at all times)
no_move_first_yr = (1)	   	                # ==0 allow movement in event of release, ==1 do not allow movement in first event
####################################################################################################################

######## Define Movement Prior/Penalty ####################################################################################################################
move_pen_switch = (0)		                    # ==0 no prior, ==1 normal prior (use T_ave and T_sigma), ==2 penalty (e.g., to avoid flat repsonse surface at small T's; use T_low, T_hi)
T_sigma = (0.5)                             # standard deviation defining the normal prior
T_ave = (0.7)                               # average value for the normal prior
T_low = (0.0000001)                         # value below which T penalty activated (typically a very small number; see notes above)
T_hi = (0.999999999)                        # value above which T penalty activated (typically a number very close to 1.0; see notes above)
Tpen_mult = (10)			                        # multiplier for Movement penalty function (to increase the penalty associated with exceeding the bounds)
###################################################################################################################################################################

############ Define How Fishing MOrtality Will Be Estimated ######################################################################################################
F_switch = (1)		
#==(-2) fix at Fixed_F (i.e., fix at assumed F),requires that a Fixed_F.csv file is in WD (same dimensions as SIM_F)
#==(-1) fix at SIM_F (i.e., fix at true F)
#==(1) F estimated for all ages/years/regions up until the  max_age_rel in the given year (or for all nages if age composition is all nonzero)
##########################################################################################################################################################################################

############ F Specifics ##############################################################################################################################################
phase_F = (2)                               # if F_switch negative this must be negative as well, F_switch=(-2) requires a Fixed_F.csv in the WD as well
initial_F = (-2)			                      # starting F estimate in log space
####################################################################################################################################################################################

########## Define Fishing Mortality Prior/Penalty #########################################################################################################################
F_pen_switch<- (1)                          # ==0 no prior, ==1 normal prior (use F_sigma, F_ave), ==2 penalty (use F_pen_hi)
F_sigma = (0.5)                             # standard deviation defining the normal prior
F_ave = (0.5)                               # average value for the normal prior
F_pen_hi = (1.5)                              # value above which F penalty activated (typically a large number representing unreasonable mortality rates, e.g., > 3.0)
F_pen_low = (.0000001)                      # value below which F penalty activated (typically a small number representing essentially 0 mortality)
F_pen_mult = (10)                            # multiplier for F penalty function (to increase the penalty associated with exceeding the bounds)
#############################################################################################################################################################################

########### Define How Reporting Rate Will Be Estimated #####################################################################################################################
#estimated reporting rate is by region and not fleet, therefore model estimated values differ from user input reporting which is by fleet then weighted by catch to get regional reporting in simulation model
report_rate_switch = (1)
#==(-2) fix at Fixed_report_rate (i.e., fix at assumed reporting rate)
#==(-1) fix at SIM_reporting rate (i.e., fix at true reporting rate by region, which is the catch weighted fleet specific reporting by region)
#==(1)  estimate time and age-invariant region-specific reporting rate
#############################################################################################################################################################################################################################

######### Reporting Rate Specifics ##################################################################################################################################
phase_report_rate = (3)                    # If report_rate_switch is negative, this must be negative as well; report_rate_switch=(-2) will use Fixed_report_rate below
initial_report_rate = (-2)                  # starting report_rate estimate in log space
report_rate_fixed = c(.25,.25,.25,.25,.25)	# if report_rate_switch=(-2) use these values for regional reporting (must be length of nstocks and values are input in numerical order by stock number)
################################################################################################################################################################################################################

####### Define Reporting Rate Prior/Penalty ##############################################################################################################################################################################
report_rate_pen_switch = (0)	              # ==0 no prior, ==1 normal prior (use Report_rate_sigma, report_rate_ave), ==2 penalty (use report_rate_pen_low, report_rate_pen_hi)
report_rate_sigma = (0.4)                   # standard deviation defining the normal prior
report_rate_ave = (0.55)                     # average value for the normal prior
report_rate_pen_hi = (0.2)                  # value below which report rate penalty activated (see notes above)
report_rate_pen_low = (0.99999999)          # value above which report rate penalty activated (typically a number very close to 1.0; see notes above)
report_rate_pen_mult = (10)                  # multiplier for reporting rate penalty function (to increase the penalty associated with exceeding the bounds)
#####################################################################################################################################################################

############# Define How Natural Mortality Will Be Estimated ##########################################################################################################################################################################
M_switch = (2)	
#==(-3) fix at region, cohort, and age-specific Fixed_M_region (i.e., fix at assumed M that varies by region, year, and age), requires that a Fixed_M_region.csv file is in WD (dimensions of (nstocks*ncohorts) X nages)
#==(-2) fix at year and region-invariant Fixed_M (i.e., fix at assumed M)
#==(-1) fix at SIM_M (i.e., fix at true M)
#==(1)  estimate one M for all years, ages, and regions (use M_est)
#==(2)  estimate time-invariant M by age that is constant by region (use M_age)
#==(3)  estimate time-invariant M by age and region (use M_region)
#==(4)  estimate M by age, region, year (use M_year)
##############################################################################################################################################################################################################################################

########### Natural Mortality Specifics ######################################################################################################################################################################
#If movement_switch is negative all following phases should be negative as well; M_switch=(-2) will use Fixed_M below; M_switch=(-3) requires Fixed_M_region.csv in WD
phase_M_CNST = (-4)
phase_M_age = (4)
phase_M_region = (-4)
phase_M_year = (-4)
initial_M = (-2)			                     # starting M estimate in log space
Fixed_M = c(0.365,0.306,0.274,0.255,0.243)  # if M_switch=(-2) use these values for age-specific M (constant by region; must be length of nages and values are input in numerical order by age)
#############################################################################################################################################################################################################

######## Natural Mortality Prior/Penalty ##################################################################################################################################################################################
M_pen_switch = (2)		                     # ==0 no prior, ==1 normal prior (use M_sigma, M_ave), ==2 penalty (use M_pen_low, M_pen_hi)
M_sigma = (0.5)                            # standard deviation defining the normal prior
M_ave = (0.6)                              # average value for the normal prior
M_pen_low = (0.1)                          # value below which M penalty activated (see notes above)
M_pen_high = (1.5)                         # value above which M penalty activated (typically a large number representing unreasonable mortality rates, e.g., > 3.0)
M_pen_mult = (10)			                   # multiplier for M penalty function (to increase the penalty associated with exceeding the bounds)
##################################################################################################################################################################################################################

############ Tag Loss and Handling Mortality Specifics ###########################################################################################################################
tag_loss_switch = (0)		                   # ==0, set tag loss==0, ==1 set tag loss to OBS_tag_loss (double-tagging empirical estimate)
handling_mortality_switch = (0)	           # ==0, set handling mortality==0, ==1 set handling mortality to OBS_Hand_Mort
#############################################################################################################################################################################################

######## Likelihood Specifics ##########################################################################################################################################################################
tag_neff_switch = (0)		
#==0 effective sample size is constant across cohorts and regions (use neff_value; typically value~100, if too small penalty/priors may be given too much influence)
#==1 set effective sample size==ntags (true multinomial, but may overemphasize larger release events)
#Probably needs further investigation for straight tagging models (setting multinomial sample size equal to # tags vs. constant had little impact in tag-integrated models; see Goethel et al., 2015)

neff_value = (1000)			                   # input constant effective sample size if tag_neff_switch==0
TagCst = (0.000001)                        # constant used for robust likelihood to avoid log(0) (see Fournier et al., 1990)
#######################################################################################################################################################################################################################
phase_dummy1 = (-1)                         # if want to simulate results with fixed inputs turn all other phases to neg values and make this positive
########## END ####################################################################################################################################################################################

