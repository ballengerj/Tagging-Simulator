#################################################################################
###		MODEL INPUTS 
#################################################################################
ntrials = 1000						# number of iterations for simulation


tstep = 365							# timestep in days
nstocks = 4						# number of regions
##################################################################
## NCOHORTS==NYRS, can't deal with issues of no releases yet
###############################################################
ncohorts = 4						# years of tag releases
nyrs  = 4							# yrs of recaps

nages = 4							    # number of age classes, the terminal age class is a plus group; each 'at-age' input below must have this many entries per region
                          # MAX(nages)<max_age_release+(nyrs-1) in order to avoid estimating parameters for ages with no data 
#################################################################################
###		TAGGING INPUTS
#################################################################################

# Age composition must have same number of columns as nages
# If releases at every age then improves estimability of mortality at older ages
# When 0 releases, F/M will be estimated up until the max age of release in the given year, in each subsequent year and additional F/M will value be estimated for the new age classes with recaptures up until the plus group (nages) is reached 
########################################################################################
#### DO NOT INPUT 0 followed by nonzero values as this will cause model stability issues because estimating F/M values with no releases/recaptures
###########################################################################################
reg_age = matrix(c(rep(
  c(0.2, 0.2, 0.2, 0.2)
  ,times=nstocks)),ncol=nages,byrow=TRUE)	   # age composition of tag releases

marks = c(2500, 2500, 2500, 2500)		# number of marked individuals per region per year
prop_double = 0						# proportion of fish double-tagged

###################################################################################################################
#####################################################################################################################


spawn_dist = 0						# define the distribution of spawning activity, 0 = uniform, 1 = normal, 2 = poisson
spawn_par = c(120, 210)					# define the spawning season distribution parameters, for uniform start and end day of year, for normal mean and standard dev, for poisson mean

#########################################################################################################################
## MUST ENSURE THAT SPAWNING SEASON AND TAGGING SEASON OVERLAP OTHERWISE IF FISH SPAWNED AFTER TAGGING
## WILL NEVER HAVE FISH RELEASED AT MAX AGE AND MODEL WILL CRASH
## May Have OPPOSITE ISSUE IF TAGGING AFTER SPAWNING
####################################################################################################################################

tag_dist = c(0, 0, 0, 0)				# define the distribution of tagging effort in region 1, 0 = uniform, 1 = normal, 2 = poisson
tag_par = matrix(c(
  150, 240,					# define the tagging season distribution parameters each region, for uniform start and end day of year, for normal mean and standard dev, for poisson mean
  150, 240, 
  150, 240,
  150, 240
),ncol=2,byrow=TRUE)
fishery_dist = c(0, 0, 0, 0)					# define the distribution of fishing effort in each region, 0 = uniform, 1 = normal, 2 = poisson
fishery_par = matrix(c(
  60, 300,				# define the fishing distribution parameters for each region , for uniform start and end day of year, for normal mean and standard dev, for poisson mean
  60, 300,
  60, 300,
  60, 300
),ncol=2,byrow=TRUE)

nfleets = 5							# number of fleets
fleet1 = 'HB'						# fleet names
fleet2 = 'GN'
fleet3 = 'HL'
fleet4 = 'LL'
fleet5 = 'PS'

catches = read.csv('catch_matrix.csv',row.names=1)

#################################################################################
###		PARAMETERS
#################################################################################

M = c(.274, .230, .206, .191)			# annual natural mortality-at-age
F = read.csv('F_matrix.csv',row.names=1)		# annual regional fishing mortality-at-age
tag_loss1 = c(0, 0, 0, 0)			# type 1 discrete tag loss from the handling process
tag_loss2 = 0						# type 2 chronic tag loss from tissue deterioration, biofouling etc.
handling_mort = c(0, 0, 0, 0)			# discrete handling mortality rates
report_rates = c(.5, .5, .5, .5, .5)			# fleet reporting rates

#################################################################################
###		MIXING ASSUMPTIONS
#################################################################################

mixing = 2							
# mixing is the probility of moving at each time step
# 0=none 
# 1=uniform 
# 2=resident with symmetric off-diagnol movement 
# 3=resident with user-defined off_diagnol movement using Move_Fract.csv 
# 4=user defined using mixing_matrix.csv
residency = c(.6, .6, .6, .6)						# rate of residency of fish in home region (in numerical order by region), only used when mixing = 2 or 3
mix_event1 = 0						# mixing on intitial sampling event, 0 = no (mixing occurs at the beginning of event 2), 1 = yes (mixing can occur instantly after release)
natal_homing = 1          # natal homing yes==1, markov movement ==0
pr_transient = 0          # probability that an individual doesn't follow natal movement (asuming natal_homing==1) (i.e., fish marked outside of natal region but assumed natal origin)
#################################################################################
###		MODEL TYPE 
#################################################################################

model_type = 1					       # define model type for data formatting, 1 = tag returns (one capture event only), 2 = capture-recapture sampling (multiple capture events)
estimation_model_switch = 1    # ==0 perform simulations only (no estimation), ==1 perform simulations and estimation
run_trial_SIM=FALSE            #TRUE/FALSE Run a trial sim to ensure all inputs are working correctly
do_parallel = TRUE             #TRUE/FALSE Use parallel computing to utilize all computer cores and speed up simulations
#################################################################################
###		END OF SIMULATION MODEL INPUTS AND PARAMETERIZATION 
#################################################################################

#################################################################################
###		Graphical Parameters
####################################################################################
###   Only change if issues with graph text sizes
#################################################################################

cex.text.sub=0.4               # text size for small font on plot
cex.text.main=1.0              # text size for large font on plot
cex.medium=0.75                # text size for inbetween font on plot
bg.text=0.2                    # graphical parameter defining the thickness of text shadow

########################################################################################################
max_age <- nages
