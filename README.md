# Tagging-Simulator
Code for Continuous Time Tagging Simulator Used in CJFAS Article (w/ Matt Lauretta)


###############################################################################################################################################################
#################################################################################################################################################################
########## REQUIRED UPDATES TO SIM CODE ############################################################
### More error traps to ensure that user inputs do not create model crashes (eg stop model if vector of dim(M)<nages, etc...)
### Error traps if Brownie inputs do not have proper dimensions (similar to above)
### Check that all user-inputs work together (e.g., if set tagging season to first day of year, but spawning later, then can never have a tagged fish released at max age and model crashes)
### Simulat High reward tags
### Simulate time-varying M
### Simulate time-varying T
##########################################################################################################################################################################################
############################################################################################################################################################################################





# Program RETURN is a spatial Brownie tag return model for simulating large-scale, multi-region tagging programs.
  Written in R and ADMB by Matthew Lauretta and Daniel Goethel at the National Marine Fisheries Service, Southeast Fisheries Science Center.
  Contact matthew.lauretta@noaa.gov for questions or feedback.

The Program consists of two major components
	1.  A simulation model for generating tagging datasets for an open population with any number of defined spatial regions, mixing rates,
		age classes, age-structure of tag releases, tag cohorts, recapture events, and duration (in days) between sampling events; as well as defined 
		spawning, tagging, and fishing seasons to create stochastic age-at-tagging and times-at-liberty of individuals.
	2.  A spatial Brownie tag return model for estimation of key life-history and fishery parameters, including mixing rates, natural mortality, 
		fishing mortality, tag loss, and reporting rates.

NOTE:  
All files to be modified by the user are contained in the main folder, and the core program files contained in the folder, "Program_files" should not be altered.

To RUN the tagging simulation model:
	1.  Open the "Simulation_control" file and define the simulation and population parameters. This file can be opened and modified using any text editor.
		A. MODEL INPUTS
		B. PARAMETERS
		c. MIXING ASSUMPTIONS
		D. TAGGING INPUTS
		E. DATA TYPE
	2.  Two data csv input matrices must be parameterized 
		A. A matrix of catches by gear type and region, titled "catch_matrix.csv"
		B. A fishing mortality matrix by age, year, and region, titled "F_matrix.csv"
	3.  A third input matrix ("mixing_matrix.csv") must be parameterized if user defined mixing estimates are to be used, mixing = 3 in "Simulation_control" file
	    Note: if mixing = 0 (no mixing), mixing = 1 (uniform mixing), or mixing =2 (residency), then the "mixing_matrix" file is not used.
	4.  Open the "Simulation_run" file in R.
		A. Set the working directory for the main folder.
		B. Install required packages.  These are only necessary for parameter estimation and graphics, and not needed for the basic simulation model.
		C. Run the tagging simulation by sourcing the "Simulation_control" and "Simulation_model" files.  This will run the intial tagging simulation to verify
			the simulation component is parameterized and running correctly.
OUTPUTS for the intitial simulation model will be written to the main folder directory, and include:
	1.  A simulated dataset titled, "simulated_data.csv," that contains the fate of each tagged individual.  The following is a description of the headers of that file:
		age_tagged:	the age of the individual at time of tagging
		cohort:		the tagging event in which the individual was released
		time:		the time-at-liberty of the individual, not applicable if the individual died of handling or natural mortality
		age:		the age of the individual at the time of recapture or end of the study duration, not applicable if the fish died of handling or natural mortality
		mark_region:	the region in which the individual was tagged and released
		region:		the region in which the individual was recaptured, or alternatively, the location of the individual at the end of the study
		hand_mort:	0 indicates the individual survived the handling process, 1 indicates mortality from handling
		doubletag:	0 indicates the individual received one mark, 1 indicates the individual was double-tagged
		tag1:		0 indicates a tag loss of the primary tag, 1 indicates the individual retained its primary tag for the duration of the study or until recapture
		tag2:		0 indicates a tag loss of the second tag, 1 indicates the individual retained its second tag for the duration of the study or until recapture
		lambda:		0 indicates that the individual either died from handling or did not retain a tag, 1 indicates the individual survived the handling 
					process and retained at least one tag, and is therefore available for tag return as part of the marked population
		survive:	0 indicates the individual died either from handling mortality or natural mortality, 1 indicates the individual is alive until recaptured or the duration of the study
		capture:	0 indicates the individual was not recaptured, either due to mortality or just not captured, 1 indicates the individual was captured by the fishery/survey
		report:		0 indicates the individual was recaptured but the tag was not reported, 1 indicates a tag return
	2.  Tag releases
	3.  Tag returns
	4.  Tag loss estimates from double tagging

NOTES:
The tagging simulation model can be ran without the Brownie estimation model component to evaluate expected tag returns for a given set of population 
and tagging parameters. A second option of the tagging simulation model is to produce a capture-recapture sampling dataset in which multiple tagging and 
recapture events can be produced.  In this case, the simulated dataset will contain capture histories by sampling event, 0 indicating the tagged individual 
was not observed, a positive integer indicating the region the individual was observed during that capture event.

To RUN the spatial Brownie tag return parameter estimation model:
	1.  Open the "Brownie_control" file and define the estimation parameters.
	2.  Source the "Brownie_control" and "Brownie_master" files from the "Simulation_run" file in R.

OUTPUTS for the Brownie tag return model will be written to the "Program_files" folder and include:
	1.  A summary matrix of tag releases for each run, titled, "releases_i.csv", where i is the iteration
	2.  A summary matrix of tag returns for each run (i), titled, "returns_i.csv"
	3.  The simulated capture history for each tagged fish for each run (i), titled, "individual_fates_i.csv"
	4.  The estimated tag loss matrix for tag returns of fish that were double-tagged and had at least one tag retained and reported, titled, "tag_loss_i.csv"
	5.  The data input file for each run (i), titled, "Brownie_i.dat" 
	6.  The Brownie model report output file for each run (i), titled, "Brownie_i.rep"
	7.  The Brownie model parameter output file for each run (i), titled, "Brownie_i.par"
	8.  The Brownie model parameter standard deviation output file for each run (i), titled, "Brownie_i.std"
	9.  The Brownie model correlation output file for each run (i), titled, "Brownie_i.cor"

Lastly, source the "Simulation_summary.r" file from the "Simulation_run" file in R to create summary plots showing the distribution of parameter estimates across
trials in terms of relative percent bias.

  

