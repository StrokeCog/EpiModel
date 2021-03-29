####Overview####

#The StrokeCog Markov probabilistic epidemiological model

#Zenodo link: https://zenodo.org/badge/latestdoi/311600790

#Model_BC_Women.r generates the base case model for women
#Model_BC_Men.r generates the base case model for men

##Set the working directory to the folder where the model inputs (parameter estimates) are stored

##Section 1 reads in relevant data, and generates a table of probabilities for initial prevalent states 
##probabilistic parameter estimates based on 1000 iterations

##Section 2 reads in population estimates for the first year of the model 

##Section 3 reads in population projections for 40 years olds for the time period of the model 

##Sections 4-14 read in relevant data, and generate a table of parameter estimates for calculation of transition probabilities
##probabilistic parameter estimates based on 1000 iterations

##Section 15 combines these parameter estimates to calculate transition probabilities for transitions between states in the model

##In Section 16, the model is calculated
##In Section 17, model results are generated
##Section 18 includes code to validate the model

#Model_BC_Results.r generates combined outputs for men and women, and the Tables and Figures in the published paper 

##This code file should be executed after Model_BC_Women.r and Model_BC_Men.r 


####Sensitivity analysis#####

#Most sensitivity analyses can be run by changing the "sheet = " part of the read.xlsx command for the relevant parameter estimate

#For some sensitivity analyses, additional code is required. These are described below. 

##For sensitivity analysis S6.3
	
	###add the following code at line 792
 
	rr_dem<-as.data.table(read.xlsx(xlsxFile="P6 Stroke Recurrence.xlsx",
                                sheet="RR_dem",
                                rows=1:2,
                                cols=1:3))


	rr_dem[,sd:=(RRhr-RRlr)/3.92]


	mean <- rr_dem[,"RR"]
	sd <- rr_dem[,"sd"]


	lnorm.pars <- lnorm_mom(mean,sd) 

	meanlog <- unlist(lnorm.pars$mu)

	sdlog <- unlist(sqrt(lnorm.pars$sigma2))

	n <- 1000
	lnorm.sample <- rlnormTrunc(n, meanlog = meanlog, sdlog = sdlog, min=0, max = Inf)


	sample.df <- data.frame(sim = c(lnorm.sample),
                        dist = rep(c("Lognormal"), each = n))

	setDT(sample.df)

	ggplot2::ggplot(sample.df, aes_string(x = "sim", fill = "dist")) +
	  geom_density(alpha = .4) + scale_fill_discrete("Distribution") +
	  xlab("Simulated values") + ylab("Density")


	#add iter variable

	sample.df[, "dist" := NULL]
	sample.df[, "iter" := 1:1000]


	#correct names
	setnames(sample.df,"sim","rec_dem_RR")

	#merge sample-df with tps_s

	tps_S <- merge(tps_S, sample.df, by = ("iter"), allow.cartesian = T )

	#calculate recurrent stroke risk specifically for dementia
	tps_S[,Rec_dem := Rec*rec_dem_RR]

	###Modify the following transition probability calculations:  


	tps_S[,"P7.7":= ((1 - Rec_dem - (Dem_mr))*Dem_stable)]
	tps_S[,"P7.16":= Rec_dem]


	tps_S[,"P21.21":= (1 - Rec_dem - (Dem_mr))*Dem_stable]


	tps_S[,"P21.29":= Rec_dem]
      

	tps_S[,"P34.34":= (1 - Rec_dem - (Dem_mr))*Dem_stable]

	tps_S[,"P34.42":= Rec_dem]

##For sensitivity analysis S9.1.2

	###replace lines 1187-1190 as follows:

	#2 times overall case fatality, based on PCSS. 

	tps_S[,"Rec_mr":=(Stroke1m*2)]

##For sensitivity analysis S9.2.2####

	###Add following code at line 688

	##dementia RR data from Oksala 2009
	##HR = 1.53, 95% CI 1.15-2.04 

	mean <- 1.53
	sd <- ((2.04-1.15)/3.92)

	lnorm.pars <- lnorm_mom(mean,sd) 

	meanlog <- unlist(lnorm.pars$mu)

	sdlog <- unlist(sqrt(lnorm.pars$sigma2))

	n <- 1000
	lnorm.sample <- rlnormTrunc(n, meanlog = meanlog, sdlog = sdlog, min=0, max = Inf)

	sample.df <- data.frame(sim = c(lnorm.sample),
                        dist = rep(c("Lognormal"), each = n))

	setDT(sample.df)

	ggplot2::ggplot(sample.df, aes_string(x = "sim", fill = "dist")) +
	  geom_density(alpha = .4) + scale_fill_discrete("Distribution") +
	  xlab("Simulated values") + ylab("Density")

	sample.df[, "dist" := NULL]
	sample.df[, "iter" := 1:1000]

	#merge value for dementia RR back into tps_S
	tps_S <- merge(tps_S, sample.df, by = ("iter"), allow.cartesian = T )

	#re-calculate Dem_RR

	tps_S[, Dem_RR := Dem_RR*sim]

	tps_S[,sim:=NULL]


####Life Expectancy (LE)####


#Model_Women_LE50.r generates the model to calculate LE for women

#The purpose of this code is to allow for generation of a life table to calculate LE at age 50
#prevalent strokes after age 50 are excluded so that the cohort who had a stroke by 50 can be followed
#incident strokes after age 50 are excluded
#the model is closed - no new 40 year old cohorts are included
#the model runs to 2064 so that the cohort age 50 in 2015 reach age 99

#The code can be adapted for men by changing the source data for 
##starting prevalences, inital pop, projections, stroke incidence, background mortality, case fatality 
##change W to M to read in the source data for men

#The code can be adapted to calculate life expectancy at ages other than 50

#When adapting the code, change the names of the output files

#Section 1 reads in relevant data, and generates a table of probabilities for initial prevalent states 
#probabilistic parameter estimates based on 1000 iterations

#Section 2 reads in population estimates for the first year of the model 

#Section 3 reads in population projections for 40 years olds for the time period of the model - this is still here but not used 

#Sections 4-14 read in relevant data, and generate a table of parameter estimates for calculation of transition probabilities
#probabilistic parameter estimates based on 1000 iterations

#Section 15 combines these parameter estimates to calculate transition probabilities for transitions between states in the model

#In Section 16, the model is calculated
#In Section 17, model results are generated
#Section 18 includes code to validate the model, based on internal checks and comparison with independent data sources

#In Section 19, the life table is generated and life expectancies are calculated
