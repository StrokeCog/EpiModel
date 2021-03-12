####Overview####

#The StrokeCog Markov probabilistic epidemiological model

#This code generates the model to calculate LE for women

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

#Sections 4-13 read in relevant data, and generate a table of parameter estimates for calculation of transition probabilities
#probabilistic parameter estimates based on 1000 iterations

#Section 15 combines these parameter estimates to calculate transition probabilities for transitions between states in the model

#In Section 16, the model is calculated
#In Section 17, model results are generated
#Section 18 includes code to validate the model, based on internal checks and comparison with independent data sources

#In Section 19, the life table is generated and life expectancies are calculated



##install and load relevant R packages

install.packages("data.table")
install.packages("openxlsx")
install.packages("ggplot2")
install.packages("EnvStats")
install.packages("hesim")
install.packages("zoo")


library(data.table)
library(openxlsx)
library(ggplot2)
library(EnvStats)
library(hesim)
library(zoo)

#set seed to ensure consistency in generation of stochastic estimates
set.seed(100)

####Additional model to run to 2055 ####

FIRST.YEAR<-2014L
LAST.YEAR<-2035L

FIRST.YEAR2<-2035L
LAST.YEAR2<-2055L

FIRST.YEAR3<-2055L
LAST.YEAR3<-2064L

##CHECKS##

##Check correct wd
##Check sex for: 
  ##starting prevalences, inital pop, projections, stroke incidence, background mortality, case fatality
  ##Also output files
  ##search for "W" or "M" and replace as appropriate 
##Check treatment effect
##Check result filenames

#NB: set working directory- local folder where model inputs and outputs are stored
setwd("FILENAME")

####1. Starting Prevalences- Community and NH ####


#stochastic values for community stroke prevalence

#read in community stroke prevalance data


 prev_c<-as.data.table(read.xlsx(xlsxFile="P2_1 Stroke Community Prevalence.xlsx",
                              sheet="W3 Stroke Prevalence C W Q",
                              rows=1:7,
                              cols=2:4))


##index for merging
prev_c[, index := 1]

##create table of iterations
iters <- data.table(index = 1, iter = 1:1000)


##merge iterations and original table by all combinations
prev_c <- merge(prev_c, iters, by = "index" , allow.cartesian = T )

#Create alpha and beta (non-events and events) 
prev_c[, c("a","b"):= list ((StrokePr_C*N), ((1-StrokePr_C)*N))]

##create stochastic values for StrokePr_C based on the beta distribution

prev_c[, live := rbeta(length(a),a,b,ncp=0)]

##change name to equivalent transition probability
setnames(prev_c,"live","StrokePr_C_st")

##delete alpha and beta columns
prev_c[,c("a","b") := NULL]


#read in proportion of total population in NH
prev_nh<-as.data.table(read.xlsx(xlsxFile="P2_2 NH residents.xlsx",
                                sheet="Pr in NH W",
                                rows=1:11,
                                cols=2:3))

#create table with age 40-100

age <- data.table(Age=40:100, index=1)

#create age group variable to match with age groups in prev_nh

age[Age<100,"Age.Group":= 10]
age[Age<95,"Age.Group":= 9]
age[Age<90,"Age.Group":= 8]
age[Age<85,"Age.Group":=7]
age[Age<80,"Age.Group":=6]
age[Age<75,"Age.Group":=5]
age[Age<70,"Age.Group":=4]
age[Age<65,"Age.Group":=3]
age[Age<55,"Age.Group":=2]
age[Age<45,"Age.Group":=1]

#merge age and nh table by all combinations
prev_nh <- merge(prev_nh, age, by = "Age.Group" , allow.cartesian = T )

#merge iterations and original table by all combinations
prev_nh <- merge(prev_nh, iters, by = "index" , allow.cartesian = T )

#delete unnecessary columns
prev_nh[,c("Age.Group","index"):= NULL]

##read in data on pr in nh with stroke by cognitive outcome
prev_nh_str<-as.data.table(read.xlsx(xlsxFile="P2_2 NH residents.xlsx",
                                 sheet="Stroke_Pr_NH",
                                 rows=1:3,
                                 cols=2:8))


#create index for merging
prev_nh_str[, index := 1]


#merge iterations and prev_nh_str by all combinations
prev_nh_str <- merge(prev_nh_str , iters, by = "index" , allow.cartesian = T )


##Beta Distribution for prevalence of stroke in NH admissions

#Create alpha and beta for stroke pr (non-events and events) 
prev_nh_str[, c("a","b"):= list ((Stroke_Pr_NH1*N_StrokePr), ((1-Stroke_Pr_NH1)*N_StrokePr))]

#generate stochastic values
prev_nh_str[, live := rbeta(length(a),a,b,ncp=0)]

#change name to equivalent transition probability
setnames(prev_nh_str,"live","Stroke_Pr_NH")

#delete alpha and beta columns
prev_nh_str[,c("a","b","index","Stroke_Pr_NH1","N_StrokePr") := NULL]


##Dirichlet distribution for cognitive profile
##Ref: https://devinincerti.com/2018/02/10/psa.html

#read in data on cognitive profile of nh stroke residents 
prev_nh_ci<-as.data.table(read.xlsx(xlsxFile="P3_2 Stroke CI Profile in NH.xlsx",
                                     sheet="Stroke_CI_NH ",
                                     rows=1:3,
                                     cols=2:6))


#sets the prior = 1
alpha <- 1


#pulls estimates from the Dirichlet distribution
prev_nh_ciD <- prev_nh_ci[,hesim::rdirichlet_mat(n = 1000, 
                                             alpha + c(NCI_Pr_NH, CIND_Pr_NH,Dementia_Pr_NH)),
                      by="Age.Group"]

prev_nh_ciD <- as.data.table(prev_nh_ciD)


#puts the stochastic estimates into a data frame with iterations and states                    
prev_nh_ci_s <- data.frame(p = c(prev_nh_ciD), 
                         iter = rep(1:1000, each=3),
                         state = rep(rep(paste0("state = ", seq(1, 3))) 
                         ))

#sets the df as a dt
setDT(prev_nh_ci_s)

#reshapes the data table
prev_nh_ci_s <- dcast(prev_nh_ci_s, p.Age.Group + iter ~ state, value.var = "p.V1")

#sets the column names back to the original
setnames(prev_nh_ci_s,c("p.Age.Group","iter","state = 1","state = 2","state = 3"), 
         c("Age.Group","iter","NCID_Pr_NH","CINDD_Pr_NH","Dementia_Pr_NH"))

##Merge estimates of stroke and health state prevalence in NH, with estimates of proportion of total population in NH##



#create age group variable to match with age groups in prev_nh_str and prev_nh_ci_s 
prev_nh[, "Age.Group":= 2]
prev_nh[Age<75,"Age.Group":=1]

#merge prev_nh and prev_nh_str 
prev_nh <- merge(prev_nh, prev_nh_str, by = c("Age.Group", "iter") , allow.cartesian = T  )

#merge prev_nh and prev_nh_ci_s 
prev_nh <- merge(prev_nh, prev_nh_ci_s, by = c("Age.Group", "iter") , allow.cartesian = T  )

#Delete unnecessary columns 

prev_nh[,"Stroke_Pr_NH":= (NH_Pr*Stroke_Pr_NH)]
prev_nh[,"NH_Pr":=NULL]
prev_nh[,"Age.Group":=NULL]


#Change age group variable in prev_nh to match up with prev_c for merging purposes
#Age groups match up with prev_c

prev_nh[, "Age.Group":= 4]
prev_nh[Age<75,"Age.Group":=3]
prev_nh[Age<65,"Age.Group":=2]
prev_nh[Age<45,"Age.Group":=1]


#merge prev_c and prev_nh by age group
prev <- merge(prev_nh, prev_c, by = c("Age.Group","iter"), allow.cartesian = T )


#delete unneccessary columns
prev[,c("index","N") := NULL]

#sort by iteration
setkey(prev,iter)


##read in in the cognitive profile estimates for community prevalent stroke
prev_cog<-as.data.table(read.xlsx(xlsxFile="P3_1 Stroke Comm Prev_Cog.xlsx",
                                 sheet="Cog Pr C",
                                 rows=1:4,
                                 cols=1:7))



#sets the prior = 1
alpha <- 1


#pulls estimates from the Dirichlet distribution
prev_cogD <- prev_cog[,hesim::rdirichlet_mat(n = 1000, 
                    alpha + c(NCInD_Pr_C, NCID_Pr_C,CINDnD_Pr_C,CINDD_Pr_C,Dementia_Pr_C)),
                    by="Age.Group"]


prev_cogD <- as.data.table(prev_cogD)


#puts the stochastic estimates into a data frame with iterations and states                    
prev_cog_s <- data.frame(p = c(prev_cogD), 
                     iter = rep(1:1000, each=5),
                     state = rep(rep(paste0("state = ", seq(1, 5)), 5) 
                     ))

#sets the df as a dt
setDT(prev_cog_s)


#reshapes the data table
prev_cog_s <- dcast(prev_cog_s, p.Age.Group + iter ~ state, value.var = "p.V1")


#sets the column names back to the original
setnames(prev_cog_s,c("p.Age.Group","iter","state = 1","state = 2","state = 3",
                "state = 4", "state = 5"), 
          c("Age.Group2","iter","NCInD_Pr_C","NCID_Pr_C","CINDnD_Pr_C","CINDD_Pr_C",
          "Dementia_Pr_C"))


#creates Age Group2 variable in prev table to allow merging of prev_cog_s
prev[, "Age.Group2":= 2]
prev[Age<75,"Age.Group2":=1]



#merge cognitive profile estimates back into prevalence table
prev <- merge(prev, prev_cog_s, by = c("Age.Group2","iter"), allow.cartesian = T )

#sort by iteration
setkey(prev,iter)



##Combine community and nh prevalence estimates to calculate overall starting prevalences

prev[,"NCInD" := (StrokePr_C_st*NCInD_Pr_C)]
prev[,"NCID" := (StrokePr_C_st*NCID_Pr_C)+(Stroke_Pr_NH*NCID_Pr_NH)]
prev[,"CINDnD" := (StrokePr_C_st*CINDnD_Pr_C)]
prev[,"CINDD" := (StrokePr_C_st*CINDD_Pr_C)+(Stroke_Pr_NH*CINDD_Pr_NH)]
prev[,"Dementia" := (StrokePr_C_st*Dementia_Pr_C)+(Stroke_Pr_NH*Dementia_Pr_NH)]


#delete unneccessary columns
prev[,c("Age.Group","Age.Group2","StrokePr_C","StrokePr_C_st",
        "NCInD_Pr_C", "NCID_Pr_C", "CINDnD_Pr_C", "CINDD_Pr_C", "Dementia_Pr_C", 
        "NCID_Pr_NH", "Stroke_Pr_NH", "CINDD_Pr_NH", "Dementia_Pr_NH") := NULL]

#change column names to reflect health states
setnames(prev,c("iter", "Age","NCInD","NCID","CINDnD","CINDD","Dementia"), 
         c("iter", "Age","s.3","s.4","s.5","s.6","s.7"))


####TAKE OUT PREVALENT STROKES AFTER AGE 50####

prev[Age>50, s.3:=0]
prev[Age>50, s.4:=0]
prev[Age>50, s.5:=0]
prev[Age>50, s.6:=0]
prev[Age>50, s.7:=0]

#outputting median prevalences to Excel for checking purposes
DMcheckW<-prev[,list(s.3=median(s.3),
                     s.4=median(s.4),
                     s.5=median(s.5),
                     s.6=median(s.6),
                     s.7=median(s.7)
                     
),
by=list(Age)]

write.xlsx(DMcheckW,"DMcheckW.xlsx", rows=1:60,cols=1:6)




####2. Initial population ####


#read in sex-specific population estimates for the first year of the model 

inipop<-as.data.table(read.xlsx(xlsxFile="P1_1 Initial Population.xlsx",
                                 sheet="Inipop W",
                                 rows=1:61,
                                 cols=1:5))

setnames(inipop,c("DiseaseFree","IncStroke"), c("s.1", "s.2"))


#create colums for remaining health states
inipop[,s.8:=0]
inipop[,s.9:=0]
inipop[,s.10:=0]
inipop[,s.11:=0]
inipop[,s.12:=0]
inipop[,s.13:=0]
inipop[,s.14:=0]
inipop[,s.15:=0]
inipop[,s.16:=0]
inipop[,s.17:=0]
inipop[,s.18:=0]
inipop[,s.19:=0]
inipop[,s.20:=0]
inipop[,s.21:=0]
inipop[,s.22:=0]
inipop[,s.23:=0]
inipop[,s.24:=0]
inipop[,s.25:=0]
inipop[,s.26:=0]
inipop[,s.27:=0]
inipop[,s.28:=0]
inipop[,s.29:=0]
inipop[,s.30:=0]
inipop[,s.31:=0]
inipop[,s.32:=0]
inipop[,s.33:=0]
inipop[,s.34:=0]
inipop[,s.35:=0]
inipop[,s.36:=0]
inipop[,s.37:=0]
inipop[,s.38:=0]
inipop[,s.39:=0]
inipop[,s.40:=0]
inipop[,s.41:=0]
inipop[,s.42:=0]
inipop[,s.43:=0]
inipop[,s.44:=0]
inipop[,s.45:=0]
inipop[,s.46:=0]
inipop[,s.47:=0]
inipop[,s.48:=0]
inipop[,s.49:=0]
inipop[,s.50:=0]
inipop[,s.51:=0]



#### 3. Projections of 40 yr old population - not used####

#popproj<-as.data.table(read.xlsx(xlsxFile="P1_2 Pop Projections LE.xlsx",
 #                             sheet="Projections W",
  #                            rows=1:45,
 #                             cols=1:2))

####4. Transition probability table#### 

tps_S<-data.table(Age = 40:100)


##index for merging
tps_S[,index := 1]

##merge in 1000 iterations
tps_S <- merge(tps_S,iters, by = "index", allow.cartesian = T)

setkey(tps_S, iter)


####5. Stroke Incidence####


##first, calculate an unadjusted incidence rate that does not take account of out of hospital
##and incident recurrent stroke


#read in stroke incidence data

stroke_inc<-as.data.table(read.xlsx(xlsxFile="P4 Stroke Incidence.xlsx",
                                    sheet="Stroke Incidence W",
                                    rows=1:11,
                                    cols=2:4))


##index for merging
stroke_inc[, index := 1]

##create table of iterations
iters <- data.table(index = 1, iter = 1:1000)


##merge iterations and original table by all combinations
stroke_inc <- merge(stroke_inc, iters, by = "index" , allow.cartesian = T )

#Create alpha and beta (non-events and events)
stroke_inc[, c("a","b"):= list ((ISR*N), ((1-ISR)*N))]

#create stochastic values for incident stroke rate (ISR) based on the beta distribution

stroke_inc[, live := rbeta(length(a),a,b,ncp=0)]

##change name 
##unadjusted = ooh and rec stroke not accounted for

setnames(stroke_inc,"live","ISR_unadj")

##delete alpha and beta columns
stroke_inc[,c("a","b") := NULL]


##creates Age Group variable in tps table to allow merging of stroke_inc
tps_S[Age<50, "Age.Group":= 1]
tps_S[Age>49,"Age.Group":=2]
tps_S[Age>54,"Age.Group":=3]
tps_S[Age>59,"Age.Group":=4]
tps_S[Age>64,"Age.Group":=5]
tps_S[Age>69,"Age.Group":=6]
tps_S[Age>74,"Age.Group":=7]
tps_S[Age>79,"Age.Group":=8]
tps_S[Age>84,"Age.Group":=9]
tps_S[Age>89,"Age.Group":=10]


##merge stroke_inc and tps_S by age group and iter
tps_S <- merge(tps_S, stroke_inc, by = c("Age.Group","iter"), allow.cartesian = T )


##delete unneccessary columns
tps_S[,c("index.x","index.y","N","ISR") := NULL]

#sort by iteration
setkey(tps_S,iter)

setcolorder(tps_S,c("Age","iter","ISR_unadj"))

##delete unneccessary columns
tps_S[,"Age.Group" := NULL]

##Estimates for proportion of out of hospital strokes

ooh<-as.data.table(read.xlsx(xlsxFile="P4 Stroke Incidence.xlsx",
                                    sheet="OOH Data",
                                    rows=1:2,
                                    cols=1:2))

#index for merging
ooh[,index:=1]

#merge iterations and original table by all combinations
ooh <- merge(ooh, iters, by = "index" , allow.cartesian = T )

#Create alpha and beta (non-events and events) 
ooh[, c("a","b"):= list ((OOH), (OOH_N-OOH))]

#create stochastic values for OOH based on the beta distribution

ooh[, live := rbeta(length(a),a,b,ncp=0)]

#change name 
setnames(ooh,"live","OOH_st")

##delete alpha and beta columns
ooh[,c("a","b","OOH","OOH_N","index") := NULL]

#merge into tps_S

##merge ooh and tps_S by iter
tps_S <- merge(tps_S, ooh, by = "iter", allow.cartesian = T )

##Estimates for proportion of incident recurrent strokes

inc_rec<-as.data.table(read.xlsx(xlsxFile="P4 Stroke Incidence.xlsx",
                             sheet="Inc Rec Data",
                             rows=1:5,
                             cols=2:4))

#index for merging
inc_rec[,index:=1]

##merge iterations and original table by all combinations
inc_rec <- merge(inc_rec, iters, by = "index" , allow.cartesian = T )

#Create alpha and beta (non-events and events) 
inc_rec[, c("a","b"):= list ((Inc_Rec), (Inc_Rec_N-Inc_Rec))]

##create stochastic values for incident recurrent stroke based on the beta distribution

inc_rec[, live := rbeta(length(a),a,b,ncp=0)]

##change name 
setnames(inc_rec,"live","inc_rec_st")

##delete alpha and beta columns
inc_rec[,c("a","b","Inc_Rec","Inc_Rec_N","index") := NULL]

##creates Age Group variable in tpr table to allow merging of inc_rec
tps_S[Age<60, "Age.Group":= 1]
tps_S[Age>59, "Age.Group":= 2]
tps_S[Age>69, "Age.Group":= 3]
tps_S[Age>79, "Age.Group":= 4]

##merge inc_rec and tps_S by age group and iter
tps_S <- merge(tps_S, inc_rec, by = c("Age.Group","iter"), allow.cartesian = T )

#Final step is to create adjusted incident stroke rate 
#by applying estimates of OOH and incident recurrent stroke

tps_S[,P1.2:=ISR_unadj*(1+OOH_st)]
tps_S[,P1.2:=P1.2*(1-inc_rec_st)]

####TAKE OUT INCIDENT STROKES FOR AGE>50####

tps_S[Age>50,P1.2:=0]

#Delete unnecessary columns

tps_S[,c("ISR_unadj","OOH_st","inc_rec_st"):=NULL]



####6. Background Mortality####


##read in sex-specific mortality rates for total population, by age

mr<-as.data.table(read.xlsx(xlsxFile="P9_2 Background mortality.xlsx",
                            sheet="mr W",
                            rows=1:7,
                            cols=2:3))

##Add in iterations to mr
#index for merging
mr[, index := 1]

#merge iterations and original table by all combinations
mr<- merge(mr, iters, by = "index" , allow.cartesian = T )


##merge mr with tps_S

#add age group columns to tps_S

#creates Age Group variable in tpr table to allow merging of mr
tps_S[Age<45, "Age.Group":=1]
tps_S[Age>44, "Age.Group":=2]
tps_S[Age>54, "Age.Group":=3]
tps_S[Age>64, "Age.Group":=4]
tps_S[Age>74, "Age.Group":=5]
tps_S[Age>84, "Age.Group":=6]




#merge mr and tps_S by age group and iter
tps_S <- merge(tps_S, mr, by = c("Age.Group","iter"), allow.cartesian = T )


#delete unneccessary columns
tps_S[,c("index","Age.Group") := NULL]

#set column order
setcolorder(tps_S,c("Age","iter","P1.2","mr"))


####7. Relative risk for Background Mortality####

#read in data

rr<-as.data.table(read.xlsx(xlsxFile="P9_2 Background Mortality.xlsx",
                            sheet="RR_PSA_5yr",
                            rows=1:2,
                            cols=3:6))


##Lognormal function for method of moments- for stochastic values for RR
##Ref https://devinincerti.com/2018/02/10/psa.html

lnorm_mom <- function(mean, sd){
  if (mean > 0){
    sigma2 <- log((sd^2 + mean^2)/mean^2)
    mu <- log(mean) - 1/2 * sigma2
  } else{
    stop("Mean must be positive")
  }
  return(list(mu = mu, sigma2 = sigma2))
}

##Calculate sd from confidence intervals

rr[,sd:=(uci-lci)/3.92]


mean <- rr[,"RR"]
sd <- rr[,"sd"]

#Apply method of moments to mean and sd of risk ratio estimate
#to estimate parameters of the lognormal distribution
lnorm.pars <- lnorm_mom(mean,sd) 

meanlog <- unlist(lnorm.pars$mu)

sdlog <- unlist(sqrt(lnorm.pars$sigma2))

#generate stochastic estimates by sampling from the lognormal distribution
#iterations = 1000

n <- 1000

lnorm.sample <- rlnormTrunc(n, meanlog = meanlog, sdlog = sdlog, min=0, max = Inf)


#puts the stochastic estimates in a data table
sample.df <- data.table(sim = c(lnorm.sample),
                        dist = rep(c("Lognormal"), each = n))


#plots estimates to check distribution
ggplot2::ggplot(sample.df, aes_string(x = "sim", fill = "dist")) +
  geom_density(alpha = .4) + scale_fill_discrete("Distribution") +
  xlab("Simulated values") + ylab("Density")


##add iter variable

sample.df[, "dist" := NULL]
sample.df[, "iter" := 1:1000]

#change var name to strRR (relative mortality risk for prevalent stroke)
setnames(sample.df,"sim","strRR")


#copy var to generate relative risks specific to health states
#In base case these are equivalent, but this is varied in sensitivity analysis 
sample.df[, "NCInD_RR" := strRR]
sample.df[, "NCID_RR" := strRR]
sample.df[, "CINDnD_RR" := strRR]
sample.df[, "CINDD_RR" := strRR]
sample.df[, "Dem_RR" := strRR]

#merge sample-df with transition probability table

tps_S <- merge(tps_S, sample.df, by = ("iter"), allow.cartesian = T )


####8. Stroke Recurrence####


#read in recurrence data 

rec12m<-as.data.table(read.xlsx(xlsxFile="P6 Stroke Recurrence.xlsx",
                                    sheet="Rec12m_all ages",
                                    rows=1:2,
                                    cols=1:2))

#index for merging
rec12m[, index := 1]


#merge iterations and original table by all combinations
rec12m <- merge(rec12m, iters, by = "index" , allow.cartesian = T )

#Create alpha and beta (non-events and events) 
rec12m[, c("a","b"):= list ((Rec12m*N), ((1-Rec12m)*N))]

#create stochastic values for rec12m based on the beta distribution

rec12m[, live := rbeta(length(a),a,b,ncp=0)]

#change var name
setnames(rec12m,"live","Rec_12m")

#delete alpha and beta columns
rec12m[,c("a","b") := NULL]



#add age group columns to tps_S 
#use this code for age-specific estimates in sensitivity analysis

#creates Age Group variable in tpr table to allow merging of rec12m
#tps_S[Age<65, "Age.Group":= 1]
#tps_S[Age>64,"Age.Group":=2]
#tps_S[Age>74,"Age.Group":=3]




#merge stroke_inc and tps_S by iter
tps_S <- merge(tps_S, rec12m, by = c("iter"), allow.cartesian = T )


#delete unneccessary columns
tps_S[,c("index","N","Rec12m") := NULL]

#sort by iteration
setkey(tps_S,iter)


#copy Rec_12m to new variable P2.11 which is the relevant tpr.  keep Rec_12m for use in calc. 
tps_S[,"P2.11" := Rec_12m ]


##Stroke recurrence - annual rate


#read in recurrence data 

rec<-as.data.table(read.xlsx(xlsxFile="P6 Stroke Recurrence.xlsx",
                                sheet="Rec_all ages",
                                rows=1:4,
                                cols=1:2))

#index for merging
rec[, index := 1]


#merge iterations and rec by all combinations
rec <- merge(rec, iters, by = "index" , allow.cartesian = T )

#Create alpha and beta (non-events and events)
rec[, c("a","b"):= list ((rec*N), ((1-rec)*N))]

#create stochastic values for rec based on the beta distribution

rec[, live := rbeta(length(a),a,b,ncp=0)]

#change name to equivalent transition probability
setnames(rec,"live","Rec")

#delete alpha and beta columns
rec[,c("a","b") := NULL]


#merge rec and tps_S by iter
tps_S <- merge(tps_S, rec, by = c("iter"), allow.cartesian = T )

#delete unneccessary columns
tps_S[,c("index","N","rec") := NULL]



#sort by iteration
setkey(tps_S,iter)

#copy Rec to new variable P3.12 which is the relevant tpr.  keep Rec for use in calc. 
tps_S[,"P3.12" := Rec ]


####9. Treatment effect####

#Hypothetical Treatment effect for cognitive rehabilitation 
#Systematic review evidence from Merriman 2019 (10.1136/bmjopen-2018-024429)and Rogers 2018 (10.1007/s11065-018-9378-4) 
#indicates a benefit for post-stroke rehab of hedges g = 0.48 (Rogers) and g = 0.38 (Merriman). 
#These are interpreted as small to medium effect sizes 
#Olivier 2017 (10.1080/03610926.2015.1134575) suggest risk ratios to correspond with established Cohen's d effect sizes
#a "small" risk ratio = 1.22, and a "medium" risk ratio = 1.86. 
#inverting these equates to a risk reduction for CIND post-stroke of 18% (small effect) and 46% (medium effect).

#bc = base case
#se = small effect
#me = medium effect

bc = 1
se1 =  0.82
me1 = 0.54

te = bc


####10. Cognitive Outcomes - 12m####

##read in cognitive profile estimates for 12m outcome

inicog<-as.data.table(read.xlsx(xlsxFile="P5 12m outcome.xlsx",
                                  sheet="Cog 12m",
                                  rows=1:3,
                                  cols=1:7))


#sets the prior = 1
alpha <- 1


#pulls estimates from the Dirichlet distribution

inicog <- inicog[,hesim::rdirichlet_mat(n = 1000, 
              alpha + c(NCInD_initial, NCID_initial,CINDnD_initial,CINDD_initial,Dementia_initial)),
                      by="Age.Group"]

inicog <- as.data.table(inicog)


#puts the stochastic estimates into a data table                   
inicog_s <- data.frame(p = c(inicog), 
                         iter = rep(1:1000, each=5),
                         state = rep(rep(paste0("state = ", seq(1, 5)), 5) 
                         ))
#sets as DT
setDT(inicog_s)

#reshapes the data table
inicog_s <- dcast(inicog_s, p.Age.Group + iter ~ state, value.var = "p.V1")


#sets the column names back to the original
setnames(inicog_s,c("p.Age.Group","iter","state = 1","state = 2","state = 3",
                      "state = 4", "state = 5"), 
         c("Age.Group","iter","NCInD_initial","NCID_initial","CINDnD_initial","CINDD_initial",
           "Dementia_initial"))

##Apply treatment effect

#create variables for total no disability and total disability (non-dementia)

inicog_s[,tot_ND_NDem:=NCInD_initial+CINDnD_initial]
inicog_s[,tot_D_NDem:=NCID_initial+CINDD_initial]

#applies treatment effect to proportion of CIND post-stroke

inicog_s[,CINDnD_initial:=CINDnD_initial*te]
inicog_s[,CINDD_initial:=CINDD_initial*te]

##create new proportions of NCI based on reduced proportions of CIND

inicog_s[,NCInD_initial:=tot_ND_NDem-CINDnD_initial]
inicog_s[,NCID_initial:=tot_D_NDem-CINDD_initial]

##check proportions add to 1

inicog_s[,check:=NCInD_initial+CINDnD_initial+NCID_initial+CINDD_initial+Dementia_initial]

#Delete total and check columns

inicog_s[,tot_ND_NDem:= NULL]
inicog_s[,tot_D_NDem:= NULL]
inicog_s[,check:= NULL]

#index for merging
tps_S[, index := 1]

#creates Age Group variable in tpr table to allow merging of ini_cog
tps_S[Age<70, "Age.Group":= 1]
tps_S[Age>69,"Age.Group":=2]

#merge inicog_s and tps_S by age group and iter
tps_S <- merge(tps_S, inicog_s, by = c("Age.Group","iter"), allow.cartesian = T )


##delete unneccessary columns
tps_S[,c("index","Age.Group") := NULL]

##sort by iteration
setkey(tps_S,iter)



####11. Health State transitions####

##read in the cognitive transition estimates

cog_tr<-as.data.table(read.xlsx(xlsxFile="P7 Cognitive Transitions.xlsx",
                                sheet="Cog Trans",
                                rows=1:11,
                                cols=2:8))

#sets the prior = 1
alpha <- 1


#pulls estimates from the Dirichlet distribution for Health State 1
cog_tr_s1 <- cog_tr[State==1,hesim::rdirichlet_mat(n = 1000, 
                                        alpha + c(NCInD, NCID, CINDnD, CINDD)), by="Age.Group"
                ]

#puts the stochastic estimates into a data table with iterations and states   

cog_tr_s1<- data.frame(p = c(cog_tr_s1), 
                       iter = rep(1:1000, each=4),
                       state = rep(rep(paste0("state = ", seq(1, 4)), 4) 
                       ))
#sets as DT
setDT(cog_tr_s1)

#reshapes the data table

cog_tr_s1 <- dcast(cog_tr_s1, p.Age.Group + iter ~ state, value.var = "p.V1")

#sets var names
setnames(cog_tr_s1,c("p.Age.Group","iter","state = 1","state = 2","state = 3",
                     "state = 4"), 
         c("Age.Group","iter","NCInD_stable","NCID","NCInD_CINDnD","NCInD_CINDD"))




#pulls estimates from the Dirichlet distribution for Health State 2
cog_tr_s2 <- cog_tr[State==2,hesim::rdirichlet_mat(n = 1000, 
                                                   alpha + c(NCID, CINDD)), by="Age.Group"]


#puts the stochastic estimates into a data table with iterations and states   

cog_tr_s2<- data.frame(p = c(cog_tr_s2), 
                       iter = rep(1:1000, each=2),
                       state = rep(rep(paste0("state = ", seq(1, 2)), 2) 
                       ))
setDT(cog_tr_s2)

#reshapes the data table
cog_tr_s2 <- dcast(cog_tr_s2, p.Age.Group + iter ~ state, value.var = "p.V1")

#sets var names
setnames(cog_tr_s2,c("p.Age.Group","iter","state = 1","state = 2"), 
         c("Age.Group","iter","NCID_stable","NCID_CINDD"))


#pulls estimates from the Dirichlet distribution for Health State 3
cog_tr_s3 <- cog_tr[State==3,hesim::rdirichlet_mat(n = 1000, 
                                                   alpha + c(CINDnD, CINDD, Dem)), by="Age.Group"]

#puts the stochastic estimates into a data table with iterations and states   

cog_tr_s3<- data.frame(p = c(cog_tr_s3), 
                       iter = rep(1:1000, each=3)
                       , state = rep(rep(paste0("state = ", seq(1, 3)))))

setDT(cog_tr_s3)

#reshapes the data table
cog_tr_s3 <- dcast(cog_tr_s3, p.Age.Group + iter ~ state, value.var = "p.V1")


#set var names
setnames(cog_tr_s3,c("p.Age.Group","iter","state = 1","state = 2","state = 3"), 
         c("Age.Group","iter","CINDnD_stable","CINDnD_CINDD","CINDnD_Dem"))


#pulls estimates from the Dirichlet distribution for Health State 4
cog_tr_s4 <- cog_tr[State==4,hesim::rdirichlet_mat(n = 1000, 
                                                   alpha + c(CINDD, Dem))
                    , by="Age.Group"]


#puts the stochastic estimates into a data table with iterations and states   

cog_tr_s4<- data.frame(p = c(cog_tr_s4), 
                       iter = rep(1:1000, each=2),
                       state = rep(rep(paste0("state = ", seq(1, 2))) 
                       ))

setDT(cog_tr_s4)

#reshapes the data table
cog_tr_s4 <- dcast(cog_tr_s4, p.Age.Group + iter ~ state, value.var = "p.V1")

#sets var names
setnames(cog_tr_s4,c("p.Age.Group","iter","state = 1","state = 2"), 
         c("Age.Group","iter","CINDD_stable","CINDD_Dem"))


##merge tps_S and each cog_tr datatable

#index for merging
tps_S[, index := 1]

#add age group columns to tps_S

#creates Age Group variable in tpr table to allow merging of cog_tr_s1
tps_S[Age<75, "Age.Group":= 1]
tps_S[Age>74,"Age.Group":=2]

#merge tps_S and each cog_tr dataframe
tps_S <- merge(tps_S, cog_tr_s1, by = c("Age.Group", "iter"), allow.cartesian = T )
tps_S <- merge(tps_S, cog_tr_s2, by = c("Age.Group", "iter"), allow.cartesian = T )
tps_S <- merge(tps_S, cog_tr_s3, by = c("Age.Group", "iter"), allow.cartesian = T )
tps_S <- merge(tps_S, cog_tr_s4, by = c("Age.Group", "iter"), allow.cartesian = T )


#delete unneccessary columns
tps_S[,c("index","Age.Group") := NULL]

#sort by iteration
setkey(tps_S,iter)


####12. Post-recurrent outcomes####

#read in post-recurrent outcome data

pr_cog_tr<-as.data.table(read.xlsx(xlsxFile="P8 Rec Stroke Cog Transitions.xlsx",
                                sheet="Rec Cog Trans",
                                rows=1:5,
                                cols=1:6))

#sets the prior = 1
alpha <- 1

#pulls estimates from the Dirichlet distribution for Health State 1
pr_cog_tr_s1 <- pr_cog_tr[State==1,hesim::rdirichlet_mat(n = 1000, 
                                                   alpha + c(NCInD, NCID, CINDnD, CINDD, Dem))
                          ]
#puts stochastic estimates into a data table
pr_cog_tr_s1 <- as.data.table(pr_cog_tr_s1)
                          
#reshapes the data table
pr_cog_tr_s1 <- dcast(pr_cog_tr_s1, V3 ~ V1 + V2, value.var = "value")
                          
#set var names                          
setnames(pr_cog_tr_s1,c("V3","1_1","1_2","1_3","1_4","1_5"), 
          c("iter","pr_NCInD_stable","pr_NCID","pr_NCInD_CINDnD","pr_NCInD_CINDD","pr_NCInD_Dem"))                          

                          
#pulls estimates from the Dirichlet distribution for Health State 2
pr_cog_tr_s2 <- pr_cog_tr[State==2,hesim::rdirichlet_mat(n = 1000, 
                                                         alpha + c(NCID, CINDD, Dem))
                          ]

#puts stochastic estimates into a data table
pr_cog_tr_s2 <- as.data.table(pr_cog_tr_s2)

#reshapes the data table
pr_cog_tr_s2 <- dcast(pr_cog_tr_s2, V3 ~ V1 + V2, value.var = "value")


setnames(pr_cog_tr_s2,c("V3","1_1","1_2","1_3"), 
         c("iter","pr_NCID_stable","pr_NCID_CINDD","pr_NCID_Dem"))                          

#pulls estimates from the Dirichlet distribution for Health State 3
pr_cog_tr_s3 <- pr_cog_tr[State==3,hesim::rdirichlet_mat(n = 1000, 
                                                         alpha + c(CINDnD, CINDD, Dem))
                          ]

#puts stochastic estimates into a data table
pr_cog_tr_s3 <- as.data.table(pr_cog_tr_s3)

#reshapes the data table
pr_cog_tr_s3 <- dcast(pr_cog_tr_s3, V3 ~ V1 + V2, value.var = "value")

#sets var names
setnames(pr_cog_tr_s3,c("V3","1_1","1_2","1_3"), 
         c("iter","pr_CINDnD_stable","pr_CINDnD_CINDD","pr_CINDnD_Dem"))     


#pulls estimates from the Dirichlet distribution for Health State 4
pr_cog_tr_s4 <- pr_cog_tr[State==4,hesim::rdirichlet_mat(n = 1000, 
                                                         alpha + c(CINDD, Dem))
                          ]

#puts stochastic estimates into a data table
pr_cog_tr_s4 <- as.data.table(pr_cog_tr_s4)

#reshapes the data table
pr_cog_tr_s4 <- dcast(pr_cog_tr_s4, V3 ~ V1 + V2, value.var = "value")

#sets var names
setnames(pr_cog_tr_s4,c("V3","1_1","1_2"), 
         c("iter","pr_CINDD_stable","pr_CINDD_Dem"))


#merge tps_S and each cog_tr dataframe
tps_S <- merge(tps_S, pr_cog_tr_s1, by = ("iter"), allow.cartesian = T )
tps_S <- merge(tps_S, pr_cog_tr_s2, by = ("iter"), allow.cartesian = T )
tps_S <- merge(tps_S, pr_cog_tr_s3, by = ("iter"), allow.cartesian = T )
tps_S <- merge(tps_S, pr_cog_tr_s4, by = ("iter"), allow.cartesian = T )


##Calculating outcomes for recurrent stroke within 12 months
#applies post-recurrent probabilities to the base line initial post-stroke probabilities

tps_S[,"NCInD_pr":= (NCInD_initial*pr_NCInD_stable)]
tps_S[,"NCID_pr":= ((NCInD_initial*pr_NCID)+(NCID_initial*pr_NCID_stable))]
tps_S[,"CINDnD_pr":=((NCInD_initial*pr_NCInD_CINDnD)+(CINDnD_initial*pr_CINDnD_stable))]
tps_S[,"CINDD_pr":=((NCInD_initial*pr_NCInD_CINDD)+(NCID_initial*pr_NCID_CINDD)+
                      (CINDnD_initial*pr_CINDnD_CINDD)+(CINDD_initial*pr_CINDD_stable))]
tps_S[,"Dem_pr":=((NCInD_initial*pr_NCInD_Dem)+(NCID_initial*pr_NCID_Dem)+(CINDnD_initial*pr_CINDnD_Dem)
                       +(CINDD_initial*pr_CINDD_Dem)+(Dementia_initial))]

####13. Case fatality - 1M####

#read in data
stroke1m<-as.data.table(read.xlsx(xlsxFile="P9_1 Case Fatality.xlsx",
                                  sheet="Case Fatality W",
                                  rows=1:6,
                                  cols=2:4))

#index for merging
stroke1m[, index := 1]

##merge iterations and stroke1m by all combinations
stroke1m <- merge(stroke1m, iters, by = "index" , allow.cartesian = T )

#Create alpha and beta (non-events and events)
stroke1m[, c("a","b"):= list ((stroke1m*n), ((1-stroke1m)*n))]

##create stochastic values for stroke1m based on the beta distribution
stroke1m[, live := rbeta(length(a),a,b,ncp=0)]

##change name to equivalent transition probability
setnames(stroke1m,"live","Stroke1m")

##delete alpha and beta columns
stroke1m[,c("a","b") := NULL]


##creates Age Group variable in tpr table to allow merging of stroke1m
tps_S[Age<55, "Age.Group":= 1]
tps_S[Age>54,"Age.Group":=2]
tps_S[Age>64,"Age.Group":=3]
tps_S[Age>74,"Age.Group":=4]
tps_S[Age>84,"Age.Group":=5]

##merge stroke1m and tps_S by age group and iter
tps_S <- merge(tps_S, stroke1m, by = c("Age.Group","iter"), allow.cartesian = T )

#sort by iteration
setkey(tps_S,iter)

#delete unneccessary columns
tps_S[,c("Age.Group","index","n","stroke1m") := NULL]

##reduce Stroke1m by 34% based on 10yr reduction in-hospital mortality (see NAHM Annual Report 2018)
##NAHM reports a 38% reduction in IS mortality 2009-2018, 17% reduction in SAH/ICH mortality over the same time period. 
##As IS accounts for 80% of total strokes in NDPSS, and SAH/ICH accounts for 20%, overall decrease = 34%.
##(((0.38*8)+(0.17*2))/10)
##=3.4% per year

#estimate reduction from 2010 (year data collected) and 2018
#reduced rate for 2018 is also applied to 2015-2017 for simplification 


reduc <- (0.034*8)

tps_S[,"Stroke1m" := Stroke1m*(1-reduc)]

##copy Stroke1m to new variable P2.10 which is the relevant tpr.  keep Stroke1m for use in calc. 
tps_S[,"P2.10" := Stroke1m ]


##Recurrent stroke mortality

#same as overall case fatality, per NDPSS
#This is varied in SA

tps_S[,"Rec_mr":=(Stroke1m)]






####15. Calculating transition probabilities####

##In this section we combine the parameter estimates to calculate the transition probabilities for the model
##Note: Health states 50 and 51 are defined as new incident cases of CIND and dementia, and these are calculated within the model 




##first calculate health-state specific mortality risk##

tps_S[,NCInDmr := (mr*NCInD_RR)]
tps_S[,NCIDmr := (mr*NCID_RR)]
tps_S[,CINDnDmr := (mr*CINDnD_RR)]
tps_S[,CINDDmr := (mr*CINDD_RR)]
tps_S[,Dem_mr := (mr*Dem_RR)]


#Delete unnecessary columns
tps_S[,"strRR":= NULL]
tps_S[,"NCInD_RR":= NULL]
tps_S[,"NCID_RR":= NULL]
tps_S[,"CINDnD_RR":= NULL]
tps_S[,"CINDD_RR":= NULL]
tps_S[,"Dem_RR":= NULL]
tps_S[,"imp":= NULL]


##See TPr Table.xlsx for transition probability table
##var names for transition probabilities indicate the health states involved
##e.g. P1.2 is the probability of moving from state 1 (disease-free) to state 2 (incident stroke)


#Probability of remaining alive and disease-free

tps_S[, "P1.1":= (1-mr-P1.2)]

#P1.2 = probability of having a first ever stroke is calculated in Stroke Incidence section above

#Probability of dying for disese-free pop
tps_S[,"P1.8" := mr ]


##Health States at 12m

tps_S[,"P2.3":= (NCInD_initial*(1-Stroke1m-Rec_12m)*(1-(NCInDmr*(11/12))))]

tps_S[,"P2.4":= (NCID_initial*(1-Stroke1m-Rec_12m)*(1-(NCIDmr*(11/12))))]

tps_S[,"P2.5":= (CINDnD_initial*(1-Stroke1m-Rec_12m)*(1-(CINDnDmr*(11/12))))]

tps_S[,"P2.6":= (CINDD_initial*(1-Stroke1m-Rec_12m)*(1-(CINDDmr*(11/12))))]

tps_S[,"P2.7":= (Dementia_initial*(1-Stroke1m-Rec_12m)*(1-(Dem_mr*(11/12))))]

tps_S[,"P2.9":= (NCInD_initial*(1-Stroke1m-Rec_12m)*(NCInDmr*(11/12)))
                 +(NCID_initial*(1-Stroke1m-Rec_12m)*(NCIDmr*(11/12)))
                 +(CINDnD_initial*(1-Stroke1m-Rec_12m)*(CINDnDmr*(11/12)))
                 +(CINDD_initial*(1-Stroke1m-Rec_12m)*(CINDDmr*(11/12)))
                 +(Dementia_initial*(1-Stroke1m-Rec_12m)*(Dem_mr*(11/12)))]



##Transitions post 12m

tps_S[,"Dem_stable":= 1]


tps_S[,"P3.3":= ((1 - Rec - (NCInDmr))*NCInD_stable)]
tps_S[,"P3.4":= ((1 - Rec - (NCInDmr))*NCID)]
tps_S[,"P3.5":= ((1  - Rec - (NCInDmr))*NCInD_CINDnD)]
tps_S[,"P3.6":= ((1 - Rec - (NCInDmr))*NCInD_CINDD)]
tps_S[,"P3.9":= (NCInDmr)]
tps_S[,"P3.12":= Rec]
tps_S[,"P4.4":= ((1 - Rec - (NCIDmr))*NCID_stable)]
tps_S[,"P4.6":= ((1 - Rec - (NCIDmr))*NCID_CINDD)]
tps_S[,"P4.9":= (NCIDmr)]
tps_S[,"P4.13":= Rec]
tps_S[,"P5.5":= ((1 - Rec - (CINDnDmr))*CINDnD_stable)]
tps_S[,"P5.6":= ((1 - Rec - (CINDnDmr))*CINDnD_CINDD)]
tps_S[,"P5.7":= ((1 - Rec - (CINDnDmr))*CINDnD_Dem)]
tps_S[,"P5.9":= (CINDnDmr)]
tps_S[,"P5.14":= Rec]
tps_S[,"P6.6":= ((1 - Rec - (CINDDmr))*CINDD_stable)]
tps_S[,"P6.7":= ((1 - Rec - (CINDDmr))*CINDD_Dem)]
tps_S[,"P6.9":= (CINDDmr)]
tps_S[,"P6.15":= Rec]
tps_S[,"P7.7":= ((1 - Rec - (Dem_mr))*Dem_stable)]
tps_S[,"P7.9":= (Dem_mr)]
tps_S[,"P7.16":= Rec]

##Recurrent outcomes if recurrence within first 12 months

tps_S[,"P11.17":= (NCInD_pr*(1-Rec_mr-Rec_12m)*(1-(NCInDmr*(11/12))))]
tps_S[,"P11.18":= (NCID_pr*(1-Rec_mr-Rec_12m)*(1-(NCIDmr*(11/12))))]
tps_S[,"P11.19":= (CINDnD_pr*(1-Rec_mr-Rec_12m)*(1-(CINDnDmr*(11/12))))]
tps_S[,"P11.20":= (CINDD_pr*(1-Rec_mr-Rec_12m)*(1-(CINDDmr*(11/12))))]
tps_S[,"P11.21":= (Dem_pr*(1-Rec_mr-Rec_12m)*(1-(Dem_mr*(11/12))))]
tps_S[,"P11.22":= Rec_mr]
tps_S[,"P11.23":= ((NCInD_pr*(1-Rec_mr-Rec_12m)*(NCInDmr*(11/12)))
                   +(NCID_pr*(1-Rec_mr-Rec_12m)*(NCIDmr*(11/12)))
                   +(CINDnD_pr*(1-Rec_mr-Rec_12m)*(CINDnDmr*(11/12)))
                   +(CINDD_pr*(1-Rec_mr-Rec_12m)*(CINDDmr*(11/12)))
                   +(Dem_pr*(1-Rec_mr-Rec_12m)*(Dem_mr*(11/12))))]

##Recurrent outcomes if recurrence is post 12 months (based on ELSA data)

tps_S[,"P12.17":= (pr_NCInD_stable*(1-Rec_mr-Rec_12m)*(1-(NCInDmr*(11/12))))]
tps_S[,"P12.18":= (pr_NCID*(1-Rec_mr-Rec_12m)*(1-(NCIDmr*(11/12))))]
tps_S[,"P12.19":= (pr_NCInD_CINDnD*(1-Rec_mr-Rec_12m)*(1-(NCInDmr*(11/12))))]
tps_S[,"P12.20":= (pr_NCInD_CINDD*(1-Rec_mr-Rec_12m)*(1-(CINDDmr*(11/12))))]
tps_S[,"P12.21":= (pr_NCInD_Dem*(1-Rec_mr-Rec_12m)*(1-(Dem_mr*(11/12))))]
tps_S[,"P12.22":= Rec_mr]
tps_S[,"P12.23":= ((pr_NCInD_stable*(1-Rec_mr-Rec_12m)*(NCInDmr*(11/12)))
                   +(pr_NCID*(1-Rec_mr-Rec_12m)*(NCIDmr*(11/12)))
                   +(pr_NCInD_CINDnD*(1-Rec_mr-Rec_12m)*(NCInDmr*(11/12)))
                   +(pr_NCInD_CINDD*(1-Rec_mr-Rec_12m)*(CINDDmr*(11/12)))
                   +(pr_NCInD_Dem*(1-Rec_mr-Rec_12m)*(Dem_mr*(11/12))))]



tps_S[,"P13.18":= (pr_NCID_stable*(1-Rec_mr-Rec_12m)*(1-(NCIDmr*(11/12))))]
tps_S[,"P13.20":= (pr_NCID_CINDD*(1-Rec_mr-Rec_12m)*(1-(CINDDmr*(11/12))))]
tps_S[,"P13.21":= ((1-(Rec_mr + (Dem_mr*(11/12)))-Rec_12m)*pr_NCID_Dem)]
tps_S[,"P13.22":= Rec_mr]
tps_S[,"P13.23":= (pr_NCID_stable*(NCIDmr*(11/12))
                   +pr_NCID_CINDD*(CINDDmr*(11/12))
                   +pr_NCID_Dem*(Dem_mr*(11/12)))]

tps_S[,"P14.19":= (pr_CINDnD_stable*(1-Rec_mr-Rec_12m)*(1-(CINDnDmr*(11/12))))]
tps_S[,"P14.20":= (pr_CINDnD_CINDD*(1-Rec_mr-Rec_12m)*(1-(CINDDmr*(11/12))))]
tps_S[,"P14.21":= (pr_CINDnD_Dem*(1-Rec_mr-Rec_12m)*(1-(Dem_mr*(11/12))))]
tps_S[,"P14.22":= Rec_mr]
tps_S[,"P14.23":=((pr_CINDnD_stable*(1-Rec_mr-Rec_12m)*(CINDnDmr*(11/12)))
                  +(pr_CINDnD_CINDD*(1-Rec_mr-Rec_12m)*(CINDDmr*(11/12)))
                  +(pr_CINDnD_Dem*(1-Rec_mr-Rec_12m)*(Dem_mr*(11/12))))]

tps_S[,"P15.20":= (pr_CINDD_stable*(1-Rec_mr-Rec_12m)*(1-(CINDDmr*(11/12))))]
tps_S[,"P15.21":= (pr_CINDD_Dem*(1-Rec_mr-Rec_12m)*(1-(Dem_mr*(11/12))))]
tps_S[,"P15.22":= Rec_mr]
tps_S[,"P15.23":= ((pr_CINDD_stable*(1-Rec_mr-Rec_12m)*(CINDDmr*(11/12)))
                   +(pr_CINDD_Dem*(1-Rec_mr-Rec_12m)*(Dem_mr*(11/12))))]

tps_S[,"P16.21":= (1-Rec_mr-Rec_12m)*(1-(Dem_mr*(11/12)))]
tps_S[,"P16.22":= Rec_mr]
tps_S[,"P16.23":=((1-Rec_mr-Rec_12m)*(Dem_mr*(11/12)))]

##Health state transitions post a recurrent stroke


tps_S[,"P17.17":= (1 - Rec - (NCInDmr))*NCInD_stable]
tps_S[,"P17.18":= (1 - Rec - (NCInDmr))*NCID]
tps_S[,"P17.19":= (1  - Rec - (NCInDmr))*NCInD_CINDnD]
tps_S[,"P17.20":= (1 - Rec - (NCInDmr))*NCInD_CINDD]
tps_S[,"P17.23":= NCInDmr]

tps_S[,"P18.18":= (1 - Rec - (NCIDmr))*NCID_stable]
tps_S[,"P18.20":= (1 - Rec - (NCIDmr))*NCID_CINDD]
tps_S[,"P18.23":= NCIDmr]


tps_S[,"P19.19":= (1 - Rec - (CINDnDmr))*CINDnD_stable]
tps_S[,"P19.20":= (1 - Rec - (CINDnDmr))*CINDnD_CINDD]
tps_S[,"P19.21":= (1 - Rec - (CINDnDmr))*CINDnD_Dem]
tps_S[,"P19.23":= CINDnDmr]

tps_S[,"P20.20":= (1 - Rec - (CINDDmr))*CINDD_stable]
tps_S[,"P20.21":= (1 - Rec - (CINDDmr))*CINDD_Dem]
tps_S[,"P20.23":= CINDDmr]

tps_S[,"P21.21":= (1 - Rec - (Dem_mr))*Dem_stable]
tps_S[,"P21.23":= Dem_mr]


##SECOND RECURRENCE##

##transitions to recurrent stroke tunnel states

#recurrence within 12 months for first recurrent strokes
#where the initial recurrence was within 12m of first stroke

tps_S[,"P11.24":= Rec_12m]


#where the initial recurrence was after 12m of first stroke 
#case had been allocated to a health state after the first stroke
#we use this information to allocate this recurrence to a tunnel health state
#based on ELSA recurrent stroke probabilities

tps_S[,"P12.25":= Rec_12m*pr_NCInD_stable]
tps_S[,"P12.26":= Rec_12m*pr_NCID]
tps_S[,"P12.27":= Rec_12m*pr_NCInD_CINDnD]
tps_S[,"P12.28":= Rec_12m*pr_NCInD_CINDD]
tps_S[,"P12.29":= Rec_12m*pr_NCInD_Dem]


tps_S[,"P13.26":= Rec_12m*pr_NCID_stable]
tps_S[,"P13.28":= Rec_12m*pr_NCID_CINDD]
tps_S[,"P13.29":= Rec_12m*pr_NCID_Dem]

tps_S[,"P14.27":= Rec_12m*pr_CINDnD_stable]
tps_S[,"P14.28":= Rec_12m*pr_CINDnD_CINDD]
tps_S[,"P14.29":= Rec_12m*pr_CINDnD_Dem]

tps_S[,"P15.28":= Rec_12m*pr_CINDD_stable]
tps_S[,"P15.29":= Rec_12m*pr_CINDD_Dem]

#all cases with dementia pre-recurrence transition to dementia
tps_S[,"P16.29":= Rec_12m]


#recurrent strokes that occur >12m after the first recurrence

tps_S[,"P17.25":= Rec]
tps_S[,"P18.26":= Rec]
tps_S[,"P19.27":= Rec]
tps_S[,"P20.28":= Rec]
tps_S[,"P21.29":= Rec]
      

##transitions from tunnel states - 2nd recurrence


##Recurrent outcomes if recurrence within first 12 months

tps_S[,"P24.30":= (NCInD_pr*(1-Rec_mr-Rec_12m)*(1-(NCInDmr*(11/12))))]
tps_S[,"P24.31":= (NCID_pr*(1-Rec_mr-Rec_12m)*(1-(NCIDmr*(11/12))))]
tps_S[,"P24.32":= (CINDnD_pr*(1-Rec_mr-Rec_12m)*(1-(CINDnDmr*(11/12))))]
tps_S[,"P24.33":= (CINDD_pr*(1-Rec_mr-Rec_12m)*(1-(CINDDmr*(11/12))))]
tps_S[,"P24.34":= (Dem_pr*(1-Rec_mr-Rec_12m)*(1-(Dem_mr*(11/12))))]
tps_S[,"P24.35":= Rec_mr]
tps_S[,"P24.36":= ((NCInD_pr*(1-Rec_mr-Rec_12m)*(NCInDmr*(11/12)))
                   +(NCID_pr*(1-Rec_mr-Rec_12m)*(NCIDmr*(11/12)))
                   +(CINDnD_pr*(1-Rec_mr-Rec_12m)*(CINDnDmr*(11/12)))
                   +(CINDD_pr*(1-Rec_mr-Rec_12m)*(CINDDmr*(11/12)))
                   +(Dem_pr*(1-Rec_mr-Rec_12m)*(Dem_mr*(11/12))))]

##Recurrent outcomes if recurrence is post 12 months (based on ELSA data)

tps_S[,"P25.30":= (pr_NCInD_stable*(1-Rec_mr-Rec_12m)*(1-(NCInDmr*(11/12))))]
tps_S[,"P25.31":= (pr_NCID*(1-Rec_mr-Rec_12m)*(1-(NCIDmr*(11/12))))]
tps_S[,"P25.32":= (pr_NCInD_CINDnD*(1-Rec_mr-Rec_12m)*(1-(CINDnDmr*(11/12))))]
tps_S[,"P25.33":= (pr_NCInD_CINDD*(1-Rec_mr-Rec_12m)*(1-(CINDDmr*(11/12))))]
tps_S[,"P25.34":= (pr_NCInD_Dem*(1-Rec_mr-Rec_12m)*(1-(Dem_mr*(11/12))))]
tps_S[,"P25.35":= Rec_mr]
tps_S[,"P25.36":= ((pr_NCInD_stable*(1-Rec_mr-Rec_12m)*(NCInDmr*(11/12)))
                   +(pr_NCID*(1-Rec_mr-Rec_12m)*(NCIDmr*(11/12)))
                   +(pr_NCInD_CINDnD*(1-Rec_mr-Rec_12m)*(NCInDmr*(11/12)))
                   +(pr_NCInD_CINDD*(1-Rec_mr-Rec_12m)*(CINDDmr*(11/12)))
                   +(pr_NCInD_Dem*(1-Rec_mr-Rec_12m)*(Dem_mr*(11/12))))]



tps_S[,"P26.31":= (pr_NCID_stable*(1-Rec_mr-Rec_12m)*(1-(NCIDmr*(11/12))))]
tps_S[,"P26.33":= (pr_NCID_CINDD*(1-Rec_mr-Rec_12m)*(1-(CINDDmr*(11/12))))]
tps_S[,"P26.34":= (pr_NCID_Dem*(1-Rec_mr-Rec_12m)*(1-(Dem_mr*(11/12))))]
tps_S[,"P26.35":= Rec_mr]
tps_S[,"P26.36":= ((pr_NCID_stable*(1-Rec_mr-Rec_12m)*(NCIDmr*(11/12)))
                   +(pr_NCID_CINDD*(1-Rec_mr-Rec_12m)*(CINDDmr*(11/12)))
                   +(pr_NCID_Dem*(1-Rec_mr-Rec_12m)*(Dem_mr*(11/12))))]




tps_S[,"P27.32":= (pr_CINDnD_stable*(1-Rec_mr-Rec_12m)*(1-(CINDnDmr*(11/12))))]
tps_S[,"P27.33":= (pr_CINDnD_CINDD*(1-Rec_mr-Rec_12m)*(1-(CINDDmr*(11/12))))]
tps_S[,"P27.34":= (pr_CINDnD_Dem*(1-Rec_mr-Rec_12m)*(1-(Dem_mr*(11/12))))]
tps_S[,"P27.35":= Rec_mr]
tps_S[,"P27.36":= ((pr_CINDnD_stable*(1-Rec_mr-Rec_12m)*(CINDnDmr*(11/12)))
                   +(pr_CINDnD_CINDD*(1-Rec_mr-Rec_12m)*(CINDDmr*(11/12)))
                   +(pr_CINDnD_Dem*(1-Rec_mr-Rec_12m)*(Dem_mr*(11/12))))]


tps_S[,"P28.33":= (pr_CINDD_stable*(1-Rec_mr-Rec_12m)*(1-(CINDDmr*(11/12))))]
tps_S[,"P28.34":= (pr_CINDD_Dem*(1-Rec_mr-Rec_12m)*(1-(Dem_mr*(11/12))))]
tps_S[,"P28.35":= Rec_mr]
tps_S[,"P28.36":= ((pr_CINDD_stable*(1-Rec_mr-Rec_12m)*(CINDDmr*(11/12)))
                   +(pr_CINDD_Dem*(1-Rec_mr-Rec_12m)*(Dem_mr*(11/12))))]

tps_S[,"P29.34":= (1-Rec_mr-Rec_12m)*(1-(Dem_mr*(11/12)))]
tps_S[,"P29.35":= Rec_mr]
tps_S[,"P29.36":=((1-Rec_mr-Rec_12m)*(Dem_mr*(11/12)))]


##Health state transitions post 2nd recurrent stroke



tps_S[,"P30.30":= (1 - Rec - (NCInDmr))*NCInD_stable]
tps_S[,"P30.31":= (1 - Rec - (NCInDmr))*NCID]
tps_S[,"P30.32":= (1  - Rec - (NCInDmr))*NCInD_CINDnD]
tps_S[,"P30.33":= (1 - Rec - (NCInDmr))*NCInD_CINDD]
tps_S[,"P30.36":= NCInDmr]



tps_S[,"P31.31":= (1 - Rec - (NCIDmr))*NCID_stable]
tps_S[,"P31.33":= (1 - Rec - (NCIDmr))*NCID_CINDD]
tps_S[,"P31.36":= NCIDmr]


tps_S[,"P32.32":= (1 - Rec - (CINDnDmr))*CINDnD_stable]
tps_S[,"P32.33":= (1 - Rec - (CINDnDmr))*CINDnD_CINDD]
tps_S[,"P32.34":= (1 - Rec - (CINDnDmr))*CINDnD_Dem]
tps_S[,"P32.36":= CINDnDmr]

tps_S[,"P33.33":= (1 - Rec - (CINDDmr))*CINDD_stable]
tps_S[,"P33.34":= (1 - Rec - (CINDDmr))*CINDD_Dem]
tps_S[,"P33.36":= CINDDmr]

tps_S[,"P34.34":= (1 - Rec - (Dem_mr))*Dem_stable]
tps_S[,"P34.36":= Dem_mr]

##THIRD RECURRENCE - TPrs copied from second recurrence


##transitions to recurrent stroke tunnel states

#recurrence within 12 months for second recurrent strokes
#where the initial recurrences were all within 12 months

tps_S[,"P24.37":= Rec_12m]


#where the initial recurrence was after 12m of 2nd recurrence
#case had been allocated to a health state after the second stroke
#we use this information to allocate this recurrence to a tunnel health state
#based on ELSA recurrent stroke probabilities

tps_S[,"P25.38":= Rec_12m*pr_NCInD_stable]
tps_S[,"P25.39":= Rec_12m*pr_NCID]
tps_S[,"P25.40":= Rec_12m*pr_NCInD_CINDnD]
tps_S[,"P25.41":= Rec_12m*pr_NCInD_CINDD]
tps_S[,"P25.42":= Rec_12m*pr_NCInD_Dem]


tps_S[,"P26.39":= Rec_12m*pr_NCID_stable]
tps_S[,"P26.41":= Rec_12m*pr_NCID_CINDD]
tps_S[,"P26.42":= Rec_12m*pr_NCID_Dem]

tps_S[,"P27.40":= Rec_12m*pr_CINDnD_stable]
tps_S[,"P27.41":= Rec_12m*pr_CINDnD_CINDD]
tps_S[,"P27.42":= Rec_12m*pr_CINDnD_Dem]

tps_S[,"P28.41":= Rec_12m*pr_CINDD_stable]
tps_S[,"P28.42":= Rec_12m*pr_CINDD_Dem]

#all cases with dementia pre-recurrence transition to dementia
tps_S[,"P29.42":= Rec_12m]


#recurrent strokes that occur >12m after the 2nd recurrence

tps_S[,"P30.38":= Rec]
tps_S[,"P31.39":= Rec]
tps_S[,"P32.40":= Rec]
tps_S[,"P33.41":= Rec]
tps_S[,"P34.42":= Rec]


##transitions from tunnel states - 3rd recurrence
##copied from P24.30 to P29.36
##NB: Last recurrence, so recurrence probabilities are not removed from these.


##Recurrent outcomes if recurrence within first 12 months

tps_S[,"P37.43":= (NCInD_pr*(1-Rec_mr)*(1-(NCInDmr*(11/12))))]
tps_S[,"P37.44":= (NCID_pr*(1-Rec_mr)*(1-(NCIDmr*(11/12))))]
tps_S[,"P37.45":= (CINDnD_pr*(1-Rec_mr)*(1-(CINDnDmr*(11/12))))]
tps_S[,"P37.46":= (CINDD_pr*(1-Rec_mr)*(1-(CINDDmr*(11/12))))]
tps_S[,"P37.47":= (Dem_pr*(1-Rec_mr)*(1-(Dem_mr*(11/12))))]
tps_S[,"P37.48":= Rec_mr]
tps_S[,"P37.49":= ((NCInD_pr*(1-Rec_mr)*(NCInDmr*(11/12)))
                   +(NCID_pr*(1-Rec_mr)*(NCIDmr*(11/12)))
                   +(CINDnD_pr*(1-Rec_mr)*(CINDnDmr*(11/12)))
                   +(CINDD_pr*(1-Rec_mr)*(CINDDmr*(11/12)))
                   +(Dem_pr*(1-Rec_mr)*(Dem_mr*(11/12))))]

##Recurrent outcomes if recurrence is post 12 months (based on ELSA data)

tps_S[,"P38.43":= (pr_NCInD_stable*(1-Rec_mr)*(1-(NCInDmr*(11/12))))]
tps_S[,"P38.44":= (pr_NCID*(1-Rec_mr)*(1-(NCIDmr*(11/12))))]
tps_S[,"P38.45":= (pr_NCInD_CINDnD*(1-Rec_mr)*(1-(NCInDmr*(11/12))))]
tps_S[,"P38.46":= (pr_NCInD_CINDD*(1-Rec_mr)*(1-(CINDDmr*(11/12))))]
tps_S[,"P38.47":= (pr_NCInD_Dem*(1-Rec_mr)*(1-(Dem_mr*(11/12))))]
tps_S[,"P38.48":= Rec_mr]
tps_S[,"P38.49":= ((pr_NCInD_stable*(1-Rec_mr)*(NCInDmr*(11/12)))
                   +(pr_NCID*(1-Rec_mr)*(NCIDmr*(11/12)))
                   +(pr_NCInD_CINDnD*(1-Rec_mr)*(NCInDmr*(11/12)))
                   +(pr_NCInD_CINDD*(1-Rec_mr)*(CINDDmr*(11/12)))
                   +(pr_NCInD_Dem*(1-Rec_mr)*(Dem_mr*(11/12))))]



tps_S[,"P39.44":= (pr_NCID_stable*(1-Rec_mr)*(1-(NCIDmr*(11/12))))]
tps_S[,"P39.46":= (pr_NCID_CINDD*(1-Rec_mr)*(1-(CINDDmr*(11/12))))]
tps_S[,"P39.47":= (pr_NCID_Dem*(1-Rec_mr)*(1-(Dem_mr*(11/12))))]
tps_S[,"P39.48":= Rec_mr]
tps_S[,"P39.49":= ((pr_NCInD_stable*(1-Rec_mr)*(NCInDmr*(11/12)))
                   +(pr_NCID*(1-Rec_mr)*(NCIDmr*(11/12)))
                   +(pr_NCInD_CINDnD*(1-Rec_mr)*(NCInDmr*(11/12)))
                   +(pr_NCInD_CINDD*(1-Rec_mr)*(CINDDmr*(11/12)))
                   +(pr_NCInD_Dem*(1-Rec_mr)*(Dem_mr*(11/12))))]

tps_S[,"P40.45":= (pr_CINDnD_stable*(1-Rec_mr)*(1-(CINDnDmr*(11/12))))]
tps_S[,"P40.46":= (pr_CINDnD_CINDD*(1-Rec_mr)*(1-(CINDDmr*(11/12))))]
tps_S[,"P40.47":= (pr_CINDnD_Dem*(1-Rec_mr)*(1-(Dem_mr*(11/12))))]
tps_S[,"P40.48":= Rec_mr]
tps_S[,"P40.49":= ((pr_CINDnD_stable*(1-Rec_mr)*(CINDnDmr*(11/12)))
                   +(pr_CINDnD_CINDD*(1-Rec_mr)*(CINDDmr*(11/12)))
                   +(pr_CINDnD_Dem*(1-Rec_mr)*(Dem_mr*(11/12))))]

tps_S[,"P41.46":= (pr_CINDD_stable*(1-Rec_mr)*(1-(CINDDmr*(11/12))))]
tps_S[,"P41.47":= (pr_CINDD_Dem*(1-Rec_mr)*(1-(Dem_mr*(11/12))))]
tps_S[,"P41.48":= Rec_mr]
tps_S[,"P41.49":= ((pr_CINDD_stable*(1-Rec_mr)*(CINDDmr*(11/12)))
                   +(pr_CINDD_Dem*(1-Rec_mr)*(Dem_mr*(11/12))))]

tps_S[,"P42.47":= (1-Rec_mr)*(1-(Dem_mr*(11/12)))]
tps_S[,"P42.48":= Rec_mr]
tps_S[,"P42.49":=((1-Rec_mr)*(Dem_mr*(11/12)))]


##Health state transitions post 3rd recurrent stroke
##copied from P30.30 to P34.36


tps_S[,"P43.43":= (1 - (NCInDmr))*NCInD_stable]
tps_S[,"P43.44":= (1 - (NCInDmr))*NCID]
tps_S[,"P43.45":= (1 - (NCInDmr))*NCInD_CINDnD]
tps_S[,"P43.46":= (1 - (NCInDmr))*NCInD_CINDD]
tps_S[,"P43.49":= NCInDmr]



tps_S[,"P44.44":= (1 - (NCIDmr))*NCID_stable]
tps_S[,"P44.46":= (1 - (NCIDmr))*NCID_CINDD]
tps_S[,"P44.49":= NCIDmr]


tps_S[,"P45.45":= (1 - (CINDnDmr))*CINDnD_stable]
tps_S[,"P45.46":= (1 - (CINDnDmr))*CINDnD_CINDD]
tps_S[,"P45.47":= (1 - (CINDnDmr))*CINDnD_Dem]
tps_S[,"P45.49":= CINDnDmr]

tps_S[,"P46.46":= (1 - (CINDDmr))*CINDD_stable]
tps_S[,"P46.47":= (1 - (CINDDmr))*CINDD_Dem]
tps_S[,"P46.49":= CINDDmr]

tps_S[,"P47.47":= (1 - (Dem_mr))*Dem_stable]
tps_S[,"P47.49":= Dem_mr]


##delete unnecessary vars

tps_S[,"Stroke1m":= NULL]
tps_S[,"mr":= NULL]
tps_S[,"Rec_12m":= NULL]
      tps_S[,"Rec":= NULL]
      tps_S[,"NCInD_initial":= NULL]
      tps_S[,"NCID_initial":= NULL]
      tps_S[,"CINDnD_initial":= NULL]
      tps_S[,"CINDD_initial":= NULL]
      tps_S[,"Dementia_initial":= NULL]
      tps_S[,NCInDmr := NULL]
      tps_S[,NCIDmr := NULL]
      tps_S[,CINDnDmr := NULL]
      tps_S[,CINDDmr := NULL]
      tps_S[,Dem_mr := NULL]
      tps_S[,"NCInD_stable":= NULL]
      tps_S[,"NCID":= NULL]
      tps_S[,"NCInD_CINDnD":= NULL]
      tps_S[,"NCInD_CINDD":= NULL]
      tps_S[,"NCID_stable":= NULL]
      tps_S[,"NCID_CINDD":= NULL]
      tps_S[,"CINDnD_stable":= NULL]
      tps_S[,"CINDnD_CINDD":= NULL]
      tps_S[,"CINDnD_Dem":= NULL]
      tps_S[,"CINDD_stable":= NULL]
      tps_S[,"CINDD_Dem":= NULL]
      tps_S[,"Dem_stable":= NULL]
      tps_S[,"NCInD_pr":= NULL]
      tps_S[,"NCID_pr":= NULL]
      tps_S[,"CINDnD_pr":= NULL]
      tps_S[,"CINDD_pr":= NULL]
      tps_S[,"Dem_pr":= NULL]
      tps_S[,"pr_NCInD_stable":= NULL]
      tps_S[,"pr_NCID":= NULL]
      tps_S[,"pr_NCInD_CINDnD":= NULL]
      tps_S[,"pr_NCInD_CINDD":= NULL]
      tps_S[,"pr_NCInD_Dem":= NULL]
      tps_S[,"pr_NCID_stable":= NULL]
      tps_S[,"pr_NCID_CINDD":= NULL]
      tps_S[,"pr_NCID_Dem":= NULL]
      tps_S[,"pr_CINDnD_stable":= NULL]
      tps_S[,"pr_CINDnD_CINDD":= NULL]
      tps_S[,"pr_CINDnD_Dem":= NULL]
      tps_S[,"pr_CINDD_stable":= NULL]
      tps_S[,"pr_CINDD_Dem":= NULL]
      tps_S[,"Rec_mr":= NULL]
      





##sort by iteration
setkey(tps_S,iter)

##set column order
setcolorder(tps_S,c("Age","iter","P1.1","P1.2","P1.8","P2.3",
                   "P2.4","P2.5","P2.6","P2.7","P2.9",
                    "P2.10","P2.11","P3.3","P3.4","P3.5","P3.6",
                   "P3.9","P3.12"))






####16. Calculating model #####

 
# calculate number of people in states in initial year
# by merging initial population estimates with table of prevalence estimates

initial.states<-merge(inipop,prev, by="Age")

initial.states[,`:=`(s.3=s.3*Pop,s.4=s.4*Pop,s.5=s.5*Pop,s.6=s.6*Pop,s.7=s.7*Pop)]
initial.states[,`:=`(s.1=Pop-s.3-s.4-s.5-s.6-s.7)]

#incident stroke is zero in first year, included within prevalent strokes for that year
initial.states[,s.2:=0]

#sets Year to first year of model 

initial.states[,Year:=FIRST.YEAR]

setcolorder(initial.states,c("iter", "Age","Year","Pop","s.1","s.2","s.3","s.4","s.5",
                             "s.6","s.7","s.8","s.9","s.10","s.11","s.12","s.13",
                             "s.14","s.15","s.16","s.17","s.18","s.19","s.20","s.21",
                             "s.22","s.23","s.24","s.25"))

#put out excel file of starting prevalences for checking
#takes a long time so only if needed
#write.xlsx(initial.states,"initial.xlsx", rows=1:250,cols=1:27)
  
#merge initial states with transition probabilities  
initial.states<-merge(initial.states, tps_S, by=c("Age","iter"))

#sort by iteration
setkey(initial.states,iter)

#saves out the results for year 1

res<-list()
res[[paste0(FIRST.YEAR)]]<-copy(initial.states)

#moves_trace <- list()

#create moves table for calculation of transitions between states

moves<-copy(initial.states)



y<-FIRST.YEAR

##Loop starts here - current year is stored in y
##loop continues until y = LAST.YEAR (defined at beginning of code)

while (y<LAST.YEAR) 
{

  
#age and year are increased by 1 year to move model along by one year  
    
moves[,Age:=Age+1]
moves[,Year:=Year+1]

#Taking out existing tprs - need to re-merge, as age not synced up

moves[,`:=`(
  P1.1=NULL,
  P1.2=NULL,
  P1.8=NULL,
  P2.3=NULL,
  P2.4=NULL,
  P2.5=NULL,
  P2.6=NULL,
  P2.7=NULL,
  P2.9=NULL,
  P2.10=NULL,
  P2.11=NULL,
  P3.3=NULL,
  P3.4=NULL,
  P3.5=NULL,
  P3.6=NULL,
  P3.9=NULL,
  P3.12=NULL,
  P4.4=NULL,
  P4.6=NULL,
  P4.9=NULL,
  P4.13=NULL,
  P5.5=NULL,
  P5.6=NULL,
  P5.7=NULL,
  P5.9=NULL,
  P5.14=NULL,
  P6.6=NULL,
  P6.7=NULL,
  P6.9=NULL,
  P6.15=NULL,
  P7.7=NULL,
  P7.9=NULL,
  P7.16=NULL,
  P11.17=NULL,
  P11.18=NULL,
  P11.19=NULL,
  P11.20=NULL,
  P11.21=NULL,
  P11.22=NULL,
  P11.23=NULL,
  P12.17=NULL,
  P12.18=NULL,
  P12.19=NULL,
  P12.20=NULL,
  P12.21=NULL,
  P12.22=NULL,
  P12.23=NULL,

  
  P13.18=NULL,
  P13.20=NULL,
  P13.21=NULL,
  P13.22=NULL,
  P13.23=NULL,
  
  P14.19=NULL,
  P14.20=NULL,
  P14.21=NULL,
  P14.22=NULL,
  P14.23=NULL,
  P15.20=NULL,
  P15.21=NULL,
  P15.22=NULL,
  P15.23=NULL,
  P16.21=NULL,
  P16.22=NULL,
  P16.23=NULL,
  P17.17=NULL,
  P17.18=NULL,
  P17.19=NULL,
  P17.20=NULL,
  P17.23=NULL,
  P18.18=NULL,
  P18.20=NULL,
  P18.23=NULL,
  P19.19=NULL,
  P19.20=NULL,
  P19.21=NULL,
  P19.23=NULL,
  P20.20=NULL,
  P20.21=NULL,
  P20.23=NULL,
  P21.21=NULL,
  P21.23=NULL,
  P11.24=NULL,
  P12.25=NULL,
  P12.26=NULL,
  P12.27=NULL,
  P12.28=NULL,
  P12.29=NULL,
  P13.26=NULL,
  P13.28=NULL,
  P13.29=NULL,
  P14.27=NULL,
  P14.28=NULL,
  P14.29=NULL,
  P15.28=NULL,
  P15.29=NULL,
  P16.29=NULL,
  P17.25=NULL,
  P18.26=NULL,
  P19.27=NULL,
  P20.28=NULL,
  P21.29=NULL,
  P24.30=NULL,
  P24.31=NULL,
  P24.32=NULL,
  P24.33=NULL,
  P24.34=NULL,
  P24.35=NULL,
  P24.36=NULL,
  P25.30=NULL,
  P25.31=NULL,
  P25.32=NULL,
  P25.33=NULL,
  P25.34=NULL,
  P25.35=NULL,
  P25.36=NULL,
  P26.31=NULL,
  P26.33=NULL,
  P26.34=NULL,
  P26.35=NULL,
  P26.36=NULL,
  P27.32=NULL,
  P27.33=NULL,
  P27.34=NULL,
  P27.35=NULL,
  P27.36=NULL,
  P28.33=NULL,
  P28.34=NULL,
  P28.35=NULL,
  P28.36=NULL,
  P29.34=NULL,
  P29.35=NULL,
  P29.36=NULL,
  P30.30=NULL,
  P30.31=NULL,
  P30.32=NULL,
  P30.33=NULL,
  P30.36=NULL,
  P31.31=NULL,
  P31.33=NULL,
  P31.36=NULL,
  P32.32=NULL,
  P32.33=NULL,
  P32.34=NULL,
  P32.36=NULL,
  P33.33=NULL,
  P33.34=NULL,
  P33.36=NULL,
  P34.34=NULL,
  P34.36=NULL,
  P24.37=NULL,
  P25.38=NULL,
  P25.39=NULL,
  P25.40=NULL,
  P25.41=NULL,
  P25.42=NULL,
  P26.39=NULL,
  P26.41=NULL,
  P26.42=NULL,
  P27.40=NULL,
  P27.41=NULL,
  P27.42=NULL,
  P28.41=NULL,
  P28.42=NULL,
  P29.42=NULL,
  P30.38=NULL,
  P31.39=NULL,
  P32.40=NULL,
  P33.41=NULL,
  P34.42=NULL,
  P37.43=NULL,
  P37.44=NULL,
  P37.45=NULL,
  P37.46=NULL,
  P37.47=NULL,
  P37.48=NULL,
  P37.49=NULL,
  P38.43=NULL,
  P38.44=NULL,
  P38.45=NULL,
  P38.46=NULL,
  P38.47=NULL,
  P38.48=NULL,
  P38.49=NULL,
  P39.44=NULL,
  P39.46=NULL,
  P39.47=NULL,
  P39.48=NULL,
  P39.49=NULL,
  P40.45=NULL,
  P40.46=NULL,
  P40.47=NULL,
  P40.48=NULL,
  P40.49=NULL,
  P41.46=NULL,
  P41.47=NULL,
  P41.48=NULL,
  P41.49=NULL,
  P42.47=NULL,
  P42.48=NULL,
  P42.49=NULL,
  P43.43=NULL,
  P43.44=NULL,
  P43.45=NULL,
  P43.46=NULL,
  P43.49=NULL,
  P44.44=NULL,
  P44.46=NULL,
  P44.49=NULL,
  P45.45=NULL,
  P45.46=NULL,
  P45.47=NULL,
  P45.49=NULL,
  P46.46=NULL,
  P46.47=NULL,
  P46.49=NULL,
  P47.47=NULL,
  P47.49=NULL

)]

#merging in tps again so that age is synced up correctly

moves<-merge(moves, tps_S, by=c("Age","iter"))

#M50 and M51 included to capture incident CIND and dementia
#Health states are multipled by the relevant transition probabilites to calculate the numbers transitioning between states
#the M var names indicate between which states the transition is occuring
#e.g. M1.2 is the numbers moving between state 1 and 2, which is calculated by applying the tpr p1.2 to state 1 (s.1)

moves[,`:=`(
  M1.1=s.1*P1.1,
  M1.2=s.1*P1.2,
  M1.8=s.1*P1.8,
  M2.3=s.2*P2.3,
  M2.4=s.2*P2.4,
  M2.5=s.2*P2.5,
  M2.6=s.2*P2.6,
  M2.7=s.2*P2.7,
  M2.9=s.2*P2.9,
  M2.10=s.2*P2.10,
  M2.11=s.2*P2.11,
  M3.3=s.3*P3.3,
  M3.4=s.3*P3.4,
  M3.5=s.3*P3.5,
  M3.6=s.3*P3.6,
  M3.9=s.3*P3.9,
  M3.12=s.3*P3.12,
  M4.4=s.4*P4.4,
  M4.6=s.4*P4.6,
  M4.9=s.4*P4.9,
  M4.13=s.4*P4.13,
  M5.5=s.5*P5.5,
  M5.6=s.5*P5.6,
  M5.7=s.5*P5.7,
  M5.9=s.5*P5.9,
  M5.14=s.5*P5.14,
  M6.6=s.6*P6.6,
  M6.7=s.6*P6.7,
  M6.9=s.6*P6.9,
  M6.15=s.6*P6.15,
  M7.7=s.7*P7.7,
  M7.9=s.7*P7.9,
  M7.16=s.7*P7.16,
  M11.17=s.11*P11.17,
  M11.18=s.11*P11.18,
  M11.19=s.11*P11.19,
  M11.20=s.11*P11.20,
  M11.21=s.11*P11.21,
  M11.22=s.11*P11.22,
  M11.23=s.11*P11.23,
  M12.17=s.12*P12.17,
  M12.18=s.12*P12.18,
  M12.19=s.12*P12.19,
  M12.20=s.12*P12.20,
  M12.21=s.12*P12.21,
  M12.22=s.12*P12.22,
  M12.23=s.12*P12.23,
  M13.18=s.13*P13.18,
  M13.20=s.13*P13.20,
  M13.21=s.13*P13.21,
  M13.22=s.13*P13.22,
  M13.23=s.13*P13.23,
  
  M14.19=s.14*P14.19,
  M14.20=s.14*P14.20,
  M14.21=s.14*P14.21,
  M14.22=s.14*P14.22,
  M14.23=s.14*P14.23,
  
  M15.20=s.15*P15.20,
  M15.21=s.15*P15.21,
  M15.22=s.15*P15.22,
  M15.23=s.15*P15.23,
  
  M16.21=s.16*P16.21,
  M16.22=s.16*P16.22,
  M16.23=s.16*P16.23,
  
  M17.17=s.17*P17.17,
  M17.18=s.17*P17.18,
  M17.19=s.17*P17.19,
  M17.20=s.17*P17.20,
  M17.23=s.17*P17.23,
  M18.18=s.18*P18.18,
  M18.20=s.18*P18.20,
  M18.23=s.18*P18.23,
  M19.19=s.19*P19.19,
  M19.20=s.19*P19.20,
  M19.21=s.19*P19.21,
  M19.23=s.19*P19.23,
  M20.20=s.20*P20.20,
  M20.21=s.20*P20.21,
  M20.23=s.20*P20.23,
  M21.21=s.21*P21.21,
  M21.23=s.21*P21.23,
  
  
  M11.24=s.11*P11.24,
  M12.25=s.12*P12.25,
  M12.26=s.12*P12.26,
  M12.27=s.12*P12.27,
  M12.28=s.12*P12.28,
  M12.29=s.12*P12.29,
  M13.26=s.13*P13.26,
  M13.28=s.13*P13.28,
  M13.29=s.13*P13.29,
  M14.27=s.14*P14.27,
  M14.28=s.14*P14.28,
  M14.29=s.14*P14.29,
  M15.28=s.15*P15.28,
  M15.29=s.15*P15.29,
  M16.29=s.16*P16.29,
  M17.25=s.17*P17.25,
  
  M18.26=s.18*P18.26,
  M19.27=s.19*P19.27,
  M20.28=s.20*P20.28,
  M21.29=s.21*P21.29,
  M24.30=s.24*P24.30,
  M24.31=s.24*P24.31,
  M24.32=s.24*P24.32,
  M24.33=s.24*P24.33,
  M24.34=s.24*P24.34,
  M24.35=s.24*P24.35,
  M24.36=s.24*P24.36,
  M25.30=s.25*P25.30,
  M25.31=s.25*P25.31,
  M25.32=s.25*P25.32,
  M25.33=s.25*P25.33,
  M25.34=s.25*P25.34,
  M25.35=s.25*P25.35,
  M25.36=s.25*P25.36,
  M26.31=s.26*P26.31,
  M26.33=s.26*P26.33,
  M26.34=s.26*P26.34,
  M26.35=s.26*P26.35,
  M26.36=s.26*P26.36,
  M27.32=s.27*P27.32,
  M27.33=s.27*P27.33,
  M27.34=s.27*P27.34,
  M27.35=s.27*P27.35,
  M27.36=s.27*P27.36,
  M28.33=s.28*P28.33,
  M28.34=s.28*P28.34,
  M28.35=s.28*P28.35,
  M28.36=s.28*P28.36,
  M29.34=s.29*P29.34,
  M29.35=s.29*P29.35,
  M29.36=s.29*P29.36,
  M30.30=s.30*P30.30,
  M30.31=s.30*P30.31,
  M30.32=s.30*P30.32,
  M30.33=s.30*P30.33,
  M30.36=s.30*P30.36,
  M31.31=s.31*P31.31,
  M31.33=s.31*P31.33,
  M31.36=s.31*P31.36,
  M32.32=s.32*P32.32,
  M32.33=s.32*P32.33,
  M32.34=s.32*P32.34,
  M32.36=s.32*P32.36,
  M33.33=s.33*P33.33,
  M33.34=s.33*P33.34,
  M33.36=s.33*P33.36,
  M34.34=s.34*P34.34,
  M34.36=s.34*P34.36,
  M24.37=s.24*P24.37,
  M25.38=s.25*P25.38,
  M25.39=s.25*P25.39,
  M25.40=s.25*P25.40,
  M25.41=s.25*P25.41,
  M25.42=s.25*P25.42,
  M26.39=s.26*P26.39,
  M26.41=s.26*P26.41,
  M26.42=s.26*P26.42,
  M27.40=s.27*P27.40,
  M27.41=s.27*P27.41,
  M27.42=s.27*P27.42,
  M28.41=s.28*P28.41,
  M28.42=s.28*P28.42,
  M29.42=s.29*P29.42,
  M30.38=s.30*P30.38,
  M31.39=s.31*P31.39,
  M32.40=s.32*P32.40,
  M33.41=s.33*P33.41,
  M34.42=s.34*P34.42,
  M37.43=s.37*P37.43,
  M37.44=s.37*P37.44,
  M37.45=s.37*P37.45,
  M37.46=s.37*P37.46,
  M37.47=s.37*P37.47,
  M37.48=s.37*P37.48,
  M37.49=s.37*P37.49,
  M38.43=s.38*P38.43,
  M38.44=s.38*P38.44,
  M38.45=s.38*P38.45,
  M38.46=s.38*P38.46,
  M38.47=s.38*P38.47,
  M38.48=s.38*P38.48,
  M38.49=s.38*P38.49,
  M39.44=s.39*P39.44,
  M39.46=s.39*P39.46,
  M39.47=s.39*P39.47,
  M39.48=s.39*P39.48,
  M39.49=s.39*P39.49,
  M40.45=s.40*P40.45,
  M40.46=s.40*P40.46,
  M40.47=s.40*P40.47,
  M40.48=s.40*P40.48,
  M40.49=s.40*P40.49,
  M41.46=s.41*P41.46,
  M41.47=s.41*P41.47,
  M41.48=s.41*P41.48,
  M41.49=s.41*P41.49,
  M42.47=s.42*P42.47,
  M42.48=s.42*P42.48,
  M42.49=s.42*P42.49,
  M43.43=s.43*P43.43,
  M43.44=s.43*P43.44,
  M43.45=s.43*P43.45,
  M43.46=s.43*P43.46,
  M43.49=s.43*P43.49,
  M44.44=s.44*P44.44,
  M44.46=s.44*P44.46,
  M44.49=s.44*P44.49,
  M45.45=s.45*P45.45,
  M45.46=s.45*P45.46,
  M45.47=s.45*P45.47,
  M45.49=s.45*P45.49,
  M46.46=s.46*P46.46,
  M46.47=s.46*P46.47,
  M46.49=s.46*P46.49,
  M47.47=s.47*P47.47,
  M47.49=s.47*P47.49,
  
  M2.50=(s.2*P2.5)+(s.2*P2.6),
  M3.50=(s.3*P3.5)+(s.3*P3.6),
  M4.50=(s.4*P4.6),
  M11.50=(s.11*P11.19)+(s.11*P11.20),
  M12.50=(s.12*P12.19)+(s.12*P12.20),
  M13.50=(s.13*P13.20),
  M17.50=(s.17*P17.19)+(s.17*P17.20),
  M18.50=(s.18*P18.20),
  M24.50=(s.24*P24.32)+(s.24*P24.33),
  M25.50=(s.25*P25.32)+(s.25*P25.33),
  M26.50=(s.26*P26.33),
  M30.50=(s.30*P30.32)+(s.30*P30.33),
  M31.50=(s.31*P31.33),
  M37.50=(s.37*P37.45)+(s.37*P37.46),
  M38.50=(s.38*P38.45)+(s.38*P38.46),
  M39.50=(s.39*P39.46),
  M43.50=(s.43*P43.45)+(s.43*P43.46),
  M44.50=(s.44*P44.46),
  
  
  M2.51=(s.2*P2.7),
  M5.51=(s.5*P5.7),
  M6.51=(s.6*P6.7),
  M11.51=(s.11*P11.21),
  M12.51=(s.12*P12.21),
  M13.51=(s.13*P13.21),
  M14.51=(s.14*P14.21),
  M15.51=(s.15*P15.21),
  M19.51=(s.19*P19.21),
  M20.51=(s.20*P20.21),
  M24.51=(s.24*P24.34),
  M25.51=(s.25*P25.34),
  M26.51=(s.26*P26.34),
  M27.51=(s.27*P27.34),
  M28.51=(s.28*P28.34),
  M32.51=(s.32*P32.34),
  M33.51=(s.33*P33.34),
  M37.51=(s.37*P37.47),
  M38.51=(s.38*P38.47),
  M39.51=(s.39*P39.47),
  M40.51=(s.40*P40.47),
  M41.51=(s.41*P41.47),
  M45.51=(s.45*P45.47),
  M46.51=(s.46*P46.47)
)]

#saves the moves dataset so that you can review the numbers moving between states each cycle (year) 
#moves_trace[[paste0(y)]]<-copy(moves)


#removes transition probabilities (these will be merged in again in the next cycle)
moves[,`:=`(
  P1.1=NULL,
  P1.2=NULL,
  P1.8=NULL,
  P2.3=NULL,
  P2.4=NULL,
  P2.5=NULL,
  P2.6=NULL,
  P2.7=NULL,
  P2.9=NULL,
  P2.10=NULL,
  P2.11=NULL,
  P3.3=NULL,
  P3.4=NULL,
  P3.5=NULL,
  P3.6=NULL,
  P3.9=NULL,
  P3.12=NULL,
  P4.4=NULL,
  P4.6=NULL,
  P4.9=NULL,
  P4.13=NULL,
  P5.5=NULL,
  P5.6=NULL,
  P5.7=NULL,
  P5.9=NULL,
  P5.14=NULL,
  P6.6=NULL,
  P6.7=NULL,
  P6.9=NULL,
  P6.15=NULL,
  P7.7=NULL,
  P7.9=NULL,
  P7.16=NULL,
  P11.17=NULL,
  P11.18=NULL,
  P11.19=NULL,
  P11.20=NULL,
  P11.21=NULL,
  P11.22=NULL,
  P11.23=NULL,
  P12.17=NULL,
  P12.18=NULL,
  P12.19=NULL,
  P12.20=NULL,
  P12.21=NULL,
  P12.22=NULL,
  P12.23=NULL,
  P13.18=NULL,
  P13.20=NULL,
  P13.21=NULL,
  P13.22=NULL,
  P13.23=NULL,
  P14.19=NULL,
  P14.20=NULL,
  P14.21=NULL,
  P14.22=NULL,
  P14.23=NULL,
  P15.20=NULL,
  P15.21=NULL,
  P15.22=NULL,
  P15.23=NULL,
  P16.21=NULL,
  P16.22=NULL,
  P16.23=NULL,
  P17.17=NULL,
  P17.18=NULL,
  P17.19=NULL,
  P17.20=NULL,
  P17.23=NULL,
  P18.18=NULL,
  P18.20=NULL,
  P18.23=NULL,
  P19.19=NULL,
  P19.20=NULL,
  P19.21=NULL,
  P19.23=NULL,
  P20.20=NULL,
  P20.21=NULL,
  P20.23=NULL,
  P21.21=NULL,
  P21.23=NULL,
  P11.24=NULL,
  P12.25=NULL,
  P12.26=NULL,
  P12.27=NULL,
  P12.28=NULL,
  P12.29=NULL,
  P13.26=NULL,
  P13.28=NULL,
  P13.29=NULL,
  P14.27=NULL,
  P14.28=NULL,
  P14.29=NULL,
  P15.28=NULL,
  P15.29=NULL,
  P16.29=NULL,
  P17.25=NULL,
  P18.26=NULL,
  P19.27=NULL,
  P20.28=NULL,
  P21.29=NULL,
  P24.30=NULL,
  P24.31=NULL,
  P24.32=NULL,
  P24.33=NULL,
  P24.34=NULL,
  P24.35=NULL,
  P24.36=NULL,
  P25.30=NULL,
  P25.31=NULL,
  P25.32=NULL,
  P25.33=NULL,
  P25.34=NULL,
  P25.35=NULL,
  P25.36=NULL,
  P26.31=NULL,
  P26.33=NULL,
  P26.34=NULL,
  P26.35=NULL,
  P26.36=NULL,
  P27.32=NULL,
  P27.33=NULL,
  P27.34=NULL,
  P27.35=NULL,
  P27.36=NULL,
  P28.33=NULL,
  P28.34=NULL,
  P28.35=NULL,
  P28.36=NULL,
  P29.34=NULL,
  P29.35=NULL,
  P29.36=NULL,
  P30.30=NULL,
  P30.31=NULL,
  P30.32=NULL,
  P30.33=NULL,
  P30.36=NULL,
  P31.31=NULL,
  P31.33=NULL,
  P31.36=NULL,
  P32.32=NULL,
  P32.33=NULL,
  P32.34=NULL,
  P32.36=NULL,
  P33.33=NULL,
  P33.34=NULL,
  P33.36=NULL,
  P34.34=NULL,
  P34.36=NULL,
  P24.37=NULL,
  P25.38=NULL,
  P25.39=NULL,
  P25.40=NULL,
  P25.41=NULL,
  P25.42=NULL,
  P26.39=NULL,
  P26.41=NULL,
  P26.42=NULL,
  P27.40=NULL,
  P27.41=NULL,
  P27.42=NULL,
  P28.41=NULL,
  P28.42=NULL,
  P29.42=NULL,
  P30.38=NULL,
  P31.39=NULL,
  P32.40=NULL,
  P33.41=NULL,
  P34.42=NULL,
  P37.43=NULL,
  P37.44=NULL,
  P37.45=NULL,
  P37.46=NULL,
  P37.47=NULL,
  P37.48=NULL,
  P37.49=NULL,
  P38.43=NULL,
  P38.44=NULL,
  P38.45=NULL,
  P38.46=NULL,
  P38.47=NULL,
  P38.48=NULL,
  P38.49=NULL,
  P39.44=NULL,
  P39.46=NULL,
  P39.47=NULL,
  P39.48=NULL,
  P39.49=NULL,
  P40.45=NULL,
  P40.46=NULL,
  P40.47=NULL,
  P40.48=NULL,
  P40.49=NULL,
  P41.46=NULL,
  P41.47=NULL,
  P41.48=NULL,
  P41.49=NULL,
  P42.47=NULL,
  P42.48=NULL,
  P42.49=NULL,
  P43.43=NULL,
  P43.44=NULL,
  P43.45=NULL,
  P43.46=NULL,
  P43.49=NULL,
  P44.44=NULL,
  P44.46=NULL,
  P44.49=NULL,
  P45.45=NULL,
  P45.46=NULL,
  P45.47=NULL,
  P45.49=NULL,
  P46.46=NULL,
  P46.47=NULL,
  P46.49=NULL,
  P47.47=NULL,
  P47.49=NULL
  
)]

#calculates health states for the cycle by summing the number of people moving to or remaining in that state 

moves[,`:=`(
  s.1=M1.1,
  s.2=M1.2,
  s.3=M2.3+M3.3,
  s.4=M2.4+M3.4+M4.4,
  s.5=M2.5+M3.5+M5.5,
  s.6=M2.6+M3.6+M4.6+M5.6+M6.6,
  s.7=M2.7+M5.7+M6.7+M7.7,
  s.8=M1.8,
  s.9=M2.9+M3.9+M4.9+M5.9+M6.9+M7.9,
  s.10=M2.10,
  s.11=M2.11,
  s.12=M3.12,
  s.13=M4.13,
  s.14=M5.14,
  s.15=M6.15,
  s.16=M7.16,
  s.17=M11.17+M12.17+M17.17,
  s.18=M11.18+M12.18+M13.18+M17.18+M18.18,
  s.19=M11.19+M12.19+M14.19+M17.19+M19.19,
  s.20=M11.20+M12.20+M13.20+M14.20+M15.20+M17.20+M18.20+M19.20+M20.20,
  s.21=M11.21+M12.21+M13.21+M14.21+M15.21+M16.21+M19.21+M20.21+M21.21,
  s.22=M11.22+M12.22+M13.22+M14.22+M15.22+M16.22,
  s.23=M11.23+M12.23+M13.23+M14.23+M15.23+M16.23+M17.23+M18.23+M19.23+M20.23+M21.23,
  s.24=M11.24,
  s.25=M12.25+M17.25,
  s.26=M12.26+M13.26+M18.26,
  s.27=M12.27+M14.27+M19.27,
  s.28=M12.28+M13.28+M14.28+M15.28+M20.28,
  s.29=M12.29+M13.29+M14.29+M15.29+M16.29+M21.29,
  s.30=M24.30+M25.30+M30.30,
  s.31=M24.31+M25.31+M26.31+M30.31+M31.31,
  s.32=M24.32+M25.32+M27.32+M30.32+M32.32,
  s.33=M24.33+M25.33+M26.33+M27.33+M28.33+M30.33+M31.33+M32.33+M33.33,
  s.34=M24.34+M25.34+M26.34+M27.34+M28.34+M29.34+M32.34+M33.34+M34.34,
  s.35=M24.35+M25.35+M26.35+M27.35+M28.35+M29.35,
  s.36=M24.36+M25.36+M26.36+M27.36+M28.36+M29.36+M30.36+M31.36+M32.36+M33.36+M34.36,
  s.37=M24.37,
  s.38=M25.38+M30.38,
  s.39=M25.39+M26.39+M31.39,
  s.40=M25.40+M27.40+M32.40,
  s.41=M25.41+M26.41+M27.41+M28.41+M33.41,
  s.42=M25.42+M26.42+M27.42+M28.42+M29.42+M34.42,
  s.43=M37.43+M38.43+M43.43,
  s.44=M37.44+M38.44+M39.44+M43.44+M44.44,
  s.45=M37.45+M38.45+M40.45+M43.45+M45.45,
  s.46=M37.46+M38.46+M39.46+M40.46+M41.46+M43.46+M44.46+M45.46+M46.46,
  s.47=M37.47+M38.47+M39.47+M40.47+M41.47+M42.47+M45.47+M46.47+M47.47,
  s.48=M37.48+M38.48+M39.48+M40.48+M41.48+M42.48,
  s.49=M37.49+M38.49+M39.49+M40.49+M41.49+M42.49+M43.49+M44.49+M45.49+M46.49+M47.49,
  s.50=M2.50+M3.50+M4.50+M11.50+M12.50+M13.50+M17.50+M18.50
        +M24.50+M25.50+M26.50+M30.50+M31.50+M37.50+M38.50+M39.50+M43.50+M44.50,
  s.51=M2.51+M5.51+M6.51+M11.51+M12.51+M13.51+M14.51+M15.51+M19.51+M20.51+M24.51+M25.51
        +M26.51+M27.51+M28.51+M32.51+M33.51+M37.51+M38.51+M39.51+M40.51+M41.51+M45.51
        +M46.51
)]  


#changes highest age in the model (100) to age 40, to allow addition of new 40 year olds into the model 
moves[Age==100,Age:=40]
setkey(moves, iter, Age, Year)

#selects value of year for first iteration, age 40 (based on order of setkey)
y<-moves[.(1,40)]$Year


#reads in the projections for age 40 for the relevant year, from the popproj data table
#newpop<-popproj[Year==y]$Pop
#sp3<-prev[Age==40 & iter==1]$s.3
#sp4<-prev[Age==40 & iter==1]$s.4
#sp5<-prev[Age==40 & iter==1]$s.5
#sp6<-prev[Age==40 & iter==1]$s.6
#sp7<-prev[Age==40 & iter==1]$s.7

####new 40 year olds not included in model####

moves[Age==40,`:=`(
  Pop=0,
  s.1=0,
  s.2=0,
  s.3=0,
  s.4=0,
  s.5=0,
  s.6=0,
  s.7=0,
  s.8=0,
  s.9=0,
  s.10=0,
  s.11=0,
  s.12=0,
  s.13=0,
  s.14=0,
  s.15=0,
  s.16=0,
  s.17=0,
  s.18=0,
  s.19=0,
  s.20=0,
  s.21=0,
  s.22=0,
  s.23=0,
  s.24=0,
  s.25=0,
  s.26=0,
  s.27=0,
  s.28=0,
  s.29=0,
  s.30=0,
  s.31=0,
  s.32=0,
  s.33=0,
  s.34=0,
  s.35=0,
  s.36=0,
  s.37=0,
  s.38=0,
  s.39=0,
  s.40=0,
  s.41=0,
  s.42=0,
  s.43=0,
  s.44=0,
  s.45=0,
  s.46=0,
  s.47=0,
  s.48=0,
  s.49=0,
  s.50=0,
  s.51=0,
  M1.1=0,
  M1.2=0,
  M1.8=0,
  M2.3=0,
  M2.4=0,
  M2.5=0,
  M2.6=0,
  M2.7=0,
  M2.9=0,
  M2.10=0,
  M2.11=0,
  M3.3=0,
  M3.4=0,
  M3.5=0,
  M3.6=0,
  M3.9=0,
  M3.12=0,
  M4.4=0,
  M4.6=0,
  M4.9=0,
  M4.13=0,
  M5.5=0,
  M5.6=0,
  M5.7=0,
  M5.9=0,
  M5.14=0,
  M6.6=0,
  M6.7=0,
  M6.9=0,
  M6.15=0,
  M7.7=0,
  M7.9=0,
  M7.16=0,
  M11.17=0,
  M11.18=0,
  M11.19=0,
  M11.20=0,
  M11.21=0,
  M11.22=0,
  M11.23=0,
  
  M12.17=0,
  M12.18=0,
  M12.19=0,
  M12.20=0,
  M12.21=0,
  M12.22=0,
  M12.23=0,
  
  M13.18=0,
  M13.20=0,
  M13.21=0,
  M13.22=0,
  M13.23=0,
  
  M14.19=0,
  M14.20=0,
  M14.21=0,
  M14.22=0,
  M14.23=0,
  
  M15.20=0,
  M15.21=0,
  M15.22=0,
  M15.23=0,
  
  M16.21=0,
  M16.22=0,
  M16.23=0,
  
  M17.17=0,
  M17.18=0,
  M17.19=0,
  M17.20=0,
  M17.23=0,
  M18.18=0,
  M18.20=0,
  M18.23=0,
  M19.19=0,
  M19.20=0,
  M19.21=0,
  M19.23=0,
  M20.20=0,
  M20.21=0,
  M20.23=0,
  M21.21=0,
  M21.23=0,
  M11.24=0,
  M12.25=0,
  M12.26=0,
  M12.27=0,
  M12.28=0,
  M12.29=0,
  M13.26=0,
  M13.28=0,
  M13.29=0,
  M14.27=0,
  M14.28=0,
  M14.29=0,
  M15.28=0,
  M15.29=0,
  M16.29=0,
  M17.25=0,
  M18.26=0,
  M19.27=0,
  M20.28=0,
  M21.29=0,
  M24.30=0,
  M24.31=0,
  M24.32=0,
  M24.33=0,
  M24.34=0,
  M24.35=0,
  M24.36=0,
  M25.30=0,
  M25.31=0,
  M25.32=0,
  M25.33=0,
  M25.34=0,
  M25.35=0,
  M25.36=0,
  M26.31=0,
  M26.33=0,
  M26.34=0,
  M26.35=0,
  M26.36=0,
  M27.32=0,
  M27.33=0,
  M27.34=0,
  M27.35=0,
  M27.36=0,
  M28.33=0,
  M28.34=0,
  M28.35=0,
  M28.36=0,
  M29.34=0,
  M29.35=0,
  M29.36=0,
  M30.30=0,
  M30.31=0,
  M30.32=0,
  M30.33=0,
  M30.36=0,
  M31.31=0,
  M31.33=0,
  M31.36=0,
  M32.32=0,
  M32.33=0,
  M32.34=0,
  M32.36=0,
  M33.33=0,
  M33.34=0,
  M33.36=0,
  M34.34=0,
  M34.36=0,
  M24.37=0,
  M25.38=0,
  M25.39=0,
  M25.40=0,
  M25.41=0,
  M25.42=0,
  M26.39=0,
  M26.41=0,
  M26.42=0,
  M27.40=0,
  M27.41=0,
  M27.42=0,
  M28.41=0,
  M28.42=0,
  M29.42=0,
  M30.38=0,
  M31.39=0,
  M32.40=0,
  M33.41=0,
  M34.42=0,
  M37.43=0,
  M37.44=0,
  M37.45=0,
  M37.46=0,
  M37.47=0,
  M37.48=0,
  M37.49=0,
  M38.43=0,
  M38.44=0,
  M38.45=0,
  M38.46=0,
  M38.47=0,
  M38.48=0,
  M38.49=0,
  M39.44=0,
  M39.46=0,
  M39.47=0,
  M39.48=0,
  M39.49=0,
  M40.45=0,
  M40.46=0,
  M40.47=0,
  M40.48=0,
  M40.49=0,
  M41.46=0,
  M41.47=0,
  M41.48=0,
  M41.49=0,
  M42.47=0,
  M42.48=0,
  M42.49=0,
  M43.43=0,
  M43.44=0,
  M43.45=0,
  M43.46=0,
  M43.49=0,
  M44.44=0,
  M44.46=0,
  M44.49=0,
  M45.45=0,
  M45.46=0,
  M45.47=0,
  M45.49=0,
  M46.46=0,
  M46.47=0,
  M46.49=0,
  M47.47=0,
  M47.49=0,
  M2.50=0,
  M3.50=0,
  M4.50=0,
  M11.50=0,
  M12.50=0,
  M13.50=0,
  M17.50=0,
  M18.50=0,
  M24.50=0,
  M25.50=0,
  M26.50=0,
  M30.50=0,
  M31.50=0,
  M37.50=0,
  M38.50=0,
  M39.50=0,
  M43.50=0,
  M44.50=0,
  M2.51=0,
  M5.51=0,
  M6.51=0,
  M11.51=0,
  M12.51=0,
  M13.51=0,
  M14.51=0,
  M15.51=0,
  M19.51=0,
  M20.51=0,
  M24.51=0,
  M25.51=0,
  M26.51=0,
  M27.51=0,
  M28.51=0,
  M32.51=0,
  M33.51=0,
  M37.51=0,
  M38.51=0,
  M39.51=0,
  M40.51=0,
  M41.51=0,
  M45.51=0,
  M46.51=0
)]

#calculates disease-free health state for 40 year old pop
#moves[Age==40,`:=`(s.1=newpop-s.3-s.4-s.5-s.6-s.7)]


#merges in transition probabilities
moves<-merge(moves, tps_S, by=c("Age","iter"))

#outputs the results for the health states for this model cycle (year)
res[[paste0(y)]]<-copy(moves)


#takes out the M vars so that these can be re-calculated in the next cycle
moves[,`:=`(M1.1=NULL,
            M1.2=NULL,
            M1.8=NULL,
            M2.3=NULL,
            M2.4=NULL,
            M2.5=NULL,
            M2.6=NULL,
            M2.7=NULL,
            M2.9=NULL,
            M2.10=NULL,
            M2.11=NULL,
            M3.3=NULL,
            M3.4=NULL,
            M3.5=NULL,
            M3.6=NULL,
            M3.9=NULL,
            M3.12=NULL,
            M4.4=NULL,
            M4.6=NULL,
            M4.9=NULL,
            M4.13=NULL,
            M5.5=NULL,
            M5.6=NULL,
            M5.7=NULL,
            M5.9=NULL,
            M5.14=NULL,
            M6.6=NULL,
            M6.7=NULL,
            M6.9=NULL,
            M6.15=NULL,
            M7.7=NULL,
            M7.9=NULL,
            M7.16=NULL,
            M11.17=NULL,
            M11.18=NULL,
            M11.19=NULL,
            M11.20=NULL,
            M11.21=NULL,
            M11.22=NULL,
            M11.23=NULL,
            
            M12.17=NULL,
            M12.18=NULL,
            M12.19=NULL,
            M12.20=NULL,
            M12.21=NULL,
            M12.22=NULL,
            M12.23=NULL,
            
            M13.18=NULL,
            M13.20=NULL,
            M13.21=NULL,
            M13.22=NULL,
            M13.23=NULL,
            
            M14.19=NULL,
            M14.20=NULL,
            M14.21=NULL,
            M14.22=NULL,
            M14.23=NULL,
            
            M15.20=NULL,
            M15.21=NULL,
            M15.22=NULL,
            M15.23=NULL,
            
            M16.21=NULL,
            M16.22=NULL,
            M16.23=NULL,
            
            M17.17=NULL,
            M17.18=NULL,
            M17.19=NULL,
            M17.20=NULL,
            M17.23=NULL,
            M18.18=NULL,
            M18.20=NULL,
            M18.23=NULL,
            M19.19=NULL,
            M19.20=NULL,
            M19.21=NULL,
            M19.23=NULL,
            M20.20=NULL,
            M20.21=NULL,
            M20.23=NULL,
            M21.21=NULL,
            M21.23=NULL,
            M11.24=NULL,
            M12.25=NULL,
            M12.26=NULL,
            M12.27=NULL,
            M12.28=NULL,
            M12.29=NULL,
            M13.26=NULL,
            M13.28=NULL,
            M13.29=NULL,
            M14.27=NULL,
            M14.28=NULL,
            M14.29=NULL,
            M15.28=NULL,
            M15.29=NULL,
            M16.29=NULL,
            M17.25=NULL,
            M18.26=NULL,
            M19.27=NULL,
            M20.28=NULL,
            M21.29=NULL,
            M24.30=NULL,
            M24.31=NULL,
            M24.32=NULL,
            M24.33=NULL,
            M24.34=NULL,
            M24.35=NULL,
            M24.36=NULL,
            M25.30=NULL,
            M25.31=NULL,
            M25.32=NULL,
            M25.33=NULL,
            M25.34=NULL,
            M25.35=NULL,
            M25.36=NULL,
            M26.31=NULL,
            M26.33=NULL,
            M26.34=NULL,
            M26.35=NULL,
            M26.36=NULL,
            M27.32=NULL,
            M27.33=NULL,
            M27.34=NULL,
            M27.35=NULL,
            M27.36=NULL,
            M28.33=NULL,
            M28.34=NULL,
            M28.35=NULL,
            M28.36=NULL,
            M29.34=NULL,
            M29.35=NULL,
            M29.36=NULL,
            M30.30=NULL,
            M30.31=NULL,
            M30.32=NULL,
            M30.33=NULL,
            M30.36=NULL,
            M31.31=NULL,
            M31.33=NULL,
            M31.36=NULL,
            M32.32=NULL,
            M32.33=NULL,
            M32.34=NULL,
            M32.36=NULL,
            M33.33=NULL,
            M33.34=NULL,
            M33.36=NULL,
            M34.34=NULL,
            M34.36=NULL,
            M24.37=NULL,
            M25.38=NULL,
            M25.39=NULL,
            M25.40=NULL,
            M25.41=NULL,
            M25.42=NULL,
            M26.39=NULL,
            M26.41=NULL,
            M26.42=NULL,
            M27.40=NULL,
            M27.41=NULL,
            M27.42=NULL,
            M28.41=NULL,
            M28.42=NULL,
            M29.42=NULL,
            M30.38=NULL,
            M31.39=NULL,
            M32.40=NULL,
            M33.41=NULL,
            M34.42=NULL,
            M37.43=NULL,
            M37.44=NULL,
            M37.45=NULL,
            M37.46=NULL,
            M37.47=NULL,
            M37.48=NULL,
            M37.49=NULL,
            M38.43=NULL,
            M38.44=NULL,
            M38.45=NULL,
            M38.46=NULL,
            M38.47=NULL,
            M38.48=NULL,
            M38.49=NULL,
            M39.44=NULL,
            M39.46=NULL,
            M39.47=NULL,
            M39.48=NULL,
            M39.49=NULL,
            M40.45=NULL,
            M40.46=NULL,
            M40.47=NULL,
            M40.48=NULL,
            M40.49=NULL,
            M41.46=NULL,
            M41.47=NULL,
            M41.48=NULL,
            M41.49=NULL,
            M42.47=NULL,
            M42.48=NULL,
            M42.49=NULL,
            M43.43=NULL,
            M43.44=NULL,
            M43.45=NULL,
            M43.46=NULL,
            M43.49=NULL,
            M44.44=NULL,
            M44.46=NULL,
            M44.49=NULL,
            M45.45=NULL,
            M45.46=NULL,
            M45.47=NULL,
            M45.49=NULL,
            M46.46=NULL,
            M46.47=NULL,
            M46.49=NULL,
            M47.47=NULL,
            M47.49=NULL,
            M2.50=NULL,
            M3.50=NULL,
            M4.50=NULL,
            M11.50=NULL,
            M12.50=NULL,
            M13.50=NULL,
            M17.50=NULL,
            M18.50=NULL,
            M24.50=NULL,
            M25.50=NULL,
            M26.50=NULL,
            M30.50=NULL,
            M31.50=NULL,
            M37.50=NULL,
            M38.50=NULL,
            M39.50=NULL,
            M43.50=NULL,
            M44.50=NULL,
            M2.51=NULL,
            M5.51=NULL,
            M6.51=NULL,
            M11.51=NULL,
            M12.51=NULL,
            M13.51=NULL,
            M14.51=NULL,
            M15.51=NULL,
            M19.51=NULL,
            M20.51=NULL,
            M24.51=NULL,
            M25.51=NULL,
            M26.51=NULL,
            M27.51=NULL,
            M28.51=NULL,
            M32.51=NULL,
            M33.51=NULL,
            M37.51=NULL,
            M38.51=NULL,
            M39.51=NULL,
            M40.51=NULL,
            M41.51=NULL,
            M45.51=NULL,
            M46.51=NULL
            
            
)]

}  

##END OF LOOP

#binds all of the results for each year together
res<-rbindlist(res, use.names=T, fill=T)

#same for the move traces table
#moves_trace<-rbindlist(moves_trace, use.names=T, fill=T)



#takes out the move values - only need health states
res[,`:=`(M1.1=NULL,
          M1.2=NULL,
          M1.8=NULL,
          M2.3=NULL,
          M2.4=NULL,
          M2.5=NULL,
          M2.6=NULL,
          M2.7=NULL,
          M2.9=NULL,
          M2.10=NULL,
          M2.11=NULL,
          M3.3=NULL,
          M3.4=NULL,
          M3.5=NULL,
          M3.6=NULL,
          M3.9=NULL,
          M3.12=NULL,
          M4.4=NULL,
          M4.6=NULL,
          M4.9=NULL,
          M4.13=NULL,
          M5.5=NULL,
          M5.6=NULL,
          M5.7=NULL,
          M5.9=NULL,
          M5.14=NULL,
          M6.6=NULL,
          M6.7=NULL,
          M6.9=NULL,
          M6.15=NULL,
          M7.7=NULL,
          M7.9=NULL,
          M7.16=NULL,
          M11.17=NULL,
          M11.18=NULL,
          M11.19=NULL,
          M11.20=NULL,
          M11.21=NULL,
          M11.22=NULL,
          M11.23=NULL,
          
          M12.17=NULL,
          M12.18=NULL,
          M12.19=NULL,
          M12.20=NULL,
          M12.21=NULL,
          M12.22=NULL,
          M12.23=NULL,
          
          M13.18=NULL,
          M13.20=NULL,
          M13.21=NULL,
          M13.22=NULL,
          M13.23=NULL,
          
          M14.19=NULL,
          M14.20=NULL,
          M14.21=NULL,
          M14.22=NULL,
          M14.23=NULL,
          
          M15.20=NULL,
          M15.21=NULL,
          M15.22=NULL,
          M15.23=NULL,
          
          M16.21=NULL,
          M16.22=NULL,
          M16.23=NULL,
          
          M17.17=NULL,
          M17.18=NULL,
          M17.19=NULL,
          M17.20=NULL,
          M17.23=NULL,
          M18.18=NULL,
          M18.20=NULL,
          M18.23=NULL,
          M19.19=NULL,
          M19.20=NULL,
          M19.21=NULL,
          M19.23=NULL,
          M20.20=NULL,
          M20.21=NULL,
          M20.23=NULL,
          M21.21=NULL,
          M21.23=NULL,
          M11.24=NULL,
          M12.25=NULL,
          M12.26=NULL,
          M12.27=NULL,
          M12.28=NULL,
          M12.29=NULL,
          M13.26=NULL,
          M13.28=NULL,
          M13.29=NULL,
          M14.27=NULL,
          M14.28=NULL,
          M14.29=NULL,
          M15.28=NULL,
          M15.29=NULL,
          M16.29=NULL,
          M17.25=NULL,
          M18.26=NULL,
          M19.27=NULL,
          M20.28=NULL,
          M21.29=NULL,
          M24.30=NULL,
          M24.31=NULL,
          M24.32=NULL,
          M24.33=NULL,
          M24.34=NULL,
          M24.35=NULL,
          M24.36=NULL,
          M25.30=NULL,
          M25.31=NULL,
          M25.32=NULL,
          M25.33=NULL,
          M25.34=NULL,
          M25.35=NULL,
          M25.36=NULL,
          M26.31=NULL,
          M26.33=NULL,
          M26.34=NULL,
          M26.35=NULL,
          M26.36=NULL,
          M27.32=NULL,
          M27.33=NULL,
          M27.34=NULL,
          M27.35=NULL,
          M27.36=NULL,
          M28.33=NULL,
          M28.34=NULL,
          M28.35=NULL,
          M28.36=NULL,
          M29.34=NULL,
          M29.35=NULL,
          M29.36=NULL,
          M30.30=NULL,
          M30.31=NULL,
          M30.32=NULL,
          M30.33=NULL,
          M30.36=NULL,
          M31.31=NULL,
          M31.33=NULL,
          M31.36=NULL,
          M32.32=NULL,
          M32.33=NULL,
          M32.34=NULL,
          M32.36=NULL,
          M33.33=NULL,
          M33.34=NULL,
          M33.36=NULL,
          M34.34=NULL,
          M34.36=NULL,
          M24.37=NULL,
          M25.38=NULL,
          M25.39=NULL,
          M25.40=NULL,
          M25.41=NULL,
          M25.42=NULL,
          M26.39=NULL,
          M26.41=NULL,
          M26.42=NULL,
          M27.40=NULL,
          M27.41=NULL,
          M27.42=NULL,
          M28.41=NULL,
          M28.42=NULL,
          M29.42=NULL,
          M30.38=NULL,
          M31.39=NULL,
          M32.40=NULL,
          M33.41=NULL,
          M34.42=NULL,
          M37.43=NULL,
          M37.44=NULL,
          M37.45=NULL,
          M37.46=NULL,
          M37.47=NULL,
          M37.48=NULL,
          M37.49=NULL,
          M38.43=NULL,
          M38.44=NULL,
          M38.45=NULL,
          M38.46=NULL,
          M38.47=NULL,
          M38.48=NULL,
          M38.49=NULL,
          M39.44=NULL,
          M39.46=NULL,
          M39.47=NULL,
          M39.48=NULL,
          M39.49=NULL,
          M40.45=NULL,
          M40.46=NULL,
          M40.47=NULL,
          M40.48=NULL,
          M40.49=NULL,
          M41.46=NULL,
          M41.47=NULL,
          M41.48=NULL,
          M41.49=NULL,
          M42.47=NULL,
          M42.48=NULL,
          M42.49=NULL,
          M43.43=NULL,
          M43.44=NULL,
          M43.45=NULL,
          M43.46=NULL,
          M43.49=NULL,
          M44.44=NULL,
          M44.46=NULL,
          M44.49=NULL,
          M45.45=NULL,
          M45.46=NULL,
          M45.47=NULL,
          M45.49=NULL,
          M46.46=NULL,
          M46.47=NULL,
          M46.49=NULL,
          M47.47=NULL,
          M47.49=NULL,
          M2.50=NULL,
          M3.50=NULL,
          M4.50=NULL,
          M11.50=NULL,
          M12.50=NULL,
          M13.50=NULL,
          M17.50=NULL,
          M18.50=NULL,
          M24.50=NULL,
          M25.50=NULL,
          M26.50=NULL,
          M30.50=NULL,
          M31.50=NULL,
          M37.50=NULL,
          M38.50=NULL,
          M39.50=NULL,
          M43.50=NULL,
          M44.50=NULL,
          M2.51=NULL,
          M5.51=NULL,
          M6.51=NULL,
          M11.51=NULL,
          M12.51=NULL,
          M13.51=NULL,
          M14.51=NULL,
          M15.51=NULL,
          M19.51=NULL,
          M20.51=NULL,
          M24.51=NULL,
          M25.51=NULL,
          M26.51=NULL,
          M27.51=NULL,
          M28.51=NULL,
          M32.51=NULL,
          M33.51=NULL,
          M37.51=NULL,
          M38.51=NULL,
          M39.51=NULL,
          M40.51=NULL,
          M41.51=NULL,
          M45.51=NULL,
          M46.51=NULL
          
)]

####SECOND MODEL CALCULATION#### 
##running series of models reduces computational load

#removes transition probabilities (these will be merged in again in the next cycle)

moves[,`:=`(
  P1.1=NULL,
  P1.2=NULL,
  P1.8=NULL,
  P2.3=NULL,
  P2.4=NULL,
  P2.5=NULL,
  P2.6=NULL,
  P2.7=NULL,
  P2.9=NULL,
  P2.10=NULL,
  P2.11=NULL,
  P3.3=NULL,
  P3.4=NULL,
  P3.5=NULL,
  P3.6=NULL,
  P3.9=NULL,
  P3.12=NULL,
  P4.4=NULL,
  P4.6=NULL,
  P4.9=NULL,
  P4.13=NULL,
  P5.5=NULL,
  P5.6=NULL,
  P5.7=NULL,
  P5.9=NULL,
  P5.14=NULL,
  P6.6=NULL,
  P6.7=NULL,
  P6.9=NULL,
  P6.15=NULL,
  P7.7=NULL,
  P7.9=NULL,
  P7.16=NULL,
  P11.17=NULL,
  P11.18=NULL,
  P11.19=NULL,
  P11.20=NULL,
  P11.21=NULL,
  P11.22=NULL,
  P11.23=NULL,
  P12.17=NULL,
  P12.18=NULL,
  P12.19=NULL,
  P12.20=NULL,
  P12.21=NULL,
  P12.22=NULL,
  P12.23=NULL,
  P13.18=NULL,
  P13.20=NULL,
  P13.21=NULL,
  P13.22=NULL,
  P13.23=NULL,
  P14.19=NULL,
  P14.20=NULL,
  P14.21=NULL,
  P14.22=NULL,
  P14.23=NULL,
  P15.20=NULL,
  P15.21=NULL,
  P15.22=NULL,
  P15.23=NULL,
  P16.21=NULL,
  P16.22=NULL,
  P16.23=NULL,
  P17.17=NULL,
  P17.18=NULL,
  P17.19=NULL,
  P17.20=NULL,
  P17.23=NULL,
  P18.18=NULL,
  P18.20=NULL,
  P18.23=NULL,
  P19.19=NULL,
  P19.20=NULL,
  P19.21=NULL,
  P19.23=NULL,
  P20.20=NULL,
  P20.21=NULL,
  P20.23=NULL,
  P21.21=NULL,
  P21.23=NULL,
  P11.24=NULL,
  P12.25=NULL,
  P12.26=NULL,
  P12.27=NULL,
  P12.28=NULL,
  P12.29=NULL,
  P13.26=NULL,
  P13.28=NULL,
  P13.29=NULL,
  P14.27=NULL,
  P14.28=NULL,
  P14.29=NULL,
  P15.28=NULL,
  P15.29=NULL,
  P16.29=NULL,
  P17.25=NULL,
  P18.26=NULL,
  P19.27=NULL,
  P20.28=NULL,
  P21.29=NULL,
  P24.30=NULL,
  P24.31=NULL,
  P24.32=NULL,
  P24.33=NULL,
  P24.34=NULL,
  P24.35=NULL,
  P24.36=NULL,
  P25.30=NULL,
  P25.31=NULL,
  P25.32=NULL,
  P25.33=NULL,
  P25.34=NULL,
  P25.35=NULL,
  P25.36=NULL,
  P26.31=NULL,
  P26.33=NULL,
  P26.34=NULL,
  P26.35=NULL,
  P26.36=NULL,
  P27.32=NULL,
  P27.33=NULL,
  P27.34=NULL,
  P27.35=NULL,
  P27.36=NULL,
  P28.33=NULL,
  P28.34=NULL,
  P28.35=NULL,
  P28.36=NULL,
  P29.34=NULL,
  P29.35=NULL,
  P29.36=NULL,
  P30.30=NULL,
  P30.31=NULL,
  P30.32=NULL,
  P30.33=NULL,
  P30.36=NULL,
  P31.31=NULL,
  P31.33=NULL,
  P31.36=NULL,
  P32.32=NULL,
  P32.33=NULL,
  P32.34=NULL,
  P32.36=NULL,
  P33.33=NULL,
  P33.34=NULL,
  P33.36=NULL,
  P34.34=NULL,
  P34.36=NULL,
  P24.37=NULL,
  P25.38=NULL,
  P25.39=NULL,
  P25.40=NULL,
  P25.41=NULL,
  P25.42=NULL,
  P26.39=NULL,
  P26.41=NULL,
  P26.42=NULL,
  P27.40=NULL,
  P27.41=NULL,
  P27.42=NULL,
  P28.41=NULL,
  P28.42=NULL,
  P29.42=NULL,
  P30.38=NULL,
  P31.39=NULL,
  P32.40=NULL,
  P33.41=NULL,
  P34.42=NULL,
  P37.43=NULL,
  P37.44=NULL,
  P37.45=NULL,
  P37.46=NULL,
  P37.47=NULL,
  P37.48=NULL,
  P37.49=NULL,
  P38.43=NULL,
  P38.44=NULL,
  P38.45=NULL,
  P38.46=NULL,
  P38.47=NULL,
  P38.48=NULL,
  P38.49=NULL,
  P39.44=NULL,
  P39.46=NULL,
  P39.47=NULL,
  P39.48=NULL,
  P39.49=NULL,
  P40.45=NULL,
  P40.46=NULL,
  P40.47=NULL,
  P40.48=NULL,
  P40.49=NULL,
  P41.46=NULL,
  P41.47=NULL,
  P41.48=NULL,
  P41.49=NULL,
  P42.47=NULL,
  P42.48=NULL,
  P42.49=NULL,
  P43.43=NULL,
  P43.44=NULL,
  P43.45=NULL,
  P43.46=NULL,
  P43.49=NULL,
  P44.44=NULL,
  P44.46=NULL,
  P44.49=NULL,
  P45.45=NULL,
  P45.46=NULL,
  P45.47=NULL,
  P45.49=NULL,
  P46.46=NULL,
  P46.47=NULL,
  P46.49=NULL,
  P47.47=NULL,
  P47.49=NULL
  
)]

#creates new initial states file based on last cycle of model 1 
initial.states2<-merge(moves, tps_S, by=c("Age","iter"))

#sort by iteration
setkey(initial.states2,iter)

#saves out the results as second results file

res2<-list()
res2[[paste0(FIRST.YEAR2)]]<-copy(initial.states2)


# creates a new moves file for the second model 
moves2<-copy(initial.states2)

#Starts at first year 2 
y<-FIRST.YEAR2

##Loop starts here - current year is stored in y
##loop continues until y = LAST.YEAR (defined at beginning of code)

while (y<LAST.YEAR2) 
{
  
  
  #age and year are increased by 1 year to move model along by one year  
  
  moves2[,Age:=Age+1]
  moves2[,Year:=Year+1]
  
  #Taking out existing tprs - need to re-merge, as age not synced up
  
  moves2[,`:=`(
    P1.1=NULL,
    P1.2=NULL,
    P1.8=NULL,
    P2.3=NULL,
    P2.4=NULL,
    P2.5=NULL,
    P2.6=NULL,
    P2.7=NULL,
    P2.9=NULL,
    P2.10=NULL,
    P2.11=NULL,
    P3.3=NULL,
    P3.4=NULL,
    P3.5=NULL,
    P3.6=NULL,
    P3.9=NULL,
    P3.12=NULL,
    P4.4=NULL,
    P4.6=NULL,
    P4.9=NULL,
    P4.13=NULL,
    P5.5=NULL,
    P5.6=NULL,
    P5.7=NULL,
    P5.9=NULL,
    P5.14=NULL,
    P6.6=NULL,
    P6.7=NULL,
    P6.9=NULL,
    P6.15=NULL,
    P7.7=NULL,
    P7.9=NULL,
    P7.16=NULL,
    P11.17=NULL,
    P11.18=NULL,
    P11.19=NULL,
    P11.20=NULL,
    P11.21=NULL,
    P11.22=NULL,
    P11.23=NULL,
    P12.17=NULL,
    P12.18=NULL,
    P12.19=NULL,
    P12.20=NULL,
    P12.21=NULL,
    P12.22=NULL,
    P12.23=NULL,
    
    
    P13.18=NULL,
    P13.20=NULL,
    P13.21=NULL,
    P13.22=NULL,
    P13.23=NULL,
    
    P14.19=NULL,
    P14.20=NULL,
    P14.21=NULL,
    P14.22=NULL,
    P14.23=NULL,
    P15.20=NULL,
    P15.21=NULL,
    P15.22=NULL,
    P15.23=NULL,
    P16.21=NULL,
    P16.22=NULL,
    P16.23=NULL,
    P17.17=NULL,
    P17.18=NULL,
    P17.19=NULL,
    P17.20=NULL,
    P17.23=NULL,
    P18.18=NULL,
    P18.20=NULL,
    P18.23=NULL,
    P19.19=NULL,
    P19.20=NULL,
    P19.21=NULL,
    P19.23=NULL,
    P20.20=NULL,
    P20.21=NULL,
    P20.23=NULL,
    P21.21=NULL,
    P21.23=NULL,
    P11.24=NULL,
    P12.25=NULL,
    P12.26=NULL,
    P12.27=NULL,
    P12.28=NULL,
    P12.29=NULL,
    P13.26=NULL,
    P13.28=NULL,
    P13.29=NULL,
    P14.27=NULL,
    P14.28=NULL,
    P14.29=NULL,
    P15.28=NULL,
    P15.29=NULL,
    P16.29=NULL,
    P17.25=NULL,
    P18.26=NULL,
    P19.27=NULL,
    P20.28=NULL,
    P21.29=NULL,
    P24.30=NULL,
    P24.31=NULL,
    P24.32=NULL,
    P24.33=NULL,
    P24.34=NULL,
    P24.35=NULL,
    P24.36=NULL,
    P25.30=NULL,
    P25.31=NULL,
    P25.32=NULL,
    P25.33=NULL,
    P25.34=NULL,
    P25.35=NULL,
    P25.36=NULL,
    P26.31=NULL,
    P26.33=NULL,
    P26.34=NULL,
    P26.35=NULL,
    P26.36=NULL,
    P27.32=NULL,
    P27.33=NULL,
    P27.34=NULL,
    P27.35=NULL,
    P27.36=NULL,
    P28.33=NULL,
    P28.34=NULL,
    P28.35=NULL,
    P28.36=NULL,
    P29.34=NULL,
    P29.35=NULL,
    P29.36=NULL,
    P30.30=NULL,
    P30.31=NULL,
    P30.32=NULL,
    P30.33=NULL,
    P30.36=NULL,
    P31.31=NULL,
    P31.33=NULL,
    P31.36=NULL,
    P32.32=NULL,
    P32.33=NULL,
    P32.34=NULL,
    P32.36=NULL,
    P33.33=NULL,
    P33.34=NULL,
    P33.36=NULL,
    P34.34=NULL,
    P34.36=NULL,
    P24.37=NULL,
    P25.38=NULL,
    P25.39=NULL,
    P25.40=NULL,
    P25.41=NULL,
    P25.42=NULL,
    P26.39=NULL,
    P26.41=NULL,
    P26.42=NULL,
    P27.40=NULL,
    P27.41=NULL,
    P27.42=NULL,
    P28.41=NULL,
    P28.42=NULL,
    P29.42=NULL,
    P30.38=NULL,
    P31.39=NULL,
    P32.40=NULL,
    P33.41=NULL,
    P34.42=NULL,
    P37.43=NULL,
    P37.44=NULL,
    P37.45=NULL,
    P37.46=NULL,
    P37.47=NULL,
    P37.48=NULL,
    P37.49=NULL,
    P38.43=NULL,
    P38.44=NULL,
    P38.45=NULL,
    P38.46=NULL,
    P38.47=NULL,
    P38.48=NULL,
    P38.49=NULL,
    P39.44=NULL,
    P39.46=NULL,
    P39.47=NULL,
    P39.48=NULL,
    P39.49=NULL,
    P40.45=NULL,
    P40.46=NULL,
    P40.47=NULL,
    P40.48=NULL,
    P40.49=NULL,
    P41.46=NULL,
    P41.47=NULL,
    P41.48=NULL,
    P41.49=NULL,
    P42.47=NULL,
    P42.48=NULL,
    P42.49=NULL,
    P43.43=NULL,
    P43.44=NULL,
    P43.45=NULL,
    P43.46=NULL,
    P43.49=NULL,
    P44.44=NULL,
    P44.46=NULL,
    P44.49=NULL,
    P45.45=NULL,
    P45.46=NULL,
    P45.47=NULL,
    P45.49=NULL,
    P46.46=NULL,
    P46.47=NULL,
    P46.49=NULL,
    P47.47=NULL,
    P47.49=NULL
    
  )]
  
  #merging in tps again so that age is synced up correctly
  
  moves2<-merge(moves2, tps_S, by=c("Age","iter"))
  
  #M50 and M51 included to capture incident CIND and dementia
  #Health states are multipled by the relevant transition probabilites to calculate the numbers transitioning between states
  #the M var names indicate between which states the transition is occuring
  #e.g. M1.2 is the numbers moving between state 1 and 2, which is calculated by applying the tpr p1.2 to state 1 (s.1)
  
  moves2[,`:=`(
    M1.1=s.1*P1.1,
    M1.2=s.1*P1.2,
    M1.8=s.1*P1.8,
    M2.3=s.2*P2.3,
    M2.4=s.2*P2.4,
    M2.5=s.2*P2.5,
    M2.6=s.2*P2.6,
    M2.7=s.2*P2.7,
    M2.9=s.2*P2.9,
    M2.10=s.2*P2.10,
    M2.11=s.2*P2.11,
    M3.3=s.3*P3.3,
    M3.4=s.3*P3.4,
    M3.5=s.3*P3.5,
    M3.6=s.3*P3.6,
    M3.9=s.3*P3.9,
    M3.12=s.3*P3.12,
    M4.4=s.4*P4.4,
    M4.6=s.4*P4.6,
    M4.9=s.4*P4.9,
    M4.13=s.4*P4.13,
    M5.5=s.5*P5.5,
    M5.6=s.5*P5.6,
    M5.7=s.5*P5.7,
    M5.9=s.5*P5.9,
    M5.14=s.5*P5.14,
    M6.6=s.6*P6.6,
    M6.7=s.6*P6.7,
    M6.9=s.6*P6.9,
    M6.15=s.6*P6.15,
    M7.7=s.7*P7.7,
    M7.9=s.7*P7.9,
    M7.16=s.7*P7.16,
    M11.17=s.11*P11.17,
    M11.18=s.11*P11.18,
    M11.19=s.11*P11.19,
    M11.20=s.11*P11.20,
    M11.21=s.11*P11.21,
    M11.22=s.11*P11.22,
    M11.23=s.11*P11.23,
    M12.17=s.12*P12.17,
    M12.18=s.12*P12.18,
    M12.19=s.12*P12.19,
    M12.20=s.12*P12.20,
    M12.21=s.12*P12.21,
    M12.22=s.12*P12.22,
    M12.23=s.12*P12.23,
    M13.18=s.13*P13.18,
    M13.20=s.13*P13.20,
    M13.21=s.13*P13.21,
    M13.22=s.13*P13.22,
    M13.23=s.13*P13.23,
    
    M14.19=s.14*P14.19,
    M14.20=s.14*P14.20,
    M14.21=s.14*P14.21,
    M14.22=s.14*P14.22,
    M14.23=s.14*P14.23,
    
    M15.20=s.15*P15.20,
    M15.21=s.15*P15.21,
    M15.22=s.15*P15.22,
    M15.23=s.15*P15.23,
    
    M16.21=s.16*P16.21,
    M16.22=s.16*P16.22,
    M16.23=s.16*P16.23,
    
    M17.17=s.17*P17.17,
    M17.18=s.17*P17.18,
    M17.19=s.17*P17.19,
    M17.20=s.17*P17.20,
    M17.23=s.17*P17.23,
    M18.18=s.18*P18.18,
    M18.20=s.18*P18.20,
    M18.23=s.18*P18.23,
    M19.19=s.19*P19.19,
    M19.20=s.19*P19.20,
    M19.21=s.19*P19.21,
    M19.23=s.19*P19.23,
    M20.20=s.20*P20.20,
    M20.21=s.20*P20.21,
    M20.23=s.20*P20.23,
    M21.21=s.21*P21.21,
    M21.23=s.21*P21.23,
    
    
    M11.24=s.11*P11.24,
    M12.25=s.12*P12.25,
    M12.26=s.12*P12.26,
    M12.27=s.12*P12.27,
    M12.28=s.12*P12.28,
    M12.29=s.12*P12.29,
    M13.26=s.13*P13.26,
    M13.28=s.13*P13.28,
    M13.29=s.13*P13.29,
    M14.27=s.14*P14.27,
    M14.28=s.14*P14.28,
    M14.29=s.14*P14.29,
    M15.28=s.15*P15.28,
    M15.29=s.15*P15.29,
    M16.29=s.16*P16.29,
    M17.25=s.17*P17.25,
    
    M18.26=s.18*P18.26,
    M19.27=s.19*P19.27,
    M20.28=s.20*P20.28,
    M21.29=s.21*P21.29,
    M24.30=s.24*P24.30,
    M24.31=s.24*P24.31,
    M24.32=s.24*P24.32,
    M24.33=s.24*P24.33,
    M24.34=s.24*P24.34,
    M24.35=s.24*P24.35,
    M24.36=s.24*P24.36,
    M25.30=s.25*P25.30,
    M25.31=s.25*P25.31,
    M25.32=s.25*P25.32,
    M25.33=s.25*P25.33,
    M25.34=s.25*P25.34,
    M25.35=s.25*P25.35,
    M25.36=s.25*P25.36,
    M26.31=s.26*P26.31,
    M26.33=s.26*P26.33,
    M26.34=s.26*P26.34,
    M26.35=s.26*P26.35,
    M26.36=s.26*P26.36,
    M27.32=s.27*P27.32,
    M27.33=s.27*P27.33,
    M27.34=s.27*P27.34,
    M27.35=s.27*P27.35,
    M27.36=s.27*P27.36,
    M28.33=s.28*P28.33,
    M28.34=s.28*P28.34,
    M28.35=s.28*P28.35,
    M28.36=s.28*P28.36,
    M29.34=s.29*P29.34,
    M29.35=s.29*P29.35,
    M29.36=s.29*P29.36,
    M30.30=s.30*P30.30,
    M30.31=s.30*P30.31,
    M30.32=s.30*P30.32,
    M30.33=s.30*P30.33,
    M30.36=s.30*P30.36,
    M31.31=s.31*P31.31,
    M31.33=s.31*P31.33,
    M31.36=s.31*P31.36,
    M32.32=s.32*P32.32,
    M32.33=s.32*P32.33,
    M32.34=s.32*P32.34,
    M32.36=s.32*P32.36,
    M33.33=s.33*P33.33,
    M33.34=s.33*P33.34,
    M33.36=s.33*P33.36,
    M34.34=s.34*P34.34,
    M34.36=s.34*P34.36,
    M24.37=s.24*P24.37,
    M25.38=s.25*P25.38,
    M25.39=s.25*P25.39,
    M25.40=s.25*P25.40,
    M25.41=s.25*P25.41,
    M25.42=s.25*P25.42,
    M26.39=s.26*P26.39,
    M26.41=s.26*P26.41,
    M26.42=s.26*P26.42,
    M27.40=s.27*P27.40,
    M27.41=s.27*P27.41,
    M27.42=s.27*P27.42,
    M28.41=s.28*P28.41,
    M28.42=s.28*P28.42,
    M29.42=s.29*P29.42,
    M30.38=s.30*P30.38,
    M31.39=s.31*P31.39,
    M32.40=s.32*P32.40,
    M33.41=s.33*P33.41,
    M34.42=s.34*P34.42,
    M37.43=s.37*P37.43,
    M37.44=s.37*P37.44,
    M37.45=s.37*P37.45,
    M37.46=s.37*P37.46,
    M37.47=s.37*P37.47,
    M37.48=s.37*P37.48,
    M37.49=s.37*P37.49,
    M38.43=s.38*P38.43,
    M38.44=s.38*P38.44,
    M38.45=s.38*P38.45,
    M38.46=s.38*P38.46,
    M38.47=s.38*P38.47,
    M38.48=s.38*P38.48,
    M38.49=s.38*P38.49,
    M39.44=s.39*P39.44,
    M39.46=s.39*P39.46,
    M39.47=s.39*P39.47,
    M39.48=s.39*P39.48,
    M39.49=s.39*P39.49,
    M40.45=s.40*P40.45,
    M40.46=s.40*P40.46,
    M40.47=s.40*P40.47,
    M40.48=s.40*P40.48,
    M40.49=s.40*P40.49,
    M41.46=s.41*P41.46,
    M41.47=s.41*P41.47,
    M41.48=s.41*P41.48,
    M41.49=s.41*P41.49,
    M42.47=s.42*P42.47,
    M42.48=s.42*P42.48,
    M42.49=s.42*P42.49,
    M43.43=s.43*P43.43,
    M43.44=s.43*P43.44,
    M43.45=s.43*P43.45,
    M43.46=s.43*P43.46,
    M43.49=s.43*P43.49,
    M44.44=s.44*P44.44,
    M44.46=s.44*P44.46,
    M44.49=s.44*P44.49,
    M45.45=s.45*P45.45,
    M45.46=s.45*P45.46,
    M45.47=s.45*P45.47,
    M45.49=s.45*P45.49,
    M46.46=s.46*P46.46,
    M46.47=s.46*P46.47,
    M46.49=s.46*P46.49,
    M47.47=s.47*P47.47,
    M47.49=s.47*P47.49,
    
    M2.50=(s.2*P2.5)+(s.2*P2.6),
    M3.50=(s.3*P3.5)+(s.3*P3.6),
    M4.50=(s.4*P4.6),
    M11.50=(s.11*P11.19)+(s.11*P11.20),
    M12.50=(s.12*P12.19)+(s.12*P12.20),
    M13.50=(s.13*P13.20),
    M17.50=(s.17*P17.19)+(s.17*P17.20),
    M18.50=(s.18*P18.20),
    M24.50=(s.24*P24.32)+(s.24*P24.33),
    M25.50=(s.25*P25.32)+(s.25*P25.33),
    M26.50=(s.26*P26.33),
    M30.50=(s.30*P30.32)+(s.30*P30.33),
    M31.50=(s.31*P31.33),
    M37.50=(s.37*P37.45)+(s.37*P37.46),
    M38.50=(s.38*P38.45)+(s.38*P38.46),
    M39.50=(s.39*P39.46),
    M43.50=(s.43*P43.45)+(s.43*P43.46),
    M44.50=(s.44*P44.46),
    
    
    M2.51=(s.2*P2.7),
    M5.51=(s.5*P5.7),
    M6.51=(s.6*P6.7),
    M11.51=(s.11*P11.21),
    M12.51=(s.12*P12.21),
    M13.51=(s.13*P13.21),
    M14.51=(s.14*P14.21),
    M15.51=(s.15*P15.21),
    M19.51=(s.19*P19.21),
    M20.51=(s.20*P20.21),
    M24.51=(s.24*P24.34),
    M25.51=(s.25*P25.34),
    M26.51=(s.26*P26.34),
    M27.51=(s.27*P27.34),
    M28.51=(s.28*P28.34),
    M32.51=(s.32*P32.34),
    M33.51=(s.33*P33.34),
    M37.51=(s.37*P37.47),
    M38.51=(s.38*P38.47),
    M39.51=(s.39*P39.47),
    M40.51=(s.40*P40.47),
    M41.51=(s.41*P41.47),
    M45.51=(s.45*P45.47),
    M46.51=(s.46*P46.47)
  )]
  
  #saves the moves dataset so that you can review the numbers moving between states each cycle (year) 
  #moves_trace[[paste0(y)]]<-copy(moves)
  
  
  #removes transition probabilities (these will be merged in again in the next cycle)
  moves2[,`:=`(
    P1.1=NULL,
    P1.2=NULL,
    P1.8=NULL,
    P2.3=NULL,
    P2.4=NULL,
    P2.5=NULL,
    P2.6=NULL,
    P2.7=NULL,
    P2.9=NULL,
    P2.10=NULL,
    P2.11=NULL,
    P3.3=NULL,
    P3.4=NULL,
    P3.5=NULL,
    P3.6=NULL,
    P3.9=NULL,
    P3.12=NULL,
    P4.4=NULL,
    P4.6=NULL,
    P4.9=NULL,
    P4.13=NULL,
    P5.5=NULL,
    P5.6=NULL,
    P5.7=NULL,
    P5.9=NULL,
    P5.14=NULL,
    P6.6=NULL,
    P6.7=NULL,
    P6.9=NULL,
    P6.15=NULL,
    P7.7=NULL,
    P7.9=NULL,
    P7.16=NULL,
    P11.17=NULL,
    P11.18=NULL,
    P11.19=NULL,
    P11.20=NULL,
    P11.21=NULL,
    P11.22=NULL,
    P11.23=NULL,
    P12.17=NULL,
    P12.18=NULL,
    P12.19=NULL,
    P12.20=NULL,
    P12.21=NULL,
    P12.22=NULL,
    P12.23=NULL,
    P13.18=NULL,
    P13.20=NULL,
    P13.21=NULL,
    P13.22=NULL,
    P13.23=NULL,
    P14.19=NULL,
    P14.20=NULL,
    P14.21=NULL,
    P14.22=NULL,
    P14.23=NULL,
    P15.20=NULL,
    P15.21=NULL,
    P15.22=NULL,
    P15.23=NULL,
    P16.21=NULL,
    P16.22=NULL,
    P16.23=NULL,
    P17.17=NULL,
    P17.18=NULL,
    P17.19=NULL,
    P17.20=NULL,
    P17.23=NULL,
    P18.18=NULL,
    P18.20=NULL,
    P18.23=NULL,
    P19.19=NULL,
    P19.20=NULL,
    P19.21=NULL,
    P19.23=NULL,
    P20.20=NULL,
    P20.21=NULL,
    P20.23=NULL,
    P21.21=NULL,
    P21.23=NULL,
    P11.24=NULL,
    P12.25=NULL,
    P12.26=NULL,
    P12.27=NULL,
    P12.28=NULL,
    P12.29=NULL,
    P13.26=NULL,
    P13.28=NULL,
    P13.29=NULL,
    P14.27=NULL,
    P14.28=NULL,
    P14.29=NULL,
    P15.28=NULL,
    P15.29=NULL,
    P16.29=NULL,
    P17.25=NULL,
    P18.26=NULL,
    P19.27=NULL,
    P20.28=NULL,
    P21.29=NULL,
    P24.30=NULL,
    P24.31=NULL,
    P24.32=NULL,
    P24.33=NULL,
    P24.34=NULL,
    P24.35=NULL,
    P24.36=NULL,
    P25.30=NULL,
    P25.31=NULL,
    P25.32=NULL,
    P25.33=NULL,
    P25.34=NULL,
    P25.35=NULL,
    P25.36=NULL,
    P26.31=NULL,
    P26.33=NULL,
    P26.34=NULL,
    P26.35=NULL,
    P26.36=NULL,
    P27.32=NULL,
    P27.33=NULL,
    P27.34=NULL,
    P27.35=NULL,
    P27.36=NULL,
    P28.33=NULL,
    P28.34=NULL,
    P28.35=NULL,
    P28.36=NULL,
    P29.34=NULL,
    P29.35=NULL,
    P29.36=NULL,
    P30.30=NULL,
    P30.31=NULL,
    P30.32=NULL,
    P30.33=NULL,
    P30.36=NULL,
    P31.31=NULL,
    P31.33=NULL,
    P31.36=NULL,
    P32.32=NULL,
    P32.33=NULL,
    P32.34=NULL,
    P32.36=NULL,
    P33.33=NULL,
    P33.34=NULL,
    P33.36=NULL,
    P34.34=NULL,
    P34.36=NULL,
    P24.37=NULL,
    P25.38=NULL,
    P25.39=NULL,
    P25.40=NULL,
    P25.41=NULL,
    P25.42=NULL,
    P26.39=NULL,
    P26.41=NULL,
    P26.42=NULL,
    P27.40=NULL,
    P27.41=NULL,
    P27.42=NULL,
    P28.41=NULL,
    P28.42=NULL,
    P29.42=NULL,
    P30.38=NULL,
    P31.39=NULL,
    P32.40=NULL,
    P33.41=NULL,
    P34.42=NULL,
    P37.43=NULL,
    P37.44=NULL,
    P37.45=NULL,
    P37.46=NULL,
    P37.47=NULL,
    P37.48=NULL,
    P37.49=NULL,
    P38.43=NULL,
    P38.44=NULL,
    P38.45=NULL,
    P38.46=NULL,
    P38.47=NULL,
    P38.48=NULL,
    P38.49=NULL,
    P39.44=NULL,
    P39.46=NULL,
    P39.47=NULL,
    P39.48=NULL,
    P39.49=NULL,
    P40.45=NULL,
    P40.46=NULL,
    P40.47=NULL,
    P40.48=NULL,
    P40.49=NULL,
    P41.46=NULL,
    P41.47=NULL,
    P41.48=NULL,
    P41.49=NULL,
    P42.47=NULL,
    P42.48=NULL,
    P42.49=NULL,
    P43.43=NULL,
    P43.44=NULL,
    P43.45=NULL,
    P43.46=NULL,
    P43.49=NULL,
    P44.44=NULL,
    P44.46=NULL,
    P44.49=NULL,
    P45.45=NULL,
    P45.46=NULL,
    P45.47=NULL,
    P45.49=NULL,
    P46.46=NULL,
    P46.47=NULL,
    P46.49=NULL,
    P47.47=NULL,
    P47.49=NULL
    
  )]
  
  #calculates health states for the cycle by summing the number of people moving to or remaining in that state 
  
  moves2[,`:=`(
    s.1=M1.1,
    s.2=M1.2,
    s.3=M2.3+M3.3,
    s.4=M2.4+M3.4+M4.4,
    s.5=M2.5+M3.5+M5.5,
    s.6=M2.6+M3.6+M4.6+M5.6+M6.6,
    s.7=M2.7+M5.7+M6.7+M7.7,
    s.8=M1.8,
    s.9=M2.9+M3.9+M4.9+M5.9+M6.9+M7.9,
    s.10=M2.10,
    s.11=M2.11,
    s.12=M3.12,
    s.13=M4.13,
    s.14=M5.14,
    s.15=M6.15,
    s.16=M7.16,
    s.17=M11.17+M12.17+M17.17,
    s.18=M11.18+M12.18+M13.18+M17.18+M18.18,
    s.19=M11.19+M12.19+M14.19+M17.19+M19.19,
    s.20=M11.20+M12.20+M13.20+M14.20+M15.20+M17.20+M18.20+M19.20+M20.20,
    s.21=M11.21+M12.21+M13.21+M14.21+M15.21+M16.21+M19.21+M20.21+M21.21,
    s.22=M11.22+M12.22+M13.22+M14.22+M15.22+M16.22,
    s.23=M11.23+M12.23+M13.23+M14.23+M15.23+M16.23+M17.23+M18.23+M19.23+M20.23+M21.23,
    s.24=M11.24,
    s.25=M12.25+M17.25,
    s.26=M12.26+M13.26+M18.26,
    s.27=M12.27+M14.27+M19.27,
    s.28=M12.28+M13.28+M14.28+M15.28+M20.28,
    s.29=M12.29+M13.29+M14.29+M15.29+M16.29+M21.29,
    s.30=M24.30+M25.30+M30.30,
    s.31=M24.31+M25.31+M26.31+M30.31+M31.31,
    s.32=M24.32+M25.32+M27.32+M30.32+M32.32,
    s.33=M24.33+M25.33+M26.33+M27.33+M28.33+M30.33+M31.33+M32.33+M33.33,
    s.34=M24.34+M25.34+M26.34+M27.34+M28.34+M29.34+M32.34+M33.34+M34.34,
    s.35=M24.35+M25.35+M26.35+M27.35+M28.35+M29.35,
    s.36=M24.36+M25.36+M26.36+M27.36+M28.36+M29.36+M30.36+M31.36+M32.36+M33.36+M34.36,
    s.37=M24.37,
    s.38=M25.38+M30.38,
    s.39=M25.39+M26.39+M31.39,
    s.40=M25.40+M27.40+M32.40,
    s.41=M25.41+M26.41+M27.41+M28.41+M33.41,
    s.42=M25.42+M26.42+M27.42+M28.42+M29.42+M34.42,
    s.43=M37.43+M38.43+M43.43,
    s.44=M37.44+M38.44+M39.44+M43.44+M44.44,
    s.45=M37.45+M38.45+M40.45+M43.45+M45.45,
    s.46=M37.46+M38.46+M39.46+M40.46+M41.46+M43.46+M44.46+M45.46+M46.46,
    s.47=M37.47+M38.47+M39.47+M40.47+M41.47+M42.47+M45.47+M46.47+M47.47,
    s.48=M37.48+M38.48+M39.48+M40.48+M41.48+M42.48,
    s.49=M37.49+M38.49+M39.49+M40.49+M41.49+M42.49+M43.49+M44.49+M45.49+M46.49+M47.49,
    s.50=M2.50+M3.50+M4.50+M11.50+M12.50+M13.50+M17.50+M18.50
    +M24.50+M25.50+M26.50+M30.50+M31.50+M37.50+M38.50+M39.50+M43.50+M44.50,
    s.51=M2.51+M5.51+M6.51+M11.51+M12.51+M13.51+M14.51+M15.51+M19.51+M20.51+M24.51+M25.51
    +M26.51+M27.51+M28.51+M32.51+M33.51+M37.51+M38.51+M39.51+M40.51+M41.51+M45.51
    +M46.51
  )]  
  
  
  #changes highest age in the model (100) to age 40, to allow addition of new 40 year olds into the model 
  moves2[Age==100,Age:=40]
  setkey(moves2, iter, Age, Year)
  
  #selects value of year for first iteration, age 40 (based on order of setkey)
  y<-moves2[.(1,40)]$Year
  
  
  #reads in the projections for age 40 for the relevant year, from the popproj data table
  #not using these
  #newpop<-popproj[Year==y]$Pop
  #sp3<-prev[Age==40 & iter==1]$s.3
  #sp4<-prev[Age==40 & iter==1]$s.4
  #sp5<-prev[Age==40 & iter==1]$s.5
  #sp6<-prev[Age==40 & iter==1]$s.6
  #sp7<-prev[Age==40 & iter==1]$s.7
  
  ####new 40 year olds not included in model####
  
  moves2[Age==40,`:=`(
    Pop=0,
    s.1=0,
    s.2=0,
    s.3=0,
    s.4=0,
    s.5=0,
    s.6=0,
    s.7=0,
    s.8=0,
    s.9=0,
    s.10=0,
    s.11=0,
    s.12=0,
    s.13=0,
    s.14=0,
    s.15=0,
    s.16=0,
    s.17=0,
    s.18=0,
    s.19=0,
    s.20=0,
    s.21=0,
    s.22=0,
    s.23=0,
    s.24=0,
    s.25=0,
    s.26=0,
    s.27=0,
    s.28=0,
    s.29=0,
    s.30=0,
    s.31=0,
    s.32=0,
    s.33=0,
    s.34=0,
    s.35=0,
    s.36=0,
    s.37=0,
    s.38=0,
    s.39=0,
    s.40=0,
    s.41=0,
    s.42=0,
    s.43=0,
    s.44=0,
    s.45=0,
    s.46=0,
    s.47=0,
    s.48=0,
    s.49=0,
    s.50=0,
    s.51=0,
    M1.1=0,
    M1.2=0,
    M1.8=0,
    M2.3=0,
    M2.4=0,
    M2.5=0,
    M2.6=0,
    M2.7=0,
    M2.9=0,
    M2.10=0,
    M2.11=0,
    M3.3=0,
    M3.4=0,
    M3.5=0,
    M3.6=0,
    M3.9=0,
    M3.12=0,
    M4.4=0,
    M4.6=0,
    M4.9=0,
    M4.13=0,
    M5.5=0,
    M5.6=0,
    M5.7=0,
    M5.9=0,
    M5.14=0,
    M6.6=0,
    M6.7=0,
    M6.9=0,
    M6.15=0,
    M7.7=0,
    M7.9=0,
    M7.16=0,
    M11.17=0,
    M11.18=0,
    M11.19=0,
    M11.20=0,
    M11.21=0,
    M11.22=0,
    M11.23=0,
    
    M12.17=0,
    M12.18=0,
    M12.19=0,
    M12.20=0,
    M12.21=0,
    M12.22=0,
    M12.23=0,
    
    M13.18=0,
    M13.20=0,
    M13.21=0,
    M13.22=0,
    M13.23=0,
    
    M14.19=0,
    M14.20=0,
    M14.21=0,
    M14.22=0,
    M14.23=0,
    
    M15.20=0,
    M15.21=0,
    M15.22=0,
    M15.23=0,
    
    M16.21=0,
    M16.22=0,
    M16.23=0,
    
    M17.17=0,
    M17.18=0,
    M17.19=0,
    M17.20=0,
    M17.23=0,
    M18.18=0,
    M18.20=0,
    M18.23=0,
    M19.19=0,
    M19.20=0,
    M19.21=0,
    M19.23=0,
    M20.20=0,
    M20.21=0,
    M20.23=0,
    M21.21=0,
    M21.23=0,
    M11.24=0,
    M12.25=0,
    M12.26=0,
    M12.27=0,
    M12.28=0,
    M12.29=0,
    M13.26=0,
    M13.28=0,
    M13.29=0,
    M14.27=0,
    M14.28=0,
    M14.29=0,
    M15.28=0,
    M15.29=0,
    M16.29=0,
    M17.25=0,
    M18.26=0,
    M19.27=0,
    M20.28=0,
    M21.29=0,
    M24.30=0,
    M24.31=0,
    M24.32=0,
    M24.33=0,
    M24.34=0,
    M24.35=0,
    M24.36=0,
    M25.30=0,
    M25.31=0,
    M25.32=0,
    M25.33=0,
    M25.34=0,
    M25.35=0,
    M25.36=0,
    M26.31=0,
    M26.33=0,
    M26.34=0,
    M26.35=0,
    M26.36=0,
    M27.32=0,
    M27.33=0,
    M27.34=0,
    M27.35=0,
    M27.36=0,
    M28.33=0,
    M28.34=0,
    M28.35=0,
    M28.36=0,
    M29.34=0,
    M29.35=0,
    M29.36=0,
    M30.30=0,
    M30.31=0,
    M30.32=0,
    M30.33=0,
    M30.36=0,
    M31.31=0,
    M31.33=0,
    M31.36=0,
    M32.32=0,
    M32.33=0,
    M32.34=0,
    M32.36=0,
    M33.33=0,
    M33.34=0,
    M33.36=0,
    M34.34=0,
    M34.36=0,
    M24.37=0,
    M25.38=0,
    M25.39=0,
    M25.40=0,
    M25.41=0,
    M25.42=0,
    M26.39=0,
    M26.41=0,
    M26.42=0,
    M27.40=0,
    M27.41=0,
    M27.42=0,
    M28.41=0,
    M28.42=0,
    M29.42=0,
    M30.38=0,
    M31.39=0,
    M32.40=0,
    M33.41=0,
    M34.42=0,
    M37.43=0,
    M37.44=0,
    M37.45=0,
    M37.46=0,
    M37.47=0,
    M37.48=0,
    M37.49=0,
    M38.43=0,
    M38.44=0,
    M38.45=0,
    M38.46=0,
    M38.47=0,
    M38.48=0,
    M38.49=0,
    M39.44=0,
    M39.46=0,
    M39.47=0,
    M39.48=0,
    M39.49=0,
    M40.45=0,
    M40.46=0,
    M40.47=0,
    M40.48=0,
    M40.49=0,
    M41.46=0,
    M41.47=0,
    M41.48=0,
    M41.49=0,
    M42.47=0,
    M42.48=0,
    M42.49=0,
    M43.43=0,
    M43.44=0,
    M43.45=0,
    M43.46=0,
    M43.49=0,
    M44.44=0,
    M44.46=0,
    M44.49=0,
    M45.45=0,
    M45.46=0,
    M45.47=0,
    M45.49=0,
    M46.46=0,
    M46.47=0,
    M46.49=0,
    M47.47=0,
    M47.49=0,
    M2.50=0,
    M3.50=0,
    M4.50=0,
    M11.50=0,
    M12.50=0,
    M13.50=0,
    M17.50=0,
    M18.50=0,
    M24.50=0,
    M25.50=0,
    M26.50=0,
    M30.50=0,
    M31.50=0,
    M37.50=0,
    M38.50=0,
    M39.50=0,
    M43.50=0,
    M44.50=0,
    M2.51=0,
    M5.51=0,
    M6.51=0,
    M11.51=0,
    M12.51=0,
    M13.51=0,
    M14.51=0,
    M15.51=0,
    M19.51=0,
    M20.51=0,
    M24.51=0,
    M25.51=0,
    M26.51=0,
    M27.51=0,
    M28.51=0,
    M32.51=0,
    M33.51=0,
    M37.51=0,
    M38.51=0,
    M39.51=0,
    M40.51=0,
    M41.51=0,
    M45.51=0,
    M46.51=0
  )]
  
  #calculates disease-free health state for 40 year old pop
  #moves2[Age==40,`:=`(s.1=newpop-s.3-s.4-s.5-s.6-s.7)]
  
  
  #merges in transition probabilities
  moves2<-merge(moves2, tps_S, by=c("Age","iter"))
  
  #outputs the results for the health states for this model cycle (year)
  res2[[paste0(y)]]<-copy(moves2)
  
  
  #takes out the M vars so that these can be re-calculated in the next cycle
  moves2[,`:=`(M1.1=NULL,
              M1.2=NULL,
              M1.8=NULL,
              M2.3=NULL,
              M2.4=NULL,
              M2.5=NULL,
              M2.6=NULL,
              M2.7=NULL,
              M2.9=NULL,
              M2.10=NULL,
              M2.11=NULL,
              M3.3=NULL,
              M3.4=NULL,
              M3.5=NULL,
              M3.6=NULL,
              M3.9=NULL,
              M3.12=NULL,
              M4.4=NULL,
              M4.6=NULL,
              M4.9=NULL,
              M4.13=NULL,
              M5.5=NULL,
              M5.6=NULL,
              M5.7=NULL,
              M5.9=NULL,
              M5.14=NULL,
              M6.6=NULL,
              M6.7=NULL,
              M6.9=NULL,
              M6.15=NULL,
              M7.7=NULL,
              M7.9=NULL,
              M7.16=NULL,
              M11.17=NULL,
              M11.18=NULL,
              M11.19=NULL,
              M11.20=NULL,
              M11.21=NULL,
              M11.22=NULL,
              M11.23=NULL,
              
              M12.17=NULL,
              M12.18=NULL,
              M12.19=NULL,
              M12.20=NULL,
              M12.21=NULL,
              M12.22=NULL,
              M12.23=NULL,
              
              M13.18=NULL,
              M13.20=NULL,
              M13.21=NULL,
              M13.22=NULL,
              M13.23=NULL,
              
              M14.19=NULL,
              M14.20=NULL,
              M14.21=NULL,
              M14.22=NULL,
              M14.23=NULL,
              
              M15.20=NULL,
              M15.21=NULL,
              M15.22=NULL,
              M15.23=NULL,
              
              M16.21=NULL,
              M16.22=NULL,
              M16.23=NULL,
              
              M17.17=NULL,
              M17.18=NULL,
              M17.19=NULL,
              M17.20=NULL,
              M17.23=NULL,
              M18.18=NULL,
              M18.20=NULL,
              M18.23=NULL,
              M19.19=NULL,
              M19.20=NULL,
              M19.21=NULL,
              M19.23=NULL,
              M20.20=NULL,
              M20.21=NULL,
              M20.23=NULL,
              M21.21=NULL,
              M21.23=NULL,
              M11.24=NULL,
              M12.25=NULL,
              M12.26=NULL,
              M12.27=NULL,
              M12.28=NULL,
              M12.29=NULL,
              M13.26=NULL,
              M13.28=NULL,
              M13.29=NULL,
              M14.27=NULL,
              M14.28=NULL,
              M14.29=NULL,
              M15.28=NULL,
              M15.29=NULL,
              M16.29=NULL,
              M17.25=NULL,
              M18.26=NULL,
              M19.27=NULL,
              M20.28=NULL,
              M21.29=NULL,
              M24.30=NULL,
              M24.31=NULL,
              M24.32=NULL,
              M24.33=NULL,
              M24.34=NULL,
              M24.35=NULL,
              M24.36=NULL,
              M25.30=NULL,
              M25.31=NULL,
              M25.32=NULL,
              M25.33=NULL,
              M25.34=NULL,
              M25.35=NULL,
              M25.36=NULL,
              M26.31=NULL,
              M26.33=NULL,
              M26.34=NULL,
              M26.35=NULL,
              M26.36=NULL,
              M27.32=NULL,
              M27.33=NULL,
              M27.34=NULL,
              M27.35=NULL,
              M27.36=NULL,
              M28.33=NULL,
              M28.34=NULL,
              M28.35=NULL,
              M28.36=NULL,
              M29.34=NULL,
              M29.35=NULL,
              M29.36=NULL,
              M30.30=NULL,
              M30.31=NULL,
              M30.32=NULL,
              M30.33=NULL,
              M30.36=NULL,
              M31.31=NULL,
              M31.33=NULL,
              M31.36=NULL,
              M32.32=NULL,
              M32.33=NULL,
              M32.34=NULL,
              M32.36=NULL,
              M33.33=NULL,
              M33.34=NULL,
              M33.36=NULL,
              M34.34=NULL,
              M34.36=NULL,
              M24.37=NULL,
              M25.38=NULL,
              M25.39=NULL,
              M25.40=NULL,
              M25.41=NULL,
              M25.42=NULL,
              M26.39=NULL,
              M26.41=NULL,
              M26.42=NULL,
              M27.40=NULL,
              M27.41=NULL,
              M27.42=NULL,
              M28.41=NULL,
              M28.42=NULL,
              M29.42=NULL,
              M30.38=NULL,
              M31.39=NULL,
              M32.40=NULL,
              M33.41=NULL,
              M34.42=NULL,
              M37.43=NULL,
              M37.44=NULL,
              M37.45=NULL,
              M37.46=NULL,
              M37.47=NULL,
              M37.48=NULL,
              M37.49=NULL,
              M38.43=NULL,
              M38.44=NULL,
              M38.45=NULL,
              M38.46=NULL,
              M38.47=NULL,
              M38.48=NULL,
              M38.49=NULL,
              M39.44=NULL,
              M39.46=NULL,
              M39.47=NULL,
              M39.48=NULL,
              M39.49=NULL,
              M40.45=NULL,
              M40.46=NULL,
              M40.47=NULL,
              M40.48=NULL,
              M40.49=NULL,
              M41.46=NULL,
              M41.47=NULL,
              M41.48=NULL,
              M41.49=NULL,
              M42.47=NULL,
              M42.48=NULL,
              M42.49=NULL,
              M43.43=NULL,
              M43.44=NULL,
              M43.45=NULL,
              M43.46=NULL,
              M43.49=NULL,
              M44.44=NULL,
              M44.46=NULL,
              M44.49=NULL,
              M45.45=NULL,
              M45.46=NULL,
              M45.47=NULL,
              M45.49=NULL,
              M46.46=NULL,
              M46.47=NULL,
              M46.49=NULL,
              M47.47=NULL,
              M47.49=NULL,
              M2.50=NULL,
              M3.50=NULL,
              M4.50=NULL,
              M11.50=NULL,
              M12.50=NULL,
              M13.50=NULL,
              M17.50=NULL,
              M18.50=NULL,
              M24.50=NULL,
              M25.50=NULL,
              M26.50=NULL,
              M30.50=NULL,
              M31.50=NULL,
              M37.50=NULL,
              M38.50=NULL,
              M39.50=NULL,
              M43.50=NULL,
              M44.50=NULL,
              M2.51=NULL,
              M5.51=NULL,
              M6.51=NULL,
              M11.51=NULL,
              M12.51=NULL,
              M13.51=NULL,
              M14.51=NULL,
              M15.51=NULL,
              M19.51=NULL,
              M20.51=NULL,
              M24.51=NULL,
              M25.51=NULL,
              M26.51=NULL,
              M27.51=NULL,
              M28.51=NULL,
              M32.51=NULL,
              M33.51=NULL,
              M37.51=NULL,
              M38.51=NULL,
              M39.51=NULL,
              M40.51=NULL,
              M41.51=NULL,
              M45.51=NULL,
              M46.51=NULL
              
              
  )]
  
}  

##END OF LOOP

#binds all of the results for each year together
res2<-rbindlist(res2, use.names=T, fill=T)

#same for the move traces table
#moves_trace<-rbindlist(moves_trace, use.names=T, fill=T)



#takes out the move values - only need health states
res2[,`:=`(M1.1=NULL,
          M1.2=NULL,
          M1.8=NULL,
          M2.3=NULL,
          M2.4=NULL,
          M2.5=NULL,
          M2.6=NULL,
          M2.7=NULL,
          M2.9=NULL,
          M2.10=NULL,
          M2.11=NULL,
          M3.3=NULL,
          M3.4=NULL,
          M3.5=NULL,
          M3.6=NULL,
          M3.9=NULL,
          M3.12=NULL,
          M4.4=NULL,
          M4.6=NULL,
          M4.9=NULL,
          M4.13=NULL,
          M5.5=NULL,
          M5.6=NULL,
          M5.7=NULL,
          M5.9=NULL,
          M5.14=NULL,
          M6.6=NULL,
          M6.7=NULL,
          M6.9=NULL,
          M6.15=NULL,
          M7.7=NULL,
          M7.9=NULL,
          M7.16=NULL,
          M11.17=NULL,
          M11.18=NULL,
          M11.19=NULL,
          M11.20=NULL,
          M11.21=NULL,
          M11.22=NULL,
          M11.23=NULL,
          
          M12.17=NULL,
          M12.18=NULL,
          M12.19=NULL,
          M12.20=NULL,
          M12.21=NULL,
          M12.22=NULL,
          M12.23=NULL,
          
          M13.18=NULL,
          M13.20=NULL,
          M13.21=NULL,
          M13.22=NULL,
          M13.23=NULL,
          
          M14.19=NULL,
          M14.20=NULL,
          M14.21=NULL,
          M14.22=NULL,
          M14.23=NULL,
          
          M15.20=NULL,
          M15.21=NULL,
          M15.22=NULL,
          M15.23=NULL,
          
          M16.21=NULL,
          M16.22=NULL,
          M16.23=NULL,
          
          M17.17=NULL,
          M17.18=NULL,
          M17.19=NULL,
          M17.20=NULL,
          M17.23=NULL,
          M18.18=NULL,
          M18.20=NULL,
          M18.23=NULL,
          M19.19=NULL,
          M19.20=NULL,
          M19.21=NULL,
          M19.23=NULL,
          M20.20=NULL,
          M20.21=NULL,
          M20.23=NULL,
          M21.21=NULL,
          M21.23=NULL,
          M11.24=NULL,
          M12.25=NULL,
          M12.26=NULL,
          M12.27=NULL,
          M12.28=NULL,
          M12.29=NULL,
          M13.26=NULL,
          M13.28=NULL,
          M13.29=NULL,
          M14.27=NULL,
          M14.28=NULL,
          M14.29=NULL,
          M15.28=NULL,
          M15.29=NULL,
          M16.29=NULL,
          M17.25=NULL,
          M18.26=NULL,
          M19.27=NULL,
          M20.28=NULL,
          M21.29=NULL,
          M24.30=NULL,
          M24.31=NULL,
          M24.32=NULL,
          M24.33=NULL,
          M24.34=NULL,
          M24.35=NULL,
          M24.36=NULL,
          M25.30=NULL,
          M25.31=NULL,
          M25.32=NULL,
          M25.33=NULL,
          M25.34=NULL,
          M25.35=NULL,
          M25.36=NULL,
          M26.31=NULL,
          M26.33=NULL,
          M26.34=NULL,
          M26.35=NULL,
          M26.36=NULL,
          M27.32=NULL,
          M27.33=NULL,
          M27.34=NULL,
          M27.35=NULL,
          M27.36=NULL,
          M28.33=NULL,
          M28.34=NULL,
          M28.35=NULL,
          M28.36=NULL,
          M29.34=NULL,
          M29.35=NULL,
          M29.36=NULL,
          M30.30=NULL,
          M30.31=NULL,
          M30.32=NULL,
          M30.33=NULL,
          M30.36=NULL,
          M31.31=NULL,
          M31.33=NULL,
          M31.36=NULL,
          M32.32=NULL,
          M32.33=NULL,
          M32.34=NULL,
          M32.36=NULL,
          M33.33=NULL,
          M33.34=NULL,
          M33.36=NULL,
          M34.34=NULL,
          M34.36=NULL,
          M24.37=NULL,
          M25.38=NULL,
          M25.39=NULL,
          M25.40=NULL,
          M25.41=NULL,
          M25.42=NULL,
          M26.39=NULL,
          M26.41=NULL,
          M26.42=NULL,
          M27.40=NULL,
          M27.41=NULL,
          M27.42=NULL,
          M28.41=NULL,
          M28.42=NULL,
          M29.42=NULL,
          M30.38=NULL,
          M31.39=NULL,
          M32.40=NULL,
          M33.41=NULL,
          M34.42=NULL,
          M37.43=NULL,
          M37.44=NULL,
          M37.45=NULL,
          M37.46=NULL,
          M37.47=NULL,
          M37.48=NULL,
          M37.49=NULL,
          M38.43=NULL,
          M38.44=NULL,
          M38.45=NULL,
          M38.46=NULL,
          M38.47=NULL,
          M38.48=NULL,
          M38.49=NULL,
          M39.44=NULL,
          M39.46=NULL,
          M39.47=NULL,
          M39.48=NULL,
          M39.49=NULL,
          M40.45=NULL,
          M40.46=NULL,
          M40.47=NULL,
          M40.48=NULL,
          M40.49=NULL,
          M41.46=NULL,
          M41.47=NULL,
          M41.48=NULL,
          M41.49=NULL,
          M42.47=NULL,
          M42.48=NULL,
          M42.49=NULL,
          M43.43=NULL,
          M43.44=NULL,
          M43.45=NULL,
          M43.46=NULL,
          M43.49=NULL,
          M44.44=NULL,
          M44.46=NULL,
          M44.49=NULL,
          M45.45=NULL,
          M45.46=NULL,
          M45.47=NULL,
          M45.49=NULL,
          M46.46=NULL,
          M46.47=NULL,
          M46.49=NULL,
          M47.47=NULL,
          M47.49=NULL,
          M2.50=NULL,
          M3.50=NULL,
          M4.50=NULL,
          M11.50=NULL,
          M12.50=NULL,
          M13.50=NULL,
          M17.50=NULL,
          M18.50=NULL,
          M24.50=NULL,
          M25.50=NULL,
          M26.50=NULL,
          M30.50=NULL,
          M31.50=NULL,
          M37.50=NULL,
          M38.50=NULL,
          M39.50=NULL,
          M43.50=NULL,
          M44.50=NULL,
          M2.51=NULL,
          M5.51=NULL,
          M6.51=NULL,
          M11.51=NULL,
          M12.51=NULL,
          M13.51=NULL,
          M14.51=NULL,
          M15.51=NULL,
          M19.51=NULL,
          M20.51=NULL,
          M24.51=NULL,
          M25.51=NULL,
          M26.51=NULL,
          M27.51=NULL,
          M28.51=NULL,
          M32.51=NULL,
          M33.51=NULL,
          M37.51=NULL,
          M38.51=NULL,
          M39.51=NULL,
          M40.51=NULL,
          M41.51=NULL,
          M45.51=NULL,
          M46.51=NULL
          
)]

####THIRD MODEL CALCULATION#### 
##running series of models reduces computational load

#removes transition probabilities (these will be merged in again in the next cycle)

moves2[,`:=`(
  P1.1=NULL,
  P1.2=NULL,
  P1.8=NULL,
  P2.3=NULL,
  P2.4=NULL,
  P2.5=NULL,
  P2.6=NULL,
  P2.7=NULL,
  P2.9=NULL,
  P2.10=NULL,
  P2.11=NULL,
  P3.3=NULL,
  P3.4=NULL,
  P3.5=NULL,
  P3.6=NULL,
  P3.9=NULL,
  P3.12=NULL,
  P4.4=NULL,
  P4.6=NULL,
  P4.9=NULL,
  P4.13=NULL,
  P5.5=NULL,
  P5.6=NULL,
  P5.7=NULL,
  P5.9=NULL,
  P5.14=NULL,
  P6.6=NULL,
  P6.7=NULL,
  P6.9=NULL,
  P6.15=NULL,
  P7.7=NULL,
  P7.9=NULL,
  P7.16=NULL,
  P11.17=NULL,
  P11.18=NULL,
  P11.19=NULL,
  P11.20=NULL,
  P11.21=NULL,
  P11.22=NULL,
  P11.23=NULL,
  P12.17=NULL,
  P12.18=NULL,
  P12.19=NULL,
  P12.20=NULL,
  P12.21=NULL,
  P12.22=NULL,
  P12.23=NULL,
  P13.18=NULL,
  P13.20=NULL,
  P13.21=NULL,
  P13.22=NULL,
  P13.23=NULL,
  P14.19=NULL,
  P14.20=NULL,
  P14.21=NULL,
  P14.22=NULL,
  P14.23=NULL,
  P15.20=NULL,
  P15.21=NULL,
  P15.22=NULL,
  P15.23=NULL,
  P16.21=NULL,
  P16.22=NULL,
  P16.23=NULL,
  P17.17=NULL,
  P17.18=NULL,
  P17.19=NULL,
  P17.20=NULL,
  P17.23=NULL,
  P18.18=NULL,
  P18.20=NULL,
  P18.23=NULL,
  P19.19=NULL,
  P19.20=NULL,
  P19.21=NULL,
  P19.23=NULL,
  P20.20=NULL,
  P20.21=NULL,
  P20.23=NULL,
  P21.21=NULL,
  P21.23=NULL,
  P11.24=NULL,
  P12.25=NULL,
  P12.26=NULL,
  P12.27=NULL,
  P12.28=NULL,
  P12.29=NULL,
  P13.26=NULL,
  P13.28=NULL,
  P13.29=NULL,
  P14.27=NULL,
  P14.28=NULL,
  P14.29=NULL,
  P15.28=NULL,
  P15.29=NULL,
  P16.29=NULL,
  P17.25=NULL,
  P18.26=NULL,
  P19.27=NULL,
  P20.28=NULL,
  P21.29=NULL,
  P24.30=NULL,
  P24.31=NULL,
  P24.32=NULL,
  P24.33=NULL,
  P24.34=NULL,
  P24.35=NULL,
  P24.36=NULL,
  P25.30=NULL,
  P25.31=NULL,
  P25.32=NULL,
  P25.33=NULL,
  P25.34=NULL,
  P25.35=NULL,
  P25.36=NULL,
  P26.31=NULL,
  P26.33=NULL,
  P26.34=NULL,
  P26.35=NULL,
  P26.36=NULL,
  P27.32=NULL,
  P27.33=NULL,
  P27.34=NULL,
  P27.35=NULL,
  P27.36=NULL,
  P28.33=NULL,
  P28.34=NULL,
  P28.35=NULL,
  P28.36=NULL,
  P29.34=NULL,
  P29.35=NULL,
  P29.36=NULL,
  P30.30=NULL,
  P30.31=NULL,
  P30.32=NULL,
  P30.33=NULL,
  P30.36=NULL,
  P31.31=NULL,
  P31.33=NULL,
  P31.36=NULL,
  P32.32=NULL,
  P32.33=NULL,
  P32.34=NULL,
  P32.36=NULL,
  P33.33=NULL,
  P33.34=NULL,
  P33.36=NULL,
  P34.34=NULL,
  P34.36=NULL,
  P24.37=NULL,
  P25.38=NULL,
  P25.39=NULL,
  P25.40=NULL,
  P25.41=NULL,
  P25.42=NULL,
  P26.39=NULL,
  P26.41=NULL,
  P26.42=NULL,
  P27.40=NULL,
  P27.41=NULL,
  P27.42=NULL,
  P28.41=NULL,
  P28.42=NULL,
  P29.42=NULL,
  P30.38=NULL,
  P31.39=NULL,
  P32.40=NULL,
  P33.41=NULL,
  P34.42=NULL,
  P37.43=NULL,
  P37.44=NULL,
  P37.45=NULL,
  P37.46=NULL,
  P37.47=NULL,
  P37.48=NULL,
  P37.49=NULL,
  P38.43=NULL,
  P38.44=NULL,
  P38.45=NULL,
  P38.46=NULL,
  P38.47=NULL,
  P38.48=NULL,
  P38.49=NULL,
  P39.44=NULL,
  P39.46=NULL,
  P39.47=NULL,
  P39.48=NULL,
  P39.49=NULL,
  P40.45=NULL,
  P40.46=NULL,
  P40.47=NULL,
  P40.48=NULL,
  P40.49=NULL,
  P41.46=NULL,
  P41.47=NULL,
  P41.48=NULL,
  P41.49=NULL,
  P42.47=NULL,
  P42.48=NULL,
  P42.49=NULL,
  P43.43=NULL,
  P43.44=NULL,
  P43.45=NULL,
  P43.46=NULL,
  P43.49=NULL,
  P44.44=NULL,
  P44.46=NULL,
  P44.49=NULL,
  P45.45=NULL,
  P45.46=NULL,
  P45.47=NULL,
  P45.49=NULL,
  P46.46=NULL,
  P46.47=NULL,
  P46.49=NULL,
  P47.47=NULL,
  P47.49=NULL
  
)]

#creates new initial states file based on last cycle of model 1 
initial.states3<-merge(moves2, tps_S, by=c("Age","iter"))

#sort by iteration
setkey(initial.states3,iter)

#saves out the results as second results file

res3<-list()
res3[[paste0(FIRST.YEAR3)]]<-copy(initial.states2)


# creates a new moves file for the second model 
moves3<-copy(initial.states3)

#Starts at first year 2 
y<-FIRST.YEAR3

##Loop starts here - current year is stored in y
##loop continues until y = LAST.YEAR (defined at beginning of code)

while (y<LAST.YEAR3) 
{
  
  
  #age and year are increased by 1 year to move model along by one year  
  
  moves3[,Age:=Age+1]
  moves3[,Year:=Year+1]
  
  #Taking out existing tprs - need to re-merge, as age not synced up
  
  moves3[,`:=`(
    P1.1=NULL,
    P1.2=NULL,
    P1.8=NULL,
    P2.3=NULL,
    P2.4=NULL,
    P2.5=NULL,
    P2.6=NULL,
    P2.7=NULL,
    P2.9=NULL,
    P2.10=NULL,
    P2.11=NULL,
    P3.3=NULL,
    P3.4=NULL,
    P3.5=NULL,
    P3.6=NULL,
    P3.9=NULL,
    P3.12=NULL,
    P4.4=NULL,
    P4.6=NULL,
    P4.9=NULL,
    P4.13=NULL,
    P5.5=NULL,
    P5.6=NULL,
    P5.7=NULL,
    P5.9=NULL,
    P5.14=NULL,
    P6.6=NULL,
    P6.7=NULL,
    P6.9=NULL,
    P6.15=NULL,
    P7.7=NULL,
    P7.9=NULL,
    P7.16=NULL,
    P11.17=NULL,
    P11.18=NULL,
    P11.19=NULL,
    P11.20=NULL,
    P11.21=NULL,
    P11.22=NULL,
    P11.23=NULL,
    P12.17=NULL,
    P12.18=NULL,
    P12.19=NULL,
    P12.20=NULL,
    P12.21=NULL,
    P12.22=NULL,
    P12.23=NULL,
    
    
    P13.18=NULL,
    P13.20=NULL,
    P13.21=NULL,
    P13.22=NULL,
    P13.23=NULL,
    
    P14.19=NULL,
    P14.20=NULL,
    P14.21=NULL,
    P14.22=NULL,
    P14.23=NULL,
    P15.20=NULL,
    P15.21=NULL,
    P15.22=NULL,
    P15.23=NULL,
    P16.21=NULL,
    P16.22=NULL,
    P16.23=NULL,
    P17.17=NULL,
    P17.18=NULL,
    P17.19=NULL,
    P17.20=NULL,
    P17.23=NULL,
    P18.18=NULL,
    P18.20=NULL,
    P18.23=NULL,
    P19.19=NULL,
    P19.20=NULL,
    P19.21=NULL,
    P19.23=NULL,
    P20.20=NULL,
    P20.21=NULL,
    P20.23=NULL,
    P21.21=NULL,
    P21.23=NULL,
    P11.24=NULL,
    P12.25=NULL,
    P12.26=NULL,
    P12.27=NULL,
    P12.28=NULL,
    P12.29=NULL,
    P13.26=NULL,
    P13.28=NULL,
    P13.29=NULL,
    P14.27=NULL,
    P14.28=NULL,
    P14.29=NULL,
    P15.28=NULL,
    P15.29=NULL,
    P16.29=NULL,
    P17.25=NULL,
    P18.26=NULL,
    P19.27=NULL,
    P20.28=NULL,
    P21.29=NULL,
    P24.30=NULL,
    P24.31=NULL,
    P24.32=NULL,
    P24.33=NULL,
    P24.34=NULL,
    P24.35=NULL,
    P24.36=NULL,
    P25.30=NULL,
    P25.31=NULL,
    P25.32=NULL,
    P25.33=NULL,
    P25.34=NULL,
    P25.35=NULL,
    P25.36=NULL,
    P26.31=NULL,
    P26.33=NULL,
    P26.34=NULL,
    P26.35=NULL,
    P26.36=NULL,
    P27.32=NULL,
    P27.33=NULL,
    P27.34=NULL,
    P27.35=NULL,
    P27.36=NULL,
    P28.33=NULL,
    P28.34=NULL,
    P28.35=NULL,
    P28.36=NULL,
    P29.34=NULL,
    P29.35=NULL,
    P29.36=NULL,
    P30.30=NULL,
    P30.31=NULL,
    P30.32=NULL,
    P30.33=NULL,
    P30.36=NULL,
    P31.31=NULL,
    P31.33=NULL,
    P31.36=NULL,
    P32.32=NULL,
    P32.33=NULL,
    P32.34=NULL,
    P32.36=NULL,
    P33.33=NULL,
    P33.34=NULL,
    P33.36=NULL,
    P34.34=NULL,
    P34.36=NULL,
    P24.37=NULL,
    P25.38=NULL,
    P25.39=NULL,
    P25.40=NULL,
    P25.41=NULL,
    P25.42=NULL,
    P26.39=NULL,
    P26.41=NULL,
    P26.42=NULL,
    P27.40=NULL,
    P27.41=NULL,
    P27.42=NULL,
    P28.41=NULL,
    P28.42=NULL,
    P29.42=NULL,
    P30.38=NULL,
    P31.39=NULL,
    P32.40=NULL,
    P33.41=NULL,
    P34.42=NULL,
    P37.43=NULL,
    P37.44=NULL,
    P37.45=NULL,
    P37.46=NULL,
    P37.47=NULL,
    P37.48=NULL,
    P37.49=NULL,
    P38.43=NULL,
    P38.44=NULL,
    P38.45=NULL,
    P38.46=NULL,
    P38.47=NULL,
    P38.48=NULL,
    P38.49=NULL,
    P39.44=NULL,
    P39.46=NULL,
    P39.47=NULL,
    P39.48=NULL,
    P39.49=NULL,
    P40.45=NULL,
    P40.46=NULL,
    P40.47=NULL,
    P40.48=NULL,
    P40.49=NULL,
    P41.46=NULL,
    P41.47=NULL,
    P41.48=NULL,
    P41.49=NULL,
    P42.47=NULL,
    P42.48=NULL,
    P42.49=NULL,
    P43.43=NULL,
    P43.44=NULL,
    P43.45=NULL,
    P43.46=NULL,
    P43.49=NULL,
    P44.44=NULL,
    P44.46=NULL,
    P44.49=NULL,
    P45.45=NULL,
    P45.46=NULL,
    P45.47=NULL,
    P45.49=NULL,
    P46.46=NULL,
    P46.47=NULL,
    P46.49=NULL,
    P47.47=NULL,
    P47.49=NULL
    
  )]
  
  #merging in tps again so that age is synced up correctly
  
  moves3<-merge(moves3, tps_S, by=c("Age","iter"))
  
  #M50 and M51 included to capture incident CIND and dementia
  #Health states are multipled by the relevant transition probabilites to calculate the numbers transitioning between states
  #the M var names indicate between which states the transition is occuring
  #e.g. M1.2 is the numbers moving between state 1 and 2, which is calculated by applying the tpr p1.2 to state 1 (s.1)
  
  moves3[,`:=`(
    M1.1=s.1*P1.1,
    M1.2=s.1*P1.2,
    M1.8=s.1*P1.8,
    M2.3=s.2*P2.3,
    M2.4=s.2*P2.4,
    M2.5=s.2*P2.5,
    M2.6=s.2*P2.6,
    M2.7=s.2*P2.7,
    M2.9=s.2*P2.9,
    M2.10=s.2*P2.10,
    M2.11=s.2*P2.11,
    M3.3=s.3*P3.3,
    M3.4=s.3*P3.4,
    M3.5=s.3*P3.5,
    M3.6=s.3*P3.6,
    M3.9=s.3*P3.9,
    M3.12=s.3*P3.12,
    M4.4=s.4*P4.4,
    M4.6=s.4*P4.6,
    M4.9=s.4*P4.9,
    M4.13=s.4*P4.13,
    M5.5=s.5*P5.5,
    M5.6=s.5*P5.6,
    M5.7=s.5*P5.7,
    M5.9=s.5*P5.9,
    M5.14=s.5*P5.14,
    M6.6=s.6*P6.6,
    M6.7=s.6*P6.7,
    M6.9=s.6*P6.9,
    M6.15=s.6*P6.15,
    M7.7=s.7*P7.7,
    M7.9=s.7*P7.9,
    M7.16=s.7*P7.16,
    M11.17=s.11*P11.17,
    M11.18=s.11*P11.18,
    M11.19=s.11*P11.19,
    M11.20=s.11*P11.20,
    M11.21=s.11*P11.21,
    M11.22=s.11*P11.22,
    M11.23=s.11*P11.23,
    M12.17=s.12*P12.17,
    M12.18=s.12*P12.18,
    M12.19=s.12*P12.19,
    M12.20=s.12*P12.20,
    M12.21=s.12*P12.21,
    M12.22=s.12*P12.22,
    M12.23=s.12*P12.23,
    M13.18=s.13*P13.18,
    M13.20=s.13*P13.20,
    M13.21=s.13*P13.21,
    M13.22=s.13*P13.22,
    M13.23=s.13*P13.23,
    
    M14.19=s.14*P14.19,
    M14.20=s.14*P14.20,
    M14.21=s.14*P14.21,
    M14.22=s.14*P14.22,
    M14.23=s.14*P14.23,
    
    M15.20=s.15*P15.20,
    M15.21=s.15*P15.21,
    M15.22=s.15*P15.22,
    M15.23=s.15*P15.23,
    
    M16.21=s.16*P16.21,
    M16.22=s.16*P16.22,
    M16.23=s.16*P16.23,
    
    M17.17=s.17*P17.17,
    M17.18=s.17*P17.18,
    M17.19=s.17*P17.19,
    M17.20=s.17*P17.20,
    M17.23=s.17*P17.23,
    M18.18=s.18*P18.18,
    M18.20=s.18*P18.20,
    M18.23=s.18*P18.23,
    M19.19=s.19*P19.19,
    M19.20=s.19*P19.20,
    M19.21=s.19*P19.21,
    M19.23=s.19*P19.23,
    M20.20=s.20*P20.20,
    M20.21=s.20*P20.21,
    M20.23=s.20*P20.23,
    M21.21=s.21*P21.21,
    M21.23=s.21*P21.23,
    
    
    M11.24=s.11*P11.24,
    M12.25=s.12*P12.25,
    M12.26=s.12*P12.26,
    M12.27=s.12*P12.27,
    M12.28=s.12*P12.28,
    M12.29=s.12*P12.29,
    M13.26=s.13*P13.26,
    M13.28=s.13*P13.28,
    M13.29=s.13*P13.29,
    M14.27=s.14*P14.27,
    M14.28=s.14*P14.28,
    M14.29=s.14*P14.29,
    M15.28=s.15*P15.28,
    M15.29=s.15*P15.29,
    M16.29=s.16*P16.29,
    M17.25=s.17*P17.25,
    
    M18.26=s.18*P18.26,
    M19.27=s.19*P19.27,
    M20.28=s.20*P20.28,
    M21.29=s.21*P21.29,
    M24.30=s.24*P24.30,
    M24.31=s.24*P24.31,
    M24.32=s.24*P24.32,
    M24.33=s.24*P24.33,
    M24.34=s.24*P24.34,
    M24.35=s.24*P24.35,
    M24.36=s.24*P24.36,
    M25.30=s.25*P25.30,
    M25.31=s.25*P25.31,
    M25.32=s.25*P25.32,
    M25.33=s.25*P25.33,
    M25.34=s.25*P25.34,
    M25.35=s.25*P25.35,
    M25.36=s.25*P25.36,
    M26.31=s.26*P26.31,
    M26.33=s.26*P26.33,
    M26.34=s.26*P26.34,
    M26.35=s.26*P26.35,
    M26.36=s.26*P26.36,
    M27.32=s.27*P27.32,
    M27.33=s.27*P27.33,
    M27.34=s.27*P27.34,
    M27.35=s.27*P27.35,
    M27.36=s.27*P27.36,
    M28.33=s.28*P28.33,
    M28.34=s.28*P28.34,
    M28.35=s.28*P28.35,
    M28.36=s.28*P28.36,
    M29.34=s.29*P29.34,
    M29.35=s.29*P29.35,
    M29.36=s.29*P29.36,
    M30.30=s.30*P30.30,
    M30.31=s.30*P30.31,
    M30.32=s.30*P30.32,
    M30.33=s.30*P30.33,
    M30.36=s.30*P30.36,
    M31.31=s.31*P31.31,
    M31.33=s.31*P31.33,
    M31.36=s.31*P31.36,
    M32.32=s.32*P32.32,
    M32.33=s.32*P32.33,
    M32.34=s.32*P32.34,
    M32.36=s.32*P32.36,
    M33.33=s.33*P33.33,
    M33.34=s.33*P33.34,
    M33.36=s.33*P33.36,
    M34.34=s.34*P34.34,
    M34.36=s.34*P34.36,
    M24.37=s.24*P24.37,
    M25.38=s.25*P25.38,
    M25.39=s.25*P25.39,
    M25.40=s.25*P25.40,
    M25.41=s.25*P25.41,
    M25.42=s.25*P25.42,
    M26.39=s.26*P26.39,
    M26.41=s.26*P26.41,
    M26.42=s.26*P26.42,
    M27.40=s.27*P27.40,
    M27.41=s.27*P27.41,
    M27.42=s.27*P27.42,
    M28.41=s.28*P28.41,
    M28.42=s.28*P28.42,
    M29.42=s.29*P29.42,
    M30.38=s.30*P30.38,
    M31.39=s.31*P31.39,
    M32.40=s.32*P32.40,
    M33.41=s.33*P33.41,
    M34.42=s.34*P34.42,
    M37.43=s.37*P37.43,
    M37.44=s.37*P37.44,
    M37.45=s.37*P37.45,
    M37.46=s.37*P37.46,
    M37.47=s.37*P37.47,
    M37.48=s.37*P37.48,
    M37.49=s.37*P37.49,
    M38.43=s.38*P38.43,
    M38.44=s.38*P38.44,
    M38.45=s.38*P38.45,
    M38.46=s.38*P38.46,
    M38.47=s.38*P38.47,
    M38.48=s.38*P38.48,
    M38.49=s.38*P38.49,
    M39.44=s.39*P39.44,
    M39.46=s.39*P39.46,
    M39.47=s.39*P39.47,
    M39.48=s.39*P39.48,
    M39.49=s.39*P39.49,
    M40.45=s.40*P40.45,
    M40.46=s.40*P40.46,
    M40.47=s.40*P40.47,
    M40.48=s.40*P40.48,
    M40.49=s.40*P40.49,
    M41.46=s.41*P41.46,
    M41.47=s.41*P41.47,
    M41.48=s.41*P41.48,
    M41.49=s.41*P41.49,
    M42.47=s.42*P42.47,
    M42.48=s.42*P42.48,
    M42.49=s.42*P42.49,
    M43.43=s.43*P43.43,
    M43.44=s.43*P43.44,
    M43.45=s.43*P43.45,
    M43.46=s.43*P43.46,
    M43.49=s.43*P43.49,
    M44.44=s.44*P44.44,
    M44.46=s.44*P44.46,
    M44.49=s.44*P44.49,
    M45.45=s.45*P45.45,
    M45.46=s.45*P45.46,
    M45.47=s.45*P45.47,
    M45.49=s.45*P45.49,
    M46.46=s.46*P46.46,
    M46.47=s.46*P46.47,
    M46.49=s.46*P46.49,
    M47.47=s.47*P47.47,
    M47.49=s.47*P47.49,
    
    M2.50=(s.2*P2.5)+(s.2*P2.6),
    M3.50=(s.3*P3.5)+(s.3*P3.6),
    M4.50=(s.4*P4.6),
    M11.50=(s.11*P11.19)+(s.11*P11.20),
    M12.50=(s.12*P12.19)+(s.12*P12.20),
    M13.50=(s.13*P13.20),
    M17.50=(s.17*P17.19)+(s.17*P17.20),
    M18.50=(s.18*P18.20),
    M24.50=(s.24*P24.32)+(s.24*P24.33),
    M25.50=(s.25*P25.32)+(s.25*P25.33),
    M26.50=(s.26*P26.33),
    M30.50=(s.30*P30.32)+(s.30*P30.33),
    M31.50=(s.31*P31.33),
    M37.50=(s.37*P37.45)+(s.37*P37.46),
    M38.50=(s.38*P38.45)+(s.38*P38.46),
    M39.50=(s.39*P39.46),
    M43.50=(s.43*P43.45)+(s.43*P43.46),
    M44.50=(s.44*P44.46),
    
    
    M2.51=(s.2*P2.7),
    M5.51=(s.5*P5.7),
    M6.51=(s.6*P6.7),
    M11.51=(s.11*P11.21),
    M12.51=(s.12*P12.21),
    M13.51=(s.13*P13.21),
    M14.51=(s.14*P14.21),
    M15.51=(s.15*P15.21),
    M19.51=(s.19*P19.21),
    M20.51=(s.20*P20.21),
    M24.51=(s.24*P24.34),
    M25.51=(s.25*P25.34),
    M26.51=(s.26*P26.34),
    M27.51=(s.27*P27.34),
    M28.51=(s.28*P28.34),
    M32.51=(s.32*P32.34),
    M33.51=(s.33*P33.34),
    M37.51=(s.37*P37.47),
    M38.51=(s.38*P38.47),
    M39.51=(s.39*P39.47),
    M40.51=(s.40*P40.47),
    M41.51=(s.41*P41.47),
    M45.51=(s.45*P45.47),
    M46.51=(s.46*P46.47)
  )]
  
  #saves the moves dataset so that you can review the numbers moving between states each cycle (year) 
  #moves_trace[[paste0(y)]]<-copy(moves)
  
  
  #removes transition probabilities (these will be merged in again in the next cycle)
  moves3[,`:=`(
    P1.1=NULL,
    P1.2=NULL,
    P1.8=NULL,
    P2.3=NULL,
    P2.4=NULL,
    P2.5=NULL,
    P2.6=NULL,
    P2.7=NULL,
    P2.9=NULL,
    P2.10=NULL,
    P2.11=NULL,
    P3.3=NULL,
    P3.4=NULL,
    P3.5=NULL,
    P3.6=NULL,
    P3.9=NULL,
    P3.12=NULL,
    P4.4=NULL,
    P4.6=NULL,
    P4.9=NULL,
    P4.13=NULL,
    P5.5=NULL,
    P5.6=NULL,
    P5.7=NULL,
    P5.9=NULL,
    P5.14=NULL,
    P6.6=NULL,
    P6.7=NULL,
    P6.9=NULL,
    P6.15=NULL,
    P7.7=NULL,
    P7.9=NULL,
    P7.16=NULL,
    P11.17=NULL,
    P11.18=NULL,
    P11.19=NULL,
    P11.20=NULL,
    P11.21=NULL,
    P11.22=NULL,
    P11.23=NULL,
    P12.17=NULL,
    P12.18=NULL,
    P12.19=NULL,
    P12.20=NULL,
    P12.21=NULL,
    P12.22=NULL,
    P12.23=NULL,
    P13.18=NULL,
    P13.20=NULL,
    P13.21=NULL,
    P13.22=NULL,
    P13.23=NULL,
    P14.19=NULL,
    P14.20=NULL,
    P14.21=NULL,
    P14.22=NULL,
    P14.23=NULL,
    P15.20=NULL,
    P15.21=NULL,
    P15.22=NULL,
    P15.23=NULL,
    P16.21=NULL,
    P16.22=NULL,
    P16.23=NULL,
    P17.17=NULL,
    P17.18=NULL,
    P17.19=NULL,
    P17.20=NULL,
    P17.23=NULL,
    P18.18=NULL,
    P18.20=NULL,
    P18.23=NULL,
    P19.19=NULL,
    P19.20=NULL,
    P19.21=NULL,
    P19.23=NULL,
    P20.20=NULL,
    P20.21=NULL,
    P20.23=NULL,
    P21.21=NULL,
    P21.23=NULL,
    P11.24=NULL,
    P12.25=NULL,
    P12.26=NULL,
    P12.27=NULL,
    P12.28=NULL,
    P12.29=NULL,
    P13.26=NULL,
    P13.28=NULL,
    P13.29=NULL,
    P14.27=NULL,
    P14.28=NULL,
    P14.29=NULL,
    P15.28=NULL,
    P15.29=NULL,
    P16.29=NULL,
    P17.25=NULL,
    P18.26=NULL,
    P19.27=NULL,
    P20.28=NULL,
    P21.29=NULL,
    P24.30=NULL,
    P24.31=NULL,
    P24.32=NULL,
    P24.33=NULL,
    P24.34=NULL,
    P24.35=NULL,
    P24.36=NULL,
    P25.30=NULL,
    P25.31=NULL,
    P25.32=NULL,
    P25.33=NULL,
    P25.34=NULL,
    P25.35=NULL,
    P25.36=NULL,
    P26.31=NULL,
    P26.33=NULL,
    P26.34=NULL,
    P26.35=NULL,
    P26.36=NULL,
    P27.32=NULL,
    P27.33=NULL,
    P27.34=NULL,
    P27.35=NULL,
    P27.36=NULL,
    P28.33=NULL,
    P28.34=NULL,
    P28.35=NULL,
    P28.36=NULL,
    P29.34=NULL,
    P29.35=NULL,
    P29.36=NULL,
    P30.30=NULL,
    P30.31=NULL,
    P30.32=NULL,
    P30.33=NULL,
    P30.36=NULL,
    P31.31=NULL,
    P31.33=NULL,
    P31.36=NULL,
    P32.32=NULL,
    P32.33=NULL,
    P32.34=NULL,
    P32.36=NULL,
    P33.33=NULL,
    P33.34=NULL,
    P33.36=NULL,
    P34.34=NULL,
    P34.36=NULL,
    P24.37=NULL,
    P25.38=NULL,
    P25.39=NULL,
    P25.40=NULL,
    P25.41=NULL,
    P25.42=NULL,
    P26.39=NULL,
    P26.41=NULL,
    P26.42=NULL,
    P27.40=NULL,
    P27.41=NULL,
    P27.42=NULL,
    P28.41=NULL,
    P28.42=NULL,
    P29.42=NULL,
    P30.38=NULL,
    P31.39=NULL,
    P32.40=NULL,
    P33.41=NULL,
    P34.42=NULL,
    P37.43=NULL,
    P37.44=NULL,
    P37.45=NULL,
    P37.46=NULL,
    P37.47=NULL,
    P37.48=NULL,
    P37.49=NULL,
    P38.43=NULL,
    P38.44=NULL,
    P38.45=NULL,
    P38.46=NULL,
    P38.47=NULL,
    P38.48=NULL,
    P38.49=NULL,
    P39.44=NULL,
    P39.46=NULL,
    P39.47=NULL,
    P39.48=NULL,
    P39.49=NULL,
    P40.45=NULL,
    P40.46=NULL,
    P40.47=NULL,
    P40.48=NULL,
    P40.49=NULL,
    P41.46=NULL,
    P41.47=NULL,
    P41.48=NULL,
    P41.49=NULL,
    P42.47=NULL,
    P42.48=NULL,
    P42.49=NULL,
    P43.43=NULL,
    P43.44=NULL,
    P43.45=NULL,
    P43.46=NULL,
    P43.49=NULL,
    P44.44=NULL,
    P44.46=NULL,
    P44.49=NULL,
    P45.45=NULL,
    P45.46=NULL,
    P45.47=NULL,
    P45.49=NULL,
    P46.46=NULL,
    P46.47=NULL,
    P46.49=NULL,
    P47.47=NULL,
    P47.49=NULL
    
  )]
  
  #calculates health states for the cycle by summing the number of people moving to or remaining in that state 
  
  moves3[,`:=`(
    s.1=M1.1,
    s.2=M1.2,
    s.3=M2.3+M3.3,
    s.4=M2.4+M3.4+M4.4,
    s.5=M2.5+M3.5+M5.5,
    s.6=M2.6+M3.6+M4.6+M5.6+M6.6,
    s.7=M2.7+M5.7+M6.7+M7.7,
    s.8=M1.8,
    s.9=M2.9+M3.9+M4.9+M5.9+M6.9+M7.9,
    s.10=M2.10,
    s.11=M2.11,
    s.12=M3.12,
    s.13=M4.13,
    s.14=M5.14,
    s.15=M6.15,
    s.16=M7.16,
    s.17=M11.17+M12.17+M17.17,
    s.18=M11.18+M12.18+M13.18+M17.18+M18.18,
    s.19=M11.19+M12.19+M14.19+M17.19+M19.19,
    s.20=M11.20+M12.20+M13.20+M14.20+M15.20+M17.20+M18.20+M19.20+M20.20,
    s.21=M11.21+M12.21+M13.21+M14.21+M15.21+M16.21+M19.21+M20.21+M21.21,
    s.22=M11.22+M12.22+M13.22+M14.22+M15.22+M16.22,
    s.23=M11.23+M12.23+M13.23+M14.23+M15.23+M16.23+M17.23+M18.23+M19.23+M20.23+M21.23,
    s.24=M11.24,
    s.25=M12.25+M17.25,
    s.26=M12.26+M13.26+M18.26,
    s.27=M12.27+M14.27+M19.27,
    s.28=M12.28+M13.28+M14.28+M15.28+M20.28,
    s.29=M12.29+M13.29+M14.29+M15.29+M16.29+M21.29,
    s.30=M24.30+M25.30+M30.30,
    s.31=M24.31+M25.31+M26.31+M30.31+M31.31,
    s.32=M24.32+M25.32+M27.32+M30.32+M32.32,
    s.33=M24.33+M25.33+M26.33+M27.33+M28.33+M30.33+M31.33+M32.33+M33.33,
    s.34=M24.34+M25.34+M26.34+M27.34+M28.34+M29.34+M32.34+M33.34+M34.34,
    s.35=M24.35+M25.35+M26.35+M27.35+M28.35+M29.35,
    s.36=M24.36+M25.36+M26.36+M27.36+M28.36+M29.36+M30.36+M31.36+M32.36+M33.36+M34.36,
    s.37=M24.37,
    s.38=M25.38+M30.38,
    s.39=M25.39+M26.39+M31.39,
    s.40=M25.40+M27.40+M32.40,
    s.41=M25.41+M26.41+M27.41+M28.41+M33.41,
    s.42=M25.42+M26.42+M27.42+M28.42+M29.42+M34.42,
    s.43=M37.43+M38.43+M43.43,
    s.44=M37.44+M38.44+M39.44+M43.44+M44.44,
    s.45=M37.45+M38.45+M40.45+M43.45+M45.45,
    s.46=M37.46+M38.46+M39.46+M40.46+M41.46+M43.46+M44.46+M45.46+M46.46,
    s.47=M37.47+M38.47+M39.47+M40.47+M41.47+M42.47+M45.47+M46.47+M47.47,
    s.48=M37.48+M38.48+M39.48+M40.48+M41.48+M42.48,
    s.49=M37.49+M38.49+M39.49+M40.49+M41.49+M42.49+M43.49+M44.49+M45.49+M46.49+M47.49,
    s.50=M2.50+M3.50+M4.50+M11.50+M12.50+M13.50+M17.50+M18.50
    +M24.50+M25.50+M26.50+M30.50+M31.50+M37.50+M38.50+M39.50+M43.50+M44.50,
    s.51=M2.51+M5.51+M6.51+M11.51+M12.51+M13.51+M14.51+M15.51+M19.51+M20.51+M24.51+M25.51
    +M26.51+M27.51+M28.51+M32.51+M33.51+M37.51+M38.51+M39.51+M40.51+M41.51+M45.51
    +M46.51
  )]  
  
  
  #changes highest age in the model (100) to age 40, to allow addition of new 40 year olds into the model 
  moves3[Age==100,Age:=40]
  setkey(moves3, iter, Age, Year)
  
  #selects value of year for first iteration, age 40 (based on order of setkey)
  y<-moves3[.(1,40)]$Year
  
  
  #reads in the projections for age 40 for the relevant year, from the popproj data table
  #not using these
  #newpop<-popproj[Year==y]$Pop
  #sp3<-prev[Age==40 & iter==1]$s.3
  #sp4<-prev[Age==40 & iter==1]$s.4
  #sp5<-prev[Age==40 & iter==1]$s.5
  #sp6<-prev[Age==40 & iter==1]$s.6
  #sp7<-prev[Age==40 & iter==1]$s.7
  
  ####new 40 year olds not included in model####
  
  moves3[Age==40,`:=`(
    Pop=0,
    s.1=0,
    s.2=0,
    s.3=0,
    s.4=0,
    s.5=0,
    s.6=0,
    s.7=0,
    s.8=0,
    s.9=0,
    s.10=0,
    s.11=0,
    s.12=0,
    s.13=0,
    s.14=0,
    s.15=0,
    s.16=0,
    s.17=0,
    s.18=0,
    s.19=0,
    s.20=0,
    s.21=0,
    s.22=0,
    s.23=0,
    s.24=0,
    s.25=0,
    s.26=0,
    s.27=0,
    s.28=0,
    s.29=0,
    s.30=0,
    s.31=0,
    s.32=0,
    s.33=0,
    s.34=0,
    s.35=0,
    s.36=0,
    s.37=0,
    s.38=0,
    s.39=0,
    s.40=0,
    s.41=0,
    s.42=0,
    s.43=0,
    s.44=0,
    s.45=0,
    s.46=0,
    s.47=0,
    s.48=0,
    s.49=0,
    s.50=0,
    s.51=0,
    M1.1=0,
    M1.2=0,
    M1.8=0,
    M2.3=0,
    M2.4=0,
    M2.5=0,
    M2.6=0,
    M2.7=0,
    M2.9=0,
    M2.10=0,
    M2.11=0,
    M3.3=0,
    M3.4=0,
    M3.5=0,
    M3.6=0,
    M3.9=0,
    M3.12=0,
    M4.4=0,
    M4.6=0,
    M4.9=0,
    M4.13=0,
    M5.5=0,
    M5.6=0,
    M5.7=0,
    M5.9=0,
    M5.14=0,
    M6.6=0,
    M6.7=0,
    M6.9=0,
    M6.15=0,
    M7.7=0,
    M7.9=0,
    M7.16=0,
    M11.17=0,
    M11.18=0,
    M11.19=0,
    M11.20=0,
    M11.21=0,
    M11.22=0,
    M11.23=0,
    
    M12.17=0,
    M12.18=0,
    M12.19=0,
    M12.20=0,
    M12.21=0,
    M12.22=0,
    M12.23=0,
    
    M13.18=0,
    M13.20=0,
    M13.21=0,
    M13.22=0,
    M13.23=0,
    
    M14.19=0,
    M14.20=0,
    M14.21=0,
    M14.22=0,
    M14.23=0,
    
    M15.20=0,
    M15.21=0,
    M15.22=0,
    M15.23=0,
    
    M16.21=0,
    M16.22=0,
    M16.23=0,
    
    M17.17=0,
    M17.18=0,
    M17.19=0,
    M17.20=0,
    M17.23=0,
    M18.18=0,
    M18.20=0,
    M18.23=0,
    M19.19=0,
    M19.20=0,
    M19.21=0,
    M19.23=0,
    M20.20=0,
    M20.21=0,
    M20.23=0,
    M21.21=0,
    M21.23=0,
    M11.24=0,
    M12.25=0,
    M12.26=0,
    M12.27=0,
    M12.28=0,
    M12.29=0,
    M13.26=0,
    M13.28=0,
    M13.29=0,
    M14.27=0,
    M14.28=0,
    M14.29=0,
    M15.28=0,
    M15.29=0,
    M16.29=0,
    M17.25=0,
    M18.26=0,
    M19.27=0,
    M20.28=0,
    M21.29=0,
    M24.30=0,
    M24.31=0,
    M24.32=0,
    M24.33=0,
    M24.34=0,
    M24.35=0,
    M24.36=0,
    M25.30=0,
    M25.31=0,
    M25.32=0,
    M25.33=0,
    M25.34=0,
    M25.35=0,
    M25.36=0,
    M26.31=0,
    M26.33=0,
    M26.34=0,
    M26.35=0,
    M26.36=0,
    M27.32=0,
    M27.33=0,
    M27.34=0,
    M27.35=0,
    M27.36=0,
    M28.33=0,
    M28.34=0,
    M28.35=0,
    M28.36=0,
    M29.34=0,
    M29.35=0,
    M29.36=0,
    M30.30=0,
    M30.31=0,
    M30.32=0,
    M30.33=0,
    M30.36=0,
    M31.31=0,
    M31.33=0,
    M31.36=0,
    M32.32=0,
    M32.33=0,
    M32.34=0,
    M32.36=0,
    M33.33=0,
    M33.34=0,
    M33.36=0,
    M34.34=0,
    M34.36=0,
    M24.37=0,
    M25.38=0,
    M25.39=0,
    M25.40=0,
    M25.41=0,
    M25.42=0,
    M26.39=0,
    M26.41=0,
    M26.42=0,
    M27.40=0,
    M27.41=0,
    M27.42=0,
    M28.41=0,
    M28.42=0,
    M29.42=0,
    M30.38=0,
    M31.39=0,
    M32.40=0,
    M33.41=0,
    M34.42=0,
    M37.43=0,
    M37.44=0,
    M37.45=0,
    M37.46=0,
    M37.47=0,
    M37.48=0,
    M37.49=0,
    M38.43=0,
    M38.44=0,
    M38.45=0,
    M38.46=0,
    M38.47=0,
    M38.48=0,
    M38.49=0,
    M39.44=0,
    M39.46=0,
    M39.47=0,
    M39.48=0,
    M39.49=0,
    M40.45=0,
    M40.46=0,
    M40.47=0,
    M40.48=0,
    M40.49=0,
    M41.46=0,
    M41.47=0,
    M41.48=0,
    M41.49=0,
    M42.47=0,
    M42.48=0,
    M42.49=0,
    M43.43=0,
    M43.44=0,
    M43.45=0,
    M43.46=0,
    M43.49=0,
    M44.44=0,
    M44.46=0,
    M44.49=0,
    M45.45=0,
    M45.46=0,
    M45.47=0,
    M45.49=0,
    M46.46=0,
    M46.47=0,
    M46.49=0,
    M47.47=0,
    M47.49=0,
    M2.50=0,
    M3.50=0,
    M4.50=0,
    M11.50=0,
    M12.50=0,
    M13.50=0,
    M17.50=0,
    M18.50=0,
    M24.50=0,
    M25.50=0,
    M26.50=0,
    M30.50=0,
    M31.50=0,
    M37.50=0,
    M38.50=0,
    M39.50=0,
    M43.50=0,
    M44.50=0,
    M2.51=0,
    M5.51=0,
    M6.51=0,
    M11.51=0,
    M12.51=0,
    M13.51=0,
    M14.51=0,
    M15.51=0,
    M19.51=0,
    M20.51=0,
    M24.51=0,
    M25.51=0,
    M26.51=0,
    M27.51=0,
    M28.51=0,
    M32.51=0,
    M33.51=0,
    M37.51=0,
    M38.51=0,
    M39.51=0,
    M40.51=0,
    M41.51=0,
    M45.51=0,
    M46.51=0
  )]
  
  #calculates disease-free health state for 40 year old pop
  #moves3[Age==40,`:=`(s.1=newpop-s.3-s.4-s.5-s.6-s.7)]
  
  
  #merges in transition probabilities
  moves3<-merge(moves3, tps_S, by=c("Age","iter"))
  
  #outputs the results for the health states for this model cycle (year)
  res3[[paste0(y)]]<-copy(moves3)
  
  
  #takes out the M vars so that these can be re-calculated in the next cycle
  moves3[,`:=`(M1.1=NULL,
               M1.2=NULL,
               M1.8=NULL,
               M2.3=NULL,
               M2.4=NULL,
               M2.5=NULL,
               M2.6=NULL,
               M2.7=NULL,
               M2.9=NULL,
               M2.10=NULL,
               M2.11=NULL,
               M3.3=NULL,
               M3.4=NULL,
               M3.5=NULL,
               M3.6=NULL,
               M3.9=NULL,
               M3.12=NULL,
               M4.4=NULL,
               M4.6=NULL,
               M4.9=NULL,
               M4.13=NULL,
               M5.5=NULL,
               M5.6=NULL,
               M5.7=NULL,
               M5.9=NULL,
               M5.14=NULL,
               M6.6=NULL,
               M6.7=NULL,
               M6.9=NULL,
               M6.15=NULL,
               M7.7=NULL,
               M7.9=NULL,
               M7.16=NULL,
               M11.17=NULL,
               M11.18=NULL,
               M11.19=NULL,
               M11.20=NULL,
               M11.21=NULL,
               M11.22=NULL,
               M11.23=NULL,
               
               M12.17=NULL,
               M12.18=NULL,
               M12.19=NULL,
               M12.20=NULL,
               M12.21=NULL,
               M12.22=NULL,
               M12.23=NULL,
               
               M13.18=NULL,
               M13.20=NULL,
               M13.21=NULL,
               M13.22=NULL,
               M13.23=NULL,
               
               M14.19=NULL,
               M14.20=NULL,
               M14.21=NULL,
               M14.22=NULL,
               M14.23=NULL,
               
               M15.20=NULL,
               M15.21=NULL,
               M15.22=NULL,
               M15.23=NULL,
               
               M16.21=NULL,
               M16.22=NULL,
               M16.23=NULL,
               
               M17.17=NULL,
               M17.18=NULL,
               M17.19=NULL,
               M17.20=NULL,
               M17.23=NULL,
               M18.18=NULL,
               M18.20=NULL,
               M18.23=NULL,
               M19.19=NULL,
               M19.20=NULL,
               M19.21=NULL,
               M19.23=NULL,
               M20.20=NULL,
               M20.21=NULL,
               M20.23=NULL,
               M21.21=NULL,
               M21.23=NULL,
               M11.24=NULL,
               M12.25=NULL,
               M12.26=NULL,
               M12.27=NULL,
               M12.28=NULL,
               M12.29=NULL,
               M13.26=NULL,
               M13.28=NULL,
               M13.29=NULL,
               M14.27=NULL,
               M14.28=NULL,
               M14.29=NULL,
               M15.28=NULL,
               M15.29=NULL,
               M16.29=NULL,
               M17.25=NULL,
               M18.26=NULL,
               M19.27=NULL,
               M20.28=NULL,
               M21.29=NULL,
               M24.30=NULL,
               M24.31=NULL,
               M24.32=NULL,
               M24.33=NULL,
               M24.34=NULL,
               M24.35=NULL,
               M24.36=NULL,
               M25.30=NULL,
               M25.31=NULL,
               M25.32=NULL,
               M25.33=NULL,
               M25.34=NULL,
               M25.35=NULL,
               M25.36=NULL,
               M26.31=NULL,
               M26.33=NULL,
               M26.34=NULL,
               M26.35=NULL,
               M26.36=NULL,
               M27.32=NULL,
               M27.33=NULL,
               M27.34=NULL,
               M27.35=NULL,
               M27.36=NULL,
               M28.33=NULL,
               M28.34=NULL,
               M28.35=NULL,
               M28.36=NULL,
               M29.34=NULL,
               M29.35=NULL,
               M29.36=NULL,
               M30.30=NULL,
               M30.31=NULL,
               M30.32=NULL,
               M30.33=NULL,
               M30.36=NULL,
               M31.31=NULL,
               M31.33=NULL,
               M31.36=NULL,
               M32.32=NULL,
               M32.33=NULL,
               M32.34=NULL,
               M32.36=NULL,
               M33.33=NULL,
               M33.34=NULL,
               M33.36=NULL,
               M34.34=NULL,
               M34.36=NULL,
               M24.37=NULL,
               M25.38=NULL,
               M25.39=NULL,
               M25.40=NULL,
               M25.41=NULL,
               M25.42=NULL,
               M26.39=NULL,
               M26.41=NULL,
               M26.42=NULL,
               M27.40=NULL,
               M27.41=NULL,
               M27.42=NULL,
               M28.41=NULL,
               M28.42=NULL,
               M29.42=NULL,
               M30.38=NULL,
               M31.39=NULL,
               M32.40=NULL,
               M33.41=NULL,
               M34.42=NULL,
               M37.43=NULL,
               M37.44=NULL,
               M37.45=NULL,
               M37.46=NULL,
               M37.47=NULL,
               M37.48=NULL,
               M37.49=NULL,
               M38.43=NULL,
               M38.44=NULL,
               M38.45=NULL,
               M38.46=NULL,
               M38.47=NULL,
               M38.48=NULL,
               M38.49=NULL,
               M39.44=NULL,
               M39.46=NULL,
               M39.47=NULL,
               M39.48=NULL,
               M39.49=NULL,
               M40.45=NULL,
               M40.46=NULL,
               M40.47=NULL,
               M40.48=NULL,
               M40.49=NULL,
               M41.46=NULL,
               M41.47=NULL,
               M41.48=NULL,
               M41.49=NULL,
               M42.47=NULL,
               M42.48=NULL,
               M42.49=NULL,
               M43.43=NULL,
               M43.44=NULL,
               M43.45=NULL,
               M43.46=NULL,
               M43.49=NULL,
               M44.44=NULL,
               M44.46=NULL,
               M44.49=NULL,
               M45.45=NULL,
               M45.46=NULL,
               M45.47=NULL,
               M45.49=NULL,
               M46.46=NULL,
               M46.47=NULL,
               M46.49=NULL,
               M47.47=NULL,
               M47.49=NULL,
               M2.50=NULL,
               M3.50=NULL,
               M4.50=NULL,
               M11.50=NULL,
               M12.50=NULL,
               M13.50=NULL,
               M17.50=NULL,
               M18.50=NULL,
               M24.50=NULL,
               M25.50=NULL,
               M26.50=NULL,
               M30.50=NULL,
               M31.50=NULL,
               M37.50=NULL,
               M38.50=NULL,
               M39.50=NULL,
               M43.50=NULL,
               M44.50=NULL,
               M2.51=NULL,
               M5.51=NULL,
               M6.51=NULL,
               M11.51=NULL,
               M12.51=NULL,
               M13.51=NULL,
               M14.51=NULL,
               M15.51=NULL,
               M19.51=NULL,
               M20.51=NULL,
               M24.51=NULL,
               M25.51=NULL,
               M26.51=NULL,
               M27.51=NULL,
               M28.51=NULL,
               M32.51=NULL,
               M33.51=NULL,
               M37.51=NULL,
               M38.51=NULL,
               M39.51=NULL,
               M40.51=NULL,
               M41.51=NULL,
               M45.51=NULL,
               M46.51=NULL
               
               
  )]
  
}  

##END OF LOOP

#binds all of the results for each year together
res3<-rbindlist(res3, use.names=T, fill=T)

#same for the move traces table
#moves_trace<-rbindlist(moves_trace, use.names=T, fill=T)



#takes out the move values - only need health states
res3[,`:=`(M1.1=NULL,
           M1.2=NULL,
           M1.8=NULL,
           M2.3=NULL,
           M2.4=NULL,
           M2.5=NULL,
           M2.6=NULL,
           M2.7=NULL,
           M2.9=NULL,
           M2.10=NULL,
           M2.11=NULL,
           M3.3=NULL,
           M3.4=NULL,
           M3.5=NULL,
           M3.6=NULL,
           M3.9=NULL,
           M3.12=NULL,
           M4.4=NULL,
           M4.6=NULL,
           M4.9=NULL,
           M4.13=NULL,
           M5.5=NULL,
           M5.6=NULL,
           M5.7=NULL,
           M5.9=NULL,
           M5.14=NULL,
           M6.6=NULL,
           M6.7=NULL,
           M6.9=NULL,
           M6.15=NULL,
           M7.7=NULL,
           M7.9=NULL,
           M7.16=NULL,
           M11.17=NULL,
           M11.18=NULL,
           M11.19=NULL,
           M11.20=NULL,
           M11.21=NULL,
           M11.22=NULL,
           M11.23=NULL,
           
           M12.17=NULL,
           M12.18=NULL,
           M12.19=NULL,
           M12.20=NULL,
           M12.21=NULL,
           M12.22=NULL,
           M12.23=NULL,
           
           M13.18=NULL,
           M13.20=NULL,
           M13.21=NULL,
           M13.22=NULL,
           M13.23=NULL,
           
           M14.19=NULL,
           M14.20=NULL,
           M14.21=NULL,
           M14.22=NULL,
           M14.23=NULL,
           
           M15.20=NULL,
           M15.21=NULL,
           M15.22=NULL,
           M15.23=NULL,
           
           M16.21=NULL,
           M16.22=NULL,
           M16.23=NULL,
           
           M17.17=NULL,
           M17.18=NULL,
           M17.19=NULL,
           M17.20=NULL,
           M17.23=NULL,
           M18.18=NULL,
           M18.20=NULL,
           M18.23=NULL,
           M19.19=NULL,
           M19.20=NULL,
           M19.21=NULL,
           M19.23=NULL,
           M20.20=NULL,
           M20.21=NULL,
           M20.23=NULL,
           M21.21=NULL,
           M21.23=NULL,
           M11.24=NULL,
           M12.25=NULL,
           M12.26=NULL,
           M12.27=NULL,
           M12.28=NULL,
           M12.29=NULL,
           M13.26=NULL,
           M13.28=NULL,
           M13.29=NULL,
           M14.27=NULL,
           M14.28=NULL,
           M14.29=NULL,
           M15.28=NULL,
           M15.29=NULL,
           M16.29=NULL,
           M17.25=NULL,
           M18.26=NULL,
           M19.27=NULL,
           M20.28=NULL,
           M21.29=NULL,
           M24.30=NULL,
           M24.31=NULL,
           M24.32=NULL,
           M24.33=NULL,
           M24.34=NULL,
           M24.35=NULL,
           M24.36=NULL,
           M25.30=NULL,
           M25.31=NULL,
           M25.32=NULL,
           M25.33=NULL,
           M25.34=NULL,
           M25.35=NULL,
           M25.36=NULL,
           M26.31=NULL,
           M26.33=NULL,
           M26.34=NULL,
           M26.35=NULL,
           M26.36=NULL,
           M27.32=NULL,
           M27.33=NULL,
           M27.34=NULL,
           M27.35=NULL,
           M27.36=NULL,
           M28.33=NULL,
           M28.34=NULL,
           M28.35=NULL,
           M28.36=NULL,
           M29.34=NULL,
           M29.35=NULL,
           M29.36=NULL,
           M30.30=NULL,
           M30.31=NULL,
           M30.32=NULL,
           M30.33=NULL,
           M30.36=NULL,
           M31.31=NULL,
           M31.33=NULL,
           M31.36=NULL,
           M32.32=NULL,
           M32.33=NULL,
           M32.34=NULL,
           M32.36=NULL,
           M33.33=NULL,
           M33.34=NULL,
           M33.36=NULL,
           M34.34=NULL,
           M34.36=NULL,
           M24.37=NULL,
           M25.38=NULL,
           M25.39=NULL,
           M25.40=NULL,
           M25.41=NULL,
           M25.42=NULL,
           M26.39=NULL,
           M26.41=NULL,
           M26.42=NULL,
           M27.40=NULL,
           M27.41=NULL,
           M27.42=NULL,
           M28.41=NULL,
           M28.42=NULL,
           M29.42=NULL,
           M30.38=NULL,
           M31.39=NULL,
           M32.40=NULL,
           M33.41=NULL,
           M34.42=NULL,
           M37.43=NULL,
           M37.44=NULL,
           M37.45=NULL,
           M37.46=NULL,
           M37.47=NULL,
           M37.48=NULL,
           M37.49=NULL,
           M38.43=NULL,
           M38.44=NULL,
           M38.45=NULL,
           M38.46=NULL,
           M38.47=NULL,
           M38.48=NULL,
           M38.49=NULL,
           M39.44=NULL,
           M39.46=NULL,
           M39.47=NULL,
           M39.48=NULL,
           M39.49=NULL,
           M40.45=NULL,
           M40.46=NULL,
           M40.47=NULL,
           M40.48=NULL,
           M40.49=NULL,
           M41.46=NULL,
           M41.47=NULL,
           M41.48=NULL,
           M41.49=NULL,
           M42.47=NULL,
           M42.48=NULL,
           M42.49=NULL,
           M43.43=NULL,
           M43.44=NULL,
           M43.45=NULL,
           M43.46=NULL,
           M43.49=NULL,
           M44.44=NULL,
           M44.46=NULL,
           M44.49=NULL,
           M45.45=NULL,
           M45.46=NULL,
           M45.47=NULL,
           M45.49=NULL,
           M46.46=NULL,
           M46.47=NULL,
           M46.49=NULL,
           M47.47=NULL,
           M47.49=NULL,
           M2.50=NULL,
           M3.50=NULL,
           M4.50=NULL,
           M11.50=NULL,
           M12.50=NULL,
           M13.50=NULL,
           M17.50=NULL,
           M18.50=NULL,
           M24.50=NULL,
           M25.50=NULL,
           M26.50=NULL,
           M30.50=NULL,
           M31.50=NULL,
           M37.50=NULL,
           M38.50=NULL,
           M39.50=NULL,
           M43.50=NULL,
           M44.50=NULL,
           M2.51=NULL,
           M5.51=NULL,
           M6.51=NULL,
           M11.51=NULL,
           M12.51=NULL,
           M13.51=NULL,
           M14.51=NULL,
           M15.51=NULL,
           M19.51=NULL,
           M20.51=NULL,
           M24.51=NULL,
           M25.51=NULL,
           M26.51=NULL,
           M27.51=NULL,
           M28.51=NULL,
           M32.51=NULL,
           M33.51=NULL,
           M37.51=NULL,
           M38.51=NULL,
           M39.51=NULL,
           M40.51=NULL,
           M41.51=NULL,
           M45.51=NULL,
           M46.51=NULL
           
)]


####17. Results ####

#sets results as data table
setDT(res)


#repeats for results 2

#sets results as data table
setDT(res2)


#sets results as data table
setDT(res3)

####calculates model outputs - take out Age<90####

res[Year>=2015,stroke_deaths:=s.10+s.22+s.35+s.48]
res[Year>=2015,nstroke_deaths:=s.9+s.23+s.36+s.49]
res[Year>=2015,NCI:=s.3+s.4+s.12+s.13+s.17+s.18+s.25+s.26+s.30+s.31
    +s.38+s.39+s.43+s.44]
res[Year>=2015,CIND:=s.5+s.6+s.14+s.15+s.19+s.20+s.27+s.28+s.32+s.33
    +s.40+s.41+s.45+s.46]
res[Year>=2015,Dem:=s.7+s.16+s.21+s.29+s.34+s.42+s.47]

res[Year>=2015,stroke_prev:=s.2+NCI+CIND+Dem+s.11+s.24+s.37]
res[Year>=2015,newCIND:=s.50]
res[Year>=2015,newDem:=s.51]

#stroke_inc = first stroke only
res[Year>=2015,stroke_inc:=s.2]

res[Year>=2015 ,LY:=stroke_prev]
res[Year>=2015,DFLY:=s.2+NCI+CIND+s.11+s.24+s.37]
res[Year>=2015,CIFLY:=s.2+NCI+s.11+s.24+s.37]
res[Year>=2015,Rec:=s.11+s.12+s.13+s.14+s.15+s.16+s.24+s.25
    +s.26+s.27+s.28+s.29+s.37+s.38+s.39+s.40+s.41+s.42]

##calculates model outputs for results 2

res2[Year>=2015,stroke_deaths:=s.10+s.22+s.35+s.48]
res2[Year>=2015,nstroke_deaths:=s.9+s.23+s.36+s.49]
res2[Year>=2015,NCI:=s.3+s.4+s.12+s.13+s.17+s.18+s.25+s.26+s.30+s.31
    +s.38+s.39+s.43+s.44]
res2[Year>=2015,CIND:=s.5+s.6+s.14+s.15+s.19+s.20+s.27+s.28+s.32+s.33
    +s.40+s.41+s.45+s.46]
res2[Year>=2015,Dem:=s.7+s.16+s.21+s.29+s.34+s.42+s.47]

res2[Year>=2015,stroke_prev:=s.2+NCI+CIND+Dem+s.11+s.24+s.37]
res2[Year>=2015,newCIND:=s.50]
res2[Year>=2015,newDem:=s.51]

#stroke_inc = first stroke only
res2[Year>=2015,stroke_inc:=s.2]

res2[Year>=2015 ,LY:=stroke_prev]
res2[Year>=2015,DFLY:=s.2+NCI+CIND+s.11+s.24+s.37]
res2[Year>=2015,CIFLY:=s.2+NCI+s.11+s.24+s.37]
res2[Year>=2015,Rec:=s.11+s.12+s.13+s.14+s.15+s.16+s.24+s.25
    +s.26+s.27+s.28+s.29+s.37+s.38+s.39+s.40+s.41+s.42]

##calculates model outputs for results 3

res3[Year>=2015,stroke_deaths:=s.10+s.22+s.35+s.48]
res3[Year>=2015,nstroke_deaths:=s.9+s.23+s.36+s.49]
res3[Year>=2015,NCI:=s.3+s.4+s.12+s.13+s.17+s.18+s.25+s.26+s.30+s.31
     +s.38+s.39+s.43+s.44]
res3[Year>=2015,CIND:=s.5+s.6+s.14+s.15+s.19+s.20+s.27+s.28+s.32+s.33
     +s.40+s.41+s.45+s.46]
res3[Year>=2015,Dem:=s.7+s.16+s.21+s.29+s.34+s.42+s.47]

res3[Year>=2015,stroke_prev:=s.2+NCI+CIND+Dem+s.11+s.24+s.37]
res3[Year>=2015,newCIND:=s.50]
res3[Year>=2015,newDem:=s.51]

#stroke_inc = first stroke only
res3[Year>=2015,stroke_inc:=s.2]

res3[Year>=2015 ,LY:=stroke_prev]
res3[Year>=2015,DFLY:=s.2+NCI+CIND+s.11+s.24+s.37]
res3[Year>=2015,CIFLY:=s.2+NCI+s.11+s.24+s.37]
res3[Year>=2015,Rec:=s.11+s.12+s.13+s.14+s.15+s.16+s.24+s.25
     +s.26+s.27+s.28+s.29+s.37+s.38+s.39+s.40+s.41+s.42]

#Age Groups

res[Age>=40 & Age<50, AgeGroup:="40-49"]
res[Age>=50 & Age<65, AgeGroup:="50-64"]
res[Age>=65 & Age<75, AgeGroup:="65-74"]
res[Age>=75 & Age<90, AgeGroup:="75-89"]
res[Age>=90, AgeGroup:="90+"]

res[Age>=40 & Age<65, AgeGroup2:="40-64"]
res[Age>=65 & Age<75, AgeGroup2:="65-74"]
res[Age>=75 & Age<90, AgeGroup2:="75-89"]
res[Age>=90, AgeGroup2:="90+"]

res2[Age>=40 & Age<50, AgeGroup:="40-49"]
res2[Age>=50 & Age<65, AgeGroup:="50-64"]
res2[Age>=65 & Age<75, AgeGroup:="65-74"]
res2[Age>=75 & Age<90, AgeGroup:="75-89"]
res2[Age>=90, AgeGroup:="90+"]

res2[Age>=40 & Age<65, AgeGroup2:="40-64"]
res2[Age>=65 & Age<75, AgeGroup2:="65-74"]
res2[Age>=75 & Age<90, AgeGroup2:="75-89"]
res2[Age>=90, AgeGroup2:="90+"]

res3[Age>=40 & Age<50, AgeGroup:="40-49"]
res3[Age>=50 & Age<65, AgeGroup:="50-64"]
res3[Age>=65 & Age<75, AgeGroup:="65-74"]
res3[Age>=75 & Age<90, AgeGroup:="75-89"]
res3[Age>=90, AgeGroup:="90+"]

res3[Age>=40 & Age<65, AgeGroup2:="40-64"]
res3[Age>=65 & Age<75, AgeGroup2:="65-74"]
res3[Age>=75 & Age<90, AgeGroup2:="75-89"]
res3[Age>=90, AgeGroup2:="90+"]


####19. Life Table####
#table following the age50 cohort from 2015-2064#

LT50<-res[Year==2015&Age==50]
LT51<-res[Year==2016&Age==51]
LT52<-res[Year==2017&Age==52]
LT53<-res[Year==2018&Age==53]
LT54<-res[Year==2019&Age==54]
LT55<-res[Year==2020&Age==55]
LT56<-res[Year==2021&Age==56]
LT57<-res[Year==2022&Age==57]
LT58<-res[Year==2023&Age==58]
LT59<-res[Year==2024&Age==59]
LT60<-res[Year==2025&Age==60]
LT61<-res[Year==2026&Age==61]
LT62<-res[Year==2027&Age==62]
LT63<-res[Year==2028&Age==63]
LT64<-res[Year==2029&Age==64]
LT65<-res[Year==2030&Age==65]
LT66<-res[Year==2031&Age==66]
LT67<-res[Year==2032&Age==67]
LT68<-res[Year==2033&Age==68]
LT69<-res[Year==2034&Age==69]
LT70<-res[Year==2035&Age==70]
LT71<-res2[Year==2036&Age==71]
LT72<-res2[Year==2037&Age==72]
LT73<-res2[Year==2038&Age==73]
LT74<-res2[Year==2039&Age==74]
LT75<-res2[Year==2040&Age==75]
LT76<-res2[Year==2041&Age==76]
LT77<-res2[Year==2042&Age==77]
LT78<-res2[Year==2043&Age==78]
LT79<-res2[Year==2044&Age==79]
LT80<-res2[Year==2045&Age==80]
LT81<-res2[Year==2046&Age==81]
LT82<-res2[Year==2047&Age==82]
LT83<-res2[Year==2048&Age==83]
LT84<-res2[Year==2049&Age==84]
LT85<-res2[Year==2050&Age==85]
LT86<-res2[Year==2051&Age==86]
LT87<-res2[Year==2052&Age==87]
LT88<-res2[Year==2053&Age==88]
LT89<-res2[Year==2054&Age==89]
LT90<-res2[Year==2055&Age==90]
LT91<-res3[Year==2056&Age==91]
LT92<-res3[Year==2057&Age==92]
LT93<-res3[Year==2058&Age==93]
LT94<-res3[Year==2059&Age==94]
LT95<-res3[Year==2060&Age==95]
LT96<-res3[Year==2061&Age==96]
LT97<-res3[Year==2062&Age==97]
LT98<-res3[Year==2063&Age==98]
LT99<-res3[Year==2064&Age==99]

#selects out from 50 to 90 to create lifetable for age 50
l = list(LT50, LT51, LT52,
         LT53, LT54, LT55,
         LT56, LT57, LT58,
         LT59, LT60, LT61,
         LT62, LT63, LT64,
          LT65, LT66, LT67,
         LT68, LT69, 
         LT70, LT71, LT72, LT73,
         LT74, LT75, LT76, LT77,
         LT78, LT79,
         LT80, LT81, LT82, LT83,
         LT84, LT85, LT86, LT87,
         LT88, LT89, LT90,
         LT91, LT92, LT93,
         LT94, LT95, LT96, LT97,
         LT98, LT99)

rbindlist(l)

lifetable50 <- rbindlist(l)



setkey(lifetable50,iter)



#take TPs out of lifetable

lifetable50[,`:=`(
  P1.1=NULL,
  P1.2=NULL,
  P1.8=NULL,
  P2.3=NULL,
  P2.4=NULL,
  P2.5=NULL,
  P2.6=NULL,
  P2.7=NULL,
  P2.9=NULL,
  P2.10=NULL,
  P2.11=NULL,
  P3.3=NULL,
  P3.4=NULL,
  P3.5=NULL,
  P3.6=NULL,
  P3.9=NULL,
  P3.12=NULL,
  P4.4=NULL,
  P4.6=NULL,
  P4.9=NULL,
  P4.13=NULL,
  P5.5=NULL,
  P5.6=NULL,
  P5.7=NULL,
  P5.9=NULL,
  P5.14=NULL,
  P6.6=NULL,
  P6.7=NULL,
  P6.9=NULL,
  P6.15=NULL,
  P7.7=NULL,
  P7.9=NULL,
  P7.16=NULL,
  P11.17=NULL,
  P11.18=NULL,
  P11.19=NULL,
  P11.20=NULL,
  P11.21=NULL,
  P11.22=NULL,
  P11.23=NULL,
  P12.17=NULL,
  P12.18=NULL,
  P12.19=NULL,
  P12.20=NULL,
  P12.21=NULL,
  P12.22=NULL,
  P12.23=NULL,
  P13.18=NULL,
  P13.20=NULL,
  P13.21=NULL,
  P13.22=NULL,
  P13.23=NULL,
  P14.19=NULL,
  P14.20=NULL,
  P14.21=NULL,
  P14.22=NULL,
  P14.23=NULL,
  P15.20=NULL,
  P15.21=NULL,
  P15.22=NULL,
  P15.23=NULL,
  P16.21=NULL,
  P16.22=NULL,
  P16.23=NULL,
  P17.17=NULL,
  P17.18=NULL,
  P17.19=NULL,
  P17.20=NULL,
  P17.23=NULL,
  P18.18=NULL,
  P18.20=NULL,
  P18.23=NULL,
  P19.19=NULL,
  P19.20=NULL,
  P19.21=NULL,
  P19.23=NULL,
  P20.20=NULL,
  P20.21=NULL,
  P20.23=NULL,
  P21.21=NULL,
  P21.23=NULL,
  P11.24=NULL,
  P12.25=NULL,
  P12.26=NULL,
  P12.27=NULL,
  P12.28=NULL,
  P12.29=NULL,
  P13.26=NULL,
  P13.28=NULL,
  P13.29=NULL,
  P14.27=NULL,
  P14.28=NULL,
  P14.29=NULL,
  P15.28=NULL,
  P15.29=NULL,
  P16.29=NULL,
  P17.25=NULL,
  P18.26=NULL,
  P19.27=NULL,
  P20.28=NULL,
  P21.29=NULL,
  P24.30=NULL,
  P24.31=NULL,
  P24.32=NULL,
  P24.33=NULL,
  P24.34=NULL,
  P24.35=NULL,
  P24.36=NULL,
  P25.30=NULL,
  P25.31=NULL,
  P25.32=NULL,
  P25.33=NULL,
  P25.34=NULL,
  P25.35=NULL,
  P25.36=NULL,
  P26.31=NULL,
  P26.33=NULL,
  P26.34=NULL,
  P26.35=NULL,
  P26.36=NULL,
  P27.32=NULL,
  P27.33=NULL,
  P27.34=NULL,
  P27.35=NULL,
  P27.36=NULL,
  P28.33=NULL,
  P28.34=NULL,
  P28.35=NULL,
  P28.36=NULL,
  P29.34=NULL,
  P29.35=NULL,
  P29.36=NULL,
  P30.30=NULL,
  P30.31=NULL,
  P30.32=NULL,
  P30.33=NULL,
  P30.36=NULL,
  P31.31=NULL,
  P31.33=NULL,
  P31.36=NULL,
  P32.32=NULL,
  P32.33=NULL,
  P32.34=NULL,
  P32.36=NULL,
  P33.33=NULL,
  P33.34=NULL,
  P33.36=NULL,
  P34.34=NULL,
  P34.36=NULL,
  P24.37=NULL,
  P25.38=NULL,
  P25.39=NULL,
  P25.40=NULL,
  P25.41=NULL,
  P25.42=NULL,
  P26.39=NULL,
  P26.41=NULL,
  P26.42=NULL,
  P27.40=NULL,
  P27.41=NULL,
  P27.42=NULL,
  P28.41=NULL,
  P28.42=NULL,
  P29.42=NULL,
  P30.38=NULL,
  P31.39=NULL,
  P32.40=NULL,
  P33.41=NULL,
  P34.42=NULL,
  P37.43=NULL,
  P37.44=NULL,
  P37.45=NULL,
  P37.46=NULL,
  P37.47=NULL,
  P37.48=NULL,
  P37.49=NULL,
  P38.43=NULL,
  P38.44=NULL,
  P38.45=NULL,
  P38.46=NULL,
  P38.47=NULL,
  P38.48=NULL,
  P38.49=NULL,
  P39.44=NULL,
  P39.46=NULL,
  P39.47=NULL,
  P39.48=NULL,
  P39.49=NULL,
  P40.45=NULL,
  P40.46=NULL,
  P40.47=NULL,
  P40.48=NULL,
  P40.49=NULL,
  P41.46=NULL,
  P41.47=NULL,
  P41.48=NULL,
  P41.49=NULL,
  P42.47=NULL,
  P42.48=NULL,
  P42.49=NULL,
  P43.43=NULL,
  P43.44=NULL,
  P43.45=NULL,
  P43.46=NULL,
  P43.49=NULL,
  P44.44=NULL,
  P44.46=NULL,
  P44.49=NULL,
  P45.45=NULL,
  P45.46=NULL,
  P45.47=NULL,
  P45.49=NULL,
  P46.46=NULL,
  P46.47=NULL,
  P46.49=NULL,
  P47.47=NULL,
  P47.49=NULL
  
)]

lifetable50_W <- lifetable50


####LE calculation####

#calc of life years lived at that age, by averaging total life years with next year
#see life table methodology https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/healthandlifeexpectancies/methodologies/guidetocalculatingnationallifetables
lifetable50[,LYlived := rollapply(LY, 2, function(x) sum(head(x,2)), fill=0),by= c("iter")]
lifetable50[,LYlived := LYlived/2]

#calc of total number of life years left to live from starting age

##total number of years lived by each year (accumulating each year)
lifetable50[,LYlivedsofar := cumsum(LYlived),by= c("iter")]

#total number of years lived by age 99
totalLY_99 <- lifetable50[Age==99,LYlivedsofar,by= c("iter")]

setnames(totalLY_99, "LYlivedsofar", "LYlivedsofar99")

lifetable50<- merge(lifetable50, totalLY_99, by = c("iter"), allow.cartesian = T )

#expected years of life up to 99 = total life years at 99 - life years lived so far
lifetable50[,LYexpected := LYlivedsofar99-LYlivedsofar]

#calculation individual LE for person in this group

lifetable50[,LE := LYexpected/LY]

#output LE for age 50
LE_50 <- lifetable50[Age==50,LE]

#check results
median(LE_50)

quantile(LE_50, 0.05)
quantile(LE_50, 0.95)

#put into a data table
LE_50 <- unlist(LE_50)
LE_50 <- data.table(iter=1:1000, LE_50=c(LE_50))

####CIFLY calculation####

#Cognitive impairment free life expectancy = expected life years free of any CI, including CIND or dementia

#calc of life years lived at that age, by averaging total life years with next year
#see life table methodology https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/healthandlifeexpectancies/methodologies/guidetocalculatingnationallifetables
lifetable50[,CIFLYlived := rollapply(CIFLY, 2, function(x) sum(head(x,2)), fill=0),by= c("iter")]
lifetable50[,CIFLYlived := CIFLYlived/2]

#calc of total number of life years left to live for this age group

lifetable50[,CIFLYlivedsofar := cumsum(CIFLYlived),by= c("iter")]

totalCIFLY_99 <- lifetable50[Age==99,CIFLYlivedsofar,by= c("iter")]

setnames(totalCIFLY_99, "CIFLYlivedsofar", "CIFLYlivedsofar99")

lifetable50<- merge(lifetable50, totalCIFLY_99, by = c("iter"), allow.cartesian = T )

lifetable50[,CIFLYexpected := CIFLYlivedsofar99-CIFLYlivedsofar]

#calculation individual CIFLE for person in this group

lifetable50[,CIFLE := CIFLYexpected/LY]

#Output CIFLE for age 50
CIFLE_50 <- lifetable50[Age==50,CIFLE]

#check results
median(CIFLE_50)

quantile(CIFLE_50, 0.05)
quantile(CIFLE_50, 0.95)

#put in datatable 
CIFLE_50 <- unlist(CIFLE_50)
CIFLE_50 <- data.table(iter=1:1000, CIFLE_50=c(CIFLE_50))

#merge with LE data
LE_50_W<- merge(LE_50, CIFLE_50, by = c("iter"), allow.cartesian = T )

#calculate percentage of expected life years free of CI
LE_50_W[,pcCIFLY := (CIFLE_50/LE_50)*100]




####DFLY calculation####

#Dementia free life expectancy

#calc of life years lived at that age, by averaging total life years with next year
#see life table methodology https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/healthandlifeexpectancies/methodologies/guidetocalculatingnationallifetables
lifetable50[,DFLYlived := rollapply(DFLY, 2, function(x) sum(head(x,2)), fill=0),by= c("iter")]
lifetable50[,DFLYlived := DFLYlived/2]

#calc of total number of life years left to live for this age group

lifetable50[,DFLYlivedsofar := cumsum(DFLYlived),by= c("iter")]

totalDFLY_99 <- lifetable50[Age==99,DFLYlivedsofar,by= c("iter")]

setnames(totalDFLY_99, "DFLYlivedsofar", "DFLYlivedsofar99")

lifetable50<- merge(lifetable50, totalDFLY_99, by = c("iter"), allow.cartesian = T )

lifetable50[,DFLYexpected := DFLYlivedsofar99-DFLYlivedsofar]

#calculation individual DFLE for person in this group
#Dementia free life years / total life years

lifetable50[,DFLE := DFLYexpected/LY]

#output DFLE for age 50

DFLE_50 <- lifetable50[Age==50,DFLE]

#check results
median(DFLE_50)

quantile(DFLE_50, 0.05)
quantile(DFLE_50, 0.95)

#put in data table
DFLE_50 <- unlist(DFLE_50)
DFLE_50 <- data.table(iter=1:1000, DFLE_50=c(DFLE_50))

#merge with other LE results
LE_50_W<- merge(LE_50_W, DFLE_50, by = c("iter"), allow.cartesian = T )

#calculate percentage of expected life years free of dementia
LE_50_W[,pcDFLY := (DFLE_50/LE_50)*100]



##generate medians and uncertainty intervals

LE_50_W_est <-  LE_50_W[,list(LE_m = median(LE_50),
                                    LE_lr = quantile(LE_50, 0.05),
                                    LE_upr = quantile(LE_50, 0.95),
                                    CIFLE_m = median(CIFLE_50),
                                    CIFLE_lr = quantile(CIFLE_50, 0.05),
                                    CIFLE_upr = quantile(CIFLE_50, 0.95),
                                    pcCIFLY_m = median(pcCIFLY),
                                    pcCIFLY_lr = quantile(pcCIFLY, 0.05),
                                    pcCIFLY_upr = quantile(pcCIFLY, 0.95),
                                    DFLE_m = median(DFLE_50),
                                    DFLE_lr = quantile(DFLE_50, 0.05),
                                    DFLE_upr = quantile(DFLE_50, 0.95),
                                    pcDFLY_m = median(pcDFLY),
                                    pcDFLY_lr = quantile(pcDFLY, 0.05),
                                    pcDFLY_upr = quantile(pcDFLY, 0.95))
                              ]
LE_50_W_est[,Sex:="Women"]
LE_50_W_est[,Age:=50]

setcolorder(LE_50_W_est,c("Sex","Age"))

write.xlsx(LE_50_W_est,"LE_50_W_est_BC.xlsx", rows=1:2,cols=1:17)
