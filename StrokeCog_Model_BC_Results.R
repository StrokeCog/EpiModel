####Overview####

##StrokeCog epidemiological model
##This code file combines the results for Men and Women
##Check that the relevant women and men .r files have been executed
##Check output file names indicating whether base case or SA are correct - e.g. files for Base Case should include _BC ("x_est_BC.xlsx")
##If treatment effect scenario, this should be added to output file name - e.g. scenario testing moderate effect should include ME ("x_est_BC ME.xlsx")


#NB: set working directory- local folder where model inputs and outputs are stored
setwd("FILENAME")


####Total Results - All Ages####


#Create table with results by sex

x_M_est[,Sex:="Men"]
x_W_est[,Sex:="Women"]

x_MF_est <- rbind(x_M_est, x_W_est)

write.xlsx(x_MF_est,"x_MF_est_BC.xlsx", rows=1:40,cols=1:45)


#Create table with results by sex and age group (40-49 and 50-64 combined)

x_M_est_age2[,Sex:="Men"]
x_W_est_age2[,Sex:="Women"]



x_MF_est_age <- rbindlist(list(x_M_est_age2, x_W_est_age2))

write.xlsx(x_MF_est_age,"x_MF_est_age_BC.xlsx", rows=1:120,cols=1:46)




##Create table with full results for total pop

#Adds sex variable to each result all ages result file 
#These files include the full PSA results

x_allages_M[,Sex:="Men"]
x_allages_W[,Sex:="Women"]

#Combines sex-specific results files, binding by row 
x_allages <- rbind(x_allages_M, x_allages_W)



#Sums together results for men and women to get results for total pop, by year and iteration

x_allages <- x_allages[, lapply(.SD, sum, na.rm = TRUE), by = c("Year","iter"), .SDcols=c("stroke_prev",
                                                                                 "stroke_deaths",
                                                                                 "nstroke_deaths",
                                                                                 "stroke_inc",
                                                                                 "NCI",
                                                                                 "CIND",
                                                                                 "Dem",
                                                                                 "newCIND",
                                                                                 "newDem",
                                                                                 "LY",
                                                                                 "DFLY",
                                                                                 "CIFLY",
                                                                                 "Rec")]
x_est<-x_allages[Year>=2016,list(
                                     stroke_deathslr=quantile(stroke_deaths, 0.05),
                                     nstroke_deathslr=quantile(nstroke_deaths, 0.05),
                                     stroke_inclr=quantile(stroke_inc, 0.05),
                                     NCIlr=quantile(NCI, 0.05),
                                     CINDlr=quantile(CIND, 0.05),
                                     Demlr=quantile(Dem, 0.05),
                                     newCINDlr=quantile(newCIND, 0.05),
                                     newDemlr=quantile(newDem, 0.05),
                                     stroke_prevlr=quantile(stroke_prev, 0.05),
                                     LYlr=quantile(LY, 0.05),
                                     DFLYlr=quantile(DFLY, 0.05),
                                     CIFLYlr=quantile(CIFLY, 0.05),
                                     Reclr=quantile(Rec, 0.05),
                                     stroke_deathsupr=quantile(stroke_deaths, 0.95),
                                     nstroke_deathsupr=quantile(nstroke_deaths, 0.95),
                                     stroke_incupr=quantile(stroke_inc, 0.95),
                                     NCIupr=quantile(NCI, 0.95),
                                     CINDupr=quantile(CIND, 0.95),
                                     Demupr=quantile(Dem, 0.95),
                                     newCINDupr=quantile(newCIND, 0.95),
                                     newDemupr=quantile(newDem, 0.95),
                                     stroke_prevupr=quantile(stroke_prev, 0.95),
                                     LYupr=quantile(LY, 0.95),
                                     DFLYupr=quantile(DFLY, 0.95),
                                     CIFLYupr=quantile(CIFLY, 0.95),
                                     Recupr=quantile(Rec, 0.95),
                                     stroke_deathsm=median(stroke_deaths),
                                     nstroke_deathsm=median(nstroke_deaths),
                                     stroke_incm=median(stroke_inc),
                                     NCIm=median(NCI),
                                     CINDm=median(CIND),
                                     Demm=median(Dem),
                                     newCINDm=median(newCIND),
                                     newDemm=median(newDem),
                                     stroke_prevm=median(stroke_prev),
                                     LYm=median(LY),
                                     DFLYm=median(DFLY),
                                     CIFLYm=median(CIFLY),
                                     Recm=median(Rec)
),
by=list(Year)]

#setting the column order
setcolorder(x_est,c("Year","stroke_deathslr","stroke_deathsm","stroke_deathsupr",
                      "nstroke_deathslr", "nstroke_deathsm", "nstroke_deathsupr",
                      "stroke_inclr","stroke_incm","stroke_incupr",
                      "NCIlr", "NCIm", "NCIupr",
                      "CINDlr","CINDm","CINDupr",
                      "Demlr","Demm","Demupr",
                      "newCINDlr","newCINDm","newCINDupr",
                      "newDemlr","newDemm","newDemupr",
                      "stroke_prevlr","stroke_prevm","stroke_prevupr",
                      "LYlr","LYm","LYupr",
                      "DFLYlr","DFLYm","DFLYupr",
                      "CIFLYlr","CIFLYm","CIFLYupr",
                      "Reclr","Recm","Recupr"
))

write.xlsx(x_est,"x_est_BC.xlsx", rows=1:20,cols=1:43)

####Total Results - By Age Group####

##Create table with full results for total pop by age group

#Adds sex variable to each by age group result file 
#These files include the full PSA results

#combine results for all sexes by age group 2 (combines 40-49 and 50-64)

x_M2[,Sex:="Men"]
x_W2[,Sex:="Women"]



#Combines sex-specific results files, binding by row 
x_MW <- rbind(x_M2, x_W2)



#Sums together results for men and women to get results for total pop, by year, age group2 and iteration

x_MW <- x_MW[, lapply(.SD, sum, na.rm = TRUE), by = c("Year","iter","AgeGroup2"), .SDcols=c("stroke_prev",
                                                                                            "stroke_deaths",
                                                                                            "nstroke_deaths",
                                                                                            "stroke_inc",
                                                                                            "NCI",
                                                                                            "CIND",
                                                                                            "Dem",
                                                                                            "newCIND",
                                                                                            "newDem",
                                                                                            "LY",
                                                                                            "DFLY",
                                                                                            "CIFLY",
                                                                                            "Rec")]


x_est_age<-x_MW[Year>=2016,list(
  stroke_deathslr=quantile(stroke_deaths, 0.05),
  nstroke_deathslr=quantile(nstroke_deaths, 0.05),
  stroke_inclr=quantile(stroke_inc, 0.05),
  NCIlr=quantile(NCI, 0.05),
  CINDlr=quantile(CIND, 0.05),
  Demlr=quantile(Dem, 0.05),
  newCINDlr=quantile(newCIND, 0.05),
  newDemlr=quantile(newDem, 0.05),
  stroke_prevlr=quantile(stroke_prev, 0.05),
  LYlr=quantile(LY, 0.05),
  DFLYlr=quantile(DFLY, 0.05),
  CIFLYlr=quantile(CIFLY, 0.05),
  Reclr=quantile(Rec, 0.05),
  stroke_deathsupr=quantile(stroke_deaths, 0.95),
  nstroke_deathsupr=quantile(nstroke_deaths, 0.95),
  stroke_incupr=quantile(stroke_inc, 0.95),
  NCIupr=quantile(NCI, 0.95),
  CINDupr=quantile(CIND, 0.95),
  Demupr=quantile(Dem, 0.95),
  newCINDupr=quantile(newCIND, 0.95),
  newDemupr=quantile(newDem, 0.95),
  stroke_prevupr=quantile(stroke_prev, 0.95),
  LYupr=quantile(LY, 0.95),
  DFLYupr=quantile(DFLY, 0.95),
  CIFLYupr=quantile(CIFLY, 0.95),
  Recupr=quantile(Rec, 0.95),
  stroke_deathsm=median(stroke_deaths),
  nstroke_deathsm=median(nstroke_deaths),
  stroke_incm=median(stroke_inc),
  NCIm=median(NCI),
  CINDm=median(CIND),
  Demm=median(Dem),
  newCINDm=median(newCIND),
  newDemm=median(newDem),
  stroke_prevm=median(stroke_prev),
  LYm=median(LY),
  DFLYm=median(DFLY),
  CIFLYm=median(CIFLY),
  Recm=median(Rec)
),
by=list(Year,AgeGroup2)]





#setting the column order
setcolorder(x_est_age,c("Year","AgeGroup2","stroke_deathslr","stroke_deathsm","stroke_deathsupr",
                        "nstroke_deathslr", "nstroke_deathsm", "nstroke_deathsupr",
                        "stroke_inclr","stroke_incm","stroke_incupr",
                        "NCIlr", "NCIm", "NCIupr",
                        "CINDlr","CINDm","CINDupr",
                        "Demlr","Demm","Demupr",
                        "newCINDlr","newCINDm","newCINDupr",
                        "newDemlr","newDemm","newDemupr",
                        "stroke_prevlr","stroke_prevm","stroke_prevupr",
                        "LYlr","LYm","LYupr",
                        "DFLYlr","DFLYm","DFLYupr",
                        "CIFLYlr","CIFLYm","CIFLYupr",
                        "Reclr","Recm","Recupr"
))

write.xlsx(x_est_age,"x_est_age_BC.xlsx", rows=1:60,cols=1:44)




####Table 2, Age 40-89####


#selection out 2016, 2017 and 2035 results for Table 2

Tab_40_89 <-  x_est[Year %in% c(2016,2017,2035)]

Tab_40_89[,AgeGroup:="40-89"]

#reshapes the data from long to wide, so that results for 2016, 2017 and 2035 are on the same row

Tab_40_89  <- dcast(Tab_40_89, AgeGroup ~ Year, value.var = c("stroke_prevlr","stroke_prevm","stroke_prevupr",
                                                                     "stroke_deathslr","stroke_deathsm","stroke_deathsupr",
                                                                "nstroke_deathslr","nstroke_deathsm","nstroke_deathsupr",
                                                                     "stroke_inclr","stroke_incm","stroke_incupr",
                                                                     "NCIlr","NCIm","NCIupr",
                                                                     "CINDlr","CINDm","CINDupr",
                                                                     "Demlr","Demm","Demupr",
                                                                     "newCINDlr","newCINDm","newCINDupr",
                                                                     "newDemlr","newDemm","newDemupr",
                                                                     "LYlr","LYm","LYupr",
                                                                     "DFLYlr","DFLYm","DFLYupr",
                                                                     "CIFLYlr","CIFLYm","CIFLYupr",
                                                                     "Reclr","Recm","Recupr"))


#setting the column order to align with Table 2

setcolorder(Tab_40_89,c("AgeGroup","stroke_prevm_2016","stroke_prevlr_2016","stroke_prevupr_2016",
                        "stroke_prevm_2035","stroke_prevlr_2035","stroke_prevupr_2035",
                        
                                            "CINDm_2016","CINDlr_2016","CINDupr_2016",
                                            "CINDm_2035","CINDlr_2035","CINDupr_2035",
                                            
                                            "Demm_2016","Demlr_2016","Demupr_2016",
                                            "Demm_2035","Demlr_2035","Demupr_2035",
                                            
                                            "newCINDm_2016","newCINDlr_2016","newCINDupr_2016",
                                            "newCINDm_2035","newCINDlr_2035","newCINDupr_2035",
                                            
                                            "newDemm_2016","newDemlr_2016","newDemupr_2016",
                                            "newDemm_2035","newDemlr_2035","newDemupr_2035",
                                            
                                            
                                            "LYm_2016","LYlr_2016","LYupr_2016",
                                            "LYm_2035","LYlr_2035","LYupr_2035",
                                            
                                            "DFLYm_2016","DFLYlr_2016","DFLYupr_2016",
                                            "DFLYm_2035","DFLYlr_2035","DFLYupr_2035",
                                            
                                            "CIFLYm_2016","CIFLYlr_2016","CIFLYupr_2016",
                                            "CIFLYm_2035","CIFLYlr_2035","CIFLYupr_2035",
                                            
                                            "Recm_2017","Reclr_2017","Recupr_2017",
                                            "Recm_2035","Reclr_2035","Recupr_2035",
                                            
                        "stroke_deathsm_2016","stroke_deathslr_2016","stroke_deathsupr_2016",
                        "stroke_deathsm_2035","stroke_deathslr_2035","stroke_deathsupr_2035",
                        "nstroke_deathsm_2016","nstroke_deathslr_2016",  "nstroke_deathsupr_2016",
                        "nstroke_deathsm_2035","nstroke_deathslr_2035",  "nstroke_deathsupr_2035",
                        "stroke_incm_2016","stroke_inclr_2016","stroke_incupr_2016",
                        "stroke_incm_2035","stroke_inclr_2035","stroke_incupr_2035",
                        "NCIm_2016","NCIlr_2016", "NCIupr_2016",
                        "NCIm_2035","NCIlr_2035", "NCIupr_2035"
                      ))

write.xlsx(Tab_40_89,"Tab_40_89_BC.xlsx", rows=1:2,cols=1:127)

####Percentage increases####


#selection out 2016, 2017 and 2035 results for % increase calculations
#need to go back to full PSA results to calculate % increase estimate for each iteration 

pc_inc <-  x_allages[Year %in% c(2016,2017,2035)]

#reshapes the data from long to wide, so that multiple years are on the same roe
#allows for calculation of % increases for each var by iteration

pc_inc  <- dcast(pc_inc, iter ~ Year, value.var = c("stroke_prev",
                                                    "stroke_deaths",
                                                    "nstroke_deaths",
                                                    "stroke_inc",
                                                    "NCI",
                                                    "CIND",
                                                    "Dem",
                                                    "newCIND",
                                                    "newDem",
                                                    "LY",
                                                    "DFLY",
                                                    "CIFLY",
                                                    "Rec"))


#calculate % increases for each iteration 


pc_inc <- pc_inc

pc_inc[,pc_inc_strokeprev:=(((stroke_prev_2035-stroke_prev_2016)/stroke_prev_2016)*100)]
pc_inc[,pc_inc_stroke_deaths:=(((stroke_deaths_2035-stroke_deaths_2016)/stroke_deaths_2016)*100)]
pc_inc[,pc_inc_nstroke_deaths:=(((nstroke_deaths_2035-nstroke_deaths_2016)/nstroke_deaths_2016)*100)]
pc_inc[,pc_inc_stroke_inc:=(((stroke_inc_2035-stroke_inc_2016)/stroke_inc_2016)*100)]
pc_inc[,pc_inc_CIND:=(((CIND_2035-CIND_2016)/CIND_2016)*100)]
pc_inc[,pc_inc_Dem:=(((Dem_2035-Dem_2016)/Dem_2016)*100)]
pc_inc[,pc_inc_newCIND:=(((newCIND_2035-newCIND_2016)/newCIND_2016)*100)]
pc_inc[,pc_inc_newDem:=(((newDem_2035-newDem_2016)/newDem_2016)*100)]
pc_inc[,pc_inc_LY:=(((LY_2035-LY_2016)/LY_2016)*100)]
pc_inc[,pc_inc_DFLY:=(((DFLY_2035-DFLY_2016)/DFLY_2016)*100)]
pc_inc[,pc_inc_CIFLY:=(((CIFLY_2035-CIFLY_2016)/CIFLY_2016)*100)]

#Recurrence is estimated from 2017
pc_inc[,pc_inc_Rec:=(((Rec_2035-Rec_2017)/Rec_2017)*100)]


pc_inc_est<-pc_inc[,list(
                        pc_inc_stroke_prevlr=quantile(pc_inc_strokeprev, 0.05),
                                  pc_inc_nstroke_deathslr=quantile(pc_inc_nstroke_deaths, 0.05),
                                  pc_inc_stroke_deathslr=quantile(pc_inc_stroke_deaths, 0.05),
                                  pc_inc_stroke_inclr=quantile(pc_inc_stroke_inc, 0.05),
                                  pc_inc_CINDlr=quantile(pc_inc_CIND, 0.05),
                                  pc_inc_Demlr=quantile(pc_inc_Dem, 0.05),
                                  pc_inc_newCINDlr=quantile(pc_inc_newCIND, 0.05),
                                  pc_inc_newDemlr=quantile(pc_inc_newDem, 0.05),
                                  pc_inc_LYlr=quantile(pc_inc_LY, 0.05),
                                  pc_inc_DFLYlr=quantile(pc_inc_DFLY, 0.05),
                                  pc_inc_CIFLYlr=quantile(pc_inc_CIFLY, 0.05),
                                  pc_inc_Reclr=quantile(pc_inc_Rec, 0.05),
                                  pc_inc_stroke_deathsupr=quantile(pc_inc_stroke_deaths, 0.95),
                                  pc_inc_nstroke_deathsupr=quantile(pc_inc_nstroke_deaths, 0.95),
                                  pc_inc_stroke_incupr=quantile(pc_inc_stroke_inc, 0.95),
                                  pc_inc_CINDupr=quantile(pc_inc_CIND, 0.95),
                                  pc_inc_Demupr=quantile(pc_inc_Dem, 0.95),
                                  pc_inc_newCINDupr=quantile(pc_inc_newCIND, 0.95),
                                  pc_inc_newDemupr=quantile(pc_inc_newDem, 0.95),
                                  pc_inc_stroke_prevupr=quantile(pc_inc_strokeprev, 0.95),
                                  pc_inc_LYupr=quantile(pc_inc_LY, 0.95),
                                  pc_inc_DFLYupr=quantile(pc_inc_DFLY, 0.95),
                                  pc_inc_CIFLYupr=quantile(pc_inc_CIFLY, 0.95),
                                  pc_inc_Recupr=quantile(pc_inc_Rec, 0.95),
                                  pc_inc_stroke_deathsm=median(pc_inc_stroke_deaths),
                                  pc_inc_nstroke_deathsm=median(pc_inc_nstroke_deaths),
                                  pc_inc_stroke_incm=median(pc_inc_stroke_inc),
                                  pc_inc_CINDm=median(pc_inc_CIND),
                                  pc_inc_Demm=median(pc_inc_Dem),
                                  pc_inc_newCINDm=median(pc_inc_newCIND),
                                  pc_inc_newDemm=median(pc_inc_newDem),
                                  pc_inc_stroke_prevm=median(pc_inc_strokeprev),
                                  pc_inc_LYm=median(pc_inc_LY),
                                  pc_inc_DFLYm=median(pc_inc_DFLY),
                                  pc_inc_CIFLYm=median(pc_inc_CIFLY),
                                  pc_inc_Recm=median(pc_inc_Rec)),by=list()
]


#Add in Age Group var

pc_inc_est[,AgeGroup:="40-89"]

#setting the column order
setcolorder(pc_inc_est,c("AgeGroup","pc_inc_stroke_prevm","pc_inc_stroke_prevlr","pc_inc_stroke_prevupr",
                         "pc_inc_CINDm","pc_inc_CINDlr","pc_inc_CINDupr",
                         "pc_inc_Demm","pc_inc_Demlr","pc_inc_Demupr",
                         "pc_inc_newCINDm","pc_inc_newCINDlr","pc_inc_newCINDupr",
                         "pc_inc_newDemm","pc_inc_newDemlr","pc_inc_newDemupr",
                         
                         "pc_inc_stroke_deathsm","pc_inc_stroke_deathslr","pc_inc_stroke_deathsupr",
                         "pc_inc_nstroke_deathsm", "pc_inc_nstroke_deathslr","pc_inc_nstroke_deathsupr",
                         "pc_inc_stroke_incm", "pc_inc_stroke_inclr","pc_inc_stroke_incupr",
                         "pc_inc_LYm","pc_inc_LYlr","pc_inc_LYupr",
                         "pc_inc_DFLYm","pc_inc_DFLYlr","pc_inc_DFLYupr",
                         "pc_inc_CIFLYm","pc_inc_CIFLYlr","pc_inc_CIFLYupr",
                         "pc_inc_Recm","pc_inc_Reclr","pc_inc_Recupr"
))



#writes out to excel
write.xlsx(pc_inc_est,"pc_inc_est_BC.xlsx", rows=1:2,cols=1:41)



####Table 2 - By Age Group####


#selection out 2016, 2017 and 2035 results for Table 2

Tab_AgeGroup <-  x_est_age[Year %in% c(2016,2017,2035)]


#reshapes the data from long to wide, so that results for 2016, 2017 and 2035 are on the same row

Tab_AgeGroup   <- dcast(Tab_AgeGroup , AgeGroup2 ~ Year, value.var = c("stroke_prevlr","stroke_prevm","stroke_prevupr",
                                                              "stroke_deathslr","stroke_deathsm","stroke_deathsupr",
                                                              "nstroke_deathslr","nstroke_deathsm","nstroke_deathsupr",
                                                              "stroke_inclr","stroke_incm","stroke_incupr",
                                                              "NCIlr","NCIm","NCIupr",
                                                              "CINDlr","CINDm","CINDupr",
                                                              "Demlr","Demm","Demupr",
                                                              "newCINDlr","newCINDm","newCINDupr",
                                                              "newDemlr","newDemm","newDemupr",
                                                              "LYlr","LYm","LYupr",
                                                              "DFLYlr","DFLYm","DFLYupr",
                                                              "CIFLYlr","CIFLYm","CIFLYupr",
                                                              "Reclr","Recm","Recupr"))


#setting the column order to align with Table 2

setcolorder(Tab_AgeGroup,c("AgeGroup2","stroke_prevm_2016","stroke_prevlr_2016","stroke_prevupr_2016",
                        "stroke_prevm_2035","stroke_prevlr_2035","stroke_prevupr_2035",
                        
                        "CINDm_2016","CINDlr_2016","CINDupr_2016",
                        "CINDm_2035","CINDlr_2035","CINDupr_2035",
                        
                        "Demm_2016","Demlr_2016","Demupr_2016",
                        "Demm_2035","Demlr_2035","Demupr_2035",
                        
                        "newCINDm_2016","newCINDlr_2016","newCINDupr_2016",
                        "newCINDm_2035","newCINDlr_2035","newCINDupr_2035",
                        
                        "newDemm_2016","newDemlr_2016","newDemupr_2016",
                        "newDemm_2035","newDemlr_2035","newDemupr_2035",
                        
                        
                        "LYm_2016","LYlr_2016","LYupr_2016",
                        "LYm_2035","LYlr_2035","LYupr_2035",
                        
                        "DFLYm_2016","DFLYlr_2016","DFLYupr_2016",
                        "DFLYm_2035","DFLYlr_2035","DFLYupr_2035",
                        
                        "CIFLYm_2016","CIFLYlr_2016","CIFLYupr_2016",
                        "CIFLYm_2035","CIFLYlr_2035","CIFLYupr_2035",
                        
                        "Recm_2017","Reclr_2017","Recupr_2017",
                        "Recm_2035","Reclr_2035","Recupr_2035",
                        
                        "stroke_deathsm_2016","stroke_deathslr_2016","stroke_deathsupr_2016",
                        "stroke_deathsm_2035","stroke_deathslr_2035","stroke_deathsupr_2035",
                        "nstroke_deathsm_2016","nstroke_deathslr_2016",  "nstroke_deathsupr_2016",
                        "nstroke_deathsm_2035","nstroke_deathslr_2035",  "nstroke_deathsupr_2035",
                        "stroke_incm_2016","stroke_inclr_2016","stroke_incupr_2016",
                        "stroke_incm_2035","stroke_inclr_2035","stroke_incupr_2035",
                        "NCIm_2016","NCIlr_2016", "NCIupr_2016",
                        "NCIm_2035","NCIlr_2035", "NCIupr_2035"))

write.xlsx(Tab_AgeGroup,"Tab_AgeGroup_BC.xlsx", rows=1:3,cols=1:127)

####Percentage increases by age group####


#selection out 2016, 2017 and 2035 results for % increase calculations
#need to go back to full PSA results to calculate % increase estimate for each iteration 

pc_inc_age <-  x_MW [Year %in% c(2016,2017,2035)]

#reshapes the data from long to wide, so that multiple years are on the same roe
#allows for calculation of % increases for each var by iteration

pc_inc_age  <- dcast(pc_inc_age, AgeGroup2 + iter ~ Year, value.var = c("stroke_prev",
                                                    "stroke_deaths",
                                                    "nstroke_deaths",
                                                    "stroke_inc",
                                                    "NCI",
                                                    "CIND",
                                                    "Dem",
                                                    "newCIND",
                                                    "newDem",
                                                    "LY",
                                                    "DFLY",
                                                    "CIFLY",
                                                    "Rec"))


#calculate % increases for each iteration 



pc_inc_age[,pc_inc_strokeprev:=(((stroke_prev_2035-stroke_prev_2016)/stroke_prev_2016)*100)]
pc_inc_age[,pc_inc_stroke_deaths:=(((stroke_deaths_2035-stroke_deaths_2016)/stroke_deaths_2016)*100)]
pc_inc_age[,pc_inc_nstroke_deaths:=(((nstroke_deaths_2035-nstroke_deaths_2016)/nstroke_deaths_2016)*100)]
pc_inc_age[,pc_inc_stroke_inc:=(((stroke_inc_2035-stroke_inc_2016)/stroke_inc_2016)*100)]
pc_inc_age[,pc_inc_CIND:=(((CIND_2035-CIND_2016)/CIND_2016)*100)]
pc_inc_age[,pc_inc_Dem:=(((Dem_2035-Dem_2016)/Dem_2016)*100)]
pc_inc_age[,pc_inc_newCIND:=(((newCIND_2035-newCIND_2016)/newCIND_2016)*100)]
pc_inc_age[,pc_inc_newDem:=(((newDem_2035-newDem_2016)/newDem_2016)*100)]
pc_inc_age[,pc_inc_LY:=(((LY_2035-LY_2016)/LY_2016)*100)]
pc_inc_age[,pc_inc_DFLY:=(((DFLY_2035-DFLY_2016)/DFLY_2016)*100)]
pc_inc_age[,pc_inc_CIFLY:=(((CIFLY_2035-CIFLY_2016)/CIFLY_2016)*100)]

#Recurrence is estimated from 2017
pc_inc_age[,pc_inc_Rec:=(((Rec_2035-Rec_2017)/Rec_2017)*100)]


pc_inc_est_age<-pc_inc_age[,list(
  pc_inc_stroke_prevlr=quantile(pc_inc_strokeprev, 0.05),
  pc_inc_nstroke_deathslr=quantile(pc_inc_nstroke_deaths, 0.05),
  pc_inc_stroke_deathslr=quantile(pc_inc_stroke_deaths, 0.05),
  pc_inc_stroke_inclr=quantile(pc_inc_stroke_inc, 0.05),
  pc_inc_CINDlr=quantile(pc_inc_CIND, 0.05),
  pc_inc_Demlr=quantile(pc_inc_Dem, 0.05),
  pc_inc_newCINDlr=quantile(pc_inc_newCIND, 0.05),
  pc_inc_newDemlr=quantile(pc_inc_newDem, 0.05),
  pc_inc_LYlr=quantile(pc_inc_LY, 0.05),
  pc_inc_DFLYlr=quantile(pc_inc_DFLY, 0.05),
  pc_inc_CIFLYlr=quantile(pc_inc_CIFLY, 0.05),
  pc_inc_Reclr=quantile(pc_inc_Rec, 0.05),
  pc_inc_stroke_deathsupr=quantile(pc_inc_stroke_deaths, 0.95),
  pc_inc_nstroke_deathsupr=quantile(pc_inc_nstroke_deaths, 0.95),
  pc_inc_stroke_incupr=quantile(pc_inc_stroke_inc, 0.95),
  pc_inc_CINDupr=quantile(pc_inc_CIND, 0.95),
  pc_inc_Demupr=quantile(pc_inc_Dem, 0.95),
  pc_inc_newCINDupr=quantile(pc_inc_newCIND, 0.95),
  pc_inc_newDemupr=quantile(pc_inc_newDem, 0.95),
  pc_inc_stroke_prevupr=quantile(pc_inc_strokeprev, 0.95),
  pc_inc_LYupr=quantile(pc_inc_LY, 0.95),
  pc_inc_DFLYupr=quantile(pc_inc_DFLY, 0.95),
  pc_inc_CIFLYupr=quantile(pc_inc_CIFLY, 0.95),
  pc_inc_Recupr=quantile(pc_inc_Rec, 0.95),
  pc_inc_stroke_deathsm=median(pc_inc_stroke_deaths),
  pc_inc_nstroke_deathsm=median(pc_inc_nstroke_deaths),
  pc_inc_stroke_incm=median(pc_inc_stroke_inc),
  pc_inc_CINDm=median(pc_inc_CIND),
  pc_inc_Demm=median(pc_inc_Dem),
  pc_inc_newCINDm=median(pc_inc_newCIND),
  pc_inc_newDemm=median(pc_inc_newDem),
  pc_inc_stroke_prevm=median(pc_inc_strokeprev),
  pc_inc_LYm=median(pc_inc_LY),
  pc_inc_DFLYm=median(pc_inc_DFLY),
  pc_inc_CIFLYm=median(pc_inc_CIFLY),
  pc_inc_Recm=median(pc_inc_Rec))
  ,by=list(AgeGroup2)
]



#setting the column order
setcolorder(pc_inc_est_age,c("AgeGroup2","pc_inc_stroke_prevm","pc_inc_stroke_prevlr","pc_inc_stroke_prevupr",
                         "pc_inc_CINDm","pc_inc_CINDlr","pc_inc_CINDupr",
                         "pc_inc_Demm","pc_inc_Demlr","pc_inc_Demupr",
                         "pc_inc_newCINDm","pc_inc_newCINDlr","pc_inc_newCINDupr",
                         "pc_inc_newDemm","pc_inc_newDemlr","pc_inc_newDemupr",
                         
                         "pc_inc_stroke_deathsm","pc_inc_stroke_deathslr","pc_inc_stroke_deathsupr",
                         "pc_inc_nstroke_deathsm", "pc_inc_nstroke_deathslr","pc_inc_nstroke_deathsupr",
                         "pc_inc_stroke_incm", "pc_inc_stroke_inclr","pc_inc_stroke_incupr",
                         "pc_inc_LYm","pc_inc_LYlr","pc_inc_LYupr",
                         "pc_inc_DFLYm","pc_inc_DFLYlr","pc_inc_DFLYupr",
                         "pc_inc_CIFLYm","pc_inc_CIFLYlr","pc_inc_CIFLYupr",
                         "pc_inc_Recm","pc_inc_Reclr","pc_inc_Recupr"
))



#writes out to excel
write.xlsx(pc_inc_est_age,"pc_inc_est_age_BC.xlsx", rows=1:3,cols=1:40)



####Graphs####

Fig1<-ggplot(x_MF_est,
             aes(x=Year, y=stroke_prevm, group=Sex))+
  geom_point()+geom_line()+
  geom_ribbon(aes(ymin=stroke_prevlr, ymax=stroke_prevupr, group=Sex, colour=Sex,fill=Sex), linetype=2, alpha=0.1)+
  scale_x_continuous(breaks = c(2016,2020,2025,2030,2035))+
  ylab("Prevalent Stroke") +  
  scale_fill_manual(values=c("#009E73", "#CC79A7")) + 
  scale_colour_manual(values=c("#009E73", "#CC79A7")) +
  expand_limits(y = 0)+
  theme_bw()

Fig1

ggsave(filename = "Fig 1 Prev.png")


Fig2<-ggplot(x_MF_est,
             aes(x=Year, y=Demm, group=Sex))+
  geom_point()+geom_line()+
  scale_x_continuous(breaks = c(2016,2020,2025,2030,2035))+
  expand_limits(y = 0)+
  theme_bw()


Fig2 + geom_ribbon(aes(ymin=x_MF_est$Demlr, ymax=x_MF_est$Demupr, group=Sex, colour=Sex,fill=Sex), linetype=2, alpha=0.1) +
  geom_ribbon(aes(ymin=x_MF_est$CINDlr, ymax=x_MF_est$CINDupr, group=Sex, colour=Sex,fill=Sex), linetype=2, alpha=0.1) +
  geom_line(aes(x=x_MF_est$Year, y=x_MF_est$CINDm, group=Sex)) +
  geom_point(aes(x=x_MF_est$Year, y=x_MF_est$CINDm, group=Sex))+ ylab("Prevalent Post-stroke Dementia and CIND")+
  scale_fill_manual(values=c("#009E73", "#CC79A7")) + scale_colour_manual(values=c("#009E73", "#CC79A7"))+
  theme(legend.position="top")

ggsave(filename = "Fig 2 Prev Cog.png")



Fig3<-ggplot(x_MF_est,
             aes(x=Year, y=newDemm, group=Sex))+
  geom_point()+geom_line()+
  scale_x_continuous(breaks = c(2016,2020,2025,2030,2035))+
  expand_limits(y = 0)+
  theme_bw()



Fig3 + geom_ribbon(aes(ymin=x_MF_est$newDemlr, ymax=x_MF_est$newDemupr, group=Sex, colour=Sex,fill=Sex), linetype=2, alpha=0.1) + 
  geom_ribbon(aes(ymin=x_MF_est$newCINDlr, ymax=x_MF_est$newCINDupr, group=Sex, colour=Sex,fill=Sex), linetype=2, alpha=0.1) +
  geom_line(aes(x=x_MF_est$Year, y=x_MF_est$newCINDm, group=Sex)) +
  geom_point(aes(x=x_MF_est$Year, y=x_MF_est$newCINDm, group=Sex)) + ylab("Incident Post-stroke Dementia and CIND") +
  scale_fill_manual(values=c("#009E73", "#CC79A7")) + scale_colour_manual(values=c("#009E73", "#CC79A7"))+
  theme(legend.position="top")

ggsave(filename = "Fig 3 Inc Cog.png")

FigS1<-ggplot(x_MF_est,
              aes(x=Year, y=stroke_deathsm, group=Sex))+
  geom_point()+geom_line()+
  scale_x_continuous(breaks = c(2016,2020,2025,2030,2035))+
  expand_limits(y = 0)+
  theme_bw()



FigS1 + geom_ribbon(aes(ymin=x_MF_est$stroke_deathslr, ymax=x_MF_est$stroke_deathsupr, group=Sex, colour=Sex,fill=Sex), linetype=2, alpha=0.1) + 
  ylab("Stroke Deaths") +
  scale_fill_manual(values=c("#009E73", "#CC79A7")) + scale_colour_manual(values=c("#009E73", "#CC79A7")) 

ggsave(filename = "Fig S1 Deaths.png")

FigS2<-ggplot(x_MF_est,
              aes(x=Year, y=nstroke_deathsm, group=Sex))+
  geom_point()+geom_line()+
  scale_x_continuous(breaks = c(2016,2020,2025,2030,2035))+
  expand_limits(y = 0)+
  theme_bw()



FigS2 + geom_ribbon(aes(ymin=x_MF_est$nstroke_deathslr, ymax=x_MF_est$nstroke_deathsupr, group=Sex, colour=Sex,fill=Sex), linetype=2, alpha=0.1) + 
  ylab("Non-Stroke Deaths (stroke survivors)") +
  scale_fill_manual(values=c("#009E73", "#CC79A7")) + scale_colour_manual(values=c("#009E73", "#CC79A7")) 

ggsave(filename = "Fig S2 Non-Stroke Deaths.png")



FigS3<-ggplot(x_MF_est,
              aes(x=Year, y=stroke_incm, group=Sex))+
  geom_point()+geom_line()+
  scale_x_continuous(breaks = c(2016,2020,2025,2030,2035))+
  expand_limits(y = 0)+
  theme_bw()



FigS3 + geom_ribbon(aes(ymin=x_MF_est$stroke_inclr, ymax=x_MF_est$stroke_incupr, group=Sex, colour=Sex,fill=Sex), linetype=2, alpha=0.1) + 
  ylab("Incident Stroke") +
  scale_fill_manual(values=c("#009E73", "#CC79A7")) + scale_colour_manual(values=c("#009E73", "#CC79A7")) 

ggsave(filename = "Fig S3 Inc Stroke.png")



FigS4<-ggplot(x_MF_est,
              aes(x=Year, y=Recm, group=Sex))+
  geom_point()+geom_line()+
  scale_x_continuous(breaks = c(2016,2020,2025,2030,2035))+
  expand_limits(y = 0)+
  theme_bw()



FigS4 + geom_ribbon(aes(ymin=x_MF_est$Reclr, ymax=x_MF_est$Recupr, group=Sex, colour=Sex,fill=Sex), linetype=2, alpha=0.1) + 
  ylab("Recurrent Stroke")+
  scale_fill_manual(values=c("#009E73", "#CC79A7")) + scale_colour_manual(values=c("#009E73", "#CC79A7")) 

ggsave(filename = "Fig S4 Rec Stroke.png")




####Validation#####

sd_W_est_age[,Sex:="Women"]
sd_W_est[,Sex:="Women"]
sd_M_est_age[,Sex:="Men"]
sd_M_est[,Sex:="Men"]

setnames(sd_W_est_age,"AgeGroup3","Age.Group")
setnames(sd_M_est_age,"AgeGroup3","Age.Group")

sd_W_est[,Age.Group:="40-84"]
sd_M_est[,Age.Group:="40-84"]


sd <- rbind(sd_W_est_age, sd_M_est_age,sd_W_est, sd_M_est)

sd[,Data:="StrokeCog"]

sd <- sd[Year==2017]

#read in CSO mortality data for 2015 - 2017
#by age
val_sd<-as.data.table(read.xlsx(xlsxFile="Val Stroke Mortality.xlsx",
                                      sheet="CSO Stroke Mortality 2015_17",
                                      rows=1:25,
                                      cols=1:4))



val_sd<-val_sd[Year==2017]
val_sd[,Data:="CSO"]

setnames(val_sd,"stroke_deaths_v","stroke_deathsm")

val_sd[,stroke_deathslr:=stroke_deathsm]
val_sd[,stroke_deathsupr:=stroke_deathsm]


val_sd <- rbind(sd, val_sd)

val_sd2017 <- val_sd[Year==2017]

write.xlsx(val_sd2017,"val_sd2017_BC.xlsx", rows=1:18,cols=1:8)

positions = c("40-64","65-74","75-84","40-84")

FigV1<-ggplot(val_sd2017,
              aes(x=Age.Group, y=stroke_deathsm, group=Data,colour=Data,shape=Data))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=stroke_deathslr, ymax=stroke_deathsupr),width=0.3)+
  expand_limits(y = 0)+
  scale_x_discrete(limits=positions)+
  ylab("Stroke Deaths")+
  scale_fill_manual(values=c("#009E73", "#CC79A7")) + scale_colour_manual(values=c("#009E73", "#CC79A7"))+ 
  theme_bw()+
  facet_wrap(~Sex)

FigV1


ggsave(filename = "Fig V1 Deaths.png")





