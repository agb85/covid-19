rm(list = ls())
options(scipen=999)
# set working directory


# file paths (Emily)
setwd("/Users/emilywu883/Documents/CU Anschutz/COVID-19")
scen <- read.csv("./Modeling Team/Current Model/model_params_with_variant.csv")
variantp <- read.csv("./Modeling Team/Data Sets/Input Data/proportionvariantovertime.csv")
epi2 <- read.csv("./Modeling Team/Data Sets/Input Data/CO_EMR_Hosp.csv")

# load packages
library(deSolve)
library(lubridate)
library(dplyr)
library(scales)
library(tidyr)
library(magrittr)
library(tidyverse)

# read in spreadsheets
setwd("C:/Users/buchwala/OneDrive - The University of Colorado Denver/Carlton Lab Folder/COVID19/Modeling_Code")

###Projection scenarios
scen <- read.csv("./model_params_spreadsheet/model_params_with_variant.csv")
###Proportion of variants - needs two columns, one for each variant - Beth to update in share drive weekly?
variantp <- read.csv("C:/Users/buchwala/OneDrive - The University of Colorado Denver/Carlton Lab Folder/COVID19/Data_Figures_Info/New_Variants/proportionvariantovertime_0412.csv", header = TRUE)
####Daily vaccination rate -  to be updated weekly by David J.
vacrate <- read.csv("C:/Users/buchwala/OneDrive - The University of Colorado Denver/Carlton Lab Folder/COVID19/Data_Figures_Info/Vaccination_Immunity/Vaccine Age Distribution Data/daily_vaccination_rates.csv", header = TRUE)

###Process vaccine data from David
vacrate1 <- vacrate %>% select(t, group, immun_gain, immun_loss, pS_mult, hosp_mult, dnh_mult) %>%
  filter(group== "0-19")
vacrate2 <- vacrate %>% select(t, group, immun_gain, immun_loss, pS_mult, hosp_mult, dnh_mult) %>%
  filter(group== "20-39")
vacrate3 <- vacrate %>% select(t, group, immun_gain, immun_loss, pS_mult, hosp_mult, dnh_mult) %>%
  filter(group== "40-64")
vacrate4 <- vacrate %>% select(t, group, immun_gain, immun_loss, pS_mult, hosp_mult, dnh_mult) %>%
  filter(group== "65+")


#####THIS IS NECESSARY TO UPDATE TO MAKE VARIANT SIMS WORK###########################
#########Not needed if not doing variant sims#################
###Update to set current "true" TC value assuming no variant
#lastmag <- 0.6836
#pvar <- 1.4  ##Proportionate estimated increase in infectiousness (weighted average using 2:1 ratio of CA to B.1.1.7)
#theta <- 0.466955447  ##Estimated prevalence of variant at fitday-13, rapid growth

# build SEIR model as an R function
seir1 <- function(t, x, parms) {
  
  with(as.list(c(parms, x)), {
    
    #proportion of variants
    theta1 <- ifelse(t < 352, 0, variantp$b117[[t]]) 
    theta2 <- ifelse(t< 352, 0, variantp$CAvariant[[t]])   
    
    ##Variant impact on hosp/death
    variant <- (pvar1*theta1) + (pvar2*theta2) + (1*(1 - (theta1+theta2)))
    vhosp   <- (v1hosp*theta1) + (v2hosp*theta2) + (1*(1 - (theta1+theta2)))
    vdh     <- (v1dh*theta1) + (v2dh*theta2) + (1*(1 - (theta1+theta2)))
    vdnh    <- (v1dnh*theta1) + (v2dnh*theta2) + (1*(1 - (theta1+theta2)))
    
    
    ef1 <- ifelse(t<t2, mag1, ifelse(t<t2a, mag2, ifelse(t<t3, mag2a, ifelse(t<t3a, mag3, ifelse(t<t4, mag3a, ifelse(t<t5,
           mag4, ifelse(t<t6, mag5, ifelse(t<t6a, mag6, ifelse (t<t6b, mag6a, ifelse(t<t7, mag6b, ifelse(t<t8, mag7,
           ifelse (t<t9, mag8, ifelse(t<t10, mag9, ifelse (t<t11, mag10, ifelse(t<t12,mag11, ifelse(t<t13,mag12,
           ifelse(t<t14, mag13, ifelse(t<t15,mag14, ifelse(t<t16, mag15, ifelse(t<t17, mag16, ifelse(t<t18, mag17,
           ifelse(t<t19, mag18, ifelse(t<t20, mag19, ifelse(t<t21, mag20, ifelse(t<t22, mag21, ifelse(t<t23, mag22,
           ifelse(t<t24, mag23, ifelse(t<t25, mag24, ifelse(t<t26, mag25, ifelse(t<ttraj, mag26, ifelse(t<tproject, ef1_1,  ef1_2)))))))))))))))))))))))))))))))
    
    #Change length of stay over time
    hlos4 <- ifelse(t<190, hlos4, hlos4a)
    hlos3 <- ifelse(t<190, hlos3,  hlos3a)
    hlos2 <- ifelse(t<190, hlos2,  hlos2a)
    hlos1 <- ifelse(t<190, hlos1,  hlos1a)
    
    #Introduce vaccination
    vac1 <- ifelse (t < 341, 0, vacrate1$immun_gain[[t]])
    vac2 <- ifelse (t < 341, 0, vacrate2$immun_gain[[t]])
    vac3 <- ifelse (t < 341, 0, vacrate3$immun_gain[[t]])
    vac4 <- ifelse (t < 341, 0, vacrate4$immun_gain[[t]])
    
    ##Loss of vaccine-derived immunity over time
    vacloss1 <- ifelse (t < 341, 0, vacrate1$immun_loss[[t]])
    vacloss2 <- ifelse (t < 341, 0, vacrate2$immun_loss[[t]])
    vacloss3 <- ifelse (t < 341, 0, vacrate3$immun_loss[[t]])
    vacloss4 <- ifelse (t < 341, 0, vacrate4$immun_loss[[t]])
    
    ##Probability of symptoms under vaccine assumptions
    pS1 <- ifelse (t < 341, pS1, pS1* vacrate1$pS_mult[[t]])
    pS2 <- ifelse (t < 341, pS2, pS2* vacrate2$pS_mult[[t]])
    pS3 <- ifelse (t < 341, pS3, pS3* vacrate3$pS_mult[[t]])
    pS4 <- ifelse (t < 341, pS4, pS4* vacrate4$pS_mult[[t]])
    
    ##Probability of hospitalizaton under vaccine and variant assumptions
    hosp1 <- ifelse (t < 341, hosp1, hosp1*vacrate1$hosp_mult[[t]])
    hosp2 <- ifelse(t < 147, hosp2, ifelse(t < 250, hosp2b, hosp2c*vacrate2$hosp_mult[[t]]))
    hosp3 <- ifelse(t < 147, hosp3, ifelse(t < 250, hosp3b, hosp3c*vacrate3$hosp_mult[[t]])*vhosp)
    hosp4 <- ifelse(t < 147, hosp4, ifelse(t < 250, hosp4b, hosp4c*vacrate4$hosp_mult[[t]])*vhosp)
    
    ##Probability of death under vaccine and variant assumptions
    dnh1 <- ifelse (t < 341, dnh1, dnh1*vacrate1$dnh_mult[[t]])
    dnh2 <- ifelse (t < 341, dnh2, dnh2*vacrate2$dnh_mult[[t]])
    dnh3 <- ifelse (t < 341, dnh3, dnh3*vacrate3$dnh_mult[[t]]*vdnh)
    dnh4 <- ifelse (t < 341, dnh4, dnh4*vacrate4$dnh_mult[[t]]*vdnh)
    
    dh3 <- ifelse (t < 341, dh3, dh3*vdh)
    dh4 <- ifelse (t < 341, dh4, dh4*vdh)
    
    
    
    dS1  <-    - (I1+I2+I3+I4)*(beta*lambda*S1*(1-ef1)*variant)/N - (beta*S1*(A1+A2+A3+A4)*(1-ef1)*variant)/N + R1*(1/dimmuneI) + RA1*(1/dimmuneA) - vac1*(S1/(n1-(V1+Ih1+D1))) + V1*vacloss1
    dE1  <-    - E1/alpha   + (I1+I2+I3+I4)*(beta*lambda*S1*(1-ef1)*variant)/N + (beta*S1*(A1+A2+A3+A4)*(1-ef1)*variant)/N 
    dI1  <- (E1*pS1)/alpha     - I1*gamma
    dIh1 <-                      I1*gamma*hosp1            - Ih1/hlos1
    dA1  <- (E1*(1-pS1))/alpha - A1*gamma
    dR1  <-                      I1*(gamma*(1-hosp1-dnh1)) +(Ih1/hlos1)*(1-dh1) - R1*(1/dimmuneI)  - vac1*(R1/(n1-(V1+Ih1+D1)))
    dRA1 <-                      A1*gamma                                       - RA1*(1/dimmuneA) - vac1*(RA1/(n1-(V1+Ih1+D1)))
    dV1  <- vac1*((S1+R1+RA1)/(n1-(V1+Ih1+D1))) - V1*vacloss1
    dD1  <-                      I1*gamma*dnh1             +(Ih1/hlos1)*dh1
    
    dS2  <-    - (I1+I2+I3+I4)*(beta*lambda*S2*(1-ef1)*variant)/N - (beta*S2*(A1+A2+A3+A4)*(1-ef1)*variant)/N + R2*(1/dimmuneI) + RA2*(1/dimmuneA) - vac2*(S2/(n2-(V2+Ih2+D2))) + V2*vacloss2
    dE2  <-    - E2/alpha   + (I1+I2+I3+I4)*(beta*lambda*S2*(1-ef1)*variant)/N + (beta*S2*(A1+A2+A3+A4)*(1-ef1)*variant)/N 
    dI2  <- (E2*pS2)/alpha     - I2*gamma
    dIh2 <-                      I2*gamma*hosp2            - Ih2/hlos2
    dA2  <- (E2*(1-pS2))/alpha - A2*gamma
    dR2  <-                      I2*(gamma*(1-hosp2-dnh2)) +(Ih2/hlos2)*(1-dh2) - R2*(1/dimmuneI)  - vac2*(R2/(n2-(V2+Ih2+D2)))
    dRA2 <-                      A2*gamma                                       - RA2*(1/dimmuneA) - vac2*(RA2/(n2-(V2+Ih2+D2)))
    dV2  <- vac2*((S2+R2+RA2)/(n2-(V2+Ih2+D2))) - V2**vacloss2
    dD2  <-                      I2*gamma*dnh2             +(Ih2/hlos2)*dh2
    
    dS3  <-    - (I1+I2+I3+I4)*(beta*lambda*S3*(1-ef1)*variant)/N - (beta*S3*(A1+A2+A3+A4)*(1-ef1)*variant)/N + R3*(1/dimmuneI) + RA3*(1/dimmuneA) - vac3*(S3/(n3-(V3+Ih3+D3))) + V3*vacloss3
    dE3  <-    - E3/alpha   + (I1+I2+I3+I4)*(beta*lambda*S3*(1-ef1)*variant)/N + (beta*S3*(A1+A2+A3+A4)*(1-ef1)*variant)/N 
    dI3  <- (E3*pS3)/alpha     - I3*gamma
    dIh3 <-                      I3*gamma*hosp3            - Ih3/hlos3
    dA3  <- (E3*(1-pS3))/alpha - A3*gamma
    dR3  <-                      I3*(gamma*(1-hosp3-dnh3)) +(Ih3/hlos3)*(1-dh3) - R3*(1/dimmuneI)  - vac3*(R3/(n3-(V3+Ih3+D3)))
    dRA3 <-                      A3*gamma                                       - RA3*(1/dimmuneA) - vac3*(RA3/(n3-(V3+Ih3+D3)))
    dV3  <- vac3*((S3+R3+RA3)/(n3-(V3+Ih3+D3))) - V3**vacloss3
    dD3  <-                      I3*gamma*dnh3             +(Ih3/hlos3)*dh3
    
    dS4  <-    - (I1+I2+I3+I4)*(beta*lambda*S4*(1-ef1)*variant)/N - (beta*S4*(A1+A2+A3+A4)*(1-ef1)*variant)/N + R4*(1/dimmuneI) + RA4*(1/dimmuneA) - vac4*(S4/(n4-(V4+Ih4+D4))) + V4*vacloss4
    dE4  <-    - E4/alpha   + (I1+I2+I3+I4)*(beta*lambda*S4*(1-ef1)*variant)/N + (beta*S4*(A1+A2+A3+A4)*(1-ef1)*variant)/N 
    dI4  <- (E4*pS4)/alpha     - I4*gamma
    dIh4 <-                      I4*gamma*hosp4            - Ih4/hlos4
    dA4  <- (E4*(1-pS4))/alpha - A4*gamma
    dR4  <-                      I4*(gamma*(1-hosp4-dnh4)) +(Ih4/hlos4)*(1-dh4) - R4*(1/dimmuneI)  - vac4*(R4/(n4-(V4+Ih4+D4)))
    dRA4 <-                      A4*gamma                                       - RA4*(1/dimmuneA) - vac4*(RA4/(n4-(V4+Ih4+D4)))
    dV4  <- vac4*((S4+R4+RA4)/(n4-(V4+Ih4+D4))) - V4**vacloss4
    dD4  <-                      I4*gamma*dnh4             +(Ih4/hlos4)*dh4

    
    
    
    return(list(c(dS1, dE1, dI1, dIh1, dA1, dR1, dRA1, dV1, dD1,
                  dS2, dE2, dI2, dIh2, dA2, dR2, dRA2, dV2, dD2,
                  dS3, dE3, dI3, dIh3, dA3, dR3, dRA3, dV3, dD3,
                  dS4, dE4, dI4, dIh4, dA4, dR4, dRA4, dV4, dD4), 
                
                incI = (I1 + I2 + I3 + I4)/9,
                incA = (A1 + A2 + A3 + A4)/9,
                Ih4n = Ih4/2,
                D4n = D4/2,
                Iht = Ih1 + Ih2 + Ih3 + Ih4,
                Rt = R1+D1+R2+D2+R3+D3+R4+D4+RA1+RA2+RA3+RA4,
                Itotal = I1+I2+I3+I4 +A1+A2+A3+A4,
                Etotal = E1 + E2 + E3 + E4,
                Einc = (E1 + E2 + E3 + E4)/4,
                Vt = V1 + V2 + V3 + V4,
                vacel1 = (S1+R1+RA1)/(n1-(V1+Ih1+D1)),
                vacel2 = (S2+R2+RA2)/(n2-(V2+Ih2+D2)),
                vacel3 = (S3+R3+RA3)/(n3-(V3+Ih3+D3)),
                vacel4 = (S4+R4+RA4)/(n4-(V4+Ih4+D4)),
                immune = R1+R2+R3+R4+RA1+RA2+RA3+RA4+V1 + V2 + V3 + V4))

  })
}

#temptheory <- read.csv('./temptheory.csv')

# rows (n) to represent scenario numbers
n <- as.numeric(nrow(scen)) 
covid_ts <- list() # empty data frame to hold the time series data

# run simulations from time 1 to 500, one simulation per scenario row for as many rows as we have

for(i in 1:n){
  # read in parameters from spreadsheet
  parms <- c(beta = scen[i, c('beta')], # transmission rate
             gamma = 1/9,
             alpha = 4,
             Cp = 5840795, 
             n1 = scen[i,c('n1')],
             n2 = scen[i, c('n2')],
             n3 = scen[i, c('n3')],
             n4 = scen[i, c('n4')],
             
             ef1_1 = scen[i,c('ef1_1')],
             ef1_2 = scen[i,c('ef1_2')],
             #ef1_3 = scen[i,c('ef1_3')],
             #ef1_4 = scen[i,c('ef1_4')],
             #ef1_5 = scen[i,c('ef1_5')],
             ef1 = 0,
             #ef2 = 0,
             #ef3 = 0,
             #ef4 = 0,
             dh1 = scen[i,c('dh1')], dh2 = scen[i,c('dh2')], dh3 = scen[i,c('dh3')],dh4 = scen[i,c('dh4')],
             dnh1 = scen[i,c('dnh1')], dnh2 = scen[i,c('dnh2')], dnh3 = scen[i,c('dnh3')],dnh4 = scen[i,c('dnh4')],
             hlos1 = scen[i,c('hlos1')],
             hlos2 = scen[i,c('hlos2')],
             hlos3 = scen[i,c('hlos3')],
             hlos4 = scen[i,c('hlos4')],
             hlos1a = scen[i,c('hlos1a')],
             hlos2a = scen[i,c('hlos2a')],
             hlos3a = scen[i,c('hlos3a')],
             hlos4a = scen[i,c('hlos4a')],
             dimmuneI = scen[i,c('dimmuneI')],
             dimmuneA = scen[i,c('dimmuneA')],
             
             pS1 = scen[i,c('pS1')], ## proportion of infectious individuals symptomatic (0-19)
             pS2 = scen[i,c('pS2')], ## proportion of infectious individuals symptomatic (20-39)
             pS3 = scen[i,c('pS3')], ## proportion of infectious individuals symptomatic (40-64)
             pS4 = scen[i,c('pS4')], ## proportion of infectious individuals symptomatic (65+)
             lambda = scen[i,c('lambda')], ##difference in infectiousness symptomatic/asymptomatic
             hosp1 = scen[i,c('hosp1')], 
             hosp2 = scen[i,c('hosp2')], 
             hosp3 = scen[i,c('hosp3')], 
             hosp4 = scen[i,c('hosp4')], 
             hosp2b = scen[i,c('hosp2b')], 
             hosp3b = scen[i,c('hosp3b')], 
             hosp4b = scen[i,c('hosp4b')], 
             hosp2c = scen[i,c('hosp2c')], 
             hosp3c = scen[i,c('hosp3c')], 
             hosp4c = scen[i,c('hosp4c')], 
             mag1 = scen[i, c('mag1')],
             mag2 = scen[i, c('mag2')],
             mag2a = scen[i, c('mag2a')],
             mag3 = scen[i, c('mag3')],
             mag3a = scen[i, c('mag3a')],
             mag4 = scen[i, c('mag4')],
             mag4a = scen[i, c('mag4a')],
             mag4b = scen[i, c('mag4b')],
             mag5 = scen[i, c('mag5')],
             mag5a = scen[i, c('mag5a')],
             mag5b = scen[i, c('mag5b')],
             mag5c = scen[i, c('mag5c')],
             mag6 = scen[i, c('mag6')],
             mag6a = scen[i, c('mag6a')],
             mag6b = scen[i, c('mag6b')],
             mag6c = scen[i, c('mag6c')],
             mag7 = scen[i, c('mag7')],
             mag8 = scen[i, c('mag8')],
             mag9 = scen[i, c('mag9')],
             mag10 = scen[i, c('mag10')],
             mag11 = scen[i, c('mag11')],
             mag12 = scen[i, c('mag12')],
             mag13 = scen[i, c('mag13')],
             mag14 = scen[i, c('mag14')],
             mag15 = scen[i, c('mag15')],
             mag16 = scen[i, c('mag16')],
             mag17 = scen[i, c('mag17')],
             mag18 = scen[i, c('mag18')],
             mag19 = scen[i, c('mag19')],
             mag20 = scen[i, c('mag20')],
             mag21 = scen[i, c('mag21')],
             mag22 = scen[i, c('mag22')],
             mag23 = scen[i, c('mag23')],
             mag24 = scen[i, c('mag24')],
             mag25 = scen[i, c('mag25')],
             mag26 = scen[i, c('mag26')],
             traj = scen[i, c("traj")],
             t1 = scen[i,c('t1')],
             t2 = scen[i,c('t2')],
             t2a = scen[i,c('t2a')],
             t3 = scen[i,c('t3')],
             t3a = scen[i,c('t3a')],
             t4 = scen[i,c('t4')],
             t4a = scen[i,c('t4a')],
             t5 = scen[i,c('t5')],
             t5a = scen[i,c('t5a')],
             t5b = scen[i,c('t5b')],
             t6 = scen[i,c('t6')],
             t6a = scen[i,c('t6a')],
             t6b = scen[i,c('t6b')],
             t7 = scen[i,c('t7')], 
             t8 = scen[i,c('t8')],
             t9 = scen[i,c('t9')],
             t10 = scen[i,c('t10')],
             t11 = scen[i,c('t11')],
             t12 = scen[i,c('t12')],
             t13 = scen[i,c('t13')],
             t14 = scen[i,c('t14')],
             t15 = scen[i,c('t15')],
             t16 = scen[i,c('t16')],
             t17 = scen[i,c('t17')],
             t18 = scen[i,c('t18')],
             t19 = scen[i,c('t19')],
             t20 = scen[i,c('t20')],
             t21 = scen[i,c('t21')],
             t22 = scen[i,c('t22')],
             t23 = scen[i,c('t23')],
             t24 = scen[i,c('t24')],
             t25 = scen[i,c('t25')],
             t26 = scen[i,c('t26')],
             ttraj = scen[i,c('ttraj')],
             tproject = scen[i,c('tproject')],
             tpa = scen[i,c('tpa')],
             tpb = scen[i,c('tpb')],
             tpc = scen[i,c('tpc')],
             ##Set key parameters for B.1.1.7 variant sims
             pvar1 = 1.5,#Proportionate increase in transmission due to variant (1.5 = B.1.1.7)
             pvar2 = 1.2,#Proportionate increase in transmission due to variant (1.2 = B.1.427/429)
             v1hosp = 1.4, #Proportionate increase in hospitalization due to variant 1 (B117)
             v2hosp = 1, #Proportionate increase in hospitalization due to variant 2 (B1427/429)
             v1dh = 1.15,#Proportionate increase in death among hospitalized individuals due to variant 1
             v1dnh = 1.6,#Proportionate increase in death among non- hospitalized individuals due to variant 1
             v2dh = 1,#Proportionate increase in death among hospitalized individuals due to variant 2
             v2dnh = 1##Proportionate increase in death among non- hospitalized individuals due to variant 2
              )
  
  dt      <- seq(0, 600, 1)
  
  n1 = scen[i,c('n1')]
  n2 = scen[i, c('n2')]
  n3 = scen[i, c('n3')]
  n4 = scen[i, c('n4')]
  
  inits      <- c(S1 = n1 - 1, E1 = 0, I1 = 1, Ih1 = 0, A1 = 0, R1 = 0, RA1 = 0, V1 = 0, D1 = 0,
                  S2 = n2,     E2 = 0, I2 = 0, Ih2 = 0, A2 = 0, R2 = 0, RA2 = 0, V2 = 0, D2 = 0,
                  S3 = n3,     E3 = 0, I3 = 0, Ih3 = 0, A3 = 0, R3 = 0, RA3 = 0, V3 = 0, D3 = 0,
                  S4 = n4,     E4 = 0, I4 = 0, Ih4 = 0, A4 = 0, R4 = 0, RA4 = 0, V4 = 0, D4 = 0)
  N  <- 5840795
  
  
  out <- lsoda(inits, dt, seir1, parms = parms)
  covid_ts[[i]] <- as.matrix(out)
}

#library(dplyr)
all <-  as.data.frame(cbind(rep(1:n, each=601), do.call("rbind", covid_ts)))
all$scenario <- all$V1
all$V1 <- NULL

all.scen <- merge(scen, all, by = "scenario")

#all.scen.temp <- merge(all.scen, temp, by = "time")

# create incrementing date vector of length 500 for all scenarios

all.scen$date <- seq(from = as.Date("2020/1/24"), to = as.Date("2020/1/24") + 600, "days")
all.scen$Ihtn <- all.scen$Ih1 + all.scen$Ih2 + all.scen$Ih3 + all.scen$Ih4n
all.scen$Dtn <- all.scen$D1 + all.scen$D2 + all.scen$D3 + all.scen$D4n
# write to csv (Andrea)
#write.csv(all.scen, './allscenarios.csv', row.names = F)

#write.csv(all.scen, './scen20210125.csv', row.names = F)

write.csv(all.scen, "/Users/emilywu883/Documents/CU Anschutz/COVID-19/Modeling Team/Data Sets/Output Data/allscenarios.csv", row.names = FALSE)

##### TABLE 2 CODE #####

##DATE of hospital Capacity reached 

#hospthresh <- all.scen %>% filter(time > max(epi2$time), Ihtn > 1) %>% select(scenario, scenalpha, date) %>% rename(capdate = date)  
#hospthresh2 <- hospthresh %>% group_by(scenario, scenalpha) %>% slice(1)

#library(plyr) ###Annoying thing, for some reason if this is loaded before it prevents rename from working 
##Date of hospital peak and need at peak
maxh <- all.scen %>% group_by(scenario, scenalpha) %>%
  filter(time > max(epi2$time), time < 526) %>%
  slice(which.max(Ihtn)) %>%
select(scenario, scenalpha, date, Ihtn)

#maxh <- ddply(all.scen, .(scenario, scenalpha), transform, maxh = max(Ihtn))
#detach(package:plyr)
#hmax <- maxh %>% 
  #select(scenario, scenalpha, date, Ihtn, maxh1)
#filter(Ihtn == maxh) %>%
  #rename(maxdate = date) %>%
  #rename(maxhosp = Ihtn)

#table <- filter(merge(hospthresh2, hmax, by = c('scenario', 'scenalpha'), all = TRUE), !grepl("week", scenalpha))
table <- as.data.frame(filter(maxh, !grepl("week", scenalpha)))
table$maxdate <- as.Date(table$date, "%m/%d/%Y")
#table$maxdate <-as.Date(table$maxdate, "%m/%d/%Y")
table$maxhosp <- label_comma(accuracy = 1)(round(table$Ihtn))

#Cumulative infections through June 1st
###THIS WILL NOT WORK ANYMORE - delete - remove from table or use the sum of Einc
#all.scen$cuminf <- all.scen$Rt + all.scen$Itotal+all.scen$Etotal
#t1 <- filter(subset(all.scen, time == 495,
#             select = c(scenario, scenalpha, cuminf)), !grepl("week", scenalpha))
#t1$cuminf <- label_comma()(signif(t1$cuminf, 3))

#Cumulative infections through July 1st
infectJ <- mutate(all.scen, cuminf = cumsum(Einc))
infectJ <- filter(subset(infectJ, time==524,
                         select = c(scenario, scenalpha, cuminf)), !grepl("week", scenalpha))
infectJ$IJ <- infectJ$cuminf

#Cumulative infections through present (used to calculate excess infections from now until July 1st)
infectP <- mutate(all.scen, cuminf = cumsum(Einc))
infectP <- filter(subset(infectP, time==as.numeric(floor_date(today(), "week", 1)-as.Date("2020-01-24")),
                         select = c(scenario, scenalpha, cuminf)), !grepl("week", scenalpha))
infectP$IP <- infectP$cuminf

#Cumulative deaths through July 1st
deathsJ <- filter(subset(all.scen, time == 524,
                           select = c(scenario, scenalpha, Dtn, D1, D2, D3, D4n)), !grepl("week", scenalpha))
deathsJ$DJ <- deathsJ$Dtn

# deaths through present (used to calculate excess deaths from now until June 1st)
deathsP <- filter(subset(all.scen, time == as.numeric(floor_date(today(), "week", 1)-as.Date("2020-01-24")),
                               select = c(scenario, scenalpha, Dtn, D1, D2, D3, D4n)), !grepl("week", scenalpha))
deathsP$DP <- deathsP$Dtn

table2 <- cbind(table, infectJ, infectP, deathsJ, deathsP)
table2$IExcess <- label_comma()(signif(table2$IJ - table2$IP, 3))
table2$DExcess <- label_comma(accuracy = 1)(round(table2$DJ - table2$DP))
table2$D1Excess <- deathsJ$D1 - deathsP$D1
table2$D2Excess <- deathsJ$D2 - deathsP$D2
table2$D12Excess <- label_comma(accuracy = 1)(round(table2$D1Excess + table2$D2Excess))
table2$D3Excess <- label_comma(accuracy = 1)(round(deathsJ$D3 - deathsP$D3))
table2$D4Excess <- label_comma(accuracy = 1)(round(deathsJ$D4n - deathsP$D4n))
table2$maxdate <- format(table2$maxdate, "%m/%d/%Y")
#table2$capdate <- ifelse(as.character(table2$capdate)=="NA", "N/A", format(table2$capdate, "%m/%d/%y"))
#table2$maxdate <- ifelse(table2$maxdate < floor_date(today(), "week", 1), "past", format(table2$maxdate, "%m/%d/%y"))
table2 <- table2[c("scenalpha", "maxdate", "maxhosp", "IExcess", "DExcess", "D12Excess", "D3Excess", "D4Excess")]
table2
write.csv(table2, "/Users/emilywu883/Documents/CU Anschutz/COVID-19/Modeling Team/Data Sets/Output Data/table2.csv", row.names = FALSE)

