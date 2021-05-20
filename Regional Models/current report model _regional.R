# set working directory

setwd("C:/Users/buchwala/OneDrive - The University of Colorado Denver/Covid-private/Rfiles/Regional-Models")

# load packages

library(deSolve)

# read in Colorado age-group population spreadsheet and age parameter spreadsheet

scen <- read.csv('./Model params_regional.csv')

# build SEIR model as an R function

seir1 <- function(t, x, parms) {
  
  with(as.list(c(parms, x)), {
    
    # change over time in efficacy of % mag SD among specific age groups
    
    ef1 <- ifelse(t<t2, mag1, ifelse(t<t2a, mag2, ifelse(t<t3, mag2a, ifelse(t<t3a, mag3, ifelse(t<t4, mag3a, ifelse(t<t5, mag4, ifelse(t<t6, mag5, ifelse(t<t6a, mag6,
           ifelse (t<t6b, mag6a, ifelse(t<t7, mag6b, ifelse(t<t8, mag7, ifelse(t<t9, mag8, ifelse(t<t10, mag9, ifelse(t<t11, mag10, ifelse(t<t12, mag11,
           ifelse(t<t13, mag12, ifelse(t<t14, mag13, ifelse(t<t15, mag14, ifelse(t<t16, mag15, ifelse(t<t17, mag16, ifelse(t<t18, mag17, ifelse(t<t19, mag18,
           ifelse(t<t20, mag19, ifelse(t<ttraj, mag20, ifelse(t<t21, mag21,traj))))))))))))))))))))))))) 
    ef2 <- ef1
    ef3 <- ef1
    ef4 <- ef1
   
    
    dS1  <-    - (I1+I2+I3+I4)*(beta*lambda*S1*(1-ef1))/N - (beta*S1*(A1+A2+A3+A4)*(1-ef1))/N 
    dE1  <-    - E1/alpha   + (I1+I2+I3+I4)*(beta*lambda*S1*(1-ef1))/N + (beta*S1*(A1+A2+A3+A4)*(1-ef1))/N 
    dI1  <- (E1*pS1)/alpha - I1*(gamma) -  I1*pID*CT*kap*pi*om
    dII1 <-                         (I1+A1)*pID*CT*kap*pi*om - II1*gamma
    dIh1 <- I1*hosp1*gamma + II1*pS1*hosp1*gamma - Ih1*1/hlos1
    dIc1 <- I1*cc1*gamma   + II1*pS1*cc1*gamma- Ic1*(1/clos1) 
    dA1  <- (E1*(1-pS1))/alpha - A1*gamma - A1*pID*CT*kap*pi*om
    dR1  <- (I1+II1*pS1)*(gamma*(1-hosp1-cc1-dnh1)) + A1*gamma 
    dRh1 <- (1-dh1)*Ih1*1/hlos1
    dRc1 <- (1-dc1)*Ic1*1/clos1
    dD1  <-     dc1*Ic1*(1/clos1) + dh1*Ih1*1/hlos1+ dnh1*I1*gamma
    
    dS2  <-    - (I1+I2+I3+I4)*(beta*lambda*S2*(1-ef2))/N - (beta*S2*(A1+A2+A3+A4)*(1-ef2))/N 
    dE2  <-    - E2/alpha   + (I1+I2+I3+I4)*(beta*lambda*S2*(1-ef2))/N + (beta*S2*(A1+A2+A3+A4)*(1-ef2))/N
    dI2  <- (E2*pS2)/alpha - I2*(gamma) -  I2*pID*CT*kap*pi*om
    dII2 <-                         (I2+A2)*pID*CT*kap*pi*om - II2*gamma
    dIh2 <- I2*hosp2*gamma + II2*pS2*hosp2*gamma - Ih2*1/hlos2
    dIc2 <- I2*cc2*gamma   + II2*pS2*cc2*gamma- Ic2*(1/clos2) 
    dA2  <- (E2*(1-pS2))/alpha - A2*gamma - A2*pID*CT*kap*pi*om
    dR2  <- (I2+II2*pS2)*(gamma*(1-hosp2-cc2-dnh2)) + A2*gamma 
    dRh2 <- (1-dh2)*Ih2*1/hlos2
    dRc2 <- (1-dc2)*Ic2*1/clos2
    dD2  <-     dc2*Ic2*(1/clos2) + dh2*Ih2*(1/hlos2)+ dnh2*I2*gamma
    
    dS3  <-    - (I1+I2+I3+I4)*(beta*lambda*S3*(1-ef3))/N - (beta*S3*(A1+A2+A3+A4)*(1-ef3))/N
    dE3  <-    - E3/alpha   + (I1+I2+I3+I4)*(beta*lambda*S3*(1-ef3))/N + (beta*S3*(A1+A2+A3+A4)*(1-ef3))/N
    dI3  <- (E3*pS3)/alpha - I3*(gamma)  - I3*pID*CT*kap*pi*om
    dII3 <-                         (I3+A3)*pID*CT*kap*pi*om - II3*gamma
    dIh3 <- I3*hosp3*gamma + II3*pS3*hosp3*gamma - Ih3*1/hlos3
    dIc3 <- I3*cc3*gamma   + II3*pS3*cc3*gamma- Ic3*(1/clos3) 
    dA3  <- (E3*(1-pS3))/alpha - A3*gamma - A3*pID*CT*kap*pi*om
    dR3  <- (I3+II3*pS3)*(gamma*(1-hosp3-cc3-dnh3)) + A3*gamma 
    dRh3 <- (1-dh3)*Ih3*1/hlos3
    dRc3 <- (1-dc3)*Ic3*(1/clos3)
    dD3  <-    dc3 *Ic3*(1/clos3) + dh3*Ih3*(1/hlos3) + dnh3*I3*gamma
    
    dS4  <-    - (I1+I2+I3+I4)*(beta*lambda*S4*(1-ef4))/N - (beta*S4*(A1+A2+A3+A4)*(1-ef4))/N
    dE4  <-    - E4/alpha   + (I1+I2+I3+I4)*(beta*lambda*S4*(1-ef4))/N + (beta*S4*(A1+A2+A3+A4)*(1-ef4))/N
    dI4  <- (E4*pS4)/alpha - I4*(gamma)  - I4*pID*CT*kap*pi*om
    dII4 <-                         (I4+A4)*pID*CT*kap*pi*om - II4*gamma
    dIh4 <- I4*hosp4*gamma + II4*pS4*hosp4*gamma - Ih4*1/hlos4
    dIc4 <- I4*cc4*gamma   + II4*pS4*cc4*gamma- Ic4*(1/clos4) 
    dA4  <- (E4*(1-pS4))/alpha - A4*gamma - A4*pID*CT*kap*pi*om
    dR4  <- (I4+II4*pS4)*(gamma*(1-hosp4-cc4-dnh4)) + A4*gamma 
    dRh4 <- (1-dh4)*Ih4*(1/hlos4)
    dRc4 <- (1-dc4)*Ic4*(1/clos4)
    dD4  <-    dc4* Ic4*(1/clos4) + dh4*Ih4*(1/hlos4) + dnh4*I4*gamma
    
    
    der <- c(dS1, dE1, dI1, dII1, dIh1, dIc1, dA1, dR1, dRh1, dRc1, dD1,
             dS2, dE2, dI2, dII2, dIh2, dIc2, dA2, dR2, dRh2, dRc2, dD2,
             dS3, dE3, dI3, dII3, dIh3, dIc3, dA3, dR3, dRh3, dRc3, dD3,
             dS4, dE4, dI4, dII4, dIh4, dIc4, dA4, dR4, dRh4, dRc4, dD4)
    
    list(der,
         incI = (I1 + I2 + I3 + I4)/9,
         incA = (A1 + A2 + A3 + A4)/9,
         Inc = (I1 + I2 + I3 + I4+A1 + A2 + A3 + A4)/9,
         IIt = II1 + II2 + II3 + II4,
         Iht =Ih1 + Ih2 + Ih3 + Ih4 + Ic1 + Ic2 + Ic3 + Ic4, 
         Iht1 = Ih1 +Ic1,
         Iht2 = Ih2 +Ic2,
         Iht3 = Ih3 +Ic3,
         Iht4 = Ih4 +Ic4,
         Ict =Ic1 + Ic2 + Ic3 + Ic4,
         Rt = R1+Rh1+Rc1+D1+R2+Rh2+Rc2+D2+R3+Rh3+Rc3+D3+R4+Rh4+Rc4+D4,
         Itotal = I1+I2+I3+I4 +A1+A2+A3+A4,
         Etotal = E1 + E2 + E3 + E4)
  })
}

#temptheory <- read.csv('./temptheory.csv')

# rows (n) to represent scenario numbers
n <- as.numeric(nrow(scen)) 
covid_ts <- list() # empty data frame to hold the time series data

# run simulations from time 1 to 500, one simulation per scenario row for as many rows as we have

for(i in 1:n){
  # define parameters that will change
  
  parms <- c(beta = scen[i, c('beta')], # transmission rate
             gamma = 1/9,
             alpha = 4,
             N = scen[i, c('Cp')], # called back from population spreadsheet
             n1 = scen[i, c('n1')],
             n2 = scen[i, c('n2')],
             n3 = scen[i, c('n3')],
             n4 = scen[i, c('n4')],
             ef1_1 = scen[i,c('ef1_1')],
             ef1_2 = scen[i,c('ef1_2')],
             ef1_3 = scen[i,c('ef1_3')],
             ef1_4 = scen[i,c('ef1_4')],
             ef1 = 0,
             ef2 = 0,
             ef3 = 0,
             ef4 = 0,
             dh1 = 0, dh2 = 0, dh3 = 0.0045, dh4 = 0.0923,
             dc1 = 0.0417, dc2 = 0.0392, dc3 = 0.1543, dc4 = 0.3956,
             dnh1 = 0.000072, dnh2 = 0.000129, dnh3 = 0.001136, dnh4 = 0.030285,
             pS1 = scen[i,c('pS1')], ## proportion of infectious individuals symptomatic (0-19)
             pS2 = scen[i,c('pS2')], ## proportion of infectious individuals symptomatic (20-39)
             pS3 = scen[i,c('pS3')], ## proportion of infectious individuals symptomatic (40-64)
             pS4 = scen[i,c('pS4')], ## proportion of infectious individuals symptomatic (65+)
             pID = scen[i,c('pID')], ## proportion of infections identified
             lambda = scen[i,c('lambda')], ##difference in infectiousness symptomatic/asymptomatic
             hosp1 = scen[i,c('hosp1')], 
             cc1 = scen[i,c('cc1')],
             hosp2 = scen[i,c('hosp2')], 
             cc2 = scen[i,c('cc2')],
             hosp3 = scen[i,c('hosp3')], 
             cc3 = scen[i,c('cc3')],
             hosp4 = scen[i,c('hosp4')], 
             cc4 = scen[i,c('cc4')],
             hlos1 = scen[i,c('hlos1')],
             hlos2 = scen[i,c('hlos2')],
             hlos3 = scen[i,c('hlos3')],
             hlos4 = scen[i,c('hlos4')],
             clos1 = scen[i,c('clos1')],
             clos2 = scen[i,c('clos2')],
             clos3 = scen[i,c('clos3')],
             clos4 = scen[i,c('clos4')],
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
             mag7a = scen[i, c('mag7a')],
             mag7b = scen[i, c('mag7b')],
             mag7c = scen[i, c('mag7c')],
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
             t6c = scen[i,c('t6c')],
             t6d = scen[i,c('t6d')],
             t7 = scen[i,c('t7')], 
             t7a = scen[i,c('t7a')],
             t7b = scen[i,c('t7b')],
             t7c = scen[i,c('t7c')],
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
             ttraj = scen[i,c('ttraj')],
             tproject = scen[i,c('tproject')],
             tpa = scen[i,c('tpa')],
             tschool = scen[i,c('tschool')],
             kap = 0, #average number of contacts traced per detected case
             pCT = 0, #proportion of identified cases with contacts traced
             pi = 0, #probability a contact traced infected individual is isolated before infecting other susceptibles 
             om = 0, #probability a contact traced individual is infected
             CT = 0
  )
  
  dt      <- seq(0, 500, 1)
  
  N <- scen[i, c('Cp')] # called back from population spreadsheet
  n1 <- scen[i, c('n1')]
  n2 <- scen[i, c('n2')]
  n3 <- scen[i, c('n3')]
  n4 <- scen[i, c('n4')] 
  
  inits      <- c(S1 = n1 - 1, E1 = 0, I1 = 1, II1 = 0, Ih1 = 0, Ic1 = 0, A1 = 0, R1 = 0, Rh1 = 0, Rc1 = 0, D1 = 0,
                  S2 = n2,     E2 = 0, I2 = 0, II2 = 0, Ih2 = 0, Ic2 = 0, A2 = 0, R2 = 0, Rh2 = 0, Rc2 = 0, D2 = 0,
                  S3 = n3,     E3 = 0, I3 = 0, II3 = 0, Ih3 = 0, Ic3 = 0, A3 = 0, R3 = 0, Rh3 = 0, Rc3 = 0, D3 = 0,
                  S4 = n4,     E4 = 0, I4 = 0, II4 = 0, Ih4 = 0, Ic4 = 0, A4 = 0, R4 = 0, Rh4 = 0, Rc4 = 0, D4 = 0)

  
  out <- lsoda(inits, dt, seir1, parms = parms)
  covid_ts[[i]] <- as.matrix(out)
}

#library(dplyr)
all <-  as.data.frame(cbind(rep(1:n, each=501), do.call("rbind", covid_ts)))
all$scenario <- all$V1
all$V1 <- NULL

all.scen <- merge(scen, all, by = "scenario")
#all.scen.temp <- merge(all.scen, temp, by = "time")

# create incrementing date vector of length 500 for all scenarios

all.scen$date <- seq(from = as.Date("2020/1/24"), to = as.Date("2020/1/24") + 500, "days")

write.csv(all.scen, './regional_0125.csv', row.names = F)



