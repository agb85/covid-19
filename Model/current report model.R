rm(list = ls())
# set working directory


# file paths (Emily)
setwd("/Users/emilywu883/covid_private")
pIDstate <- read.csv('/Users/emilywu883/Documents/CU Anschutz/COVID-19/Modeling Team/Data Sets/Input Data/pIDcomparison 10302020.csv', header = TRUE, na.strings =".")
scen <- read.csv('./Model params.csv')


# load packages
library(deSolve)
library(lubridate)

# read in spreadsheet
setwd("C:/Users/buchwala/OneDrive - The University of Colorado Denver/Covid-private/Rfiles/covid_private")
scen <- read.csv('./Model params.csv')

# build SEIR model as an R function
seir1 <- function(t, x, parms) {
  
  with(as.list(c(parms, x)), {
    
   
    ef1 <- ifelse(t<t2, mag1, ifelse(t<t2a, mag2, ifelse(t<t3, mag2a, ifelse(t<t3a, mag3, ifelse(t<t4, mag3a, ifelse(t<t5, mag4, ifelse(t<t6, mag5, ifelse(t<t6a, mag6, ifelse (t<t6b, mag6a, ifelse(t<t7, mag6b, 
           ifelse(t<t8, mag7, ifelse (t<t9, mag8, ifelse(t<t10, mag9, ifelse (t<t11, mag10, ifelse(t<t12,mag11, ifelse(t<t13,mag12, ifelse(t<t14, mag13, ifelse(t<t15,mag14, ifelse(t<t16, mag15, ifelse(t<t17, mag16, 
           ifelse(t<t18, mag17, ifelse(t<t19, mag18, ifelse(t<t20, mag19, ifelse(t<t21, mag20, ifelse(t<t22, mag21, ifelse(t<t23, mag22, ifelse(t<ttraj, mag23, ifelse(t<tproject,traj, ifelse(t<tpa, ef1_2, ifelse (t<tpb, ef1_3, ifelse(t<tpc, ef1_4, ef1_5)))))))))))))))))))))))))))))))
    
    
    hlos4 <- ifelse(t<190, hlos4, hlos4a)
    hlos3 <- ifelse(t<190, hlos3,  hlos3a)
    hlos2 <- ifelse(t<190, hlos2,  hlos2a)
    hlos1 <- ifelse(t<190, hlos1,  hlos1a)
    
    hosp2 <- ifelse(t < 147, hosp2, ifelse(t < 250, hosp2b, hosp2c))
    hosp3 <- ifelse(t < 147, hosp3, ifelse(t < 250, hosp3b, ifelse(t<tproject, hosp3c, ifelse(t<tpa, hosp3v1, ifelse(t<tpb, hosp3v2, ifelse(t<tpc, hosp3v3, hosp3v4))))))
    hosp4 <- ifelse(t < 147, hosp4, ifelse(t < 250, hosp4b, ifelse(t<tproject, hosp4c, ifelse(t<tpa, hosp4v1, ifelse(t<tpb, hosp4v2, ifelse(t<tpc, hosp4v3, hosp4v4))))))
    
    
    dh3 <- ifelse(t<tproject, dh3, ifelse(t<tpa, dh3v1, ifelse(t<tpb, dh3v2, ifelse(t<tpc, dh3v3, dh3v4))))
    dh4 <- ifelse(t<tproject, dh4, ifelse(t<tpa, dh4v1, ifelse(t<tpb, dh4v2, ifelse(t<tpc, dh4v3, dh4v4))))
    
    dnh3 <- ifelse(t<tproject, dnh3, ifelse(t<tpa, dnh3v1, ifelse(t<tpb, dnh3v2, ifelse(t<tpc, dnh3v3, dnh3v4))))
    dnh4 <- ifelse(t<tproject, dnh4, ifelse(t<tpa, dnh4v1, ifelse(t<tpb, dnh4v2, ifelse(t<tpc, dnh4v3, dnh4v4))))
    

    vac1 <- ifelse(t<tv1+14, 0, ifelse(t< tv2+14, vac1_*0.52, ifelse(t<tv3+14, vac1a*0.52, ifelse(t<tv4+14, vac1b*0.52,ifelse(t<(tvacend+14), vac1c*0.52,0)))))
    vac2 <- ifelse(t<tv1+14, 0, ifelse(t< tv2+14, vac2_*0.52, ifelse(t<tv3+14, vac2a*0.52, ifelse(t<tv4+14, vac2b*0.52,ifelse(t<(tvacend+14), vac2c*0.52,0)))))
    vac3 <- ifelse(t<tv1+14, 0, ifelse(t< tv2+14, vac3_*0.52, ifelse(t<tv3+14, vac3a*0.52, ifelse(t<tv4+14, vac3b*0.52,ifelse(t<(tvacend+14), vac3c*0.52,0)))))
    vac4 <- ifelse(t<tv1+14, 0, ifelse(t< tv2+14, vac4_*0.52, ifelse(t<tv3+14, vac4a*0.52, ifelse(t<tv4+14, vac4b*0.52,ifelse(t<(tvacend+14), vac4c*0.52,0)))))
    
    vac1_2 <- ifelse(t<tv1+32, 0, ifelse(t< tv2+32, vac1_*0.38, ifelse(t<tv3+32, vac1a*0.38, ifelse(t<tv4+32, vac1b*0.38, ifelse(t<tvacend+32, vac1c*0.38,0)))))
    vac2_2 <- ifelse(t<tv1+32, 0, ifelse(t< tv2+32, vac2_*0.38, ifelse(t<tv3+32, vac2a*0.38, ifelse(t<tv4+32, vac2b*0.38,ifelse(t<tvacend+32, vac2c*0.38,0)))))
    vac3_2 <- ifelse(t<tv1+32, 0, ifelse(t< tv2+32, vac3_*0.38, ifelse(t<tv3+32, vac3a*0.38, ifelse(t<tv4+32, vac3b*0.38,ifelse(t<tvacend+32, vac3c*0.38,0)))))
    vac4_2 <- ifelse(t<tv1+32, 0, ifelse(t< tv2+32, vac4_*0.38, ifelse(t<tv3+32, vac4a*0.38, ifelse(t<tv4+32, vac4b*0.38,ifelse(t<tvacend+32, vac4c*0.38,0)))))
    
    vj1 <- ifelse(t<tvja+28, 0, ifelse(t<tvjb, vj1*0.72, vj1a*0.72))
    vj2 <- ifelse(t<tvja+28, 0, ifelse(t<tvjb, vj2*0.72, vj2a*0.72))
    vj3 <- ifelse(t<tvja+28, 0, ifelse(t<tvjb, vj3*0.72, vj3a*0.72))
    vj4 <- ifelse(t<tvja+28, 0, vj4*0.72)
    
    
    
    dS1  <-    - (I1+I2+I3+I4)*(beta*lambda*S1*(1-ef1))/N - (beta*S1*(A1+A2+A3+A4)*(1-ef1))/N + R1*(1/dimmuneI) + RA1*(1/dimmuneA) - (vac1+vac1_2+vj1)*(S1/(n1-(V1+Ih1+D1))) + V1*(1/vd)
    dE1  <-    - E1/alpha   + (I1+I2+I3+I4)*(beta*lambda*S1*(1-ef1))/N + (beta*S1*(A1+A2+A3+A4)*(1-ef1))/N 
    dI1  <- (E1*pS1)/alpha     - I1*gamma
    dIh1 <-                      I1*gamma*hosp1            - Ih1/hlos1
    dA1  <- (E1*(1-pS1))/alpha - A1*gamma
    dR1  <-                      I1*(gamma*(1-hosp1-dnh1)) +(Ih1/hlos1)*(1-dh1) - R1*(1/dimmuneI)  - (vac1+vac1_2+vj1)*(R1/(n1-(V1+Ih1+D1)))
    dRA1 <-                      A1*gamma                                       - RA1*(1/dimmuneA) - (vac1+vac1_2+vj1)*(RA1/(n1-(V1+Ih1+D1)))
    dV1  <- (vac1+vac1_2+vj1)*((S1+R1+RA1)/(n1-(V1+Ih1+D1))) - V1*(1/vd)
    dD1  <-                      I1*gamma*dnh1             +(Ih1/hlos1)*dh1
    
    dS2  <-    - (I1+I2+I3+I4)*(beta*lambda*S2*(1-ef1))/N - (beta*S2*(A1+A2+A3+A4)*(1-ef1))/N + R2*(1/dimmuneI) + RA2*(1/dimmuneA) - (vac2+vac2_2+vj2)*(S2/(n2-(V2+Ih2+D2))) + V2*(1/vd)
    dE2  <-    - E2/alpha   + (I1+I2+I3+I4)*(beta*lambda*S2*(1-ef1))/N + (beta*S2*(A1+A2+A3+A4)*(1-ef1))/N 
    dI2  <- (E2*pS2)/alpha     - I2*gamma
    dIh2 <-                      I2*gamma*hosp2            - Ih2/hlos2
    dA2  <- (E2*(1-pS2))/alpha - A2*gamma
    dR2  <-                      I2*(gamma*(1-hosp2-dnh2)) +(Ih2/hlos2)*(1-dh2) - R2*(1/dimmuneI)  - (vac2+vac2_2+vj2)*(R2/(n2-(V2+Ih2+D2)))
    dRA2 <-                      A2*gamma                                       - RA2*(1/dimmuneA) - (vac2+vac2_2+vj2)*(RA2/(n2-(V2+Ih2+D2)))
    dV2  <- (vac2+vac2_2+vj2)*((S2+R2+RA2)/(n2-(V2+Ih2+D2))) - V2*(1/vd)
    dD2  <-                      I2*gamma*dnh2             +(Ih2/hlos2)*dh2
    
    dS3  <-    - (I1+I2+I3+I4)*(beta*lambda*S3*(1-ef1))/N - (beta*S3*(A1+A2+A3+A4)*(1-ef1))/N + R3*(1/dimmuneI) + RA3*(1/dimmuneA) - (vac3+vac3_2+vj3)*(S3/(n3-(V3+Ih3+D3))) + V3*(1/vd)
    dE3  <-    - E3/alpha   + (I1+I2+I3+I4)*(beta*lambda*S3*(1-ef1))/N + (beta*S3*(A1+A2+A3+A4)*(1-ef1))/N 
    dI3  <- (E3*pS3)/alpha     - I3*gamma
    dIh3 <-                      I3*gamma*hosp3            - Ih3/hlos3
    dA3  <- (E3*(1-pS3))/alpha - A3*gamma
    dR3  <-                      I3*(gamma*(1-hosp3-dnh3)) +(Ih3/hlos3)*(1-dh3) - R3*(1/dimmuneI)  - (vac3+vac3_2+vj3)*(R3/(n3-(V3+Ih3+D3)))
    dRA3 <-                      A3*gamma                                       - RA3*(1/dimmuneA) - (vac3+vac3_2+vj3)*(RA3/(n3-(V3+Ih3+D3)))
    dV3  <- (vac3+vac3_2+vj3)*((S3+R3+RA3)/(n3-(V3+Ih3+D3))) - V3*(1/vd)
    dD3  <-                      I3*gamma*dnh3             +(Ih3/hlos3)*dh3
    
    dS4  <-    - (I1+I2+I3+I4)*(beta*lambda*S4*(1-ef1))/N - (beta*S4*(A1+A2+A3+A4)*(1-ef1))/N + R4*(1/dimmuneI) + RA4*(1/dimmuneA) - (vac4+vac4_2+vj4)*(S4/(n4-(V4+Ih4+D4))) + V4*(1/vd)
    dE4  <-    - E4/alpha   + (I1+I2+I3+I4)*(beta*lambda*S4*(1-ef1))/N + (beta*S4*(A1+A2+A3+A4)*(1-ef1))/N 
    dI4  <- (E4*pS4)/alpha     - I4*gamma
    dIh4 <-                      I4*gamma*hosp4            - Ih4/hlos4
    dA4  <- (E4*(1-pS4))/alpha - A4*gamma
    dR4  <-                      I4*(gamma*(1-hosp4-dnh4)) +(Ih4/hlos4)*(1-dh4) - R4*(1/dimmuneI)  - (vac4+vac4_2+vj4)*(R4/(n4-(V4+Ih4+D4)))
    dRA4 <-                      A4*gamma                                       - RA4*(1/dimmuneA) - (vac4+vac4_2+vj4)*(RA4/(n4-(V4+Ih4+D4)))
    dV4  <- (vac4+vac4_2+vj4)*((S4+R4+RA4)/(n4-(V4+Ih4+D4))) - V4*(1/vd)
    dD4  <-                      I4*gamma*dnh4             +(Ih4/hlos4)*dh4

    
    
    
    return(list(c(dS1, dE1, dI1, dIh1, dA1, dR1, dRA1, dV1, dD1,
                  dS2, dE2, dI2, dIh2, dA2, dR2, dRA2, dV2, dD2,
                  dS3, dE3, dI3, dIh3, dA3, dR3, dRA3, dV3, dD3,
                  dS4, dE4, dI4, dIh4, dA4, dR4, dRA4, dV4, dD4), 
                
                incI = (I1 + I2 + I3 + I4)/9,
                incA = (A1 + A2 + A3 + A4)/9,
                Iht = Ih1+Ih2+Ih4+Ih3,
                Dt = D1 + D2 + D3 + D4,
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
             ef1_3 = scen[i,c('ef1_3')],
             ef1_4 = scen[i,c('ef1_4')],
             ef1_5 = scen[i,c('ef1_5')],
             ef1 = 0,
             ef2 = 0,
             ef3 = 0,
             ef4 = 0,
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
             vac1_ = scen[i,c('vac1')], 
             vac2_ = scen[i,c('vac2')],
             vac3_ = scen[i,c('vac3')],
             vac4_ = scen[i,c('vac4')],
             vac1a = scen[i,c('vac1a')], 
             vac2a = scen[i,c('vac2a')],
             vac3a = scen[i,c('vac3a')],
             vac4a = scen[i,c('vac4a')],
             vac1b = scen[i,c('vac1b')], 
             vac2b = scen[i,c('vac2b')],
             vac3b = scen[i,c('vac3b')],
             vac4b = scen[i,c('vac4b')],
             vac1c = scen[i,c('vac1c')], 
             vac2c = scen[i,c('vac2c')],
             vac3c = scen[i,c('vac3c')],
             vac4c = scen[i,c('vac4c')],
             vac1d = scen[i,c('vac1d')], 
             vac2d = scen[i,c('vac2d')],
             vac3d = scen[i,c('vac3d')],
             vac1e = scen[i,c('vac1e')], 
             vac2e = scen[i,c('vac2e')],
             vac3e = scen[i,c('vac3e')],##Vaccinated number by age group 
             vj1= scen[i,c('vj1')],
             vj2= scen[i,c('vj2')],
             vj3= scen[i,c('vj3')],
             vj4= scen[i,c('vj4')],##J&J vaccine by age group
             vj1a= scen[i,c('vj1a')],
             vj2a= scen[i,c('vj2a')],
             vj3a= scen[i,c('vj3a')],
             vj1b= scen[i,c('vj1b')],
             vj2b= scen[i,c('vj2b')],
             vj3b= scen[i,c('vj3b')],
             vd = scen[i,c('vd')], #Duration of immunity from vaccination
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
             ttraj = scen[i,c('ttraj')],
             tvacend = scen[i,c('tvacend')],
             tv1 = scen[i,c('tv1')],
             tv2 = scen[i,c('tv2')],
             tv3 = scen[i,c('tv3')],
             tv4 = scen[i,c('tv4')],
             tv5 = scen[i,c('tv5')],
             tv6 = scen[i,c('tv6')],
             tvja = scen[i,c('tvja')],
             tvjb = scen[i,c('tvjb')],
             tproject = scen[i,c('tproject')],
             tpa = scen[i,c('tpa')],
             tpb = scen[i,c('tpb')],
             tpc = scen[i,c('tpc')],
  
             ####Add in variant variables
             hosp3v1 = scen[i,c('hosp3v1')],hosp3v2 = scen[i,c('hosp3v2')],
             hosp3v3 = scen[i,c('hosp3v3')],hosp3v4 = scen[i,c('hosp3v4')],
             hosp4v1 = scen[i,c('hosp4v1')],hosp4v2 = scen[i,c('hosp4v2')],
             hosp4v3 = scen[i,c('hosp4v3')],hosp4v4 = scen[i,c('hosp4v4')],
             dh3v1 = scen[i,c('dh3v1')],dh3v2 = scen[i,c('dh3v2')],dh3v3 = scen[i,c('dh3v3')],dh3v4 = scen[i,c('dh3v4')],
             dh4v1 = scen[i,c('dh4v1')],dh4v2 = scen[i,c('dh4v2')],dh4v3 = scen[i,c('dh4v3')],dh4v4 = scen[i,c('dh4v4')],
             dnh3v1 = scen[i,c('dnh3v1')],dnh3v2 = scen[i,c('dnh3v2')],dnh3v3 = scen[i,c('dnh3v3')],dnh3v4 = scen[i,c('dnh3v4')],
             dnh4v1 = scen[i,c('dnh4v1')],dnh4v2 = scen[i,c('dnh4v2')],dnh4v3 = scen[i,c('dnh4v3')],dnh4v4 = scen[i,c('dnh4v4')]
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

# write to csv (Andrea)
write.csv(all.scen, './allscenarios.csv', row.names = F)

#write.csv(all.scen, './scen20210125.csv', row.names = F)

write.csv(all.scen, "/Users/emilywu883/Documents/CU Anschutz/COVID-19/Modeling Team/Data Sets/Output Data/allscenarios.csv", row.names = FALSE)

##### TABLE 2 CODE #####

##DATE of ICU Capacity reached 

#ICUthresh <- all.scen %>% filter(Ict > 1325) %>% select(scenario, scenalpha, date) %>% rename(capdate = date)  
#ICUthresh2 <- ICUthresh %>% group_by(scenario, scenalpha) %>% slice(1)

#library(plyr) ###Annoying thing, for some reason if this is loaded before it prevents rename from working 
##Date of ICU peak and need at peak
#maxc <- ddply(all.scen, .(scenario, scenalpha), transform, maxc = max(Ict))
#detach(package:plyr)
#cmax <- maxc %>% filter(Ict == maxc) %>%
  #select(scenario, scenalpha, date, Ict) %>%
  #rename(maxdate = date) %>%
  #rename(maxICU = Ict)

#table <- filter(merge(ICUthresh2, cmax, by = c('scenario', 'scenalpha'), all = TRUE), !grepl("week", scenalpha))
#table$capdate <- as.Date(table$capdate, "%m/%d/%Y")
#table$maxdate <-as.Date(table$maxdate, "%m/%d/%Y")
#table$maxICU <- label_comma()(ceiling(table$maxICU))

#Cumulative infections through June 1st
###THIS WILL NOT WORK ANYMORE - delete - remove from table or use the sum of Einc
#all.scen$cuminf <- all.scen$Rt + all.scen$Itotal+all.scen$Etotal
#t1 <- filter(subset(all.scen, time == 495,
#             select = c(scenario, scenalpha, cuminf)), !grepl("week", scenalpha))
#t1$cuminf <- label_comma()(signif(t1$cuminf, 3))

#Cumulative infections through June 1st
infectJun <- mutate(all.scen, cuminf = cumsum(Einc))
infectJun <- filter(subset(infectJun, time==494,
                select = c(scenario, scenalpha, cuminf)), !grepl("week", scenalpha))
infectJun$IJun <- infectJun$cuminf

#Cumulative infections through present (used to calculate excess infections from now until June 1st)
infectPresent <- mutate(all.scen, cuminf = cumsum(Einc))
infectPresent <- filter(subset(infectPresent, time==as.numeric(floor_date(today(), "week", 1)-as.Date("2020-01-24")),
                               select = c(scenario, scenalpha, cuminf)), !grepl("week", scenalpha))
infectPresent$IPresent <- infectPresent$cuminf

#Cumulative deaths through June 1st
deathsJun <- filter(subset(all.scen, time == 494,
                           select = c(scenario, scenalpha, Dt)), !grepl("week", scenalpha))
deathsJun$DJun <- deathsJun$Dt

# deaths through present (used to calculate excess deaths from now until June 1st)
deathsPresent <- filter(subset(all.scen, time == as.numeric(floor_date(today(), "week", 1)-as.Date("2020-01-24")),
                               select = c(scenario, scenalpha, Dt)), !grepl("week", scenalpha))
deathsPresent$DPresent <- deathsPresent$Dt

table2 <- cbind(infectJun, infectPresent, deathsJun, deathsPresent)
table2$IExcess <- label_comma()(signif(table2$IJun - table2$IPresent, 3))
table2$DExcess <- label_comma(accuracy = 1)(ceiling(table2$DJun - table2$DPresent))
#table2$capdate <- ifelse(as.character(table2$capdate)=="NA", "N/A", format(table2$capdate, "%m/%d/%y"))
#table2$maxdate <- ifelse(table2$maxdate < floor_date(today(), "week", 1), "past", format(table2$maxdate, "%m/%d/%y"))
table2 <- table2[c("scenalpha", "IExcess", "DExcess")]
table2

write.csv(table2, "/Users/emilywu883/Documents/CU Anschutz/COVID-19/Modeling Team/Data Sets/Output Data/table2.csv", row.names = FALSE)

