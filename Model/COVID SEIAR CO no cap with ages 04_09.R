## Load deSolve package
library(deSolve)

#start: january 21, 2020


seir1 <- function(t, x, parms) {
  
  with(as.list(c(parms, x)), {
    
    
    ef <- ifelse(t < t2, 0, ifelse( t< t3, 0.45, ifelse(t<t4, mag, mag1))) ## Change social distancing efficacy
    siI <- ifelse (t < t1, 0, siI) ##Change proportion of symptomatics that self-isolate
    I <- ifelse(t<t0, 0, 1) ##Timing of infection introduction
    
    dS  <-    - (I+I2+I3)*(beta*lambda*S*(1-siI)*(1-ef))/N - (beta*S*(A+A2+A3)*(1-ef))/N 
    dE  <-    - E/alpha   + (I+I2+I3)*(beta*lambda*S*(1-siI)*(1-ef))/N + (beta*S*(A+A2+A3)*(1-ef))/N 
    dI  <- (E*pS)/alpha - I*(gam) 
    dIh <- I*hosp*gam - Ih*1/8
    dIc <- I*cc*gam - dc*Ic*(1/10) 
    dA  <- (E*(1-pS))/alpha - A*gam
    dR  <- I*(gam*(1-hosp+cc)) + A*gam 
    dRh <- Ih*1/8
    dRc <- Ic*1/10
    dD  <- dc*Ic*(1/10) 
    
    
    dS2  <-    - (I+I2+I3)*(beta*lambda*S2*(1-siI)*(1-ef))/N - (beta*S2*(A+A2+A3)*(1-ef))/N 
    dE2  <-    - E2/alpha   + (I+I2+I3)*(beta*lambda*S2*(1-siI)*(1-ef))/N + (beta*S2*(A+A2+A3)*(1-ef))/N 
    dI2  <- (E2*pS2)/alpha - I2*(gam) 
    dIh2 <- I2*hosp2*gam - Ih2*1/8
    dIc2 <- I2*cc2*gam - dc*Ic2*(1/10) 
    dA2  <- (E2*(1-pS2))/alpha - A2*gam
    dR2  <- I2*(gam*(1-hosp2+cc2)) + A2*gam 
    dRh2 <- Ih2*1/8
    dRc2 <- Ic2*1/10
    dD2  <- dc*Ic2*(1/10) 
    
    dS3  <-    - (I+I2+I3)*(beta*lambda*S3*(1-siI)*(1-ef))/N - (beta*S3*(A+A2+A3)*(1-ef))/N 
    dE3  <-    - E3/alpha   + (I+I2+I3)*(beta*lambda*S3*(1-siI)*(1-ef))/N + (beta*S3*(A+A2+A3)*(1-ef))/N 
    dI3  <- (E3*pS3)/alpha - I3*(gam) 
    dIh3 <- I3*hosp3*gam - Ih3*1/8
    dIc3 <- I3*cc3*gam - dc*Ic3*(1/10) 
    dA3  <- (E3*(1-pS3))/alpha - A3*gam
    dR3  <- I3*(gam*(1-hosp3+cc3)) + A3*gam 
    dRh3 <- Ih3*1/8
    dRc3 <- Ic3*1/10
    dD3  <- dc*Ic3*(1/10) 
    
    der <- c(dS,  dE,  dI,  dIh,  dIc,  dA,  dR,  dRh,  dRc,  dD,
             dS2, dE2, dI2, dIh2, dIc2, dA2, dR2, dRh2, dRc2, dD2,
             dS3, dE3, dI3, dIh3, dIc3, dA3, dR3, dRh3, dRc3, dD3)
    
    list(der,Id = (I+Ih+Ic+I2+Ih2+Ic2+I3+Ih3+Ic3)*pID)
  })
}

scen <- read.csv('data/model 3ages.csv')

n <- as.numeric(nrow(scen)) 
covid_ts <- list() # empty data frame to hold the time series data

gam <- 1/8
alpha <- 5.1 ## incubation period
Cp <- 5840795
n2 <- 2332422
n3 <- 1227124



for(i in 1:n){
## Define parameters that will change

parms <- c(beta = scen[i, c('beta')], ## Transmission rate
           Cp = Cp,
           n2 = n2,
           n3 = n3,
           gam = gam, ## recovery rate (1/average length of infection)
           alpha = alpha, ##duration of latency period
           dc = 0.5,
           ef = 0, ## effectiveness of SD (vary from 0.5 - 1)
           pS = scen[i,c('pS')], ## proportion of infectious individuals symptomatic
           pS2 = scen[i,c('pS2')], ## proportion of infectious individuals symptomatic
           pS3 = scen[i,c('pS3')], ## proportion of infectious individuals symptomatic
           pID = scen[i,c('pID')], ## proportion of symptomatic individuals Identified
           siI = scen[i,c('siI')],## Proportion of symptomatic individuals self isolate
           lambda = scen[i,c('lambda')], ##difference in infectiousness symptomatic/asymptomatic
           hosp = scen[i,c('hosp')], 
           cc = scen[i,c('cc')],
           hosp2 = scen[i,c('hosp2')], 
           cc2 = scen[i,c('cc2')],
           hosp3 = scen[i,c('hosp3')], 
           cc3 = scen[i,c('cc3')],
           mag = scen[i, c('mag')],
           mag1 = scen[i, c('mag1')],
           t0 = scen[i,c('t0')],
           t1 = scen[i,c('t1')],
           t2 = scen[i,c('t2')],
           t3 = scen[i,c('t3')],
           t4 = scen[i,c('t4')])

dt      <- seq(0, 500, 1)

inits      <- c(S = Cp - n2 - n3 -1, E = 0, I = 1, Ih = 0, Ic = 0, A = 0, R = 0, Rh = 0, Rc = 0, D = 0,
                S2 = n2, E2 = 0, I2 = 0, Ih2 = 0, Ic2 = 0, A2 = 0, R2 = 0, Rh2 = 0, Rc2 = 0, D2 = 0,
                S3 = n3, E3 = 0, I3 = 0, Ih3 = 0, Ic3 = 0, A3 = 0, R3 = 0, Rh3 = 0, Rc3 = 0, D3 = 0)

N  <- Cp


out <- lsoda(inits, dt, seir1, parms = parms)
covid_ts[[i]] <- as.matrix(out)
}

#library(dplyr)
all <-  as.data.frame(cbind(rep(1:36, each=501), do.call("rbind", covid_ts)))
all$scenario <- all$V1
all$V1 <- NULL

all.scen <- merge(scen, all, by = "scenario")

write.csv(all.scen, 'allscenarios5.csv', row.names = F)

