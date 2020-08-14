## Load deSolve package
library(deSolve)
#R0: 2.5, 3, 4
#start: january 21, 2020
# 25% or 37.5 mag

seir1 <- function(t, x, parms) {
  
  with(as.list(c(parms, x)), {
    
    ef <- ifelse(t <= 56, 0, mag)##ifelse(t < dStart & t >56, phaseI, mag))
    siI <- ifelse (t < 44, 0, siI)
    cap <- ifelse (t < 91, 709*0.8, ifelse( t < 152, 1709*0.8, 4291*0.8) )
    ##cc <- ifelse (Ic < 708, cc, 0)
    ##dcapac <- ifelse (Ic < 708, 0, 0.009166)
    
    
    dS  <-    - (I)*(beta*lambda*S*(1-siI)*(1-ef))/N - (beta*S*A*(1-ef))/N 
    dE  <-    - E/alpha   + (I)*(beta*lambda*S*(1-siI)*(1-ef))/N + (beta*S*A*(1-ef))/N 
    dI  <- (E*pS)/alpha - I*(gam) 
    dIh <- I*hosp*gam - Ih*1/8
    dIc <- I*cc*gam - dc*min(Ic,cap)*(1/10) - max(((Ic+I*cc*gam)-cap),0)
    dA  <- (E*(1-pS))/alpha - A*gam
    dR  <- I*(gam*(1-hosp+cc)) + A*gam 
    dRh <- Ih*1/8
    dRc <- Ic*1/10
    dD  <- dc*min(Ic,cap)*(1/10) + max(((Ic+(I*cc*gam))-cap),0) 

    der <- c(dS, dE, dI, dIh, dIc, dA, dR, dRh, dRc, dD)
    
    list(der)
  })
}
setwd('/Users/Andrea/Documents/COVID19')
scen <- read.csv('./parms3.csv')

n <- as.numeric(nrow(scen)) 
covid_ts <- list() # empty data frame to hold the time series data

Cp <- 4705120 ## state population
##cases <- 29 ## county cases

cap <- 709
pS <- 0.73426 ## probability of symptoms from (0.179 asymp Mizumoto, 30.8% Nishiura)
pID <- 0.67 ##probability of identifying if symptomatic
gam <- 1/8
alpha <- 5.1 ## incubation period
#dc <- 0.5 ## probability of those in ICU bed that die
siI <- 0.26322 ## Proportion of symptomatic cases that self isolate (ferguson est 2/3, assuming 50% efficacy here)
lambda <- 1.04 ## Difference in infectiousness symptomatic to asymptomatic
sigma <- 0
hosp <- 0.069346 ## percent of symptomatic cases hospitalized
cc <- 0.009166 ## percent of cases requiring CC


for(i in 1:n){
## Define parameters that will change

parms <- c(R0 = scen[i, c('R0')],
           beta = scen[i, c('R0')]*gam, ## Transmission rate
           gam = gam, ## recovery rate (1/average length of infection)
           alpha = alpha, ##duration of latency period
           dc = 0.5,
           ef = 0, ## effectiveness of SD (vary from 0.5 - 1)
           pS = pS, ## proportion of infectious individuals symptomatic
           pID = pID, ## proportion of symptomatic individuals Identified
           ##cases = cases, 
           siI = siI,## Proportion of symptomatic individuals self isolate
           lambda = lambda, ##difference in infectiousness symptomatic/asymptomatic
           hosp = hosp, 
           cc = cc, dcapac = 0,
           mag = scen[i, c('mag')]
           #dStart = scen[i, c('dStart')],
           #phaseI = scen[i, c("phaseI")]
           )

dt      <- seq(0, 500, .1)

inits      <- c(S = Cp-(2), E = 0, I = 1, Ih = 0, Ic = 0, A = 1, R = 0, Rh = 0, Rc = 0, D = 0)

N <- Cp

out <- lsoda(inits, dt, seir1, parms = parms)
covid_ts[[i]] <- as.matrix(out)
}

#library(dplyr)
all <-  as.data.frame(cbind(rep(1:6, each=5001), do.call("rbind", covid_ts)))
all$scenario <- all$V1
all$V1 <- NULL

all.scen <- merge(scen, all, by = "scenario")

write.csv(all.scen, './allscenarios3.csv', row.names = F)

