## Load deSolve package
library(deSolve)

#start: january 21, 2020

seir1 <- function(t, x, parms) {
    
    with(as.list(c(parms, x)), {
        
        N <- Cp
        
        # change over time in efficacy of % mag SD among specific age groups
        ef1 <- ifelse(t<t2, 0, ifelse(t<t3, mag1, ifelse(t<t4, mag2, ifelse(t<t4a, mag3, ifelse(t<t5, mag3a, ifelse(t<t5a, mag4, ifelse(t<t5b, mag4a, ifelse(t<t6,  mag4b, ifelse( t< t7, ef1_1,ifelse(t<t8,  ef1_2, ef1_3))))))))))
        #ef2 <- ifelse(t<t2, 0, ifelse(t<t3, mag1, ifelse(t<t4, mag2, ifelse(t<t5, mag3, ifelse(t<t6,  mag4, ifelse(t< t7, ef2_1,ifelse(t<t8,  ef2_2, ef2_3)))))))
        #ef3 <- ifelse(t<t2, 0, ifelse(t<t3, mag1, ifelse(t<t4, mag2, ifelse(t<t5, mag3, ifelse(t<t6,  mag4, ifelse(t< t7, ef3_1,ifelse(t<t8,  ef3_2, ef3_3)))))))
        #ef4 <- ifelse(t<t2, 0, ifelse(t<t3, mag1, ifelse(t<t4, mag2, ifelse(t<t5, mag3, ifelse(t<t6,  mag4, ifelse(t< t7, ef4_1,ifelse(t<t8,  ef4_2, ef4_3)))))))
        ef4 <- ifelse(t<t7, ef1, (ef1*(1-ef4p)) + (ef4p*0.8)) #ef4p is the proportion of adults over 65 practicing high (80%) social distancing
        
        siI <- ifelse (t < t1, 0, siI) ##Turn on symptomatics that self-isolate after 03/05
        ramp <-ifelse(t < 129, 0, ifelse(t<168,(t-129)*ramp, 39*ramp)) #For ramp up in case isolation : increases proportion of symptomatic case isoaltion over time
        maska <- ifelse(t< 73, 0, ifelse(t< t4,maska,maskb))
        CT  <- ifelse(t < t6, 0, pCT)
        #temp <- ifelse (t > 1, ifelse(temp_on == 1, temptheory$temp.param[[t]],1), 1)
        temp <-ifelse(temp_on == 1, 0.5*cos((t+45)*0.017)+1.5, 1)
        
        dS1  <-    - (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S1*(1-(siI+ramp))*(1-ef1))/N - (beta*temp*S1*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
        dE1  <-    - E1/alpha   + (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S1*(1-(siI+ramp))*(1-ef1))/N + (beta*temp*S1*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
        dI1  <- (E1*pS1)/alpha - I1*(gamma) -  I1*pID*CT*kap*pi*om
        dII1 <-                         (I1+A1)*pID*CT*kap*pi*om - II1*gamma
        dIh1 <- I1*hosp1*gamma + II1*pS1*hosp1*gamma - Ih1*1/3.6
        dIc1 <- I1*cc1*gamma   + II1*pS1*cc1*gamma- Ic1*(1/4.4) 
        dA1  <- (E1*(1-pS1))/alpha - A1*gamma - A1*pID*CT*kap*pi*om
        dR1  <- (I1+II1*pS1)*(gamma*(1-hosp1-cc1)) + A1*gamma 
        dRh1 <- (1-dh1)*Ih1*1/3.6
        dRc1 <- (1-dc1)*Ic1*1/4.4
        dD1  <-     dc1*Ic1*(1/4.4) + dh1*Ih1*1/3.6
        
        dS2  <-    - (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S2*(1-(siI+ramp))*(1-ef1))/N - (beta*temp*S2*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
        dE2  <-    - E2/alpha   + (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S2*(1-(siI+ramp))*(1-ef1))/N + (beta*temp*S2*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
        dI2  <- (E2*pS2)/alpha - I2*(gamma) -  I2*pID*CT*kap*pi*om
        dII2 <-                         (I2+A2)*pID*CT*kap*pi*om - II2*gamma
        dIh2 <- I2*hosp2*gamma + II2*pS2*hosp2*gamma - Ih2*1/4.6
        dIc2 <- I2*cc2*gamma   + II2*pS2*cc2*gamma- Ic2*(1/7.3) 
        dA2  <- (E2*(1-pS2))/alpha - A2*gamma - A2*pID*CT*kap*pi*om
        dR2  <- (I2+II2*pS2)*(gamma*(1-hosp2-cc2)) + A2*gamma 
        dRh2 <- (1-dh2)*Ih2*1/4.6
        dRc2 <- (1-dc2)*Ic2*1/7.3
        dD2  <-     dc*Ic2*(1/7.3) + dh2*Ih2*(1/4.6)
        
        dS3  <-    - (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S3*(1-(siI+ramp))*(1-ef1))/N - (beta*temp*S3*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
        dE3  <-    - E3/alpha   + (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S3*(1-(siI+ramp))*(1-ef1))/N + (beta*temp*S3*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
        dI3  <- (E3*pS3)/alpha - I3*(gamma)  - I3*pID*CT*kap*pi*om
        dII3 <-                         (I3+A3)*pID*CT*kap*pi*om - II3*gamma
        dIh3 <- I3*hosp3*gamma + II3*pS3*hosp3*gamma - Ih3*1/6.8
        dIc3 <- I3*cc3*gamma   + II3*pS3*cc3*gamma- Ic3*(1/11.4) 
        dA3  <- (E3*(1-pS3))/alpha - A3*gamma - A3*pID*CT*kap*pi*om
        dR3  <- (I3+II3*pS3)*(gamma*(1-hosp3-cc3)) + A3*gamma 
        dRh3 <- (1-dh3)*Ih3*1/6.8
        dRc3 <- (1-dc3)*Ic3*(1/11.4)
        dD3  <-    dc3 *Ic3*(1/11.4) + dh3*Ih3*(1/6.8) 
        
        dS4  <-    - (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S4*(1-(siI+ramp))*(1-ef4))/N - (beta*temp*S4*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef4))/N 
        dE4  <-    - E4/alpha   + (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S4*(1-(siI+ramp))*(1-ef4))/N + (beta*temp*S4*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef4))/N 
        dI4  <- (E4*pS4)/alpha - I4*(gamma)  - I4*pID*CT*kap*pi*om
        dII4 <-                         (I4+A4)*pID*CT*kap*pi*om - II4*gamma
        dIh4 <- I4*hosp4*gamma + II4*pS4*hosp4*gamma - Ih4*1/9.7
        dIc4 <- I4*cc4*gamma   + II4*pS4*cc4*gamma- Ic4*(1/9.9) 
        dA4  <- (E4*(1-pS4))/alpha - A4*gamma - A4*pID*CT*kap*pi*om
        dR4  <- (I4+II4*pS4)*(gamma*(1-hosp4-cc4)) + A4*gamma 
        dRh4 <- (1-dh4)*Ih4*(1/9.7)
        dRc4 <- (1-dc4)*Ic4*(1/9.9)
        dD4  <-    dc4* Ic4*(1/9.9) + dh4*Ih4*(1/9.7) 
        
        
        der <- c(dS1, dE1, dI1, dII1, dIh1, dIc1, dA1, dR1, dRh1, dRc1, dD1,
                 dS2, dE2, dI2, dII2, dIh2, dIc2, dA2, dR2, dRh2, dRc2, dD2,
                 dS3, dE3, dI3, dII3, dIh3, dIc3, dA3, dR3, dRh3, dRc3, dD3,
                 dS4, dE4, dI4, dII4, dIh4, dIc4, dA4, dR4, dRh4, dRc4, dD4)
        
        list(der,
             It = I1 + I2 + I3 + I4,
             IIt = II1 + II2 + II3 + II4,
             Iht =Ih1 + Ih2 + Ih3 + Ih4 + Ic1 + Ic2 + Ic3 + Ic4, 
             Iht1 = Ih1 +Ic1,
             Iht2 = Ih2 +Ic2,
             Iht3 = Ih3 +Ic3,
             Iht4 = Ih4 +Ic4,
             Ict =Ic1 + Ic2 + Ic3 + Ic4)
    })
}

## Model for Deaths
seir1D <- function(t, x, parms) {
    
    with(as.list(c(parms, x)), {
        
        N  <- Cp
        
        # change over time in efficacy of % mag SD among specific age groups
        ef1 <- ifelse(t<t2, 0, ifelse(t<t3, mag1, ifelse(t<t4, mag2, ifelse(t<t5, mag3, ifelse(t<t6,  mag4, ifelse(t< t7, ef1_1,ifelse(t<t8,  ef1_2, ef1_3)))))))
        #ef2 <- ifelse(t<t2, 0, ifelse(t<t3, mag1, ifelse(t<t4, mag2, ifelse(t<t5, mag3, ifelse(t<t6,  mag4, ifelse(t< t7, ef2_1,ifelse(t<t8,  ef2_2, ef2_3)))))))
        #ef3 <- ifelse(t<t2, 0, ifelse(t<t3, mag1, ifelse(t<t4, mag2, ifelse(t<t5, mag3, ifelse(t<t6,  mag4, ifelse(t< t7, ef3_1,ifelse(t<t8,  ef3_2, ef3_3)))))))
        ef4 <- ifelse(t<t7, ef1, (ef1*(1-ef4p)) + (ef4p*0.8)) 
        
        siI <- ifelse (t < t1, 0, siI) ##Change proportion of symptomatics that self-isolate after 4/27
        ramp <-ifelse(t < 129, 0, ifelse(t<168,(t-129)*ramp, 39*ramp)) #For ramp up in case isolation : increases proportion of symptomatic case isoaltion over time
        maska <- ifelse(t< 73, 0, ifelse(t< t4,maska,maskb))
        CT  <- ifelse(t < t6, 0, pCT)
        temp <-ifelse(temp_on == 1, 0.5*cos((t+45)*0.017)+1.5, 1)
        
        
        dS1  <-    - (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S1*(1-(siI+ramp))*(1-ef1))/N - (beta*temp*S1*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
        dE1  <-    - E1/alpha   + (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S1*(1-(siI+ramp))*(1-ef1))/N + (beta*temp*S1*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
        dI1  <- (E1*pS1)/alpha - I1*(gamma) -  I1*pID*CT*kap*pi*om
        dII1 <-                         (I1+A1)*pID*CT*kap*pi*om - II1*gamma
        dIh1 <- I1*hosp1*gamma + II1*pS1*hosp1*gamma - Ih1*1/3.6
        dA1  <- (E1*(1-pS1))/alpha - A1*gamma - A1*pID*CT*kap*pi*om
        dR1  <- (I1+II1*pS1)*(gamma*(1-hosp1-cc1)) + A1*gamma 
        dRh1 <- (1-dh1)*Ih1*1/3.6

        
        dS2  <-    - (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S2*(1-(siI+ramp))*(1-ef1))/N - (beta*temp*S2*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
        dE2  <-    - E2/alpha   + (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S2*(1-(siI+ramp))*(1-ef1))/N + (beta*temp*S2*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
        dI2  <- (E2*pS2)/alpha - I2*(gamma) -  I2*pID*CT*kap*pi*om
        dII2 <-                         (I2+A2)*pID*CT*kap*pi*om - II2*gamma
        dIh2 <- I2*hosp2*gamma + II2*pS2*hosp2*gamma - Ih2*1/4.6
        dA2  <- (E2*(1-pS2))/alpha - A2*gamma - A2*pID*CT*kap*pi*om
        dR2  <- (I2+II2*pS2)*(gamma*(1-hosp2-cc2)) + A2*gamma 
        dRh2 <- (1-dh2)*Ih2*1/4.6
      
        
        dS3  <-    - (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S3*(1-(siI+ramp))*(1-ef1))/N - (beta*temp*S3*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
        dE3  <-    - E3/alpha   + (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S3*(1-(siI+ramp))*(1-ef1))/N + (beta*temp*S3*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
        dI3  <- (E3*pS3)/alpha - I3*(gamma)  - I3*pID*CT*kap*pi*om
        dII3 <-                         (I3+A3)*pID*CT*kap*pi*om - II3*gamma
        dIh3 <- I3*hosp3*gamma + II3*pS3*hosp3*gamma - Ih3*1/6.8
        dA3  <- (E3*(1-pS3))/alpha - A3*gamma - A3*pID*CT*kap*pi*om
        dR3  <- (I3+II3*pS3)*(gamma*(1-hosp3-cc3)) + A3*gamma 
        dRh3 <- (1-dh3)*Ih3*1/6.8
  
        
        dS4  <-    - (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S4*(1-(siI+ramp))*(1-ef4))/N - (beta*temp*S4*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef4))/N 
        dE4  <-    - E4/alpha   + (I1+I2+I3+I4)*(beta*temp*(1-(maska*0.03))*lambda*S4*(1-(siI+ramp))*(1-ef4))/N + (beta*temp*S4*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef4))/N 
        dI4  <- (E4*pS4)/alpha - I4*(gamma)  - I4*pID*CT*kap*pi*om
        dII4 <-                         (I4+A4)*pID*CT*kap*pi*om - II4*gamma
        dIh4 <- I4*hosp4*gamma + II4*pS4*hosp4*gamma - Ih4*1/9.7
        dA4  <- (E4*(1-pS4))/alpha - A4*gamma - A4*pID*CT*kap*pi*om
        dR4  <- (I4+II4*pS4)*(gamma*(1-hosp4-cc4)) + A4*gamma 
        dRh4 <- (1-dh4)*Ih4*(1/9.7)
        
        
        
        
        dIc <- ((I1+II1*pS1)*cc1 + (I2+II2*pS2)*cc2 + (I3+II3*pS3)*cc3 + (I4+II4*pS4)*cc4)*gamma - min(Ic,cap)*(1/10) - max(((Ic + ((I1+II1*pS1)*cc1 + (I2+II2*pS2)*cc2 + (I3+II3*pS3)*cc3 + (I4+II4*pS4)*cc4)*gamma)-cap),0)    
        dRc <- (1 - dc)*min(Ic,cap)*(1/10)
        dD  <-      dc *min(Ic,cap)*(1/10)                       + max(((Ic + I1*cc1*gamma + I2*cc2*gamma + I3*cc3*gamma + I4*cc4*gamma)-cap),0) + 0.16*(Ih1+Ih2+Ih3+Ih4)*(1/8)
        
              
        der <- c(dS1, dE1, dI1, dII1, dIh1, dA1, dR1, dRh1, dIc, dRc, dD,
                 dS2, dE2, dI2, dII2, dIh2, dA2, dR2, dRh2, 
                 dS3, dE3, dI3, dII3, dIh3, dA3, dR3, dRh3, 
                 dS4, dE4, dI4, dII4, dIh4, dA4, dR4, dRh4)
        
        list(der,
             It = I1 + I2 + I3 + I4,
             IIt = II1 + II2 + II3 + II4)
    })
}