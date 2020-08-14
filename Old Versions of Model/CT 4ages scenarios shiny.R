## load deSolve package
library(deSolve)

# change the populations as needed according to grouping scheme

Cp <- 5840795 #5359295 # Colorado population = 5840795, 9% higher for previous runs, multiply all by 1.08984
n1 <- 1513005 #1388277 # population age group 1
n2 <- 1685869 #1546890 # population age group 2
n3 <- 1902963 #1746088 # population age group 3
n4 <- 738958 #678040 # population age group 4

# build SEIR model as an R function

seir1 <- function(t, x, parms) {
  
  with(as.list(c(parms, x)), {
    
    # change over time in efficacy of % mag SD among specific age groups
    ef1 <- ifelse(t<t2, 0, ifelse(t<t3, mag1, ifelse(t<t4, mag2, ifelse(t<t5, mag3, ifelse(t<t6,  mag4, ifelse(t< t7, ef1_1,ifelse(t<t8,  ef1_2, ef1_3)))))))
    #ef2 <- ifelse(t<t2, 0, ifelse(t<t3, mag1, ifelse(t<t4, mag2, ifelse(t<t5, mag3, ifelse(t<t6,  mag4, ifelse(t< t7, ef2_1,ifelse(t<t8,  ef2_2, ef2_3)))))))
    #ef3 <- ifelse(t<t2, 0, ifelse(t<t3, mag1, ifelse(t<t4, mag2, ifelse(t<t5, mag3, ifelse(t<t6,  mag4, ifelse(t< t7, ef3_1,ifelse(t<t8,  ef3_2, ef3_3)))))))
    ef4 <- ifelse(t<t2, 0, ifelse(t<t3, mag1, ifelse(t<t4, mag2, ifelse(t<t5, mag3, ifelse(t<t6,  mag4, ifelse(t< t7, ef4_1,ifelse(t<t8,  ef4_2, ef4_3)))))))
    
    siI <- ifelse (t < t1, 0, siI) ##Change proportion of symptomatics that self-isolate after 4/27
    ramp <-ifelse(t < 129, 0, ifelse(t<168,(t-129)*ramp, 39*ramp)) #For ramp up in case isolation : increases proportion of symptomatic case isoaltion over time
    maska <- ifelse(t< 73, 0, ifelse(t< t4,maska,maskb))
    CT  <- ifelse(t < t6, 0, pCT)
    
    dS1  <-    - (I1+I2+I3+I4)*(beta*(1-(maska*0.03))*lambda*S1*(1-(siI+ramp))*(1-ef1))/N - (beta*S1*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
    dE1  <-    - E1/alpha   + (I1+I2+I3+I4)*(beta*(1-(maska*0.03))*lambda*S1*(1-(siI+ramp))*(1-ef1))/N + (beta*S1*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
    dI1  <- (E1*pS1)/alpha - I1*(gamma) -  I1*pID*CT*kap*pi*om
    dII1 <-                         (I1+A1)*pID*CT*kap*pi*om - II1*gamma
    dIh1 <- I1*hosp1*gamma + II1*pS1*hosp1*gamma - Ih1*1/8
    dIc1 <- I1*cc1*gamma   + II1*pS1*cc1*gamma- Ic1*(1/10) 
    dA1  <- (E1*(1-pS1))/alpha - A1*gamma - A1*pID*CT*kap*pi*om
    dR1  <- (I1+II1*pS1)*(gamma*(1-hosp1-cc1)) + A1*gamma 
    dRh1 <- Ih1*1/8
    dRc1 <- (1-dc)*Ic1*1/10
    dD1  <-     dc*Ic1*(1/10) 
    
    dS2  <-    - (I1+I2+I3+I4)*(beta*(1-(maska*0.03))*lambda*S2*(1-(siI+ramp))*(1-ef1))/N - (beta*S2*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
    dE2  <-    - E2/alpha   + (I1+I2+I3+I4)*(beta*(1-(maska*0.03))*lambda*S2*(1-(siI+ramp))*(1-ef1))/N + (beta*S2*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
    dI2  <- (E2*pS2)/alpha - I2*(gamma) -  I2*pID*CT*kap*pi*om
    dII2 <-                         (I2+A2)*pID*CT*kap*pi*om - II2*gamma
    dIh2 <- I2*hosp2*gamma + II2*pS2*hosp2*gamma - Ih2*1/8
    dIc2 <- I2*cc2*gamma   + II2*pS2*cc2*gamma- Ic2*(1/10) 
    dA2  <- (E2*(1-pS2))/alpha - A2*gamma - A2*pID*CT*kap*pi*om
    dR2  <- (I2+II2*pS2)*(gamma*(1-hosp2-cc2)) + A2*gamma 
    dRh2 <- Ih2*1/8
    dRc2 <- (1-dc)*Ic2*1/10
    dD2  <-     dc*Ic2*(1/10) 
    
    dS3  <-    - (I1+I2+I3+I4)*(beta*(1-(maska*0.03))*lambda*S3*(1-(siI+ramp))*(1-ef1))/N - (beta*S3*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
    dE3  <-    - E3/alpha   + (I1+I2+I3+I4)*(beta*(1-(maska*0.03))*lambda*S3*(1-(siI+ramp))*(1-ef1))/N + (beta*S3*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef1))/N 
    dI3  <- (E3*pS3)/alpha - I3*(gamma)  - I3*pID*CT*kap*pi*om
    dII3 <-                         (I3+A3)*pID*CT*kap*pi*om - II3*gamma
    dIh3 <- I3*hosp3*gamma + II3*pS3*hosp3*gamma - Ih3*1/8
    dIc3 <- I3*cc3*gamma   + II3*pS3*cc3*gamma- Ic3*(1/10) 
    dA3  <- (E3*(1-pS3))/alpha - A3*gamma - A3*pID*CT*kap*pi*om
    dR3  <- (I3+II3*pS3)*(gamma*(1-hosp3-cc3)) + A3*gamma 
    dRh3 <- Ih3*1/8
    dRc3 <- (1-dc)*Ic3*1/10
    dD3  <-    dc*Ic3*(1/10) 
    
    dS4  <-    - (I1+I2+I3+I4)*(beta*(1-(maska*0.03))*lambda*S4*(1-(siI+ramp))*(1-ef4))/N - (beta*S4*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef4))/N 
    dE4  <-    - E4/alpha   + (I1+I2+I3+I4)*(beta*(1-(maska*0.03))*lambda*S4*(1-(siI+ramp))*(1-ef4))/N + (beta*S4*(1-(maska*0.2667))*(A1+A2+A3+A4)*(1-ef4))/N 
    dI4  <- (E4*pS4)/alpha - I4*(gamma)  - I4*pID*CT*kap*pi*om
    dII4 <-                         (I4+A4)*pID*CT*kap*pi*om - II4*gamma
    dIh4 <- I4*hosp4*gamma + II4*pS4*hosp4*gamma - Ih4*1/8
    dIc4 <- I4*cc4*gamma   + II4*pS4*cc4*gamma- Ic4*(1/10) 
    dA4  <- (E4*(1-pS4))/alpha - A4*gamma - A4*pID*CT*kap*pi*om
    dR4  <- (I4+II4*pS4)*(gamma*(1-hosp4-cc4)) + A4*gamma 
    dRh4 <- Ih4*1/8
    dRc4 <- (1-dc)*Ic4*1/10
    dD4  <-    dc*Ic4*(1/10) 
    
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

scen <- read.csv('CT 4ages params 0609.csv')

# rows (n) to represent scenario numbers
n <- as.numeric(nrow(scen)) 
covid_ts <- list() # empty data frame to hold the time series data

# run simulations from time 1 to 500, one simulation per scenario row for as many rows as we have

for(i in 1:n){
  # define parameters that will change
  
  parms <- c(beta = scen[i, c('beta')], # transmission rate
             gamma = 1/8,
             alpha = 5.1,
             Cp = Cp,
             n1 = n1,
             n2 = n2,
             n3 = n3,
             n4 = n4,
             dc = 0.5, # death rate from ICU
             ef1_1 = scen[i,c('ef1_1')],
             ef1_2 = scen[i,c('ef1_2')],
             ef1_3 = scen[i,c('ef1_3')],
             ef2_1 = scen[i,c('ef2_1')],
             ef2_2 = scen[i,c('ef2_2')],
             ef2_3 = scen[i,c('ef2_3')],
             ef3_1 = scen[i,c('ef3_1')],
             ef3_2 = scen[i,c('ef3_2')],
             ef3_3 = scen[i,c('ef3_3')],
             ef4_1 = scen[i,c('ef4_1')],
             ef4_2 = scen[i,c('ef4_2')],
             ef4_3 = scen[i,c('ef4_3')],
             ef1 = 0,
             ef2 = 0,
             ef3 = 0,
             ef4 = 0,
             pS1 = scen[i,c('pS1')], ## proportion of infectious individuals symptomatic (0-19)
             pS2 = scen[i,c('pS2')], ## proportion of infectious individuals symptomatic (20-39)
             pS3 = scen[i,c('pS3')], ## proportion of infectious individuals symptomatic (40-59)
             pS4 = scen[i,c('pS4')], ## proportion of infectious individuals symptomatic (60-79)
             pID = scen[i,c('pID')], ## proportion of symptomatic individuals Identified
             siI = scen[i,c('siI')],## Proportion of symptomatic individuals self isolate
             lambda = scen[i,c('lambda')], ##difference in infectiousness symptomatic/asymptomatic
             hosp1 = scen[i,c('hosp1')], 
             cc1 = scen[i,c('cc1')],
             hosp2 = scen[i,c('hosp2')], 
             cc2 = scen[i,c('cc2')],
             hosp3 = scen[i,c('hosp3')], 
             cc3 = scen[i,c('cc3')],
             hosp4 = scen[i,c('hosp4')], 
             cc4 = scen[i,c('cc4')],
             mag1 = scen[i, c('mag1')],
             mag2 = scen[i, c('mag2')],
             mag3 = scen[i, c('mag3')],
             mag4 = scen[i, c('mag4')],
             t1 = scen[i,c('t1')],
             t2 = scen[i,c('t2')],
             t3 = scen[i,c('t3')],
             t4 = scen[i,c('t4')],
             t5 = scen[i,c('t5')],
             t6 = scen[i,c('t6')],
             t7 = scen[i,c('t7')],
             t8 = scen[i,c('t8')],
             ramp = scen[i,c('ramp')],
             maska = scen[i,c('maska')],
             maskb = scen[i,c('maskb')],
             kap = scen[i,c("kap")], #average number of contacts traced per detected case
             pCT = scen[i,c("pCT")], #proportion of identified cases with contacts traced
             pi = scen[i,c("pi")], #probability a contact traced infected individual is isolated before infecting other susceptibles 
             om = scen[i,c("om")] #probability a contact traced individual is infected
  )
  
  dt      <- seq(0, 500, 1)
  
  inits      <- c(S1 = n1 - 1, E1 = 0, I1 = 1, II1 = 0, Ih1 = 0, Ic1 = 0, A1 = 0, R1 = 0, Rh1 = 0, Rc1 = 0, D1 = 0,
                  S2 = n2,     E2 = 0, I2 = 0, II2 = 0, Ih2 = 0, Ic2 = 0, A2 = 0, R2 = 0, Rh2 = 0, Rc2 = 0, D2 = 0,
                  S3 = n3,     E3 = 0, I3 = 0, II3 = 0, Ih3 = 0, Ic3 = 0, A3 = 0, R3 = 0, Rh3 = 0, Rc3 = 0, D3 = 0,
                  S4 = n4,     E4 = 0, I4 = 0, II4 = 0, Ih4 = 0, Ic4 = 0, A4 = 0, R4 = 0, Rh4 = 0, Rc4 = 0, D4 = 0)
  N  <- Cp
  
  
  out <- lsoda(inits, dt, seir1, parms = parms)
  covid_ts[[i]] <- as.matrix(out)
}

#library(dplyr)
all <-  as.data.frame(cbind(rep(1:n, each=501), do.call("rbind", covid_ts)))
all$scenario <- all$V1
all$V1 <- NULL

all.scen <- merge(scen, all, by = "scenario")

write.csv(all.scen, './allscenarios_4ageCT.csv', row.names = F)

# plot total daily hospitalizations for all the proposed scenarios for each of the age groups

library(ggplot2)

# Multiple plot function

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# subset into Scenario 1
# add an incrementing date vector with January 24th, 2020 as the starting value

scen1 <- subset(all.scen, all.scen$scenalpha == "1A" | all.scen$scenalpha == "1B" |
                  all.scen$scenalpha == "1C" | all.scen$scenalpha == "1D")

scen1$date <- seq(from = as.Date("2020/1/24"), to = as.Date("2020/1/24") + 500, "days")

# generate the total hospitalization and ICU plots for Scenario 1

p1hosp <- ggplot(data = scen1, aes(x = date, y = Iht, group = scenalpha, color = scenalpha)) +
  geom_line() + ylim(0, 17500) + labs(x = "Date", y = "Count",
                                     title = "Total Hospitalizations, Scenario 1(A-D)",
                                     caption = "*May 27th   **May 27th -> June 27th",
                                     col = "Scenario") + scale_x_date(date_labels = "%d%b%y", date_breaks = "3 months") +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(name = "Scenario",
                     labels = c("65%*", "55%*", "45%*", "55% -> 45%**"),
                     values = c("green", "turquoise", "gold", "deeppink")) +
  theme(plot.background = element_rect(fill = "darkblue"),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.background = element_rect(fill = "darkblue"),
        legend.background = element_rect(fill = "darkblue"),
        legend.key = element_rect(color = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 24, color = "white", family = "serif", hjust = 0.5),
        plot.caption = element_text(size = 12, color = "white", family = "serif", hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "white", family = "serif"),
        axis.line = element_line(color = "white"),
        axis.title.x = element_text(size = 16, color = "white", family = "serif", angle = 30),
        axis.title.y = element_text(size = 16, color = "white", family = "serif", angle = 30),
        axis.text.x = element_text(size = 14, color = "white", family = "serif", angle = 30, margin = margin(15, 5, 5, 5, "pt")),
        axis.text.y = element_text(size = 14, color = "white", family = "serif", angle = 30, margin = margin(0, 5, 0, 10, "pt")))

p1hosp

p1icu <- ggplot(data = scen1, aes(x = date, y = Ict, group = scenalpha, color = scenalpha)) +
  geom_line() + labs(x = "Date", y = "Count",
                     title = "Total ICU Need, Scenario 1(A-D)",
                     caption = "*May 27th   **May 27th -> June 27th",
                     col = "Scenario") + scale_x_date(date_labels = "%d%b%y", date_breaks = "3 months") +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(name = "Scenario",
                     labels = c("65%*", "55%*", "45%*", "55% -> 45%**"),
                     values = c("green", "turquoise", "gold", "deeppink")) +
  theme(plot.background = element_rect(fill = "darkblue"),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.background = element_rect(fill = "darkblue"),
        legend.background = element_rect(fill = "darkblue"),
        legend.key = element_rect(color = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 24, color = "white", family = "serif", hjust = 0.5),
        plot.caption = element_text(size = 12, color = "white", family = "serif", hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "white", family = "serif"),
        axis.line = element_line(color = "white"),
        axis.title.x = element_text(size = 16, color = "white", family = "serif", angle = 30),
        axis.title.y = element_text(size = 16, color = "white", family = "serif", angle = 30),
        axis.text.x = element_text(size = 14, color = "white", family = "serif", angle = 30, margin = margin(15, 5, 5, 5, "pt")),
        axis.text.y = element_text(size = 14, color = "white", family = "serif", angle = 30, margin = margin(0, 5, 0, 10, "pt")))

p1icu

# Repeat for Scenarios 2 and 3

# subset into Scenario 2
# add an incrementing date vector with January 24th, 2020 as the starting value

scen2 <- subset(all.scen, all.scen$scenalpha == "2A" | all.scen$scenalpha == "2B" |
                  all.scen$scenalpha == "2C" | all.scen$scenalpha == "2D")

scen2$date <- seq(from = as.Date("2020/1/24"), to = as.Date("2020/1/24") + 500, "days")

# generate the total hospitalization and ICU plots for Scenario 2

p2hosp <- ggplot(data = scen2, aes(x = date, y = Iht, group = scenalpha, color = scenalpha)) +
  geom_line() + ylim(0, 2500) + labs(x = "Date", y = "Count",
                                     title = "Total Hospitalizations, Scenario 2(A-D)",
                                     caption = "*May 27th   **May 27th -> June 27th",
                                     col = "Scenario") + scale_x_date(date_labels = "%d%b%y", date_breaks = "3 months") +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(name = "Scenario",
                     labels = c("65%*", "55%*", "45%*", "55% -> 45%**"),
                     values = c("green", "turquoise", "gold", "deeppink")) +
  theme(plot.background = element_rect(fill = "darkblue"),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.background = element_rect(fill = "darkblue"),
        legend.background = element_rect(fill = "darkblue"),
        legend.key = element_rect(color = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 24, color = "white", family = "serif", hjust = 0.5),
        plot.caption = element_text(size = 12, color = "white", family = "serif", hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "white", family = "serif"),
        axis.line = element_line(color = "white"),
        axis.title.x = element_text(size = 16, color = "white", family = "serif", angle = 30),
        axis.title.y = element_text(size = 16, color = "white", family = "serif", angle = 30),
        axis.text.x = element_text(size = 14, color = "white", family = "serif", angle = 30, margin = margin(15, 5, 5, 5, "pt")),
        axis.text.y = element_text(size = 14, color = "white", family = "serif", angle = 30, margin = margin(0, 5, 0, 10, "pt")))

p2hosp

p2icu <- ggplot(data = scen2, aes(x = date, y = Ict, group = scenalpha, color = scenalpha)) +
  geom_line() + ylim(0, 800) + labs(x = "Date", y = "Count",
                                    title = "Total ICU Need, Scenario 2(A-D)",
                                    caption = "*May 27th   **May 27th -> June 27th",
                                    col = "Scenario") + scale_x_date(date_labels = "%d%b%y", date_breaks = "3 months") +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(name = "Scenario",
                     labels = c("65%*", "55%*", "45%*", "55% -> 45%**"),
                     values = c("green", "turquoise", "gold", "deeppink")) +
  theme(plot.background = element_rect(fill = "darkblue"),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.background = element_rect(fill = "darkblue"),
        legend.background = element_rect(fill = "darkblue"),
        legend.key = element_rect(color = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 24, color = "white", family = "serif", hjust = 0.5),
        plot.caption = element_text(size = 12, color = "white", family = "serif", hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "white", family = "serif"),
        axis.line = element_line(color = "white"),
        axis.title.x = element_text(size = 16, color = "white", family = "serif", angle = 30),
        axis.title.y = element_text(size = 16, color = "white", family = "serif", angle = 30),
        axis.text.x = element_text(size = 14, color = "white", family = "serif", angle = 30, margin = margin(15, 5, 5, 5, "pt")),
        axis.text.y = element_text(size = 14, color = "white", family = "serif", angle = 30, margin = margin(0, 5, 0, 10, "pt")))

p2icu

# subset into Scenario 3
# add an incrementing date vector with January 24th, 2020 as the starting value

scen3 <- subset(all.scen, all.scen$scenalpha == "3A" | all.scen$scenalpha == "3B" |
                  all.scen$scenalpha == "3C" | all.scen$scenalpha == "3D")

scen3$date <- seq(from = as.Date("2020/1/24"), to = as.Date("2020/1/24") + 500, "days")

# generate the total hospitalization and ICU plots for Scenario 3

p3hosp <- ggplot(data = scen3, aes(x = date, y = Iht, group = scenalpha, color = scenalpha)) +
  geom_line() + ylim(0, 4000) + labs(x = "Date", y = "Count",
                                     title = "Total Hospitalizations, Scenario 3(A-D)",
                                     caption = "*May 27th   **May 27th -> June 27th",
                                     col = "Scenario") + scale_x_date(date_labels = "%d%b%y", date_breaks = "3 months") +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(name = "Scenario",
                     labels = c("65%*", "55%*", "45%*", "55% -> 45%**"),
                     values = c("green", "turquoise", "gold", "deeppink")) +
  theme(plot.background = element_rect(fill = "darkblue"),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.background = element_rect(fill = "darkblue"),
        legend.background = element_rect(fill = "darkblue"),
        legend.key = element_rect(color = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 24, color = "white", family = "serif", hjust = 0.5),
        plot.caption = element_text(size = 12, color = "white", family = "serif", hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "white", family = "serif"),
        axis.line = element_line(color = "white"),
        axis.title.x = element_text(size = 16, color = "white", family = "serif", angle = 30),
        axis.title.y = element_text(size = 16, color = "white", family = "serif", angle = 30),
        axis.text.x = element_text(size = 14, color = "white", family = "serif", angle = 30, margin = margin(15, 5, 5, 5, "pt")),
        axis.text.y = element_text(size = 14, color = "white", family = "serif", angle = 30, margin = margin(0, 5, 0, 10, "pt")))

p3hosp

p3icu <- ggplot(data = scen3, aes(x = date, y = Ict, group = scenalpha, color = scenalpha)) +
  geom_line() + ylim(0, 1200) + labs(x = "Date", y = "Count",
                                     title = "Total ICU Need, Scenario 3(A-D)",
                                     caption = "*May 27th   **May 27th -> June 27th",
                                     col = "Scenario") + scale_x_date(date_labels = "%d%b%y", date_breaks = "3 months") +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(name = "Scenario",
                     labels = c("65%*", "55%*", "45%*", "55% -> 45%**"),
                     values = c("green", "turquoise", "gold", "deeppink")) +
  theme(plot.background = element_rect(fill = "darkblue"),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.background = element_rect(fill = "darkblue"),
        legend.background = element_rect(fill = "darkblue"),
        legend.key = element_rect(color = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 24, color = "white", family = "serif", hjust = 0.5),
        plot.caption = element_text(size = 12, color = "white", family = "serif", hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "white", family = "serif"),
        axis.line = element_line(color = "white"),
        axis.title.x = element_text(size = 16, color = "white", family = "serif", angle = 30),
        axis.title.y = element_text(size = 16, color = "white", family = "serif", angle = 30),
        axis.text.x = element_text(size = 14, color = "white", family = "serif", angle = 30, margin = margin(15, 5, 5, 5, "pt")),
        axis.text.y = element_text(size = 14, color = "white", family = "serif", angle = 30, margin = margin(0, 5, 0, 10, "pt")))

p3icu
