library(shiny)
library(shinyWidgets)
library(ggplot2)
library(gridExtra)
library(deSolve)
library(dplyr)
library(tidyr)
library(plotly)
library(hrbrthemes)
require(scales)
require(gridExtra)
source("models.R")
source("figures.R")


# Define UI
ui <- navbarPage("Modeling COVID-19 in Colorado",
                 
                 
                 tabPanel("Model",
                          # Application title
                          fluidRow(column(1), 
                                   column(10, h3("Modeling COVID-19 in Colorado")),
                                   column(1)),
                          
                          # App introductory text
                          fluidRow(
                            column(1),
                            column(10, 
                                   h3("How to use this app"),
                                   p("This app allows the user to generate projections of COVID-19 cases, hospital needs
                                     and deaths in Colorado under scenarios that can be altered interactively. It is 
                                     intended to show the potential impact of the four key control measures: social 
                                     distancing, mask wearing, identification and isolation of cases, and contact 
                                     tracing—the key tools used to limit the spread of the SARS-CoV-2 virus, which 
                                     causes the illness COVID-19."),
                                   p(HTML(paste0("For details regarding the model and underlying assumptions, see the 
                                                  documentation tab. For a brief video tutorial on how to use this app, ", 
                                                 a(href="zoom_3.mp4", "click here.", target = "_blank"))))),
                            column(1)
                          ),
                          
                          br(),
                          
                          # Plot(s) and checkboxes
                          fluidRow(
                            column(1),
                            column(10, 
                                   sidebarLayout(
                                     # Checkboxes
                                     sidebarPanel(
                                       prettyCheckbox(inputId = "Infections", label = "Symptomatic Infections",
                                                      outline = TRUE, status = "success", value = FALSE),
                                       prettyCheckbox(inputId = "Hospitalizations", label = "Non-ICU Hospitalizations",
                                                      outline = TRUE, status = "info", value = TRUE),
                                       prettyCheckbox(inputId = "CriticalCare", label = "Critical Care",
                                                      outline = TRUE, status = "warning", value = TRUE),
                                       prettyCheckbox(inputId = "Deaths", label = "Deaths",
                                                      outline = TRUE, status = "danger", value = TRUE),
                                       prettyCheckbox(inputId = "ICUCapLine", label = "ICU Capacity",
                                                      outline = TRUE, status = "default", value = TRUE),
                                       prettyCheckbox(inputId = "CTCapLine", label = "Contact Tracing Capacity",
                                                      outline = TRUE, value = FALSE),
                                       prettyRadioButtons(inputId = "Daily", label = " ",
                                                          choices = list("Daily" = TRUE, "Cumulative" = FALSE),
                                                          shape = "square", status = "default"),
                                       br(),
                                       p("Contact:"),
                                       p(HTML(paste0(a(href = 'mailto:co.sph.covid@ucdenver.edu',
                                                       'CU COVID-19 Modeling Team')))),
                                       p(HTML(paste0(a(href = 'https://github.com/agb85/covid-19',
                                                       'Github Link')))),
                                       p("Updated 09/07/20"),
                                       width = 3
                                     ),
                                     
                                     
                                     # Plot(s)
                                     mainPanel(plotlyOutput("p5", height = "500px"), width = 9)
                                   ),
                            ),
                            column(1)
                          ),
                          
                          br(),
                          
                          fluidRow(
                            column(1),
                            column(3, chooseSliderSkin("Nice"),
                                   sliderInput(inputId="ef1_2",  
                                               label="What is the level of social distancing among those age 65 and under? (from present onward)",
                                               value=0.65, min=0, max=1, width='100%', step = .01)),
                            column(3, 
                                   sliderInput(inputId="ef4p", 
                                               label="What proportion of adults age 65+ practice high social distancing?",
                                               value=0.65, min=0, max=1, width='100%', step = .01)),
                            
                            column(3,
                                   sliderInput(inputId="ef1_3",  
                                               label="What is the level of social distancing among those age 65 and under for the time around and during the winter holidays?",
                                               value=0.4, min=0, max=1, width='100%', step = .01)),
                            column(1)
                            
                          ),
                          
                          fluidRow(
                            column(1),
                            
                            column(3, sliderInput(inputId="maskc", label="What proportion of the population wears masks 
                                                  in public spaces? ",
                                                  value=0.7, min=0, max=1, width='100%', step = .01)),
                          
                        
                            
            
                            
                            # column(3, sliderInput(inputId="maskc", label="What proportion of the population wears masks 
                            #                       in public spaces? ",
                            #                       value=0.7, min=0, max=1, width='100%', step = .01)),
                            # column(3, 
                            #        sliderTextInput(inputId="pi",
                            #                        label = "How quickly are contacts successfully traced after case report?",
                            #                        grid = TRUE,
                            #                        choices = c("24 Hours", "48 Hours", "72 Hours"),
                            #                        width='100%')),
                            # column(3, 
                            #        sliderInput(inputId="kap", 
                            #                    label = "How many contacts are successfully traced per case?",
                            #                    value=0, min=0, max=5, width='100%', step = 1)),
                            column(1)
                            
                          ),
                          
                          fluidRow(
                            column(1),
                            column(3, p(tags$b("How quickly are contacts successfully traced after case report?"))),
                            column(3, p(tags$b("How many contacts are successfully traced per case?"))),
                            column(2, p(tags$b("Improve case detection and isolation"))),
                          ),
                          
                          fluidRow(
                            column(1),
                            column(3,
                                   sliderTextInput(inputId="pi",
                                                   label = "", # How quickly are contacts successfully traced after case report?
                                                   grid = TRUE,
                                                   choices = c("24 Hours", "48 Hours", "72 Hours"),
                                                   width='100%')),
                            column(3,
                                   sliderInput(inputId="kap",
                                               label = "", # How many contacts are successfully traced per case?
                                               value=0, min=0, max=5, width='100%', step = 1)),
                            
                            column(3,
                                   materialSwitch(inputId="ramp",
                                                  status = "success",
                                                  value = FALSE)),
                          ),
                          
                          br(),
                          
                          fluidRow(
                            column(1),
                            column(10,
                                   h4("Adjustable Model Parameters"), 
                                   p("This app allows you to generate “what if” scenarios and generate projections of 
                                   the course of COVID-19 in the coming months based on the particular scenario 
                                   selected. The different sliders include: "),
                                   p(HTML(paste0("<b>", "Social Distancing.", "</b>",
                                                 " Social distancing refers to the steps taken to reduce infectious 
                                                 contacts between people. We know the virus spreads from person to 
                                                 person when people are physically close to each other. Social 
                                                 distancing can include both policy measures to prevent contact and 
                                                 personal behaviors such as maintaining physical distance from people 
                                                 outside of one’s household and handwashing. In this model social 
                                                 distancing is modeled as a percent reduction in contacts. There are 
                                                 three sliders that address social distancing – you can adjust social 
                                                 distancing levels for the future, before and after the winter holiday season, 
                                                 and you can adjust the social distancing level around the winter holidays.
                                                 The winter holidays are assumed to begin the week before Thanksgiving and extend
                                                 until January 3rd. 
                                                 Based on CDC recommendations that older adults take extra 
                                                 precautions, you can adjust the proportion of adults age 65+ that 
                                                 maintain high levels of social distancing (High social distancing is 
                                                 equivalent to an 80% reduction in contacts in our model)."))),
                                  
                                   p(HTML(paste0("<b>", "Masks.", "</b>",
                                                 " A mask can prevent the spread of infections by containing droplets 
                                                 from an individual’s mouth or nose when they cough, sneeze, or talk. 
                                                 The proportion of the population wearing masks can be adjusted in the 
                                                 model."))),
                                   p(HTML(paste0("<b>", "Case detection and isolation.", "</b>",
                                                 " A key strategy for controlling an epidemic is identifying those who 
                                                 are infected and isolating them, usually at home, so that they do not 
                                                 infect others. You can turn on or off a scenario that involves 
                                                 increasing testing capacity and outreach to cases so that the percent 
                                                 of symptomatic cases successfully isolated increases by 5% per week up 
                                                 to a maximum of 80% isolated."))),
                                   p(HTML(paste0("<b>", "Contact tracing.", "</b>",
                                                 " Contact tracing involves identifying those people who have been close
                                                 enough to a contagious person to have possibly become infected. Such 
                                                 individuals are identified, contacted, and encouraged to isolate to 
                                                 prevent further spread. The success of contact tracing depends on its 
                                                 timing—how quickly contacts can be reached—and on its effectiveness—how
                                                 many contacts can be traced and agree to isolate. To introduce contact 
                                                 tracing into the model, increase the number of contacts successfully 
                                                 traced to greater than 0. You can vary the number of contacts 
                                                 successfully traced per case, and the time between case report and 
                                                 successful contact tracing.")))
                            ),
                            column(1)
                          ),
                 ),
                 
                 tabPanel("Documentation",
                          fluidRow(column(1),
                                   column(10,
                                          h4("Introduction"),
                                          p("This app provides a tool to model the COVID-19 pandemic as it affects the 
                                          population of the State of Colorado.  It is based on a mathematical 
                                          representation of the pandemic in Colorado and it is intended to provide an 
                                          understanding of the potential impact of the four key control measures: social
                                          distancing, mask wearing, identification and isolation of cases, and contact 
                                          tracing—the key tools used to limit the spread of the SARS-CoV-2 virus, which 
                                          causes the illness COVID-19.  It allows users to look into the future of the 
                                          pandemic in Colorado as it will be affected by how these tools are used. "),
                                          p(HTML(paste0("The underlying model is termed an SEIR model where S refers to
                                                        susceptible, E to exposure, I to infected, and R to recovered. 
                                                        The model that underlies the app has been developed by the 
                                                        Colorado COVID-19 Modeling Group, which includes public health 
                                                        scientists from the University of Colorado and Colorado State 
                                                        University. The technical details of the model are provided in 
                                                        Further Documentation, below, and recent reports by the Colorado
                                                        COVID-19 Modeling Group are available ",
                                                        a(href = "http://www.ucdenver.edu/academics/colleges/PublicHealth/coronavirus/Pages/Modeling-Results.aspx",
                                                          "here.",
                                                          target = "blank_"),
                                                        "  The model is periodically updated to continue to reflect the
                                                          situation in the state.  "))),
                                          p("This app allows you to generate “what if” scenarios and generate 
                                          projections of the course of COVID-19 in the coming months based on the 
                                          particular scenario selected.  We caution that this model, like all models, 
                                          incorporates assumptions and its outputs are subject to uncertainty.  We 
                                          recommend that in interpreting the numbers it produces consideration needs to 
                                          be given to various sources of uncertainty that are described in the 
                                          documentation.  Nonetheless, the general patterns described in the app’s 
                                          outputs should be useful for planning. "),
                                          h4("What is modeled"),
                                          p("The output of the model is the epidemic curve, which describes the course 
                                          of the pandemic over time.  The app allows you to see the epidemic curve with 
                                          four different measures of COVID-19: 1) symptomatic infections; 2) 
                                          hospitalizations that don’t require a stay in an intensive care unit (ICU); 3)
                                          hospitalizations that do require ICU care; and 4) deaths.  Each of these 
                                          measures is important for tracking the pandemic.  Of course, we want to avoid
                                          deaths and to do so we need enough regular hospital beds and also critical 
                                          care capacity.  For hospital planning, it is important to know how much bed 
                                          capacity may be needed.  For public health, we need to know the number of 
                                          people with symptomatic infections that will need to be identified and 
                                          isolated, and their contacts need to be traced.  "),
                                          p("To control the epidemic, a set of approaches has been used since March 2020
                                          and will continue to be used. Each of these can be changed in the app to 
                                          explore the effects of different approaches on the epidemic curve. These 
                                            measures include:"),
                                          tags$ul(
                                            tags$li("Social distancing: This term refers to the steps taken to
                                                    reduce infectious contacts between people. We know the virus spreads
                                                    from person to person when people are physically close to each 
                                                    other. Social distancing can include both policy measures to prevent
                                                    contact (e.g., closing workplaces and schools) and personal 
                                                    behaviors such as maintaining physical distance from people outside 
                                                    of one’s household and handwashing. In this model social distancing 
                                                    is modeled as a percent reduction in contacts."),
                                            tags$li("Identification and isolation of infected individuals: A key 
                                                    strategy for controlling an epidemic is identifying those who are 
                                                    infected and isolating them, usually at home, so that they do not 
                                                    infect others."),
                                            tags$li("Contact tracing:  This has long been a public health 
                                                    strategy for controlling outbreaks of infection and epidemics.  
                                                    Contact tracing involves identifying those people who have been 
                                                    close enough to a contagious person to have possibly become 
                                                    infected.  Such individuals are identified, contacted, isolated, 
                                                    and tested for infection. The success of contact tracing depends 
                                                    on its timing—how quickly contacts can be reached—and on its 
                                                    effectiveness—how many contacts can be traced and agree to 
                                                    isolate."),
                                            tags$li("Masks: A mask can prevent the spread of infections by containing 
                                                    droplets from an individual’s mouth or nose when they cough, sneeze,
                                                    or talk. The proportion of the population wearing masks can be 
                                                    adjusted in the model.")
                                          ),
                                          h4("Further documentation"),
                                          p(HTML(paste0("Documentation of the modeling approaches are available here (", 
                                                        a(href = "SEIR Documentation_20200821_updating.pdf", "LINK", target = "_blank"), 
                                                        ") and our model equations are available here (", 
                                                        a(href = "covid19model.pdf", 
                                                          "LINK", target = "_blank"),
                                                        ")."))
                                          ),
                                          p("Our model is re-fit and time-varying model parameters re-estimated weekly."),
                                          tags$ul(
                                            tags$li(HTML(paste0("Model fit and parameter estimates June 30, 2020 (", 
                                                                a(href = "ParameterEstimatesAndModelFit_20200630.pdf", 
                                                                  "LINK", target = "_blank"),
                                                                ")"))),
                                            tags$li(HTML(paste0("Model fit and parameter estimates July 10, 2020 (", 
                                                                a(href = "ParameterEstimatesAndModelFit_20200710.pdf", 
                                                                  "LINK", target = "_blank"),
                                                                ")"))),
                                            tags$li(HTML(paste0("Model fit and parameter estimates July 14, 2020 (", 
                                                                a(href = "ParameterEstimatesAndModelFit_20200714.pdf", 
                                                                  "LINK", target = "_blank"),
                                                                ")"))),
                                            tags$li(HTML(paste0("Model fit and parameter estimates July 20, 2020 (", 
                                                                a(href = "ParameterEstimatesAndModelFit_20200720.pdf", 
                                                                  "LINK", target = "_blank"),
                                                                ")"))),
                                            tags$li(HTML(paste0("Model fit and parameter estimates July 27, 2020 (", 
                                                                a(href = "ParameterEstimatesAndModelFit_20200727.pdf", 
                                                                  "LINK", target = "_blank"),
                                                                ")"))),
                                            tags$li(HTML(paste0("Model fit and parameter estimates August 03, 2020 (", 
                                                                a(href = "ParameterEstimatesAndModelFit_20200803.pdf", 
                                                                  "LINK", target = "_blank"),
                                                                ")"))),
                                            tags$li(HTML(paste0("Model fit and parameter estimates August 10, 2020 (", 
                                                                a(href = "ParameterEstimatesAndModelFit_20200810.pdf", 
                                                                  "LINK", target = "_blank"),
                                                                ")"))),
                                            tags$li(HTML(paste0("Model fit and parameter estimates August 17, 2020 (", 
                                                                a(href = "ParameterEstimatesAndModelFit_20200817.pdf", 
                                                                  "LINK", target = "_blank"),
                                                                ")"))),
                                            tags$li(HTML(paste0("Model fit and parameter estimates August 24, 2020 (", 
                                                                a(href = "ParameterEstimatesAndModelFit_20200824.pdf", 
                                                                  "LINK", target = "_blank"),
                                                                ")"))),
                                            tags$li(HTML(paste0("Model fit and parameter estimates August 31, 2020 (", 
                                                                a(href = "ParameterEstimatesAndModelFit_20200831.pdf", 
                                                                  "LINK", target = "_blank"),
                                                                ")"))),
                                            tags$li(HTML(paste0("Model fit and parameter estimates September 07, 2020 (", 
                                                                a(href = "ParameterEstimatesAndModelFit_20200907.pdf", 
                                                                  "LINK", target = "_blank"),
                                                                ")")))
                                          ),
                                          p(HTML(paste0("Code for our models is posted on Github (", 
                                                        a(href = "https://github.com/agb85/covid-19", 
                                                          "LINK", target = "_blank"),
                                                        ")"))),
                                   ),
                                   column(1)
                          )
                 )
)

server <- function(input, output) {
  
  output$p5 <- renderPlotly({
    
    Cp <- 5840795
    n1 <- 1513005
    n2 <- 1685869
    n3 <- 1902963
    n4 <- 738958
    
    parms <- c(beta = 0.4793, # transmission rate
               gamma = 1/9,
               alpha = 4,
               Cp = 5840795, # called back from population spreadsheet
               n1 = 1513005,
               n2 = 1685869,
               n3 = 1902963,
               n4 = 738958,
               ef1_2 = input$ef1_2,
               ef1_3 = input$ef1_3,
               ef1_4 = input$ef1_2,
               ef4p =  input$ef4p, #proportion of adults over 65 social distancing at 80%
               ef2_2 = input$ef1_2,
               ef2_3 = input$ef1_3,
               ef2_4 = input$ef1_2,
               ef3_2 = input$ef1_2,
               ef3_3 = input$ef1_3,
               ef3_4 = input$ef1_2,
               ef4_2 = input$ef1_2,
               ef4_3 = input$ef1_3,
               ef4_4 = input$ef1_2,
               ef1 = 0,
               ef2 = 0,
               ef3 = 0,
               ef4 = 0,
               dnh1 = 0.000080592, dnh2 = 0.000144394, dnh3 = 0.001271003, dnh4 = 0.03389901,
               hlos1=4.05, hlos2 = 5.49, hlos3=7.36, hlos4=10.02,
               clos1=6.35, clos2=9.91, clos3=13.47, clos4=10.82,
               hosp1 = 0.01643, ## percent of symptomatic cases hospitalized
               hosp2= 0.02909127, hosp3= 0.043322, hosp4= 0.072384,
               cc1= 0.003964, ## percent of cases requiring CC
               cc2= 0.00751899, cc3= 0.0181659, cc4= 0.0303687,
               dh1 = 0, dh2 = 0, dh3 = 0.0086, dh4 = 0.144959,
               dc1 = 0.05714, dc2 = 0.06154, dc3 = 0.1674, dc4 = 0.39197,
               cap = 1800, #Colorado ICU capacity
               pS1 = 0.110023, ## proportion of infectious individuals symptomatic (0-19)
               pS2 = 0.35705, ## proportion of infectious individuals symptomatic (20-39)
               pS3 = 0.561205, ## proportion of infectious individuals symptomatic (40-64)
               pS4 = 0.774879, ## proportion of infectious individuals symptomatic (65+)
               siI = 0.438,## Proportion of symptomatic individuals self isolate
               siI1 = 0.3,
               lambda = 1.395, ##difference in infectiousness symptomatic/asymptomatic
               mag1 = 0,
               mag2 = 0.3,
               mag2a = 0.7738,
               mag3 = 0.7792,
               mag3a = 0.7700,
               mag4 = 0.8533,
               mag5 = 0.8126,
               mag6 = 0.9395,
               mag6a = 0.5997,
               mag6b = 0.3771,
               mag7 = 0.6899,
               mag8 = 0.8375,
               mag9 = 0.7155,
               mag10 = 0.6977,
               mag11 = 0.5905,
               traj = 0.5905,
               t1  = 41,
               t2  = 52,
               t2a = 59,
               t3 = 66,
               t3a = 80,
               t4 = 94,
               t5 = 108,
               t6 = 122,
               t6a = 136,
               t6b = 150,
               t7 = 164,
               t8 = 178,
               t9 = 192,
               t10 = 206,
               t11 = 220,
               ttraj = 237,
               tproject = 253,
               tschool = 222,
               tpa = 302,
               tpb = 346,
               ramp = ifelse(input$ramp, .00407, 0),
               maska = 0.5,
               maskb = 0.7,
               maskc = input$maskc, #proportion wearing masks for projections
               kap = input$kap, #average number of contacts traced per detected case
               pCT = 0.4, #proportion of identified cases with contacts traced
               pi = case_when(input$pi == "72 Hours" ~ 0.364, #probability a contact traced infected individual is isolated before infecting other susceptibles
                              input$pi == "48 Hours" ~ 0.4,
                              input$pi == "24 Hours" ~ 0.455),
               om = 0.061, #probability a contact traced individual is infected
               pID = 0.38, ## proportion of infections identified
               temp_on = 0
    )
    
    ## Run model for CC, H, and I
    dt      <- seq(0, 500, 1)
    
    inits      <- c(S1 = n1 - 1, E1 = 0, I1 = 1, II1 = 0, Ih1 = 0, Ic1 = 0, A1 = 0, R1 = 0, Rh1 = 0, Rc1 = 0, D1 = 0,
                    S2 = n2,     E2 = 0, I2 = 0, II2 = 0, Ih2 = 0, Ic2 = 0, A2 = 0, R2 = 0, Rh2 = 0, Rc2 = 0, D2 = 0,
                    S3 = n3,     E3 = 0, I3 = 0, II3 = 0, Ih3 = 0, Ic3 = 0, A3 = 0, R3 = 0, Rh3 = 0, Rc3 = 0, D3 = 0,
                    S4 = n4,     E4 = 0, I4 = 0, II4 = 0, Ih4 = 0, Ic4 = 0, A4 = 0, R4 = 0, Rh4 = 0, Rc4 = 0, D4 = 0)
    
    out1 <- lsoda(inits, dt, seir1, parms = parms)
    out <- as.data.frame(out1)
    
    ## Calculate daily totals for CC, H, I
    out$dailyInfections <-  round(out$I1 + out$I2 + out$I3 + out$I4, 0)
    out$dailyHospitalizations <- round(out$Ih1+out$Ih2+out$Ih3 + out$Ih4, 0)
    out$dailyCriticalCare <- round(out$Ic1 + out$Ic2 + out$Ic3 + out$Ic4, 0)
    out$cumulativeInfections <- round(out$R1 + out$R2 + out$R3 + out$R4 +
                                        out$I1 + out$I2 + out$I3 + out$I4)
    out$cumulativeHospitalizations <- round(out$Rh1 + out$Rh2 + out$Rh3 + out$Rh4 +
                                              out$Ih1 + out$Ih2 + out$Ih3 + out$Ih4)
    out$cumulativeCriticalCare <- round(out$Rc1 + out$Rc2 + out$Rc3 + out$Rc4 +
                                          out$Ic1 + out$Ic2 + out$Ic3 + out$Ic4)
    out$date <- as.Date(out$time, format="%m/%d/%Y", origin="01/24/2020")
    
    out <- out %>%
      select(date, dailyInfections, dailyHospitalizations, dailyCriticalCare,
             cumulativeInfections, cumulativeHospitalizations, cumulativeCriticalCare)
    
    ## Run model for deaths
    initsD    <- c(S1 = n1 - 1, E1 = 0, I1 = 1, II1 = 0, Ih1 = 0, A1 = 0, R1 = 0, Rh1 = 0, Ic = 0, Rc = 0, D = 0,
                   S2 = n2,     E2 = 0, I2 = 0, II2 = 0, Ih2 = 0, A2 = 0, R2 = 0, Rh2 = 0,
                   S3 = n3,     E3 = 0, I3 = 0, II3 = 0, Ih3 = 0, A3 = 0, R3 = 0, Rh3 = 0,
                   S4 = n4,     E4 = 0, I4 = 0, II4 = 0, Ih4 = 0, A4 = 0, R4 = 0, Rh4 = 0 )
    
    outD <- lsoda(initsD, dt, seir1D, parms = parms) %>% as.data.frame()
    
    ## Calculate cumulative and daily deaths
    outD$cumulativeDeaths <- round(outD$D, 0)
    outD$dailyDeaths <- c(outD$cumulativeDeaths[1], diff(outD$cumulativeDeaths))
    outD$date <- as.Date(outD$time, format="%m/%d/%Y", origin="01/24/2020")
    
    outD <- outD %>%
      select(date, cumulativeDeaths, dailyDeaths)
    
    ## Combine CC, H, I and Deaths outputs to pass to function
    out <- full_join(out, outD, by = "date") 
    
    ## Plot
    p4 <- mainBarPlot(out, input$Infections, input$Hospitalizations,
                      input$CriticalCare, input$Deaths, input$Daily, input$ICUCapLine,
                      input$CTCapLine)
    
    print(p4)
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)


