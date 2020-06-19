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
ui <- fluidPage(
  
  # Application title
  titlePanel("Modeling COVID-19 in Colorado"),
  
  # App introductory text
  fluidRow(
    column(12, 
           h3("How to use this app."),
           p("This app allows the user to estimate epidemic trajectories 
             under scenarios that can be altered interactively. In this 
             version of the app, the only aspects that are changeable are 
             those concerning social distancing, mask wearing, and contact tracing. With the sliders at the 
             bottom of the page, you can alter the efficacy of certain social 
             distancing levels, when those take effect, and whether they 
             differentially apply to those over/under 65 years of age."),
           p(HTML(paste0('For full details of the implementation of the other 
                         aspects of the model see the Colorado COVID-19 
                         Modeling Report ',
                         a(href = 'http://www.ucdenver.edu/academics/colleges/PublicHealth/coronavirus/Pages/Modeling-Results.aspx', 
                           'available here'),'.')))
           # h4("Cases and Hospitalizations"),
           # p(HTML("The first set of outputs indicate the projected number of 
           #        new <i>daily</i> cases and hospitalizations.")),
           # h4("Deaths"),
           # p(HTML("Additionally, the model also provides projections of either the daily or  
           #        <i>cumulative</i> number of deaths up to the specified date."))
           )
  ),
  
  br(),
  
  # Plot(s) and checkboxes
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
      p(HTML(paste0(a(href = 'mailto:kathryn.colborn@cuanschutz.edu',
                      'Dr. Kathryn Colborn')))),
      p(HTML(paste0(a(href = 'https://github.com/agb85/covid-19',
                       'Github Link')))),
      p("Updated 6/17/20"),
      width = 3
    ),
    
    
    #position = "right",
    
    # Plot(s)
    mainPanel(plotlyOutput("p5", height = "500px"), width = 9)
  ),
  
  fluidRow(
    column(12,
           h4("Adjustable Model Parameters"), 
           p("In the first row of sliders below, you can adjust the level of
             social distancing that applies to the general population (all 
             residents under age 65).
             This model uses a single parameter to summarize social distancing 
             as the percent decrease in effective contacts between susceptible 
             and infectious individuals. This parameter accounts for social 
             distancing policies intended to avoid contact altogether (e.g., 
             through workplace and school closures) as well as policies and 
             individual behaviors to reduce potential contact with the virus 
             (e.g., maintaining at least 6 feet of distance between people 
             outside of one's household, and handwashing)."),
           p(HTML(paste0('The second row of sliders allows the level of 
                         social distancing to be independently adjusted for 
                         residents aged 65 and older. Since the expiration 
                         of the "Stay at Home" orders (on April 27) (',
                         a(href = 'https://covid19.colorado.gov/safer-at-home', 
                           'details here'),
                         ') recommendations for reductions in social 
                         contacts have varied by age.'))),
           p('The final row of sliders allows you to adjust what proportion 
             of all residents wear masks in public, and the nature of contact tracing implementation.
             To introduce contact tracing into the model, increase the number of contacts successfully traced
             and quarantined to greater than 0.'), 
    )
  ),
  
  br(),
  
  fluidRow(
    column(4, chooseSliderSkin("Nice"),
           sliderInput(inputId="ef1_2",
                       label="Current social distancing level among under 65",
                       value=0.65, min=0, max=1, width='100%', step = .01)),
    column(4, 
           sliderInput(inputId="ef4_2", 
                       label="Current social distancing level among over 65",
                       value=0.65, min=0, max=1, width='100%', step = .01)),
    
    column(4, sliderInput(inputId="maskb", label="Proportion using masks, starting 4/27",
                          value=0.5, min=0, max=1, width='100%', step = .01))

  ),
  
  fluidRow(
    column(4, 
           sliderInput(inputId="ef1_3", 
                       label="Social distancing level among under 65, starting 8/15",
                       value=0.65, min=0, max=1, width='100%', step = .01)),
    column(4, 
           sliderInput(inputId="ef4_3", 
                       label="Social distancing level among over 65, starting 8/15",
                       value=0.65, min=0, max=1, width='100%', step = .01)),
    column(4, 
           sliderTextInput(inputId="ramp",
                           label="Ramp-up of case detection and isolation",
                           grid = TRUE,
                           choices = c("No", "Yes"),
                           width='100%'))
  ),
  
  fluidRow(
    column(4, 
           sliderInput(inputId="kap", 
                       label="Average number of contacts successfully traced and quarantined per case, starting 5/27",
                       value=0, min=0, max=5, width='100%', step = 1)),
    column(4, 
           sliderTextInput(inputId="pi",
                           label="Average time lag between report and contact tracing",
                           grid = TRUE,
                           choices = c("24 Hours", "48 Hours", "72 Hours"),
                           width='100%'))
  )
)

server <- function(input, output) {
  
  output$p5 <- renderPlotly({
    
    Cp = 5840795
    n1 = 1513005
    n2 = 1685869
    n3 = 1902963
    n4 = 738958
    
    parms <- c(beta = 0.48, # transmission rate
               cap = 1800,
               gamma = 1/8,
               alpha = 5.1,
               Cp = Cp,
               n1 = n1,
               n2 = n2,
               n3 = n3,
               n4 = n4,
               dc = 0.5, # death rate from ICU
               ef1_1 = 0.7544,
               ef1_2 = input$ef1_2,
               ef1_3 = input$ef1_3,
               # ef2_1 = scen[i,c('ef2_1')],
               # ef2_2 = scen[i,c('ef2_2')],
               # ef2_3 = scen[i,c('ef2_3')],
               # ef3_1 = scen[i,c('ef3_1')],
               # ef3_2 = scen[i,c('ef3_2')],
               # ef3_3 = scen[i,c('ef3_3')],
               ef4_1 = 0.7544,
               ef4_2 = input$ef4_2,
               ef4_3 = input$ef4_3,
               ef1 = 0,
               ef2 = 0,
               ef3 = 0,
               ef4 = 0,
               pS1 = 0.110023, ## proportion of infectious individuals symptomatic (0-19)
               pS2 = 0.35705, ## proportion of infectious individuals symptomatic (20-39)
               pS3 = 0.561205, ## proportion of infectious individuals symptomatic (40-59)
               pS4 = 0.774879, ## proportion of infectious individuals symptomatic (60-79)
               pID = 0.25, ## proportion of symptomatic individuals Identified
               siI = 0.3,## Proportion of symptomatic individuals self isolate
               lambda = 1.65, ##difference in infectiousness symptomatic/asymptomatic
               hosp1 = 0.006522, 
               cc1 = 0.0016447,
               hosp2 = 0.031819196, 
               cc2 = 0.011989122,
               hosp3 = 0.06362, 
               cc3 = 0.02884,
               hosp4 = 0.101899, 
               cc4 = 0.064882,
               mag1 = 0.5205,
               mag2 = 0.8044,
               mag3 = 0.8150,
               mag4 = 0.9030,
               t1 = 41,
               t2 = 53,
               t3 = 63,
               t4 = 94,
               t5 = 107,
               t6 = 125,
               t7 = difftime(as.Date("2020-06-05"), as.Date("2020-01-23")),
               t8 = 205,
               ramp = ifelse(input$ramp == "Yes", .00407, 0),
               maska = 0.4,
               maskb = input$maskb,
               kap = input$kap, #average number of contacts traced per detected case
               pCT = 0.66, #proportion of identified cases with contacts traced
               pi = case_when(input$pi == "72 Hours" ~ 0.41667, #probability a contact traced infected individual is isolated before infecting other susceptibles
                              input$pi == "48 Hours" ~ 0.4545,
                              input$pi == "24 Hours" ~ 0.5),  
               om = 0.046 #probability a contact traced individual is infected
    )
    
    ## Run model for CC, H, and I
    dt      <- seq(0, 500, 1)
    inits      <- c(S1 = n1 - 1, E1 = 0, I1 = 1, II1 = 0, Ih1 = 0, Ic1 = 0, A1 = 0, R1 = 0, Rh1 = 0, Rc1 = 0, D1 = 0,
                    S2 = n2,     E2 = 0, I2 = 0, II2 = 0, Ih2 = 0, Ic2 = 0, A2 = 0, R2 = 0, Rh2 = 0, Rc2 = 0, D2 = 0,
                    S3 = n3,     E3 = 0, I3 = 0, II3 = 0, Ih3 = 0, Ic3 = 0, A3 = 0, R3 = 0, Rh3 = 0, Rc3 = 0, D3 = 0,
                    S4 = n4,     E4 = 0, I4 = 0, II4 = 0, Ih4 = 0, Ic4 = 0, A4 = 0, R4 = 0, Rh4 = 0, Rc4 = 0, D4 = 0)
    N  <- Cp
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
    initsD      <- c(S1 = n1 - 1, E1 = 0, I1 = 1, II1 = 0, Ih1 = 0, A1 = 0, R1 = 0, Rh1 = 0, Ic = 0, Rc = 0, D = 0,
                     S2 = n2,     E2 = 0, I2 = 0, II2 = 0, Ih2 = 0, A2 = 0, R2 = 0, Rh2 = 0,
                     S3 = n3,     E3 = 0, I3 = 0, II3 = 0, Ih3 = 0, A3 = 0, R3 = 0, Rh3 = 0,
                     S4 = n4,     E4 = 0, I4 = 0, II4 = 0, Ih4 = 0, A4 = 0, R4 = 0, Rh4 = 0)
    N  <- Cp
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
