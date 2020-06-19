mainBarPlot <- function(out,
                        Infections,
                        Hospitalizations,
                        CriticalCare,
                        Deaths,
                        Daily,
                        ICUCapLine,
                        CTCapLine) {
    
    ## Handle case where no variables are selected
    if (all(!Deaths, !Hospitalizations, !CriticalCare, !Infections)) {
        CriticalCare <- TRUE
        Hospitalizations <- TRUE
        Deaths <- TRUE
    }
    
    ## Select daily or cumulative data
    columns <- c()
    if(Daily) {
        if (Infections) { columns <- c(columns, "dailyInfections") }
        if (Hospitalizations) { columns <- c(columns, "dailyHospitalizations") }
        if (CriticalCare) { columns <- c(columns, "dailyCriticalCare") }
        if (Deaths) { columns <- c(columns, "dailyDeaths") }
    } else {
        if (Infections) { columns <- c(columns, "cumulativeInfections") }
        if (Hospitalizations) { columns <- c(columns, "cumulativeHospitalizations") }
        if (CriticalCare) { columns <- c(columns, "cumulativeCriticalCare") }
        if (Deaths) { columns <- c(columns, "cumulativeDeaths") }
    }
    
    selectColumns <- c("date", columns)
    
    ## Prepare data for plot
    plotData <- out %>%
        select(all_of(selectColumns)) %>%
        gather(key = "Type", value = "Count", columns) %>%
        arrange(date) %>%
        group_by(date) %>%
        mutate("Sum" = ifelse(sum(Count) >= 0, sum(Count), 0),
               "Count" = ifelse(Count >= 0, Count, 0)) %>%
        mutate("Type" = case_when(.data$Type == "dailyInfections" ~ "Symptomatic Infections",
                                  .data$Type == "dailyCriticalCare" ~ "Critical Care",
                                  .data$Type == "dailyHospitalizations" ~ "Non-ICU Hospitalizations",
                                  .data$Type == "dailyDeaths" ~ "Deaths",
                                  .data$Type == "cumulativeInfections" ~ "Symptomatic Infections",
                                  .data$Type == "cumulativeCriticalCare" ~ "Critical Care",
                                  .data$Type == "cumulativeHospitalizations" ~ "Non-ICU Hospitalizations",
                                  .data$Type == "cumulativeDeaths" ~ "Deaths"))
    
    ## Update factor level's order to control order data 'stacks' in plot
    if(CriticalCare) {plotData$Type <- relevel(as.factor(plotData$Type), "Critical Care")}
    if(Hospitalizations) {plotData$Type <- relevel(as.factor(plotData$Type), "Non-ICU Hospitalizations")}
    if(Infections) {plotData$Type <- relevel(as.factor(plotData$Type), "Symptomatic Infections")}
    
    ## Create ggplot
    plot <- ggplot(data = plotData, aes(fill = Type, x = date, y = Count)) +
        geom_bar(position = "stack", stat = "identity") +  
        scale_y_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10)) +
        coord_cartesian(xlim=c(as.Date("2020-01-24"), as.Date("2021-6-07"))) +
        ylab("Total Individuals") +
        xlab("") +
        scale_x_date(date_labels="%b %Y", date_breaks  ="1 month") +
        theme(axis.text.x = element_text(angle = 90)) + 
        scale_fill_manual(" ", 
                          values=c("Critical Care" = "#E7B800",
                                   "Non-ICU Hospitalizations" = "#00AFBB",
                                   "Symptomatic Infections" = "#5CB85C",
                                   "Deaths" = "#D9534F"))
    
    ## Add line for ICU Capacity if selected (line label added differently when Infections included in model)
    ## Adjust label position based on maximum values being plotted
    offset <- max(max(plotData$Sum) / 22, 200)
    offsetCTC <- ifelse(offset > 800 & ICUCapLine, -offset, offset)
    
    ## Add lines
    if (Daily) {
        if (ICUCapLine & Infections){
            plot <- plot + geom_hline(aes(yintercept = 1800), linetype = "dashed", size = .4, color = "black") +
                geom_text(aes(as.Date("2020-04-01"), 1800, label = "ICU Capacity"), size = 3, nudge_y = offset, color = "black")
        } else if (ICUCapLine) {
            plot <- plot + geom_hline(aes(yintercept = 1800), linetype = "dashed", size = .4, color = "black") +
                geom_text(aes(as.Date("2020-04-01"), 1800, label = "ICU Capacity"), size = 4, nudge_y = offset, color = "black")
        }

        # Add line for Contact Tracing Capacity if selected
        if (CTCapLine & Infections){
            plot <- plot + geom_hline(aes(yintercept = 500), linetype = "dashed", size = .4, color = "blue") +
                geom_text(aes(as.Date("2021-01-01"), 500, label = "Contact Tracing Capacity"), size = 3, nudge_y = offsetCTC, color = "blue")
        } else if (CTCapLine) {
            plot <- plot + geom_hline(aes(yintercept = 500), linetype = "dashed", size = .4, color = "blue") +
                geom_text(aes(as.Date("2021-01-01"), 500, label = "Contact Tracing Capacity"), size = 4, nudge_y = offsetCTC, color = "blue")
        }
    }
    
    ## Add appropriate title
    if (Daily) {
        plot <- plot + ggtitle("Daily Totals")
    } else {
        plot <- plot + ggtitle("Cumulative Totals")
    }
    
    ## Convert ggplot to Plotly
    ggplotly(plot, dynamicTicks = TRUE) %>%
        hide_legend()
}