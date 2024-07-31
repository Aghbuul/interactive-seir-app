library(shiny)
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(plotly)
library(shinythemes)
library(shinyBS)
library(shinyjs)


# SEIRS model function
seirs_comprehensive_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    total_pop <- S_young + E_young + I_young + R_young + Q_young + C_young +
      S_old + E_old + I_old + R_old + Q_old + C_old
    
    if (time >= vacc_start & time <= vacc_end) {
      vacc_young <- vacc_rate * S_young
      vacc_old <- vacc_rate * S_old
    } else {
      vacc_young <- 0
      vacc_old <- 0
    }
    
    if (time >= quarantine_start & time <= quarantine_end) {
      new_quarantined_young <- quarantine_rate * I_young
      new_quarantined_old <- quarantine_rate * I_old
    } else {
      new_quarantined_young <- 0
      new_quarantined_old <- 0
    }
    
    # ODEs for young compartment
    dS_young <- birth_rate * total_pop - beta_young * S_young * (I_young + I_old + kappa * C_young + kappa * C_old) - death_rate * S_young + delta * R_young - vacc_young
    dE_young <- beta_young * S_young * (I_young + I_old + kappa * C_young + kappa * C_old) - sigma * E_young - death_rate * E_young
    dI_young <- (1-p) * sigma * E_young - gamma * I_young - death_rate * I_young - new_quarantined_young
    dR_young <- gamma * I_young - delta * R_young - death_rate * R_young + vacc_young + gamma * Q_young + gamma * C_young
    dQ_young <- new_quarantined_young - gamma * Q_young - death_rate * Q_young
    dC_young <- p * sigma * E_young - gamma * C_young - death_rate * C_young
    
    # ODEs for old compartment
    dS_old <- - beta_old * S_old * (I_young + I_old + kappa * C_young + kappa * C_old) - death_rate * S_old + delta * R_old - vacc_old
    dE_old <- beta_old * S_old * (I_young + I_old + kappa * C_young + kappa * C_old) - sigma * E_old - death_rate * E_old
    dI_old <- (1-p) * sigma * E_old - gamma * I_old - death_rate * I_old - new_quarantined_old
    dR_old <- gamma * I_old - delta * R_old - death_rate * R_old + vacc_old + gamma * Q_old + gamma * C_old
    dQ_old <- new_quarantined_old - gamma * Q_old - death_rate * Q_old
    dC_old <- p * sigma * E_old - gamma * C_old - death_rate * C_old
    
    return(list(c(dS_young, dE_young, dI_young, dR_young, dQ_young, dC_young,
                  dS_old, dE_old, dI_old, dR_old, dQ_old, dC_old)))
  })
}


# UI
ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  tags$head(
    tags$style(HTML("
    .dark-mode {
      background-color: #2d2d2d;
      color: white;
    }
    .light-mode {
      background-color: #ffffff;
      color: black;
    }
    .dark-mode .shiny-input-container {
      color: white;
    }
    .dark-mode .form-control {
      background-color: #2d2d2d;
      color: white;
      border: 1px solid #555;
    }
    .dark-mode .selectize-input {
      background-color: #2d2d2d;
      color: white;
    }
    .dark-mode label, .dark-mode h3 {
      color: white;
    }
    .dark-mode .well {
      background-color: #3b3b3b;
    }
  "))
  ),
  titlePanel("SEIRS Model Simulation"),
  sidebarLayout(
    sidebarPanel(
      h3("Model Parameters"),
      checkboxInput("dark_mode", "Dark Mode", value = FALSE),
      sliderInput("beta_young", "Transmission rate (young):", min = 0, max = 1, value = 0.1, step = 0.01),
      bsTooltip("beta_young", "Transmission rate for the young population", "right"),
      sliderInput("beta_old", "Transmission rate (old):", min = 0, max = 1, value = 0.1, step = 0.01),
      bsTooltip("beta_old", "Transmission rate for the old population", "right"),
      sliderInput("sigma", "Rate of progression from exposed to infected:", min = 0, max = 1, value = 0.05, step = 0.01),
      bsTooltip("sigma", "Rate at which exposed individuals become infected", "right"),
      sliderInput("gamma", "Recovery rate:", min = 0, max = 1, value = 0.01, step = 0.01),
      bsTooltip("gamma", "Rate at which infected individuals recover", "right"),
      sliderInput("delta", "Loss of immunity rate:", min = 0, max = 1, value = 0.01, step = 0.01),
      bsTooltip("delta", "Rate at which recovered individuals lose immunity", "right"),
      sliderInput("birth_rate", "Birth rate:", min = 0, max = 0.1, value = 0.01, step = 0.001),
      bsTooltip("birth_rate", "Birth rate in the population", "right"),
      sliderInput("death_rate", "Natural death rate:", min = 0, max = 0.1, value = 0.01, step = 0.001),
      bsTooltip("death_rate", "Natural death rate in the population", "right"),
      sliderInput("vacc_rate", "Vaccination rate:", min = 0, max = 1, value = 0.05, step = 0.01),
      bsTooltip("vacc_rate", "Rate of vaccination in the population", "right"),
      sliderInput("vacc_start", "Vaccination start time:", min = 0, max = 100, value = 5, step = 1),
      bsTooltip("vacc_start", "Time when vaccination starts", "right"),
      sliderInput("vacc_end", "Vaccination end time:", min = 0, max = 100, value = 20, step = 1),
      bsTooltip("vacc_end", "Time when vaccination ends", "right"),
      sliderInput("quarantine_rate", "Quarantine rate:", min = 0, max = 1, value = 0.2, step = 0.01),
      bsTooltip("quarantine_rate", "Rate of quarantining infected individuals", "right"),
      sliderInput("quarantine_start", "Quarantine start time:", min = 0, max = 100, value = 20, step = 1),
      bsTooltip("quarantine_start", "Time when quarantine starts", "right"),
      sliderInput("quarantine_end", "Quarantine end time:", min = 0, max = 100, value = 40, step = 1),
      bsTooltip("quarantine_end", "Time when quarantine ends", "right"),
      sliderInput("kappa", "Relative transmissibility of carriers:", min = 0, max = 1, value = 0.5, step = 0.01),
      bsTooltip("kappa", "Relative transmissibility of carriers compared to infected individuals", "right"),
      sliderInput("p", "Proportion progressing to carrier state:", min = 0, max = 1, value = 0.3, step = 0.01),
      bsTooltip("p", "Proportion of exposed individuals progressing to the carrier state", "right"),
      selectInput("plot_choice", "Select plot type:", choices = c("Plots by Age" = "age", "Plots by Compartment" = "compartment", "Interventions Comparison" = "intervention"))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plots by Age", plotlyOutput("age_plot", width = "100%", height = "600px")),
        tabPanel("Plots by Compartment", plotlyOutput("compartment_plot", width = "100%", height = "600px")),
        tabPanel("Interventions Comparison", plotlyOutput("intervention_plot", width = "100%", height = "600px"))
      )
    )
  )
)


# Server
server <- function(input, output, session) {
  observe({
    if (input$dark_mode) {
      runjs("document.body.classList.add('dark-mode');")
      runjs("document.body.classList.remove('light-mode');")
    } else {
      runjs("document.body.classList.add('light-mode');")
      runjs("document.body.classList.remove('dark-mode');")
    }
  })
  
  # Parameters
  parameters <- reactive({
    validate(
      need(input$vacc_start < input$vacc_end, "Vaccination start time must be less than end time"),
      need(input$quarantine_start < input$quarantine_end, "Quarantine start time must be less than end time")
    )
    
    c(
      beta_young = input$beta_young, beta_old = input$beta_old, sigma = input$sigma, gamma = input$gamma, delta = input$delta,
      birth_rate = input$birth_rate, death_rate = input$death_rate, vacc_rate = input$vacc_rate, vacc_start = input$vacc_start,
      vacc_end = input$vacc_end, quarantine_rate = input$quarantine_rate, quarantine_start = input$quarantine_start,
      quarantine_end = input$quarantine_end, kappa = input$kappa, p = input$p
    )
  })
  
  output_data <- reactive({
    initial_state <- c(S_young = 50, E_young = 1, I_young = 0, R_young = 0, Q_young = 0, C_young = 0,
                       S_old = 49, E_old = 0, I_old = 0, R_old = 0, Q_old = 0, C_old = 0)
    
    times <- seq(0, 100, by = 1)
    
    tryCatch({
      result <- ode(y = initial_state, times = times, func = seirs_comprehensive_model, parms = parameters())
      as.data.frame(result)
    }, error = function(e) {
      showNotification(paste("Error in ODE solver:", e$message), type = "error")
      NULL
    })
  })
  
  # Theme
  theme_dynamic <- reactive({
    if (input$dark_mode) {
      theme_minimal(base_size = 15) +
        theme(
          plot.background = element_rect(fill = "#2d2d2d", color = NA),
          panel.background = element_rect(fill = "#2d2d2d", color = NA),
          text = element_text(color = "white"),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white"),
          legend.background = element_rect(fill = "#2d2d2d", color = NA),
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white")
        )
    } else {
      theme_minimal(base_size = 15)
    }
  })
  
  output$simple_plot <- renderPlotly({
    output_df <- as.data.frame(output_data())
    output_long <- reshape2::melt(output_df, id.vars = "time")
    
    output_long$age_group <- factor(ifelse(grepl("_young", output_long$variable), "Young", "Old"), levels = c("Young", "Old"))
    output_long$compartment <- factor(sub("_young|_old", "", output_long$variable), 
                                      levels = c("S", "E", "I", "C", "Q", "R"),
                                      labels = c("Susceptible", "Exposed", "Infected", "Carrier", "Quarantined", "Recovered"))
    
    p <- ggplot(output_long, aes(x = time, y = value, color = variable)) +
      geom_line(size = 1) +
      theme_dynamic() +
      labs(x = "Time", y = "Number of Individuals", color = "Compartment") +
      ggtitle("Comprehensive SEIRS Model Dynamics") +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggplotly(p)
  })
  
  # Age plot
  output$age_plot <- renderPlotly({
    output_df <- as.data.frame(output_data())
    output_long <- reshape2::melt(output_df, id.vars = "time")
    
    output_long$age_group <- factor(ifelse(grepl("_young", output_long$variable), "Young", "Old"), levels = c("Young", "Old"))
    output_long$compartment <- factor(sub("_young|_old", "", output_long$variable), 
                                      levels = c("S", "E", "I", "C", "Q", "R"),
                                      labels = c("Susceptible", "Exposed", "Infected", "Carrier", "Quarantined", "Recovered"))
    
    p <- ggplot(output_long, aes(x = time, y = value, color = compartment)) +
      geom_line(size = 1, linetype = "solid") + # Set linetype to solid
      facet_wrap(~age_group, scales = "free_y") +
      theme_dynamic() +
      labs(x = "Time", y = "Number of Individuals", color = "Compartment") +
      ggtitle("Comprehensive SEIRS Model Dynamics") +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") # Adjust legend position here
    
    ggplotly(p)
  })
  
  # Compartment plot
  output$compartment_plot <- renderPlotly({
    output_df <- as.data.frame(output_data())
    output_long <- reshape2::melt(output_df, id.vars = "time")
    
    output_long$age_group <- factor(ifelse(grepl("_young", output_long$variable), "Young", "Old"), levels = c("Young", "Old"))
    output_long$compartment <- factor(sub("_young|_old", "", output_long$variable), 
                                      levels = c("S", "E", "I", "C", "Q", "R"),
                                      labels = c("Susceptible", "Exposed", "Infected", "Carrier", "Quarantined", "Recovered"))
    
    p <- ggplot(output_long, aes(x = time, y = value, color = age_group, linetype = age_group)) +
      geom_line(size = 1) +
      facet_wrap(~compartment, scales = "free_y") +
      theme_dynamic() +
      labs(x = "Time", y = "Number of Individuals", color = "Age Group", linetype = "Age Group") +
      ggtitle("Comprehensive SEIRS Model Dynamics") +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
    
    ggplotly(p)
  })
  
  # Intervention plot
  output$intervention_plot <- renderPlotly({
    initial_state <- c(S_young = 50, E_young = 1, I_young = 0, R_young = 0, Q_young = 0, C_young = 0,
                       S_old = 49, E_old = 0, I_old = 0, R_old = 0, Q_old = 0, C_old = 0)
    times <- seq(0, 100, by = 1)
    

    # Calculate the epidemic dynamics without vaccination
    parameters_no_vaccination <- parameters()
    parameters_no_vaccination["vacc_rate"] <- 0
    output_no_vaccination <- ode(y = initial_state, times = times, func = seirs_comprehensive_model, parms = parameters_no_vaccination)
    
    # Calculate the epidemic dynamics with vaccination
    output_vaccination <- ode(y = initial_state, times = times, func = seirs_comprehensive_model, parms = parameters())
    
    # Convert to data frames
    output_no_vaccination_df <- as.data.frame(output_no_vaccination)
    output_vaccination_df <- as.data.frame(output_vaccination)
    
    # Calculate cumulative infections for both scenarios
    output_no_vaccination_df$Cumulative_Infected_no_vaccination <- cumsum(rowSums(output_no_vaccination_df[, c("I_young", "I_old")]))
    output_vaccination_df$Cumulative_Infected_vaccination <- cumsum(rowSums(output_vaccination_df[, c("I_young", "I_old")]))
    
    # Calculate cumulative percentage decrease in infections
    percentage_prevented_vaccination <- (output_no_vaccination_df$Cumulative_Infected_no_vaccination - output_vaccination_df$Cumulative_Infected_vaccination) / output_no_vaccination_df$Cumulative_Infected_no_vaccination * 100
    
    # Combine data for plotting
    df_prevented_vaccination <- data.frame(
      time = times,
      percentage_prevented = percentage_prevented_vaccination,
      Scenario = "Vaccination"
    )
    
    
    # Calculate the epidemic dynamics without quarantine
    parameters_no_quarantine <- parameters()
    parameters_no_quarantine["quarantine_rate"] <- 0
    output_no_quarantine <- ode(y = initial_state, times = times, func = seirs_comprehensive_model, parms = parameters_no_quarantine)
    
    # Calculate the epidemic dynamics with quarantine
    output_quarantine <- ode(y = initial_state, times = times, func = seirs_comprehensive_model, parms = parameters())
    
    # Convert to data frames
    output_no_quarantine_df <- as.data.frame(output_no_quarantine)
    output_quarantine_df <- as.data.frame(output_quarantine)
    
    # Calculate cumulative infections for both scenarios
    output_no_quarantine_df$Cumulative_Infected_no_quarantine <- cumsum(rowSums(output_no_quarantine_df[, c("I_young", "I_old")]))
    output_quarantine_df$Cumulative_Infected_quarantine <- cumsum(rowSums(output_quarantine_df[, c("I_young", "I_old")]))
    
    # Calculate cumulative percentage decrease in infections
    percentage_prevented_quarantine <- (output_no_quarantine_df$Cumulative_Infected_no_quarantine - output_quarantine_df$Cumulative_Infected_quarantine) / output_no_quarantine_df$Cumulative_Infected_no_quarantine * 100
    
    # Combine data for plotting
    df_prevented_quarantine <- data.frame(
      time = times,
      percentage_prevented = percentage_prevented_quarantine,
      Scenario = "Quarantine"
    )
    
    # Combine both datasets
    df_prevented_combined <- rbind(df_prevented_vaccination, df_prevented_quarantine)
    
    # Plot cumulative percentage decrease in infections over time for both scenarios
    p <- ggplot(df_prevented_combined, aes(x = time, y = percentage_prevented, color = Scenario)) +
      geom_line(size = 1.5) +
      theme_dynamic() +
      labs(x = "Time", y = "Cumulative Percentage Decrease in Infected Individuals", title = "Cumulative Percentage Decrease in Infected Individuals Over Time") +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)) # Add margin to y-axis title
      )
    
    ggplotly(p)
  })
}


shinyApp(ui = ui, server = server)