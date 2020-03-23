## app.R ##
library(shinydashboard)
library(learnr)
knitr::opts_chunk$set(echo = FALSE)
require(socialmixr)
require(deSolve)
require(dplyr)
require(forcats)
require(ggplot2)
require(plotly)
require(tibble)
require(tidyr)
data(polymod)


ui <- dashboardPage(
  dashboardHeader(title = "Covid-19 Age-structured Model", titleWidth = 450),
  dashboardSidebar(disable=TRUE),
  dashboardBody(
    
    fluidRow(
      column(width = 3,
             box(width=NULL,title = "Model Inputs", status = "warning", "R0 controls the average infectiousness of a case:",
                 sliderInput("R0", "R0:", min = 1.0, max = 5.0, value = 2.8, step = 0.1)
             )
      ),
    # Boxes need to be put in a row (or column)
    fluidRow(
      column(width = 4,
             box( width=NULL, "The colors in these figures show the expected number of cases and deaths generated in the age group on the y-axis (assuming everyone is fully susceptible) by a single case in the age group on the x-axis. These are not weighted by the proportion of the population in each age group, so the rows sum to > R0."),
             box(title = "Age Specific Infectiousness", width=NULL,  plotOutput("agePlot")),
             box(title = "Age Specific Mortality Impact", width=NULL, plotOutput("ageMortPlot"))
         
      ),
      column(width = 4,
             box(width = NULL, "The lines in these figures show the age-specific COVID-19 prevalence and mortality over a 1-year period of uncontrolled transmission." 
             ),
             box(title = "Age Specific Prevalence", width=NULL, plotOutput("prevPlot")),
             box(title = "Age Specific Mortality ", width=NULL, plotOutput("deathPlot"))
      )
    )
  )
))

server <- function(input, output) {
  
  
  contact_data <- contact_matrix(polymod, 
                                 countries = "United Kingdom", 
                                 symmetric = TRUE, 
                                 age.limits = c(0, 5, 18, 50, 65,75)
  )
  
  cmat <- contact_data$matrix
  proportions <- contact_data$demography$proportion
  
  
  N <-100000*proportions
  S0 <- N
  S0[1] <- N[1] - 1
  avg_c <- sum(colSums(cmat)*proportions)
  
  times <- seq(from=0,to=1,by=1/365)
  gamma <- 365/7
  epsilon <- 365/5
  
  
  age_order <- c("[0,5)", "[5,18)", "[18,50)", "[50,65)", "[65,75)", "75+")
  
  age_mortality <- c(0.0001, 0.0003, 0.0025, 0.025, 0.1, 0.22)
  
  num_ag <- length(age_mortality)
  
  xstart <- c(S0, rep(0, num_ag), c(1,rep(0,num_ag-1)), rep(0,num_ag))
  
  mort_data <- data.frame(age_group = age_order, mortality = age_mortality)
  
  
  SEIR0 <- function(t,x,parameters){
    xmat <- matrix(x, nrow = 6, ncol = 4)
    S <- xmat[,1]
    E <- xmat[,2]
    I <- xmat[,3]
    R <- xmat[,4]
    
    N <- rowSums(xmat)
    group_exposure <- as.vector(parameters$beta %*% I)
    
    foi <- S*group_exposure/N
    dS <- -foi  
    dE <- foi - (E*parameters$epsilon)
    dI <- E*parameters$epsilon - (I*parameters$gamma)
    dR <- parameters$gamma*I
    
    res <- list(c(dS,dE,dI,dR))
    
    return(res)
    
  }
  
  modelResults <- reactive({
    beta <- input$R0*gamma/avg_c
    parameters <- list(beta = beta*cmat,gamma=gamma, epsilon = epsilon)
    
    out <- ode(xstart,times,SEIR0,parameters)
    ag_labels <- contact_data$demography$age.group
    
    low_index <- (ncol(cmat)*2)+2
    high_index <- (ncol(cmat)*3)+1
    df <- data.frame(out[,low_index:high_index])
    colnames(df) <- ag_labels
    df$t <- times*365
    d1 <- tidyr::pivot_longer(df, -t, names_to = "age_group", values_to = "I")
    
    low_index <- (ncol(cmat)*3)+2
    high_index <- (ncol(cmat)*4)+1
    df <- data.frame(out[,low_index:high_index])
    colnames(df) <- ag_labels
    df$t <- times*365
    d2 <- tidyr::pivot_longer(df, -t, names_to = "age_group", values_to = "R")
    
    
    dm <- data.frame(age_group = ag_labels, mortality = age_mortality, N = N)
    
    df <- inner_join(d1, d2) %>% inner_join(dm)
    
    
    df$age_group <- fct_relevel(df$age_group, c("[0,5)", "[5,18)", "[18,50)", "[50,65)", "[65,75)", "75+"))
    df <- df %>% group_by(age_group)
    
    return(df)
    
  })
  
  
  output$prevPlot <- renderPlot({
    g <-ggplot(modelResults(), aes(x = t, y = I)) + geom_line(aes(colour=age_group)) + xlab("days")
    
    return(g)
  })
  
  output$deathPlot <- renderPlot({
    g <- ggplot(modelResults(), aes(x = t, y = 100*mortality*R/N)) + geom_line(aes(colour=age_group)) + xlab("days") + ylab("Mortality Rate (%)")
    
    return(g)
  })
  
  cmat_to_df <- function(x) {
    
    df <- x %>%
      as_tibble() %>%
      mutate(to_group = colnames(cmat)) %>% 
      gather(key = "age_group", value="r", -to_group)
    
    df$age_group <- fct_relevel(df$age_group, age_order)
    df$to_group <- fct_relevel(df$to_group, age_order)  
    
    df
    
    return(df)
    
  }
  
  ageData <- reactive({
    df <- cmat_to_df(input$R0*cmat/avg_c) %>%
      inner_join(mort_data) %>%
      mutate(deaths = r * mortality)
    
    return(df)
  })
  
  output$agePlot <- renderPlot({
    g <- ggplot(ageData(), aes(x=to_group, y=age_group, fill = r)) + 
      geom_tile() + 
      xlab("From age group") +
      ylab("To age group") +
      coord_equal()
    
    return(g)
  })
  
  output$ageMortPlot <- renderPlot({
    g <- ggplot(ageData(), aes(x=to_group, y=age_group, fill = deaths)) + 
      geom_tile() + 
      xlab("From age group") +
      ylab("To age group") +
      coord_equal()
    
    return(g)
  })
}

shinyApp(ui, server)