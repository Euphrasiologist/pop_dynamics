#' @title SIS Model (Susceptible, Infectious, Susceptible)
#' @description The simplest model from a class of 'SIR' models.
#' 
#' 
#' @param state initial state of the population classes, susceptible, S, and Infectious, I.
#' @param parameters vector of named values, B - birth rate per individual, d is the natural death rate, g is the recovery rate and beta is the contact rate.
#' @param times vector of time period (or generations) for SIS equation to be solved over, best generated using the seq() function.
#' @param method see ?deSolve::ode for methods and algorithms with which to solve differential equations
#' 
#' @import deSolve
#' @importFrom data.table as.data.table
#' @export
#' @examples
#' # library(ggplot2); library(dplyr)
#' 
#' # parameters
#' parameters <- c(B = 0.21,
#'                 d = 0.2,
#'                 g = 1,
#'                 beta = 3)
#'                 
#' # states of the initial populations
#' state = c(S = 1,
#'           I = 1)
#' # time in generations
#' times <- seq(0, 100, by = 0.01)
#' 
#' # run the model
#' SIS(state = state, parameters = parameters, times = times) %>%
#' melt.data.table(measure.vars = c("S", "I")) %>% 
#' ggplot(aes(x = time, y = value)) + geom_line(aes(colour = variable))



SIS <- function(state = c(S = NULL,
                          I = NULL), 
                   parameters = c(B = NULL,
                                  d = NULL,
                                  g = NULL,
                                  beta = NULL), 
                   times, 
                   method = NULL){
  
  # define the SIS function
  susceptible_infectious_susceptible <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      # SIS equations
      N <- S + I
      dS <- B*N - ((beta*S*I)/N) + g*I - d*S
      dI <- ((beta*S*I)/N) - g*I - d*I
      
      # return the rate of change of growth
      res <- c(dS, dI)
      list(res)
    })
  }
  # solve using ode
  sol <- ode(y = state, 
             times = times, 
             func = susceptible_infectious_susceptible, 
             parms = parameters, 
             # user can specify solving method here
             method = if(!is.null(method)){
               as.character(method)
             } else NULL
  )
  # return as a data table object
  return(as.data.table(sol))
}
