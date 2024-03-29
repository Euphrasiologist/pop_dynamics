#' @title Logistic growth
#' @description Simple differential equation solution to a population undergoing logistic growth. 
#' It is essentially a wrapper for the deSolve package function ode().
#' 
#' @param state initial state of the population, X
#' @param parameters vector of named values, r (pop growth rate), and K (carrying capacity of pop)
#' @param times vector of time period (or generations) for logistic equation to be solved over, best generated using the seq() function.
#' @param method see ?deSolve::ode for methods and algorithms with which to solve differential equations
#' 
#' @import deSolve
#' @importFrom data.table as.data.table
#' @export
#' @examples
#' # library(ggplot2); library(dplyr)
#' 
#' # set up the parameters r, population growth rate
#' # and K, the carrying capacity of the population
#' parameters <- c(r = 0.1,
#'                 K = 50)
#' # set the initial parameter for X (pop size)
#' state <- c(X = 1)
#' # time in generations
#' times <- seq(0, 100, by = 0.01)
#' 
#' # run the model
#' loggr(state = state, parameters = parameters, times = times) %>%
#'    ggplot(aes(x = time, y = X)) + geom_line()


loggr <- function(state = c(X = NULL), 
                  parameters = c(r = NULL,
                                 K = NULL), 
                  times, 
                  method = NULL){
  
  # define the logistic growth equation function
  logistic_growth <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      # logistic growth equation
      dX <- r * X * (1 - (X / K))
      # return the rate of change of growth
      list(dX)
    })
  }
  # solve using ode
  sol <- ode(y = state, 
             times = times, 
             func = logistic_growth, 
             parms = parameters, 
             # user can specify solving method here
             method = if(!is.null(method)){
               as.character(method)
             } else NULL
             )
  # return as a data table object
  return(as.data.table(sol))
}
