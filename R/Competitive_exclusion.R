#' @title Competitive exclusion (logistic growth) of two species
#' @description Modifying the logistic regression to two species, adding an extra parameter, alpha
#' which describes the competitive ability of the population.
#' 
#' @param state initial state of the populations, X1 and X2
#' @param parameters vector of named values, r(1,2) (pop growth rates), K(1,2) (carrying capacities of pops) and A(12, 21) competitive advantages of each species.
#' @param times vector of time period (or generations) for logistic equation to be solved over, best generated using the seq() function.
#' @param method see ?deSolve::ode for methods and algorithms with which to solve differential equations
#' 
#' @import deSolve
#' @importFrom data.table as.data.table
#' @export
#' @examples
#' # library(ggplot2); library(dplyr)
#' 
#' # lots of parameters
#' parameters <- c(r1 = 0.15,
#'                 r2 = 0.1,
#'                 K1 = 100,
#'                 K2 = 100,
#'                 A12 = 0.2,
#'                 A21 = 0.3)
#' # states of the initial populations
#' state = c(X1 = 1,
#'           X2 = 1)
#' # time in generations
#' times <- seq(0, 100, by = 0.01)
#' 
#' # run the model
#' compex(state = state, parameters = parameters, times = times) %>%
#' melt.data.table(measure.vars = c("X1", "X2")) %>% 
#' ggplot(aes(x = time, y = value)) + geom_line(aes(colour = variable))



compex <- function(state = c(X1 = NULL,
                             X2 = NULL), 
                   parameters = c(r1 = NULL,
                                  r2 = NULL,
                                  K1 = NULL,
                                  K2 = NULL,
                                  A12 = NULL,
                                  A21 = NULL), 
                   times, 
                   method = NULL){
  
  # define the logistic growth equation function
  competitive_exclusion <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      # logistic growth equations
      dX1 <- r1 * X1 * (1 - (X1 + A12*X2)/ K1)
      dX2 <- r2 * X2 * (1 - (X2 + A21*X1)/ K2)
      # return the rate of change of growth
      res <- c(dX1, dX2)
      list(res)
    })
  }
  # solve using ode
  sol <- ode(y = state, 
             times = times, 
             func = competitive_exclusion, 
             parms = parameters, 
             # user can specify solving method here
             method = if(!is.null(method)){
               as.character(method)
             } else NULL
  )
  # return as a data table object
  return(as.data.table(sol))
}
