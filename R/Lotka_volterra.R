#' @title Lotka-Volterra model
#' @description The classic Lotka-Volterra model of predator-victim dynamics
#' 
#' @param state initial state of the populations, the victim, V and the predator, P.
#' @param parameters B, births of victims, D are deaths of predators. A and C are interactions between victims and predators.
#' @param times vector of time period (or generations) to be solved over, best generated using the seq() function.
#' @param method see ?deSolve::ode for methods and algorithms with which to solve differential equations
#' 
#' @import deSolve
#' @importFrom data.table as.data.table
#' @export
#' @examples
#' # library(ggplot2); library(dplyr)
#' 
#' # lots of parameters
#' parameters <- c(A = 4,
#'                 B = 5,
#'                 C = 6,
#'                 D = 7)
#' # states of the initial populations
#' state = c(V = 10,
#'           P = 10)
#' # time in generations
#' times <- seq(0, 25, by = 0.01)
#' 
#' # run the model
#' lotvol(state = state, parameters = parameters, times = times) %>% 
#'   melt.data.table(measure.vars = c("V", "P")) %>% 
#'   ggplot(aes(x = time, y = value)) + geom_line(aes(colour = variable))
#'   
#' @references 
#' Lotka, A.J. 1925 The Elements of Physical Biology. Williams and Williams Co., Baltimore.
#' Volterra, V. 1926 Fluctuations in the Abundance of a Species Considered Mathematically. Nature 118 558-560

lotvol <- function(state = c(V = NULL,
                             P = NULL), 
                   parameters = c(A = NULL,
                                  B = NULL,
                                  C = NULL,
                                  D = NULL), 
                   times, 
                   method = NULL){
  
  # define the logistic growth equation function
  lotka_volterra <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      # logistic growth equations
      dV <- B*V - C*V*P
      dP <- -D*P + A*V*P
      # return the rate of change of growth
      res <- c(dV, dP)
      list(res)
    })
  }
  # solve using ode
  sol <- ode(y = state, 
             times = times, 
             func = lotka_volterra, 
             parms = parameters, 
             # user can specify solving method here
             method = if(!is.null(method)){
               as.character(method)
             } else NULL
  )
  # return as a data table object
  return(as.data.table(sol))
}
