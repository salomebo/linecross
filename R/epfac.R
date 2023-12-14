#' Calculate epistasis factor
#'
#'\code{epfac} calculates the epistasis factor of a genetic background given a model

#' @param object fitted model object of the format linecross
#' @param background the genetic bakground in which the epistasis factor should be calculated
#' @param n number of samples in the Montecarlo simulation
#' @author Geir H. Bolstad & Salomé Bourg
#' @return The epistasis factor with uncertainty
#' @details
#' The epistasis factor measures how a genetic background modifies the reference effects. The 
#' uncertainty in the epistasis factor is estimated using Monte carlo simulation assuming
#' that the estimated genetic parameters have an uncertainty that is multivariate normal.  
#'
#' @examples
#' ## Run the model function linecross_models
#' mod=lcross(model="general", ref="P1", data=ipomopsis, maxeval=3)
#'
#' ## Calculate the epistasis factor
#' epfac(object = mod, background = "P2", n = 10)
#' # NB you should use n = 1000 or more to get an accurate error estimate
#'
#' @importFrom MASS mvrnorm
#'
#' @export
#'
epfac <- function(object, background = "P2", n = 1000){
  model <- object$model.information["model"]
  ref <- object$model.information["reference"]
  if(model %in% c("multilinear", "canalization","multilinear_add") == TRUE){
    read.model <-function(S,H,par) eval(parse(text=Mlist[[ref]][[model]][[1]]))
  } else {
    read.model <-function(S,H,par) eval(parse(text=Mlist[[ref]][[model]][[3]]))
  }
  param1 <- object$parameters[-1]
  param2 <- param1
  param2[!names(param2) %in% c("Yh", "Y1", "Y2")] <- 0
  with_ep <- read.model(S[background],H[background], param1)
  no_ep   <- read.model(S[background],H[background], param2)
  epfac <- with_ep/no_ep
  param1_dist <- MASS::mvrnorm(n, param1, object$vcov[-1,-1])
  param2_dist <- param1_dist
  param2_dist[, !colnames(param2_dist) %in% c("Yh", "Y1", "Y2")] <- 0
  with_ep_dist <- apply(param1_dist, 1, function(x) read.model(S[background],H[background], x))
  no_ep_dist <- apply(param2_dist, 1, function(x) read.model(S[background],H[background], x))
  epfac_dist <- with_ep_dist/no_ep_dist
  return(c(epfac = epfac, SE = sd(epfac_dist), quantile(epfac_dist, c(0.025, 0.975))))
}
