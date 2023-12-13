#' Campbell et al. (2008) data
#' 
#' The data is from a cross between two species of *Ipomopsis*, the scarlet gilia *I. aggregata* 
#' and the slendertube skyrocket *I. tenuituba*. The trait is proportion of survived seeds 
#' surviving until reproduction planted in the environment native to *I. aggregata*. The values
#' were obtained from Figure 3A in Campbell et al. (2008).
#'
#' @format ##
#' A data frame with 9 rows and 3 columns:
#' \describe{
#'   \item{pop}{name of each experimental population}
#'   \item{mean}{the poroporiton of survived seeds}
#'   \item{se}{standard error of each mean}
#'
#' }
#' @references Campbell, D. R., Waser, N. M., Aldridge, G. & Wu, C. A. (2008) Lifetime fitness in 
#' two generations of *Ipomopsis* hybrids Evolution 62:2616-2627.
#'
#' @source <https://doi.org/10.1111/j.1558-5646.2008.00460.x>
"ipomopsis"


# d <- as.data.frame(cbind(
#   pop = c("P1", "P2", "F1_12", "F1_21",  "F2_12",  "B1_12",  "B1_11b", "B2_21", "B2_22b"),
#   mean = c(0.1012, 0.0376, 0.0964, 0.0601, 0.0832, 0.0577, 0.0421, 0.0448, 0.0262),
#   se = c(0.0253, 0.0069, 0.0215, 0.0248, 0.0188, 0.0157, 0.0138, 0.0218, 0.0088)
# ))
# 
# 
# d$mean<-as.numeric(d$mean)
# d$se<-as.numeric(d$se)
# ipomopsis <- d
# 
# save(ipomopsis, file= "/data/Egenutvikling/13099000_egenutvikling_geir_bolstad/linecross/data/ipomopsis.rda")
# 

