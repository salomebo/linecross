#' 3D plot of line-cross GP-map
#'
#'\code{plot_linecross} draws a 3D plot of linecross GP-map. Black dots represent the observed mean for
#'the populations and derivatives and grey dots are the predicted values.

#' @param object fitted model object of the format linecross
#' @param theta as in function \code{\link{persp}}
#' @param phi as in function \code{\link{persp}}
#' @param margins as in function \code{mar} in \code{\link{par}}
#' @param col.triangle colors for the lines of the triangle formed by linecross analysis. See details.
#' @param main as in function \code{\link{plot}}. The default value is the model name.
#' @author Geir H. Bolstad & Salomé Bourg
#' @return A 3D plot of line-cross GP-map
#' @details
#' The axes are genotypes (proportion of P2 alleles, denoted S), heterozygosity (the probability that
#' an individual has one P1 and one P2 allele at a locus, denoted H) and phenotypes. The two variables
#' S and H define a two-dimensional genetic-content space, with P1, F1 and P2 at the edges. Using these
#' variables in a regression with the phenotype as the response variable, this triangle becomes a surface
#' where the height gives the genotypic value (i.e. the trait value of each genotype). This surface is a
#' description of a two-dimensional genotype-phenotype map representing the composite (genome wide sum)
#' genetic effects.
#'
#' If the surface is flat and tilted in the S-direction this indicates additive effects. A flat surface
#' with tilt in the H-direction indicates dominance. Any non-linearities in the surface would indicate
#' epistasis, and different types of epistasis will generate different non-linear surfaces.
#'
#' @examples
#' ## Run the model function linecross_models
#' mod=lcross(model="general", reference="F1", data=tribolium, maxeval=3)
#'
#' ## Plot the corresponding linecross GPmap
#' plot_linecross(object=mod, data=tribolium, theta=-60, phi=20,
#' margins=c(2,2,2,2), col.triangle="gray80")
#'
#' @seealso for X-Y-Z plotting see ‘persp’
#'
#' @importFrom graphics persp par segments lines points mtext
#' @importFrom grDevices trans3d
#'
#' @export
#'
plot_linecross <- function(object, theta=-40, phi=20, margins=c(2,2,2,2), col.triangle="gray80", main = "default", 
                           cex.text = 1.2, ticktype = "detailed", zlab = "z"){
  par(mar=margins)
  model=object$model.information[[1]]
  reference=object$model.information[[2]]
  data <- object$data
  AA <- seq(0,1,0.1)
  ss <- c(seq(0,1,0.02), seq(1,0,-0.01))
  hh <- c(rep(0, 50), seq(0,1,0.02), seq(1,0,-0.02))
  l1 <- seq(0.1, 0.9, 0.1)
  l2 <- 1-abs(1-2*l1)
  l3 <- seq(0.05,0.45, 0.05)

  linear_mod=c("additive", "dominance", "add_dom", "general", "general_dom", "generalW1B", "generalW1B_dom", 
               "generalW1", "generalW1_dom", "generalW2B", "generalW2B_dom", 
               "generalW2", "generalW2_dom", "generalB", "generalB_dom", "classic")
  if(model %in% linear_mod) {
    read.model <-function(S,H,par) eval(parse(text=Mlist[[reference]][[model]][[3]]))
    no.parameters=length(Mlist[[reference]][[model]][[2]])+1
  }
  non_linear_mod= c("directional", "canalization","directional_add")
  if(model %in% non_linear_mod){
    read.model <-function(S,H,par) eval(parse(text=Mlist[[reference]][[model]][[1]]))
    no.parameters=length(Mlist[[reference]][[model]][[2]])+1
  }

  FUNK <- function(X,Y) as.numeric(object$parameters[1]) + read.model(X, Y, as.numeric(object$parameters[2:no.parameters]))

  original_mat= outer(AA,AA,FUNK)
  reverse_mat =original_mat[nrow(original_mat):1,]
  minimum=min(data$mean,reverse_mat) - 0.01*abs(max(data$mean,reverse_mat))
  maximum=max(data$mean,reverse_mat) + 0.02*abs(max(data$mean,reverse_mat))
  persp(x=AA, y=AA, reverse_mat, theta = theta, phi = phi, xlab="1-S",
        ylab="H",zlab=zlab, cex.lab=1, border="gray88", zlim=c(minimum, maximum), 
        scale = TRUE,
        ticktype = ticktype
        #xlim = c(-0.01, 1.01), ylim = c(-0.01, 1.01)
        ) -> res
  # data.model_2D = trans3d(1-data$S,data$H, FUNK(data$S, data$H), pmat=res)
  data.model_2D = trans3d(1-data$S, data$H, data$est,  pmat=res)#, continuous = FALSE)
  data_2D =       trans3d(1-data$S, data$H, data$mean, pmat=res)#, continuous = FALSE)
  for(j in 1:9){
    x<-rep(l1[j], 50)
    y<-seq(0,l2[j], length.out=50)
    lines(trans3d(1-x,y,FUNK(x,y), pmat=res), col=col.triangle)
  }
  for(k in 1:9){
    y<-rep(l1[k], 50)
    x<-seq(l3[k],1-l3[k], length.out=50)
    lines(trans3d(1-x,y,FUNK(x,y), pmat=res), col=col.triangle)
  }
  lines(trans3d(1-ss, hh, FUNK(ss, hh), pmat=res), lwd=2, col = col.triangle)
  segments(x0=data.model_2D$x, x1=data_2D$x, y0=data.model_2D$y, y1=data_2D$y, col="black",lty=2, lwd=0.8)
  points(data.model_2D, pch=19, cex=1.2, col="gray68")
  points(data_2D, pch=19, cex=1.2, col="black")
  text(trans3d(0,0,FUNK(1,0), pmat=res), "P1", col = "blue", pos = 3)
  text(trans3d(1,0,FUNK(0,0), pmat=res), "P2", col = "blue", pos = 3)
  text(trans3d(0.5,1,FUNK(0.5,1), pmat=res), "F1", col = "blue", pos = 3)
  if(main == "default") main = model
  mtext(main, font=2, cex=1.2, side=3, adj=0.5, padj=0)
  #sub.title=expression(R^2 ~ "=" ~ as.numeric(round(object$Rsquare,3)))
  mtext(bquote(R^2 == .(as.numeric(round(object$Rsquare,3))) ~ " AIC" == .(as.numeric(round(object$AIC,1)))), font=2, cex=cex.text, side=1, adj=0.75, padj=1)
  par(mar=c(5, 4, 4, 2) + 0.1)
}

