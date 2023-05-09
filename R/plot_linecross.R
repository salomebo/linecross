#' 3D plot of line-cross GP-map
#'
#'\code{plot_linecross} draws a 3D plot of linecross GPmap. Black dots represent the observed mean for
#'the populations and derivatives and grey dots are the predicted values.

#' @param output.model fitted model object of the format linecross
#' @param data data.frame with three columns (derivatives, means, standard errors). See data_tribolium for an example.
#' @param theta same as function \code{\link{persp}}
#' @param phi same as function \code{\link{persp}}
#' @param margins same as function \code{mar} in \code{\link{par}}
#' @param col.triangle colors for the lines of the triangle formed by linecross analysis. See details.
#' @author Geir H. Bolstad & Salom`e Bourg
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
#' mod.epistasis=linecross_models(data=data_tribolium, model="general",reference="F1", maxeval=3)
#'
#' ## Plot the corresponding linecross GPmap
#' plot_linecross(output.model=mod.epistasis, data=data_tribolium, theta=-60, phi=20,
#' margins=c(2,2,2,2), col.triangle="gray80")
#'
#' @seealso for X-Y-Z plotting see ‘persp’
#'
#' @importFrom graphics persp par segments lines points mtext
#' @importFrom grDevices trans3d
#'
#' @export
#'
plot_linecross <- function(output.model,data, theta, phi, margins,col.triangle){
  par(mar=margins)
  model=output.model$model.information[[1]]
  reference=output.model$model.information[[2]]
  S_complete= c(1, 0, 1/2,  1/2, 1/2, 1/2, 1/2, 1/2, 3/4, 3/4, 3/4, 3/4, 1/4, 1/4, 1/4, 1/4)
  H_complete= c(0, 0, 1, 1, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2)
  data=cbind(data,S_complete, H_complete)
  sub_data=data[!is.na(data$means),]
  S = sub_data$S_complete
  names(S)=sub_data$derivatives
  H = sub_data$H_complete
  names(H)=sub_data$derivatives
  AA <- seq(0,1,0.1)
  ss <- c(seq(0,1,0.02), seq(1,0,-0.01))
  hh <- c(rep(0, 50), seq(0,1,0.02), seq(1,0,-0.02))
  l1 <- seq(0.1, 0.9, 0.1)
  l2 <- 1-abs(1-2*l1)
  l3 <- seq(0.05,0.45, 0.05)

  linears=c("additive", "dominance", "add_dom", "general", "general_dom", "generalWB", "generalWB_dom", "generalW", "generalW_dom", "generalB", "generalB_dom", "classic")
  non_linears= c("multilinear", "canalization","multilinear_add")

  if(model %in% linears == TRUE) {
    read.model <-function(S,H,par) eval(parse(text=Mlist[[reference]][[model]][[3]]))
    no.parameters=length(Mlist[[reference]][[model]][[2]])+1
  }
  if(model %in% non_linears == TRUE){
    read.model <-function(S,H,par) eval(parse(text=Mlist[[reference]][[model]][[1]]))
    no.parameters=length(Mlist[[reference]][[model]][[2]])+1
  }

  FUNK <- function(X,Y) as.numeric(output.model$parameters[1]) +read.model(X,Y, as.numeric(output.model$parameters[2:no.parameters]))

  original_mat= outer(AA,AA,FUNK)
  reverse_mat =original_mat[nrow(original_mat):1,]
  minimum=min(sub_data$means,reverse_mat)-0.01*min(sub_data$means,reverse_mat)
  maximum=max(sub_data$means,reverse_mat)+0.01*max(sub_data$means,reverse_mat)
  persp(x=AA, y=AA, reverse_mat, theta =theta, phi = phi, xlab="1-S" ,ylab="H",zlab="z", cex.lab=1, border="gray88",zlim=c(minimum, maximum)) -> res

  data.model_2D = trans3d(1-S,H, FUNK(S,H), pmat=res)
  data_2D =trans3d(1-S,H, sub_data$means, pmat=res)
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
  lines(trans3d(1-ss,hh, FUNK(ss,hh), pmat=res), lwd=2, col=col.triangle)
  segments(x0=data.model_2D$x, x1=data_2D$x, y0=data.model_2D$y, y1=data_2D$y, col="black",lty=2,lwd=0.8)
  points(data.model_2D, pch=19, cex=1.2, col="gray68")
  points(data_2D, pch=19, cex=1.2, col="black")
  mtext(model, font=2, cex=1.2, side=3, adj=0.5, padj=0)
  sub.title=expression(R^2 ~ "=" ~ as.numeric(round(output.model$Rsquare,3)))
  mtext(bquote(R^2 == .(as.numeric(round(output.model$Rsquare,3))) ~ "AIC" == .(as.numeric(round(output.model$AIC,3)))), font=2, cex=1.2, side=1, adj=0.75, padj=0)
}
