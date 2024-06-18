#' linecross model fitting
#'
#'\code{lcross} fits different models derived from the multilinear model of gene
#'interaction adapted to linecross data.
#'
#' @param model a model to fit is required. See ‘Details’. Currently, '"additive"',  '"dominance"',  '"add_dom"',
#' '"general"', '"general_dom"', '"generalWB"', '"generalWB_dom"', '"generalW"', '"generalW_dom"', '"generalB"',
#'  '"generalB_dom"', '"classic"','"multilinear"', '"canalization"' and '"multilinear_add"'.
#' @param ref the reference population (can be P1, P2, F1 or F2) as a character
#' @param data data.frame with three columns (name of population, mean, standard error). See ipomopsis for an example.
#' @param maxeval for the non-linear models ('"multilinear"', '"canalization"' and '"multilinear_add"'), an optimization algorithm is used
#' to estimate the parameters requiring starting values. A grid search is performed on the values Yh, Y1 and Y2, their respective log
#' likelihood is computed for each point in the grid. If the grid-search has highest likelihood for extreme values of Yh, Y1 and/or Y2, the search
#' grid is expanded in this direction(s). The maximum number of new search grids produced is defined by 'maxeval'.
#' @author Geir H. Bolstad & Salomé Bourg
#' @return fitted model object of the format linecross_models
#' @references Bolstad, G. H., Bourg, S., Griffin, D., Pélabon, C. & Hansen, T. F. (2023) Quantifying
#' genome-wide functional epistasis by line-cross analysis. in prep.
#'
#' Hansen, T. F., & Wagner, G. P. (2001)
#' Modeling genetic architecture: a multilinear theory of gene interaction.
#' Theoretical population biology. 59:61-86.
#'
#' Lynch, M., & Walsh, B. (1998) Genetics and analysis of quantitative traits
#' (Vol. 1, pp. 535-557). Sunderland, MA: Sinauer.
#'
#' @details
#' The data need to be in the required format for the function to run. See the tribolium vignette.
#'
#' Model '"multilinear"' is useful for studying directional epistasis in the divergence of one population
#' from another. For this question it makes most sense to have the ancestral population as P1 or P2 and use
#' this as the reference. This model is non-linear.
#'
#' Model '"multilinear_add"' is a modified version of the multilinear model without dominance effect. This
#' model is non-linear.
#'
#' Model '"canalization"' is a modified version of the multilinear model can be used to study whether the
#' parental populations are canalized (i.e. have a flatter GP-map) compared to the hybrid populations. This
#' model makes most sense with F2 as the reference background in which to compare the parental populations.
#' This model is non-linear.
#'
#' Models '"additive"', "dominance"', and the combined version'"add_dom" are simple models that respectively
#' only consider additive effect, dominance effect or both.
#'
#' Model '"classic"' is the classic line cross model parameterized with additive (alpha) and dominance (delta)
#' effects rather than effects of substitutions (Lynch and Walsh 1998), where Eaa, Ead, and Edd are AxA, AxD and
#' DxD epsitasis components respectively.
#'
#'Model '"general"' is the most general model in the family. When P1 is the reference, Ehh describes epistatic
#'interactions between heterozygote genotypes in the background of the P1 , Eh2 describes epistatic interactions
#'between P2 homozygotes and heterozygotes , and E22 describes epistatic interactions between P2 homozygotes.
#'This model is most useful as general model to study deviances from multilinearity.
#'
#'Model '"general_dom"' is a modified version of the general epistasis model without additive effect.
#'
#'Model '"generalWB"' is a modified version of the general epistasis model that specifically tests the
#'difference between within and between population epistasis in the background of the reference hybrid population.
#'This model makes most sense with an hybrid population (F1 or F2) as the reference background.
#'
#'Model '"generalWB_dom"' is a modified version of the '"generalWB"' model without additive effect.
#'
#'Model '"generalW"'is a modified version of the '"generalWB"' model that only estimates the within population
#'epistasis in the background of the reference population.
#'
#'Model '"generalW_dom"' is a modified version of the '"generalW"' model without additive effect.
#'
#'Model '"generalB"' is a modified version of the '"generalWB"' model that only estimates the between population
#'epistasis in the background of the reference population.
#'
#'Model '"generalB_dom"' is a modified version of the '"generalB"' model without additive effect.
#'
#'
#' @examples
#' lcross(model="general", ref="F1", data=ipomopsis, maxeval=3)
#'
#' @importFrom HelpersMG SEfromHessian
#' @importFrom stats optim weighted.mean
#'
#'
#' @export
#'


lcross = function(model, ref, data, maxeval = 3){
  d <- subset(data, !is.na(mean) & !is.na(se))
  data_order <- match(d$pop, names(S))
  d$S = S[data_order]
  d$H = H[data_order]
  error=d$se^2 # assuming measures are independent
  V<-diag(error)
  inv_V=solve(V)

  # AIC function
  logLikelihood<-function(resid, se){
    (-1/2)*sum(resid^2/se^2 + log(se^2) + log(2*pi))
    }
  AICfunc<-function(logLik, k){-2*logLik+2*k}

  # Linear models
  linear = c("additive", "dominance", "add_dom", "general", "general_dom", "generalWB",
             "generalWB_dom", "generalW", "generalW_dom", "generalB", "generalB_dom", "classic")
  if(model %in% linear) {
    M=Mlist[[ref]][[model]][[1]]
    M=as.matrix(cbind(1, M[data_order,])) # the model matrix
    C<- solve(t(M) %*% inv_V %*% M) 	    # the sampling covariance matrix
    parameters<- C %*% t(M) %*% inv_V%*% d$mean
    parameters.se<- sqrt(diag(C))
    d$est=as.vector(M%*%parameters)
    d$resid=as.vector(d$mean - d$est)
    AIC.est = AICfunc(logLikelihood(d$resid, d$se), length(parameters))
    vcov <- C
  }

  # Non-linear models
  non_linear = c("multilinear", "canalization", "multilinear_add")
  if(model %in% non_linear){
    func.grid.ref1=function(S,H,par) eval(parse(text=Mlist[[ref]][[model]][[3]]))
    func.grid.ref2=function(S,H,par) eval(parse(text=Mlist[[ref]][[model]][[4]]))
    func.grid.complete=function(S,H,par) eval(parse(text=Mlist[[ref]][[model]][[1]]))
    func.grid.complete2=function(S,H,par) eval(parse(text=Mlist[[ref]][[model]][[5]]))
    func.pred.non.linear <- function(par){
      pred <- func.grid.complete2(d$S,d$H,par)
      # t(d$mean-pred) %*% inv_V %*% (d$mean-pred)
      -logLikelihood(d$mean-pred, d$se)
    }
    d_sub = subset(d, pop%in%Mlist[[ref]]$ref[[1]])
    ref.mean = weighted.mean(d_sub$mean, 1/d_sub$se^2)
    grid.Y1 = seq((d$mean[1]-ref.mean)*(-10),(d$mean[1]-ref.mean)*10, length.out=100)
    grid.Y2 = seq((d$mean[2]-ref.mean)*(-10),(d$mean[2]-ref.mean)*10, length.out=100)
    d_sub = subset(d, pop%in%c("F1", "F1_12","F1_21"))
    wm.Yh = weighted.mean(d_sub$mean, 1/d_sub$se^2)
    if(is.na(wm.Yh)){
      d_sub = subset(d, pop%in%c("F2", "F2_12","F2_21","F2_11","F2_22"))
      wm.Yh = weighted.mean(d_sub$mean, 1/d_sub$se^2)
    }
    if(wm.Yh==ref.mean) wm.Yh <- ref.mean + ref.mean/100 
    grid.Yh = seq((wm.Yh-ref.mean)*(-10),(wm.Yh-ref.mean)*10, length.out=100)
    grid = data.frame(cbind(grid.Yh,grid.Y1, grid.Y2))
    names(grid)=c("Yh","Y1","Y2") ## it's actually only Y for pure_epi when ref=F1/F2...
    param.names = Mlist[[ref]][[model]][[2]]
    grid.original = as.matrix(grid[,names(grid)%in%param.names])
    grid = grid.original
    loglik.grid = matrix(NA, ncol=dim(grid)[2]*dim(grid)[1]-(98+dim(grid)[2]), nrow=length(grid[,1]))
    epsilon.grid = matrix(NA, ncol=dim(grid)[2]*dim(grid)[1]-(98+dim(grid)[2]), nrow=length(grid[,1]))
    ref.grid = matrix(NA, ncol=dim(grid)[2]*dim(grid)[1]-(98+dim(grid)[2]), nrow=length(grid[,1]))
    edges_grid = c(1,100)

    if(ncol(grid)==1){
      for(n_it in 1:maxeval){
        grid[abs(grid)<0.001] = 0 # to avoid very little variation in the model matrix
        for(i in 1:dim(grid)[1]){
          Y=d$mean-func.grid.ref1(d$S,d$H, grid[i,1])
          if(grid[i,1]==0){
            M=rep(1,nrow(d))
            C<- solve(t(M) %*% inv_V %*% M) 	#sampling covariance matrix
            parameters<- rbind(C %*% t(M) %*% inv_V%*% Y,0)
          } else{
            M=as.matrix(cbind(1,func.grid.ref2(d$S,d$H,grid[i,1])))
            C<- solve(t(M) %*% inv_V %*% M) 	#sampling covariance matrix
            parameters<- C %*% t(M) %*% inv_V%*% Y
          }
          estimates=func.grid.complete(d$S,d$H,c(grid[i,1],parameters[2,1])) + parameters[1,1]
          residuals=d$mean-estimates
          loglik.grid[i,1]= logLikelihood(residuals, d$se)
          epsilon.grid[i,1]= parameters[2,1]
          ref.grid[i,1]= parameters[1,1]
        }
        i.j= which(loglik.grid==max(loglik.grid),arr.ind=TRUE)
        if(i.j[1] %in% edges_grid){
          grid[,1] = grid.original[,1]+grid.original[i.j[1],1]
        }
        else break
      }
      parameters.optim = optim(par=c(grid[i.j[1],1],epsilon.grid[i.j[1],1], ref.grid[i.j[1],1]), fn=func.pred.non.linear, hessian=TRUE,method = "BFGS")
      estimates.optim = func.grid.complete2(d$S,d$H,c(parameters.optim$par[1],parameters.optim$par[2],parameters.optim$par[3]))
    }

    if(ncol(grid)==2){
      for(n_it in 1:maxeval){
        grid[abs(grid)<0.001]=0 # to avoid very little variation in the model matrix
        for(i in 1:dim(grid)[1]){
          for(j in 1:dim(grid)[1]){
            Y = d$mean - func.grid.ref1(d$S,d$H,c(grid[i,1],grid[j,2]))
            if(grid[i,1]==0 & grid[j,2] == 0){
              M=rep(1, nrow(d))
              C<- solve(t(M) %*% inv_V %*% M) 	#sampling covariance matrix
              parameters<- rbind(C %*% t(M) %*% inv_V%*% Y,0)
            }
            else {
              M=as.matrix(cbind(1,func.grid.ref2(d$S,d$H,c(grid[i,1],grid[j,2]))))
              C<- solve(t(M) %*% inv_V %*% M)
              parameters<- C %*% t(M) %*% inv_V%*% Y
          }

          estimates=func.grid.complete(d$S, d$H, c(grid[i,1], grid[j,2], parameters[2,1])) + parameters[1,1]
          residuals=d$mean-estimates
          loglik.grid[i,j]= logLikelihood(residuals, d$se)
          epsilon.grid[i,j]= parameters[2,1]
          ref.grid[i,j]= parameters[1,1]
        }
      }
      i.j= which(loglik.grid==max(loglik.grid),arr.ind=TRUE)
      if(i.j[1] %in% edges_grid | i.j[2] %in% edges_grid){
        grid[,1] = grid.original[,1]+grid.original[i.j[1],1]
        grid[,2] = grid.original[,2]+grid.original[i.j[2],2]

      }
      else break
      #print(n_it)
    }
    parameters.optim<-optim(par=c(grid[i.j[1],1], grid[i.j[2],2], epsilon.grid[i.j[1],i.j[2]], ref.grid[i.j[1],i.j[2]]), fn=func.pred.non.linear, hessian=TRUE,method = "BFGS")
    estimates.optim=func.grid.complete2(d$S,d$H,c(parameters.optim$par[1],parameters.optim$par[2],parameters.optim$par[3],parameters.optim$par[4]))
  }
  no.param=length(parameters.optim$par)
  Order <- c(no.param,1:(no.param-1))
  Hessian_list <- HelpersMG::SEfromHessian(parameters.optim$hessian, hessian=TRUE, silent=FALSE)
  vcov <- solve(Hessian_list$hessian)[Order, Order]
  parameters.se=Hessian_list$SE[Order]
  parameters= parameters.optim$par[Order]
  #loglik.optim= logLikelihood(residuals, SE)
  d$est = estimates.optim
  d$resid = d$mean - d$est
  AIC.est = AICfunc(logLikelihood(d$resid, d$se), length(parameters))

}

# All models
G_SSe = t(d$resid) %*% inv_V %*% d$resid
global.mean = weighted.mean(d$mean, 1/d$se^2)
G_SSt = t(d$mean-global.mean) %*% inv_V %*% (d$mean - global.mean)
Rsquare = 1 - G_SSe/G_SSt

L <- t(chol(V))
d$transformed.est <- forwardsolve(L, d$est)
d$transformed.resid <- forwardsolve(L, d$mean) - d$transformed.est

if (model =="additive"){
  dom = NA ; dom.se = NA
} else {
    dom = eval(parse(text=Mlist[[ref]]$dominance.estimate[[1]]))
    dom.se = eval(parse(text=Mlist[[ref]]$dominance.estimate[[2]]))
}

output=list(c(model, ref), as.vector(parameters), as.vector(parameters.se), c(dom,dom.se), as.vector(Rsquare), AIC.est, vcov, d)
names(output) = c("model.information","parameters", "parameters.se", "dominance.estimate", "Rsquare", "AIC", "vcov", "data")
names(output$model.information)=c("model", "reference")
names(output$parameters)=c(ref,Mlist[[ref]][[model]][[2]])
names(output$parameters.se)=c(ref, Mlist[[ref]][[model]][[2]])
names(output$dominance.estimate)=c("dom", "se")
class(output) <- "lcross"
return(output)
}


#' @export
print.lcross <- function(object, ...) {
  cat("Model: "); cat(object$model.information[1])
  cat("\n")
  cat(paste("\nReference:", object$model.information[2]))
  cat("\n\nParameters:\n")
  print(cbind(Estimate = object$parameters, Standard_error = object$parameters.se))
  cat("\nR-squared:\n")
  print(object$Rsquare)
  cat("\nAIC:\n")
  print(object$AIC)
}


