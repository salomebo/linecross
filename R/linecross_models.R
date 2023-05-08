#' linecross model fitting
#'
#'\code{linecross_models} fits different models derived from the multilinear model of gene
#'interaction adapted to linecross data.
#'
#' @param data data.frame with three columns (derivatives, means, standard errors). See data_tribolium for an example.
#' @param model a model to fit is required. See ‘Details’. Currently, '"additive"',  '"dominance"',  '"add_dom"',
#' '"general"', '"general_dom"', '"generalWB"', '"generalWB_dom"', '"generalW"', '"generalW_dom"', '"generalB"',
#'  '"generalB_dom"', '"classic"','"multilinear"', '"canalization"' and '"multilinear_add"'.
#' @param reference the reference population (can be P1, P2, F1 or F2) as a character
#' @param maxeval for the non-linear models ('"multilinear"', '"canalization"' and '"multilinear_add"'), an optimization algorithm is used
#' to estimate the parameters requiring starting values. A grid search is performed on the values Yh, Y1 and Y2, their respective log
#' likelihood is computed for each point in the grid. If the grid-search has highest likelihood for extreme values of Yh, Y1 and/or Y2, the search
#' grid is expanded in this direction(s). The maximum number of new search grids produced is defined by 'maxeval'.
#' @author Geir H. Bolstad & Salomé Bourg
#' @return fitted model object of the format linecross_models
#' @references Hansen, T. F., & Wagner, G. P. (2001).
#' Modeling genetic architecture: a multilinear theory of gene interaction.
#' Theoretical population biology, 59(1), 61-86.
#'
#' Lynch, M., & Walsh, B. (1998).
#' Genetics and analysis of quantitative traits
#' (Vol. 1, pp. 535-557). Sunderland, MA: Sinauer.
#'
#' @details
#'Model '"multilinear"' is useful for studying directional epistasis in the divergence of one population
#'from another. For this question it makes most sense to have the ancestral population as P1 or P2 and use
#'this as the reference. This model is non-linear.
#'
#'Model '"multilinear_add"' is a modified version of the multilinear model without dominance effect. This
#'model is non-linear.
#'
#'Model '"canalization"' is a modified version of the multilinear model can be used to study whether the
#'parental populations are canalized (i.e. have a flatter GP-map) compared to the hybrid populations. This
#'model makes most sense with F2 as the reference background in which to compare the parental populations.
#'This model is non-linear.
#'
#'Models '"additive"', "dominance"', and the combined version'"add_dom" are simple models that respectively
#'only consider additive effect, dominance effect or both.
#'
#'Model '"classic"' is the classic line cross model parameterized with additive (\alpha) and dominance (\delta)
#'effects rather than effects of substitutions (Lynch and Walsh 1998), where Eaa, Ead, and Edd are AxA, AxD and
#'DxD epsitasis components respectively.
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
#' mod.epistasis=linecross_models(data=data_tribolium, model="general",reference="F1", maxeval=3)
#'
#' @importFrom HelpersMG SEfromHessian
#'
#' @export
#'
#'
#'
linecross_models=function(data, model, reference, maxeval){
  linears=c("additive", "dominance", "add_dom", "general", "general_dom", "generalWB", "generalWB_dom", "generalW", "generalW_dom", "generalB", "generalB_dom", "classic")
  non_linears= c("multilinear", "canalization","multilinear_add")
  obs=data$means[!is.na(data$means)]
  SE=data$se[!is.na(data$means)]
  S=propP1[!is.na(data$means)]
  H=hybrid.index[!is.na(data$means)]
  ref.mean=weighted.mean(data$means[data$derivatives==Mlist[[reference]]$reference[[1]]], 1/(data$se[data$derivatives==Mlist[[reference]]$reference[[1]]])^2, na.rm=TRUE)
  se.mean=mean(data$se[data$derivatives==Mlist[[reference]]$reference[[1]]], na.rm=TRUE) # standard error of the reference
  if (is.na(min(SE))) lacksSE<-TRUE else lacksSE<-FALSE #To remove AIC if SE is lacking
  SE[is.na(SE)]<-1
  z=obs
  error=SE^2 # assuming measures are independent
  w<-1/error
  V<-diag(error)
  inv_V=solve(V)

  # AIC function
  logLikelihood<-function(DerivateResiduals, DerivateStdErr){(-1/2)*sum(DerivateResiduals^2 / DerivateStdErr^2 + log(DerivateStdErr^2)+6*log(2*pi))} ## Thomas's version
  AICfunc<-function(logLikelihood, k){-2*logLikelihood+2*k}

  if(model %in% linears == TRUE) {
    M=Mlist[[reference]][[model]][[1]] ## selection of the right model matrix
    M=as.matrix(cbind(1,M[!is.na(data$means),])) # 1st column intercept (reference mean), 2nd column S index (slope)
    C<- solve(t(M) %*% inv_V %*% M) 	#sampling covariance matrix
    parameters<- C %*% t(M) %*% inv_V%*% z
    parameters.se<- sqrt(diag(C))
    estimates=as.vector(M%*%parameters)
    residuals=as.vector(t(z)- estimates)
    AIC.estimation = AICfunc(logLikelihood(residuals, SE), dim(M)[2])
  }

  if(model %in% non_linears == TRUE){
    func.grid.ref1=function(S,H,par) eval(parse(text=Mlist[[reference]][[model]][[3]]))
    func.grid.ref2=function(S,H,par) eval(parse(text=Mlist[[reference]][[model]][[4]]))
    func.grid.complete=function(S,H,par) eval(parse(text=Mlist[[reference]][[model]][[1]]))
    func.grid.complete2=function(S,H,par) eval(parse(text=Mlist[[reference]][[model]][[5]]))
    func.pred.non.linears <- function(par){
      pred <- func.grid.complete2(S,H,par)
      t(z-pred) %*% inv_V %*% (z-pred)
    }
    grid.Y1=seq((obs[1]-ref.mean)*(-10),(obs[1]-ref.mean)*10, length.out=100)
    grid.Y2=seq((obs[2]-ref.mean)*(-10),(obs[2]-ref.mean)*10, length.out=100)
    wm.Yh=weighted.mean(data$means[data$derivatives%in%c("F1_12","F1_21")], 1/(data$se[data$derivatives%in%c("F1_12","F1_21")])^2, na.rm=TRUE)
    if(is.na(wm.Yh)) wm.Yh=weighted.mean(data$means[data$derivatives%in%c("F2_12","F2_21","F2_11","F2_22")], 1/(data$se[data$derivatives%in%c("F2_12","F2_21","F2_11","F2_22")])^2, na.rm=TRUE)
    grid.Yh=seq((wm.Yh-ref.mean)*(-10),(wm.Yh-ref.mean)*10, length.out=100)
    grid=data.frame(cbind(grid.Yh,grid.Y1, grid.Y2))
    names(grid)=c("Yh","Y1","Y2") ## it's actually only Y for pure_epi when ref=F1/F2...
    param.names=Mlist[[reference]][[model]][[2]]
    grid.original=as.matrix(grid[,names(grid)%in%param.names])
    grid=grid.original
    loglik.grid=matrix(NA, ncol=dim(grid)[2]*dim(grid)[1]-(98+dim(grid)[2]), nrow=length(grid[,1]))
    epsilon.grid=matrix(NA, ncol=dim(grid)[2]*dim(grid)[1]-(98+dim(grid)[2]), nrow=length(grid[,1]))
    ref.grid=matrix(NA, ncol=dim(grid)[2]*dim(grid)[1]-(98+dim(grid)[2]), nrow=length(grid[,1]))
    edges_grid=c(1,100)

    if(ncol(grid)==1){	## ncol(grid)==1 when model is pure_epi	and ref = P1/P2
      for(n_it in 1:maxeval){
        grid[abs(grid)<0.001]=0 #to avoid very litlle variation in the model matrix
        for(i in 1:dim(grid)[1]){
          Y=z-func.grid.ref1(S,H, grid[i,1])
          if(grid[i,1]==0){
            M=rep(1,length(obs))
            C<- solve(t(M) %*% inv_V %*% M) 	#sampling covariance matrix
            parameters<- rbind(C %*% t(M) %*% inv_V%*% Y,0)
          }
          else{
            M=as.matrix(cbind(1,func.grid.ref2(S,H,grid[i,1])))
            C<- solve(t(M) %*% inv_V %*% M) 	#sampling covariance matrix
            parameters<- C %*% t(M) %*% inv_V%*% Y
          }
          estimates=func.grid.complete(S,H,c(grid[i,1],parameters[2,1])) + parameters[1,1]
          residuals=z-estimates
          loglik.grid[i,1]= logLikelihood(residuals, SE)
          epsilon.grid[i,1]= parameters[2,1]
          ref.grid[i,1]= parameters[1,1]
        }
        i.j= which(loglik.grid==max(loglik.grid),arr.ind=TRUE)
        if(i.j[1] %in% edges_grid){
          grid[,1] = grid.original[,1]+grid.original[i.j[1],1]
        }
        else break
      }
      parameters.optim<-optim(par=c(grid[i.j[1],1],epsilon.grid[i.j[1],1], ref.grid[i.j[1],1]), fn=func.pred.non.linears, hessian=TRUE,method = "BFGS")
      estimates.optim=func.grid.complete2(S,H,c(parameters.optim$par[1],parameters.optim$par[2],parameters.optim$par[3]))
    }

    if(ncol(grid)==2){	## ncol(grid)==2 when model is multilinear or canalization
      for(n_it in 1:maxeval){
        grid[abs(grid)<0.001]=0 #to avoid very litlle variation in the model matrix
        for(i in 1:dim(grid)[1]){
          for(j in 1:dim(grid)[1]){
            Y=z-func.grid.ref1(S,H,c(grid[i,1],grid[j,2]))
            if(grid[i,1]==0 & grid[j,2] == 0){
              M=rep(1,length(obs))
              C<- solve(t(M) %*% inv_V %*% M) 	#sampling covariance matrix
              parameters<- rbind(C %*% t(M) %*% inv_V%*% Y,0)
            }
            else {
              M=as.matrix(cbind(1,func.grid.ref2(S,H,c(grid[i,1],grid[j,2]))))
              C<- solve(t(M) %*% inv_V %*% M)
              parameters<- C %*% t(M) %*% inv_V%*% Y
          }

          estimates=func.grid.complete(S,H,c(grid[i,1],grid[j,2],parameters[2,1])) + parameters[1,1]
          residuals=z-estimates
          loglik.grid[i,j]= logLikelihood(residuals, SE)
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
    parameters.optim<-optim(par=c(grid[i.j[1],1], grid[i.j[2],2],epsilon.grid[i.j[1],i.j[2]], ref.grid[i.j[1],i.j[2]]), fn=func.pred.non.linears, hessian=TRUE,method = "BFGS")
    estimates.optim=func.grid.complete2(S,H,c(parameters.optim$par[1],parameters.optim$par[2],parameters.optim$par[3],parameters.optim$par[4]))
  }
  no.param=length(parameters.optim$par)
  parameters.se=HelpersMG::SEfromHessian(parameters.optim$hessian, hessian=FALSE, silent=FALSE)[c(no.param,1:(no.param-1))]
  parameters=parameters.optim$par[c(no.param,1:(no.param-1))]
  #loglik.optim= logLikelihood(residuals, SE)
  residuals=z-estimates.optim
  AIC.estimation = AICfunc(logLikelihood(residuals, SE), no.param) ## not sure for the no of parameters

}



G_SSe = t(residuals) %*% inv_V %*% residuals
M.int=matrix(rep(1,length(residuals)), ncol=1)
C<- solve(t(M.int) %*% inv_V %*% M.int) 	#sampling covariance matrix
mean.z<- c(C %*% t(M.int) %*% inv_V%*% z)
#mean.z=weighted.mean(z, error)
G_SSt = t(z-mean.z) %*% inv_V %*% (z - mean.z)
Rsquare = 1 - G_SSe/G_SSt

if (model =="additive"){ dom = NA ; dom.se = NA}

else { dom = eval(parse(text=Mlist[[reference]]$dominance.estimate[[1]])) ; dom.se = eval(parse(text=Mlist[[reference]]$dominance.estimate[[2]]))}

#output.model=list(c(model, reference), as.vector(parameters), as.vector(parameters.se), c(dom,dom.se),as.vector(Rsquare), AIC.estimation)
output.model=list(c(model, reference), as.vector(parameters), as.vector(parameters.se), c(dom,dom.se),as.vector(Rsquare), AIC.estimation)
names(output.model) = c("model.information","parameters", "parameters.se", "dominance.estimate", "Rsquare", "AIC")
names(output.model$model.information)=c("model", "reference")
names(output.model$parameters)=c(reference,Mlist[[reference]][[model]][[2]])
names(output.model$parameters.se)=c(reference, Mlist[[reference]][[model]][[2]])
names(output.model$dominance.estimate)=c("dom", "se")
return(output.model)
}
