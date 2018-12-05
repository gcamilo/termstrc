#############################################################################
### Nelson/Siegel-type yield curve estimation method for 'dyncouponbonds' ###
#############################################################################

estim_nss.dyncouponbonds <- function(dataset, group, matrange="all",method="ns",
                                     lambda=0.0609*12,          # yearly lambda-value for "Diebold/Li" estimation
                                     tauconstr = NULL,              # tau constraints
                                     optimtype = "firstglobal", # 'firstglobal' of 'allglobal'
                                     constrOptimOptions = list(control = list(maxit = 2000), 
                                                               outer.iterations = 200, outer.eps = 1e-04),
                                     parallel = FALSE, ...
){
  
  res <- list()
  
  # If not in parallel or if firstglobal
  if(parallel == FALSE || optimtype == "firstglobal"){
    
    for (i in seq(length(dataset))) {
      if(i>1 && optimtype == "firstglobal"){
        ## use optimal parameters from previous period as start parameters
        b <- t(mapply(function(j) res[[i-1]]$opt_result[[j]]$par,  seq_along(group)))
        rownames(b) <- group                               
      } else b <- NULL
      
      ## perform sequence of term structure estimations
      bonddata <- dataset[[i]]
      res[[i]] <- estim_nss(bonddata, group, matrange, method=method, startparam=b, lambda=lambda,tauconstr,constrOptimOptions)
      res[[i]]$ddate <- bonddata$RD$TODAY
    }
  } else {
    
    cl <- parallel::makeCluster(detectCores())
    doParallel::registerDoParallel(cl)
    
    ## perform sequence of term structure estimations
    res <- foreach(i = seq(length(dataset)),.packages = c("termstrc","RQuantLib")) %dopar% {
      
      ## static estimation
      bonddata <- dataset[[i]]
      res_here <- estim_nss(bonddata, group, matrange, method=method, startparam=NULL, lambda=lambda,tauconstr,constrOptimOptions)
      res_here$ddate <- bonddata$RD$TODAY
      return(res_here)
    }
    parallel::stopCluster(cl)
  }
  
  
  class(res) <- "dyntermstrc_nss"
  
  res
}