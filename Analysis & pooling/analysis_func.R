library(MLmetrics)

# outList <- out[[1]]
# what <- "mi"

calcOutcomes <- function(outList, what, betaTrue = NULL){
  #Get values from the complete dataset analysis
  indices <- 1:nrow(outList$comp[[1]]$PI) ### KML
  
  preds_comp  <- lapply(lapply(outList[["comp"]], "[[", x = "PI"), '[', indices)
  
  y_val_rbind <- do.call("cbind", lapply(outList[["comp"]], "[[", x = "y"))
  colnames(y_val_rbind) <- paste("val_y", indices, sep = "_")
  
  #Get average true y value per prediction over 500 replications
  
  y_true_avg <- colMeans(y_val_rbind)
  
  #Get the statistics for the dataset you want to look at : "si", "mi", or "ld"
  
  preds   <- do.call("rbind",
                     lapply(
                       lapply(outList[[what]], "[[", x = "PI"),
                       '[',
                       indices
                     )
  )
  colnames(preds) <- paste("pred", indices, sep = "_")
  
  lwr_PI  <- do.call("rbind",
                     lapply(
                       lapply(outList[[what]], "[[", x = "PI"),
                       '[',
                       indices + length(indices)
                     )
  )
  lowerNames_PI <- colnames(lwr_PI) <- paste("lwr_PI", indices, sep = "_")
  
  upr_PI  <- do.call("rbind",
                     lapply(
                       lapply(outList[[what]], "[[", x = "PI"),
                       '[',
                       indices + length(indices) * 2
                     )
  )
  upperNames_PI <- colnames(upr_PI) <- paste("upr_PI", indices, sep = "_")
  
  lwr_CI  <- do.call("rbind",
                     lapply(lapply(outList[[what]], "[[", x = "CI"),
                            '[',
                            indices
                     )
  )
  lowerNames_CI <- colnames(lwr_CI) <- paste("lwr_CI", indices, sep = "_")
  
  upr_CI  <- do.call("rbind",
                     lapply(
                       lapply(outList[[what]], "[[", x = "CI"),
                       '[',
                       indices + length(indices)
                     )
  )
  upperNames_CI <- colnames(upr_CI) <- paste("upr_CI", indices, sep = "_")
  
  ##Put them all into one list
  
  complete_miss    <- cbind(preds, lwr_PI, upr_PI, lwr_CI, upr_CI, y_val_rbind)
  
  ##CI coverage
  
  tmp_ci <- apply(X      = complete_miss,
                  MARGIN = 1,
                  FUN    = function(x, y, u, l) x[u] > y & x[l] < y,
                  y      = y_true_avg,
                  u      = upperNames_CI,
                  l      = lowerNames_CI)
  cic <- rowMeans(tmp_ci)
  names(cic) <- paste("CI_cov", indices, sep = "_")
  
  ciw <-colMeans(complete_miss[ , upperNames_CI] - complete_miss[ , lowerNames_CI])
  MCSD_ciw <- apply(complete_miss[ , upperNames_CI] - complete_miss[ , lowerNames_CI], 2, sd)
  
  ##PI coverage
  
  ### KML: You don't need the 'val_y_names' subsetting
  
  pic <- colMeans(complete_miss[, upperNames_PI] > y_val_rbind & complete_miss[, lowerNames_PI] < y_val_rbind)
  names(pic) <- paste("PI_cov", indices, sep = "_")
  
  piw <-colMeans(complete_miss[ , upperNames_PI] - complete_miss[ , lowerNames_PI])
  MCSD_piw <- apply(complete_miss[ , upperNames_PI] - complete_miss[ , lowerNames_PI], 2, sd)
  
  
  #Compute bias/variance related outcome measures
  y_bar <- colMeans(preds)
  
  bias_mcsd <- apply(preds - y_val_rbind,2,sd)
  
  bias  <- y_bar - y_true_avg
  names(bias) <- paste("bias_pred", indices, sep = "_")
  mcsd_pred    <- apply(preds, 2, sd)
  
  MCSE_bias <- sqrt(mean(colSums((100*(preds-colMeans(preds))/y_true_avg)^2/(500*499))))
  MCSE_prb  <- MCSE_bias *100/mean(y_true_avg)
  
  
  prb <- 100 * bias / y_true_avg
  
  ##MCSE calcs
  MCSE_cic <- sqrt(cic * (1-cic)/500) 
  MCSE_pic <- sqrt(pic * (1-pic)/500) 
  
  ##MSE calculations
  
  mse_list <- c()
  for (i in 1:nrow(preds)){
    mse_list <- append(mse_list, MSE(y_pred = preds[i, ], y_true = y_val_rbind[i, ]))}
  
  MCSE_MSE <- sqrt(sum(((rowMeans(preds) - rowMeans(y_val_rbind))^2-mean(mse_list))^2)/(500*(500-1)))
  MCSD_MSE <- sd(mse_list)
  
  out <- list(bias = bias,
              y_true_avg = y_true_avg,
              prb = prb,
              MCSE_prb = MCSE_bias,
              MCSD_bias = bias_mcsd,
              mcsd_pred = mcsd_pred,
              cic = cic,
              ciw = ciw,
              MCSD_ciw = MCSD_ciw,
              MCSE_cic = MCSE_cic,
              pic = pic,
              piw = piw,
              MCSD_piw = MCSD_piw,
              MCSE_pic = MCSE_pic,
              MSE = mse_list,
              MCSD_MSE = MCSD_MSE,
              MCSE_MSE = MCSE_MSE)
  out
}
