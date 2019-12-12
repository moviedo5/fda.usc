globalVariables('icnt')

#' @keywords internal
 

#' @import fda 
#' @import splines
#' @import MASS
#' @import mgcv
#' @import nlme
#' @import methods
#' @import utils
#' @import grDevices
# @import graphics
#' @import stats
  
  # @import rpart deleted
  
  #import(foreach,"getDoParWorkers","getDoParRegistered","getDoParName","getDoParVersion","foreach","registerDoSEQ","setDoSeq","setDoPar")
  #importFrom(parallel, "makeCluster", "stopCluster", "detectCores", "clusterExport", "clusterEvalQ")
#' @importFrom doParallel registerDoParallel
#' @importFrom iterators  icount

#' @importFrom foreach %dopar% getDoParWorkers getDoParRegistered getDoParName getDoParVersion foreach registerDoSEQ setDoSeq setDoPar
#' @import parallel
# @importFrom grDevices adjustcolor colorRampPalette dev.cur dev.interactive dev.new devAskNewPage gray heat.colors palette
#' @importFrom graphics abline boxplot contour curve filled.contour image legend lines matlines pairs par persp plot points polygon rect rug stars text title
# @importFrom methods callGeneric
# @importFrom utils combn installed.packages modifyList setTxtProgressBar txtProgressBar globalVariables
  
  

  .onAttach <- function(lib, pkg,...){
    pkg.info <- drop(read.dcf(file=system.file("DESCRIPTION", package="fda.usc"),
                              fields=c("Title","Version","Date")))
    foo <- suppressWarnings(foreach::"%dopar%"(foreach::foreach(i=1), {}))
     packageStartupMessage(
      paste("-----------------------------------------------------------------\n",pkg.info["Title"]),"\n",
      #	 "Functional Data Analysis in R \n",
      paste(" fda.usc version ", pkg.info["Version"]," (built on ", pkg.info["Date"], ") is now loaded\n", sep=""),
      paste(" fda.usc is running sequentially usign foreach package\n"),
      paste(" Please, execute ops.fda.usc() once to run in local parallel mode\n"),
      paste(" Please, execute ops.fda.usc() once to run in local parallel mode\n"),
      paste("Deprecated functions: min.basis, min.np, anova.hetero, anova.onefactor, anova.RPm\n"),
      paste("New functions: optim.basis, optim.np, fanova.hetero, fanova.onefactor, fanova.RPm\n"),
            #" ops.fda.usc() changes the parameters of the package\n",
      "-----------------------------------------------------------------\n"
    )
  }
  
