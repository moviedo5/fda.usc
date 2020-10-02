################################################################################
# Set or query graphical and prompt output parameters
# Allow the user to set and examine a variety of global or local options which affect
# the way in which fda.usc functions computes and displays its results.
# internal parameter
par.fda.usc<-list()
par.fda.usc$verbose <- FALSE
par.fda.usc$trace <- FALSE
par.fda.usc$warning <- FALSE
par.fda.usc$ncores <- 1 
par.fda.usc$int.method <- "TRAPZ"
par.fda.usc$eps <- as.double(.Machine[[1]]*10)
###############################################################################

#' ops.fda.usc  Options Settings
#' 
#' Set or query graphical and prompt output parameters. Allow the user to set and examine a variety of global or local options which affect the way in which fda.usc functions computes and displays its results.
#' @aliases ops.fda.usc
#' @param verbose \code{logical}. Should R report extra information on progress? Set to \code{TRUE} by the command-line option --verbose.
#' @param trace \code{logical}. Show internal information of procedure.
#' @param warning \code{logical}: If true, warnings are shown.
# @param lty a vector of line types, see \code{\link[graphics]{par}}.
#' @param int.method see \code{method} argument in \code{\link{int.simpson}} function.
#' @param eps epsilon parameter.
#' @param ncores integer. Number of CPU cores on the current host.
#' @param reset \code{logical}.  If \code{TRUE} creates a new Parallel Socket Cluster (ncores>1) or a sequential parallel backend (ncores=1). It is useful when worker initialization failed or after a crush.
#' @author Manuel Oviedo de la Fuente (\email{manuel.oviedo@@usc.es}).
#' @examples
#' \dontrun{
#' # If worker initialization failed, please execute this code
#'  ncores <- max(parallel::detectCores() -1,1)
#'  if (ncores==1) {
#'      foreach::registerDoSEQ()
#'  }  else{
#'  cl <-  suppressWarnings(parallel::makePSOCKcluster(ncores ))
#'  doParallel::registerDoParallel(cl)
#'  }
#'  ops.fda.usc()
#' }
#'
#' @export 
ops.fda.usc = function(verbose = FALSE,trace = FALSE,warning = FALSE,
                       ncores = NULL,
                       int.method = "TRAPZ",reset = FALSE,
                       eps = as.double(.Machine[[1]]*10)){
  
  if (reset)   foreach::registerDoSEQ()
  if (is.null(ncores)) ncores = max(parallel::detectCores() -1,1)
  # If worker initialization failed, please execute this code
 
   #print("entra ops.fda.usc")
  .par.fda.usc = list()
  .par.fda.usc$verbose = verbose
  .par.fda.usc$trace = trace
  .par.fda.usc$warning = warning 
  .par.fda.usc$ncores = ncores
  .par.fda.usc$int.method = int.method
  .par.fda.usc$eps =  eps
 # stp <- FALSE
  #par.fda.usc <- eval(parse(text="par.fda.usc"), envir=.GlobalEnv)
 #cat(" ncores:",ncores)
  #if (foreach:::getDoParRegistered())
  ## Use a dummy loop to suppress possible (non-)warning from
  ## initial call to %dopar% with a sequential backend...
  
  if (ncores==1) {
     foreach::registerDoSEQ()
  }  else{
    if (foreach::getDoParWorkers()!=ncores){
      # cat("getDoParWorkers != ncores")
      cl <-  suppressWarnings(parallel::makePSOCKcluster(ncores ))
      doParallel::registerDoParallel(cl)
    }
  }
  
# foo <- suppressWarnings(foreach::"%dopar%"(foreach::foreach(i=1), {}))
#  e <- new.env()
  #unlockEnvironment <- function (env) {
   # return (new.env(parent=env))
  #}
  #e <- unlockEnvironment(e)
  e<-environment(ops.fda.usc)
  #unlockBinding("par.fda.usc", e)
  par.unlock<-list("sym"="par.fda.usc","env"=e)
  do.call("unlockBinding",par.unlock)
  assign("par.fda.usc", .par.fda.usc, envir = e)
  get("par.fda.usc", envir = e)
  #globalVariables(names=c("i'), package="fda.usc", add=F)
  return(.par.fda.usc)
}

# max(parallel::detectCores() -1,1)
# stopCluster(cl)
# ops.fda.usc(ncores=1)