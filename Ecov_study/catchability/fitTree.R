#' @title Regression tree analysis
#' @description Run regression tree analysis for catchability simulation
#' 
#' @param data A data frame containing columns of performance metrics for which regression trees should be fit and explanatory variables. No default.
#' @param perfMet A vector of performance metrics for which regression trees should be fit, names must match column in data. No default.
#' @param cp A number for the rpart.control complexity parameter, splits that don't decrease lack of fit by this margin are not attempted. Default = 0.001
#' @param outdir A file path to the directory where results should be stored, default = here::here().
#' 
#' @return A time-stamped Rdata file titled "treeResults" in outdir and a list object containing:
#' \itemize{
#'   \item{initialTree - A list of initial fitted trees for each perfMet}
#'   \item{cp - A list of complexity parameter profiles for initial tree}
#'   \item{optimalSplits - A list of optimal splits determined by the pruning protocol for each tree}
#'   \item{optimalCP - A list of optimal complexity parameter determined by the pruning protocol for each tree}
#'   \item{optimalTree - A list of optimal fitted trees, after pruning}
#'   \item{importance - A list of variable importance for each optimal tree}
#'   \item{optimalVals - A list of optimal values used in each optimal tree}
#' }
#' NOTE: Regression trees fit using rpart
#' NOTE: Categorical and oolean performance metrics and response variables are treated as factors
#' NOTE: Code assumes all columns that are not performance metrics contain explanatory variables

data = cbind(perf1 = c(1,2,3,4), perf2 = c("TRUE", "FALSE", "TRUE", "TRUE"), var1 = c("this", "that", "that", "this"), var2 = c(1,2,3,4))

fitTree <- function(data = NULL,
                    perfMet = NULL,
                    cp = 0.001,
                    outdir = here::here()){
  # Storage objects
  store_initialTree <- list()
  store_cp <- list()
  store_optimalSplits <- list()
  store_optimalCP <- list()
  store_optimalTree <- list()
  store_importance <- list()
  store_optimalVals <- list()
  
  # ID perfMet columns
  perfMetCol <- which(colnames(data) %in% perfMet)
  
  # Check if data is categorical (string) or boolean and update name so treated as factor in formula
  colnames(data)[which(is.character(data[1,]))] <- paste0(paste0("as.factor(", colnames(data)[which(is.character(data[1,]))]), ")")
  
  # Update perfMet names so categorical and boolean metrics treated as factor in formula
  perfMet <- colnames(data)[perfMetCol]
  
  # Determine method to use for tree fitting
  treeMethod = rep(NA, length(perfMet))
  treeMethod[grep("as.factor(", perfMet)] <- "class" # Use class method when perfMet treated as a factor
  treeMethod[is.na(perfMet)] <- "anova" # If continuous perfMet, use anova method
  
  perfMet_factor <- grep("as.factor(", perfMet) # Need list indicating if perfMet is treated as a factor (used to set tree method)
  
  # Pull out explanatory variables (i.e. columns that are not performance metrics)
  explanatory <- data[,which(colnames(data) %in% perfMet)] %>% colnames()
  
  # Fit regression tree
  for(iperf in perfMet){
    ##### Set up formula & set seed so same tree plot each time
    formula <- as.formula(paste(perfMet[iperf], "~", paste(explanatory, sep=" + "), sep = " "))
    set.seed(12345)
    
    ##### Fit initial tree
    initialTree <- rpart(formula, 
                         data = data, 
                         method = treeMethod[iperf],
                         control = rpart.control(cp = cp))
    
    store_initialTree[[iperf]] <- initialTree # Store initial tree
    
    ##### Prune tree
    store_cp[[iperf]] <- initialTree$cptable
    # pick optimal complexity
    MinError <-min(rowSums(store_cp[[iperf]][,c("xerror","xstd")])) # minimum error represents the line drawn when plotcp() is run
    PickCP <- min(which(store_cp[[iperf]][,"xerror"] < MinError)) # gives position of xerrors that meet are less than MinError, pick the minimum position (the leftmost cp on graph that is below dotted line(MinError))
    store_optimalSplits[[iperf]] <- store_cp[[iperf]][PickCP,"nsplit"]
    store_optimalCP[[iperf]] <- store_cp[[iperf]][PickCP,"CP"]
    
    ##### Produce optimal tree based on pruning
    optimalTree <- rpart(formula, 
                         data = data, 
                         method = treeMethod[iperf],
                         control = rpart.control(cp = store_cp[[iperf]]))
    
    store_optimalTree[[iperf]] <- optimalTree
    store_importance[[iperf]] <- optimalTree$variable.importance
    frameVars <- optimalTree$frame[,"var"]
    leaves <- frameVars=="<leaf>"
    store_optimalVars[[iperf]] <- unique(frameVars[!leaves]) # Store variable used in splits (i.e. those that aren't leaves)
    
    # Plot optimal tree
    png(paste0(outdir, "/optimal_", perfMet[iperf], ".png"))
    par(xpd = NA, mar = c(2.5, 5, 2.5, 5))
    plot(optimalTree, main=paste("Optimal", perfMet[iperf], sep=""))
    text(optimalTree, cex = 1)
    par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))
    dev.off()
    
  } # End loop over performance metrics
  
  # Save results
  store <- list(initialTree = store_initialTree,
                cp = store_cp,
                optimalSplits = store_optimalSplits,
                optimalCP = store_optimalCP,
                optimalTree = store_optimalTree,
                importance = store_importance,
                optimalVals = store_optimalVals)
  timeStamp <- Sys.time() %>% gsub(" ", "_",.) %>% gsub(":", "-", .) # Change millisecond half of Sys.time() output to avoid having spaces/weird characters in filenames
  
  write_rds(store, file = paste(outdir, paste0("treeResults", timeStamp, ".RDS"), sep="/"))
  
  # Return
  return(store)
}














