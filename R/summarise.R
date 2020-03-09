#' Summary statistics for grouped data
#' 
#' Easily get summary statistics for each group present in the data
#' 
#' @aliases summarize
#' @param x a matrix of multivariate observations, a list of summary statistics from
#' multivariate observations, a data.frame of multivariate observations, or a formula
#' with a multivariate response on the right hand side, and a grouping variable/factor
#' on the left hand side.
#' @param y a matrix of multivariate observations, a list of summary statistics from
#' multivariate observations, OR a data.frame of multivariate observations
#' @param data a data.frame containing the variables used in a formula
#' @param stats a named list of summary statistics to compute on each variable in each 
#' group. Note 1: Quantiles are not supported yet because I can't think of a good way to
#' handle the extra arguments. Help welcome. Note 2: The names of the elements in the list 
#' are used to label the columns of the output. They probably should be unique.
#' @param \ldots other arguments such as another matrix of multivariate observations:
#' see \code{summarise.default}, or a data to be used with a formula: see
#' \code{summarise.formula}
#' @examples 
#' 
#' data(container.df)
#' split.data = split(container.df[,-1],container.df$gp)
#' x = split.data[[1]]
#' y = split.data[[2]]
#' summarise(x, y)
#' 
#' ## Using the formula interface
#' data(container.df)
#' summarise(gp~., data = container.df)
#' 
#' summarise(gp~Al+Ti, data = container.df)
#' 
#' @importFrom stats median sd
#' @export
summarise = function(x, ...){
  UseMethod("summarise")
}

calcStats = function(x, y, stats){
  if((is.matrix(x) && is.matrix(y)) || (is.data.frame(x) && is.data.frame(y))){
    
    colStats = function(datam){
      results = lapply(stats, function(statFun)apply(datam, 2, statFun))
      
      results = do.call(rbind, results)
      
      return(results)
    }
    
    results = list(x = colStats(x),
                   y = colStats(y))
    
    return(results)
    
  }else if(is.list(x) && is.list(y)){
    
  }else{
    stop("x and y must either both matrices, data.frames, or lists")
  }
}

#' @describeIn summarise Summary statistics for grouped data
#' @export
summarise.default = function(x, 
                             y, 
                             stats = list(Mean = mean, 
                                                Median = median, 
                                                'Std. Dev.' = sd, 
                                                N = length),
                             ...){
  stats = calcStats(x, y, stats)
  cat("Group 1\n")
  cat("=======\n")
  print(stats$x)
  cat("\nGroup 2\n")
  cat("=======\n")
  print(stats$y)
  
  invisible(stats)
}

#' @describeIn summarise Summary statistics for grouped data
#' @export
summarise.formula = function(x, 
                             data = NULL, 
                             stats = list(Mean = mean, 
                                                          Median = median, 
                                                          'Std. Dev.' = sd, 
                                                          N = length),
                             ...){
  form = x
  
  if(is.null(data)){
    data = parent.frame()
  }
  
  mf = model.frame(form, data)
  
  group = factor(model.response(mf))
  lvls = levels(group)
  numLevels = length(lvls)
  groupName = colnames(mf)[1]
  
  if(numLevels < 2){
    msg = "Although this function will work for one group I am stopping because
    this package implements the two-sample version of Hotelling's T^2"
    stop(msg)
  }
  
  if(numLevels > 2){
    msg = "There are more than two groups in this data. Only the first two will be used"  
    warning(msg)
  }
  
  variables = mf[,-1, drop = FALSE]
  
  split.data = split(variables, group)
  
  stats = calcStats(split.data[[1]], split.data[[2]], stats)
  
  cat(paste0("Group 1: (", groupName, "==",  lvls[1], ")\n"))
  cat("=======\n")
  print(stats$x)
  cat(paste0("\nGroup 2: (", groupName, "==",  lvls[2], ")\n"))
  cat("=======\n")
  print(stats$y)
  
  invisible(stats)
  
}

#' @describeIn summarise Summary statistics for grouped data
#' @export
summarise.data.frame = function(x, y, ...){
  summarise.default(x, y, ...)
}

summarize = function(x, ...){
  summarise(x, ...)
}