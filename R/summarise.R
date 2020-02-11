#' Summary statistics for grouped data
#' 
#' Easily get summary statistics for each group present in the data
#' 
#' @aliases summarize
#' @param x a matrix of multivariate observations, a list of summary statistics from
#' multivariate observations, a data.frame of multivariate observations, or a formula
#' with a multivariate response on the left hand side, and a grouping variable/factor
#' on the right hand side.
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
#' @export
summarise = function(x, ...){
  UseMethod("summarise")
}

#' @describeIn summarise Summary statistics for grouped data
#' @export
summarise.default = function(x, y, stats = c(mean, median, sd, nrow), ...){
  if((is.matrix(x) && is.matrix(y)) || (is.data.frame(x) && is.data.frame(y))){
    
    colStats = function(datam){
      results = vector(length = length(stats), mode = "list")
      names(results) = as.character(stats)
      
      for(statFun in stats){
        results[[as.character(statFun)]] = apply(datam, 2, statFun)
      }
      
      results = do.call(rbind, )
      
      return(results)
    }
    
    print(colStats(x))
    print(colStats(y))
    
  }else if(is.list(x) && is.list(y)){
    
  }else{
    stop("x and y must either both matrices, data.frames, or lists")
  }
}

summarise.data.frame = function(x, y, ...){
  summarise.default(x, y, ...)
}

summarize = function(x, ...){
  summarise(x, ...)
}