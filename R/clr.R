#' Centered log ratio transformation
#' 
#' Aitchison's centered log ratio tranformation for compositional data
#' 
#' This function will give a warning if zeros are present because the
#' transformed data will have -Infs.
#' 
#' @param data a data frame in which the data is stored
#' @param group if not NULL then a character string specifying the name of the
#' grouping variable
#' @return a data frame with the CLR transformation applied to data. Each row
#' in the data frame is standardized by dividing by the geometric mean of that
#' row. The logarithms of the resulting ratios are returned. If a grouping
#' variable is specified, then this is preserved.
#' @author James M. Curran
#' @seealso alr
#' @references Aitchison, J. (1986). ``The Statistical Analysis of
#' Compositional Data'', Chapman and Hall, reprinted in 2003 with additional
#' material by The Blackburn Press
#' @keywords transformation
#' @examples
#' 
#' data(bottle.df)
#' 
#' ## transform preserving grouping
#' clr(bottle.df, "Number")
#' 
#' ## transform the data but remove the
#' ## grouping in column 1
#' clr(bottle.df[,-1])
#' 
#' @export
clr = function(data, group = NULL){
    geoMean = function(x){
        p = length(x)

        return(prod(x)^(1/p))
    }

    if(!is.null(group)){
        f = formula(paste(group,"~."),sep="")
        mf = model.frame(f, data)
        g = model.response(mf)
        names(g) = group
        data = mf[,-1]
    }

    gms = apply(data, 1, geoMean)

    if(any(gms==0))
        warning("Data has zeros which means the transformation will be unusable")

    data = log(data/gms)

    if(!is.null(group))
        return(data.frame(group, data))

    return(data)
}
