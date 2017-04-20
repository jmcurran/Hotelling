#' Additive log ratio transformation
#' 
#' Aitchison's additive log ratio tranformation for compositional data
#' 
#' This function will give a warning if zeros are present because the
#' transformed data will have -Infs.
#' 
#' @param form a formula which specifies the denominator variable as the
#' response
#' @param data a data frame in which the data is stored
#' @param group if not NULL then a character string specifying the name of the
#' grouping variable
#' @return a data frame with the ALR transformation applied to data. Each row
#' in the data frame is standardized with respect to a specific variable by
#' dividing by that variable. The logarithms of the resulting ratios are
#' returned. If a grouping variable is specified, then this is preserved.
#' @author James M. Curran
#' @seealso clr
#' @references Aitchison, J. (1986). ``The Statistical Analysis of
#' Compositional Data'', Chapman and Hall, reprinted in 2003 with additional
#' material by The Blackburn Press
#' @keywords transformation
#' @examples
#' 
#' data(bottle.df)
#' 
#' ## transform with respect to manganese
#' alr(Mn~., bottle.df, "Number")
#' 
#' ## transform the data with respect to barium, but removing the
#' ## grouping in column 1
#' alr(Ba~., bottle.df[,-1])
#' @export
alr = function(form, data, group = NULL){

    if(!is.null(group)){
        f = formula(paste(group,"~."),sep="")
        mf = model.frame(f, data)
        g = model.response(mf)
        names(g) = group
    }
    mf = model.frame(form, data)
    denom = model.response(mf)
    variables = mf[,-1]

    if(sum(variables==0) > 0)
        warning("Data has zeros which will make transformed data unusable")

    variables = log(variables/denom)

    if(!is.null(group))
        return(data.frame(g, variables))

    return(variables)
}
