#' Calculate Hotelling's two sample T-squared test statistic
#' 
#' Calculate Hotelling's T-squared test statistic for the difference in two
#' multivariate means.
#' 
#' Note, the sample size requirements are that nx + ny - 1 > p. The procedure
#' will stop if this is not met and the shrinkage estimator is not being used.
#' The shrinkage estimator has not been rigorously tested for this application
#' (small p, smaller n).
#' 
#' @aliases hotelling.stat hotel.stat
#' @param x a nx by p matrix containing the data points from sample 1
#' @param y a ny by p matrix containg the data points from sample 2
#' @param shrinkage set to \code{TRUE} if the covariance matrices are to be estimated
#' using Schaefer and Strimmer's James-Stein shrinkage estimator
#' @return A list containing the following components:
#' \item{statistic}{Hotelling's (unscaled) T-squared statistic} \item{m}{The
#' scaling factor - this can be used by by multiplying it with the test
#' statistic, or dividing the critical F value} \item{df}{a vector of length
#' containing the numerator and denominator degrees of freedom} \item{nx}{The
#' sample size of sample 1} \item{ny}{The sample size of sample 2} \item{p}{The
#' number of variables to be used in the comparison}
#' @author James M. Curran
#' @references Hotelling, H. (1931). ``The generalization of Student's ratio.''
#' Annals of Mathematical Statistics 2 (3): 360--378.
#' 
#' Schaefer, J., and K. Strimmer (2005). ``A shrinkage approach to large-scale
#' covariance matrix estimation and implications for functional genomics.''
#' Statist. Appl. Genet. Mol. Biol. 4: 32.
#' 
#' Opgen-Rhein, R., and K. Strimmer (2007). ``Accurate ranking of
#' differentially expressed genes by a distribution-free shrinkage approach.''
#' Statist. Appl. Genet. Mol. Biol. 6: 9.
#' @keywords htest
#' @examples
#' 
#' data(container.df)
#' split.data = split(container.df[,-1],container.df$gp)
#' x = split.data[[1]]
#' y = split.data[[2]]
#' hotelling.stat(x, y)
#' hotelling.stat(x, y, TRUE)
#' 
#' @export
hotelling.stat = function(x, y, shrinkage = FALSE){
    ## get the sample sizes for each sample
    nx = nrow(x)
    ny = nrow(y)

    px = ncol(x)
    py = ncol(y)

    if(px!=py)
        stop("Both samples must have the same number of variables (columns)")

    p = px

    if(nx + ny < p - 1 & !shrinkage)
        stop("The sample sizes (nx + ny) must be 1 greater than the number of columns")

    mx = apply(x, 2, mean)
    my = apply(y, 2, mean)

    sx = NULL
    sy = NULL

    if(!shrinkage){
        sx = cov(x)
        sy = cov(y)
    }else{
        sx = cov.shrink(x, verbose = FALSE)
        sy = cov.shrink(y, verbose = FALSE)
    }

    sPooled = ((nx - 1) * sx + (ny - 1) * sy) / (nx + ny - 2)
    sPooledInv = solve(sPooled)

    T2 = t(mx - my) %*% sPooledInv %*% (mx - my) * nx * ny/(nx + ny)
    m = (nx + ny - p - 1)/(p * (nx + ny - 2))

    invisible(list(statistic = as.vector(T2), m = m, df = c(p, nx + ny - p - 1),
                   nx = nx, ny = ny, p = p))
}



#' Two-sample Hotelling's T-squared test
#' 
#' Performs a two-sample Hotelling's T-squared test for the difference in two
#' multivariate means
#' 
#' 
#' @aliases hotelling.test hotelling.test.default hotelling.test.formula
#' hotel.test
#' @param x a matrix containing the data points from sample 1 or a formula
#' specifying the elements to be used as a response and the grouping variable
#' as a predictor
#' @param y a matrix containing the data points from sample 2
#' @param shrinkage if \code{TRUE} then Shaefer and Strimmer's James-Stein shrinkage
#' estimator is used to calculate the sample covariance matrices
#' @param perm if \code{TRUE} then permutation testing is used to estimate the
#' non-parametric P-value for the hypothesis test
#' @param B if perm is TRUE, then B is the number of permutations to perform
#' @param progBar if \code{TRUE} and \code{perm} is TRUE then a progress bar will be
#' displayed whilst the permutation procedure is carried out
#' @param data a data frame needs to be specified if a formula is to be used to
#' perform the test
#' @param pair a vector of length two which can be used when the grouping
#' factor has more than two levels to select different pairs of groups. For
#' example for a 3-level factor, pairs could be set to \code{c(1,3)} to perform
#' Hotelling's test between groups 1 an 3
#' @param \dots any additional arguments. This is useful to pass the optional
#' arguments for the default call from the formula version
#' @return A list (which is also of class 'hotelling.test') with the following
#' elements:
#' 
#' \item{stats}{a list containing all of the output from \code{hotelling.stat}}
#' \item{pval}{the P-value from the test} \item{results}{if \code{perm == TRUE}, then
#' all of the permuation test statisics are stored in results}
#' @author James M. Curran
#' @seealso hotelling.stat
#' @references Hotelling, H. (1931). ``The generalization of Student's ratio.''
#' Annals of Mathematical Statistics 2 (3): 360--378.
#' 
#' Schaefer, J., and K. Strimmer (2005). ``A shrinkage approach to large-scale
#' covariance matrix estimation and implications for functional genomics.''
#' Statist. Appl. Genet. Mol. Biol. 4: 32.
#' 
#' Opgen-Rhein, R., and K. Strimmer (2007). ``Accurate ranking of
#' differentially expressed genes by a distribution-free shrinkage approach.''
#' Statist. Appl. Genet. Mol. Biol. 6: 9.
#' 
#' Campbell, G.P. and J. M. Curran (2009). ``The interpretation of elemental
#' composition measurements from forensic glass evidence III.'' Science and
#' Justice, 49(1),2-7.
#' @keywords htest
#' @examples
#' 
#' data(container.df)
#' fit = hotelling.test(.~gp, data = container.df)
#' fit
#' 
#' subs.df = container.df[1:10,]
#' subs.df$gp = rep(1:2, c(5,5))
#' fitPerm = hotelling.test(Al+Fe~gp, data  = subs.df, perm =  TRUE)
#' fitPerm
#' plot(fitPerm)
#' 
#' data(bottle.df)
#' fit12 = hotelling.test(.~Number, data = bottle.df)
#' fit12
#' 
#' fit23 = hotelling.test(.~Number, data = bottle.df, pair = c(2,3))
#' fit23
#' 
#' @export
hotelling.test = function(x, ...){
    UseMethod("hotelling.test")
}

#' @describeIn hotelling.test Two-sample Hotelling's T-squared test
#' @export
hotelling.test.default = function(x, y, shrinkage = FALSE, perm = FALSE,
                                  B = 10000, progBar = (perm && TRUE), ...){
    if(!perm){
        stats = hotelling.stat(x, y, shrinkage)
        pVal = with(stats, 1 - pf(m*statistic, df[1], df[2]))
        output = list(stats = stats, pval = pVal)
        class(output) = "hotelling.test"
        invisible(output)
    }else{
        stats = hotelling.stat(x, y, shrinkage)
        res = rep(0, B)

        nx = stats$nx
        ny = stats$ny
        N = nx + ny
        T0 = stats$statistic

        idx = 1:N
        X = rbind(x, y)
        
        onePercent = floor(B / 100)
        pb = NULL
        if(progBar){
          pb = txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
        }
        j = 0
        k = 0

        for(i in 1:B){
            i1 = sample(idx, nx)
            x1 = X[i1,]
            x2 = X[-i1,]

            res[i] = hotelling.stat(x1, x2, shrinkage)$statistic
            j = j + 1
            if(j == onePercent && progBar){
              k = k + 1
              j = 0
              setTxtProgressBar(pb, k)
            }
        }

        pVal = sum(res > T0)/B
        output = list(stats = stats, pval = pVal , results = res)
        class(output) = "hotelling.test"
        invisible(output)
    }
}

#' @describeIn hotelling.test Two-sample Hotelling's T-squared test
#' @export
hotelling.test.formula = function(x, data = NULL, pair = c(1,2), ...){
    if(missing(x) || class(x) != "formula")
        stop("missing or incorrect formula")

    form = x
    form[[3]] = x[[2]]
    form[[2]] = x[[3]]

    mf = model.frame(form, data)

    group = model.response(mf)
    variables = mf[,-1, drop = FALSE]

    split.data = split(variables,group)


    x1 = as.matrix(split.data[[pair[1]]])
    x2 = as.matrix(split.data[[pair[2]]])

    hotelling.test(x1, x2, ...)
}



#' Plots the results from a permutation based version of Hotelling's T-squared
#' test for the difference in two multivariate sample means
#' 
#' Plots a histogram of the distribution of the permuted test statistics for a
#' permutation version of Hotelling's T-squared
#' 
#' This function only works if you have performed a permutation test. It will
#' return an error message if not. It could be programmed to draw the relevant
#' F distribution in the standard case, but this seems rather pointless.
#' 
#' @param x an object of type hotelling.test
#' @param \dots any additional arguments to be passed to the hist command
#' @author James M. Curran
#' @keywords plot
#' @examples
#' 
#' data(bottle.df)
#' bottle.df = subset(bottle.df, Number == 1)
#' bottle.df$Number = rep(1:2,c(10,10))
#' fit = hotelling.test(.~Number, bottle.df, perm = TRUE)
#' plot(fit)
#' plot(fit, col = "lightblue")
#' 
#' @export
plot.hotelling.test = function(x,...){
    if(is.na(match("results",names(x))))
        stop("Plotting only works if you have used the permutation test")

    dotArgNames = names(list(...))
    
    if("xlim" %in% dotArgNames){
      with(x,{
         h = hist(stats$m*results, main = "Distribution of permuted test stats",
                  xlab = expression(T^2),...);
         T0 = with(stats, m*statistic)
         lines(rep(T0, 2), c(0, max(h$counts)), lwd = 2)
         }
      )
    }else{
      ## Try to do something sensible to make sure the test statistic is visible
      ## in the plotting window, ONLY if the user isn't over-riding xlim
      with(x,{
        r = range(c(stats$m * results, stats$m * stats$statistic));
        pr = pretty(r);
        h = hist(stats$m * results, main = "Distribution of permuted test stats",
             xlab = expression(T^2), xlim = c(min(pr), max(pr)), ...);
        #abline(v = with(stats, m*statistic, lwd = 2))
        T0 = with(stats, m*statistic)
        lines(rep(T0, 2), c(0, max(h$counts)), lwd = 2)
        }
      )
    }
}



#' Prints the results from a Hotelling's T-squared test for the difference in
#' two multivariate sample means
#' 
#' Prints the test stastic, degrees of freedom and P-value from Hotelling's
#' T-squared test for the difference in two multivariate sample means
#' 
#' 
#' @param x an object of type hotelling.test
#' @param \dots any additional arguments to be passed to the hist command
#' @author James M. Curran
#' @keywords print
#' @examples
#' 
#' data(bottle.df)
#' bottle.df = subset(bottle.df, Number == 1)
#' bottle.df$Number = rep(1:2,c(10,10))
#' fit = hotelling.test(.~Number, bottle.df, perm = TRUE)
#' fit
#' fit = hotelling.test(.~Number, bottle.df)
#' fit
#' 
#' ## an explict call
#' print(fit)
#' 
#' @export
print.hotelling.test = function(x, ...){
    if(is.na(match("results",names(x)))){
        with(x,{
            with(stats,{
                cat(paste("Test stat: ", signif(m*statistic, 5), "\n"));
                cat(paste("Numerator df: ", stats$df[1], "\n"))});
            cat(paste("Denominator df: ", x$stats$df[2], "\n"));
            cat(paste("P-value: ", signif(x$pval, 4), "\n"))})
    }else{
        with(x,{
            with(stats,{
                cat(paste("Test stat: ", signif(m*statistic, 5), "\n"));
                cat(paste("Numerator df: ", stats$df[1], "\n"))});
            cat(paste("Denominator df: ", x$stats$df[2], "\n"));
            cat(paste("Permutation P-value: ", signif(x$pval, 4), "\n"));
            cat(paste("Number of permutations :", length(x$results), "\n"))          })
    }
}



hotel.stat = function(x, y, shrinkage = FALSE){
    hotelling.stat(x, y, shrinkage)
}

hotel.test = function(x, ...){
    hotelling.test(x, ...)
}
