#' Hotelling
#' @description A set of R functions and data sets which implements Hotelling's T^2 test, and some variants of it. 
#' Functions are also included for Aitchison's additive log ratio and centred log ratio transformations.
#' 
#' @name Hotelling-package
#' @aliases Hotelling-package Hotelling
#' @docType package
#' @author James Curran Maintainer: James Curran <j.curran@@auckland.ac.nz> ~~
#' The author and/or maintainer of the package ~~
#' @importFrom stats cov formula model.frame model.response
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom corpcor cov.shrink
#' @keywords package
NULL

#' Bottle data
#' 
#' This data contains the elemental concentration of five different elements
#' (Manganese, Barium, Strontium, Zirconium, and Titanium) in samples of glass
#' taken from six different Heineken beer bottles. 20 measurements were taken
#' from each bottle.
#' 
#' 
#' @name bottle.df
#' @docType data
#' @references R. L. Bennett. \emph{Aspects of the analysis and interpretation
#' of glass trace evidence}. Master's thesis, Department of Chemistry,
#' University of Waikato, 2002.
#' @keywords datasets
NULL


#' Container data
#' 
#' This data contains the elemental concentration of nine different elements
#' (Titanium, Aluminium, Iron, Manganese, Magnesium, Calcium, Barium,
#' Strontium, and Zirconium) in specimens of glass taken from two different
#' containers. Ten measurements were taken from each container.
#' 
#' 
#' @name container.df
#' @docType data
#' @references Jose R. Almirall. Discrimination of glass samples by solution
#' based ICP-OES PhD thesis, Department of Chemistry, Florida International
#' University, 1998.
#' @keywords datasets
NULL

#' manova1 data
#' 
#' The data contains example data for testing the unequal variance option in the 
#' package. The dataset has four varibles, \code{wratr}, \code{wrata}, \code{treatment}, 
#' and \code{disability}. \code{treatment is the grouping variable} and \code{wratr} and \code{wrata}
#' are the responses. \code{disability} is not used.
#' 
#' @name manova1.df
#' @docType data
#' @references \url{https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Hotellings_Two-Sample_T2.pdf}
#' @keywords datasets
NULL
