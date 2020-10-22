#' Title
#'
#' Find the posterior distribution of the value of a predictor variable at a
#' particular point along the derivative. For instance, one might be interested
#' in knowing the age value at which a slope is steepest, or closest to zero.
#'
#' @param object a \code{brmsfit} object
#' @param term character; the smooth term for which the derivative is required.
#' @param n numeric; the number of points to evaluate the derivative at.
#' @param eps numeric; the finite difference.
#' @param find the spot along the derivative to find; can be 'min', 'max', or 'value', in which case the \code{value} parameter must be set.
#' @param value if \code{find} is set to 'value', this must be set to a specific value, like 0.
#' @param multimodal if the expected posterior might be multimodal, set this to \code{TRUE} so that the summary is a true, possibly discontinuous, highest posterior density interval.
#' @param prob the proportion of mass within the HPDI
#' @param summary_only set this to \code{FALSE} if you want all samples returned.
#'
#' @return either a summary of the posterior density, or a list containing the summary and the posterior samples.
#' @export
#'
#' @examples
get_X_at_slope_point <- function(object, term, n, eps, find = c('min', 'max', 'value'), value = NULL, multimodal = FALSE, prob = .95, summary_only = TRUE){

  if(!find %in% c('max', 'min', 'value')){
    stop('Argument `find` must be "max", "min", or "value"')
  }

  if(find %in% c('max', 'min')){
      find_spot <- get(find)
  } else {
    if(is.null(value))
      stop('`find` is set to "value", but `value` is not provided')
    find_spot <- function(x) value
  }

  if(!grepl('s\\(.*\\)', term)){
    data_term <- term
    smooth_term <- sprintf('s(%s)', term)
  } else {
    data_term <- gsub('s\\((.*)\\)', '\\1', term)
    smooth_term <- term
  }

  p <- curvish::derivatives(object, term = term, n = n, eps = eps, summary_only = FALSE)
  x <- apply(p$posterior, 1, function(a_sample) {
    the_spot <- find_spot(a_sample)
    closest_to_index <- abs(a_sample - the_spot) == min(abs(a_sample - the_spot))
    deriv_p$interval[[data_term]][which(closest_to_index)]
  })
  x <- array(data = x, dim = list(length(x), 1), dimnames = list(NULL, data_term))
  if(!multimodal){
    summary_value <- HDInterval::hdi(x, credMass = prob)
  } else {
    summary_value <- HDInterval::hdi(density(x), credMass = prob, allowSplit = TRUE)
  }
  if(summary_only){
    r <- summary_value
  } else {
    r <- list(summary = summary_value, posterior = x)
  }
  return(r)
}
