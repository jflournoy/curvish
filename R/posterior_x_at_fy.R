#' posterior_x_at_fy
#'
#' Find the posterior distribution of the value of a predictor variable at a
#' particular point along the derivative. For instance, one might be interested
#' in knowing the age value at which a slope is steepest, or closest to zero.
#'
#' \code{posterior_x_at_miny}, \code{posterior_x_at_maxy}, and
#' \code{posterior_x_at_yequalto} are aliases of \code{posterior_x_at_fy} where
#' \code{find} is set to 'min', 'max', or 'value', respectively and are provided
#' in case they make for more readable code.
#'
#' @param object a \code{brmsfit} object
#' @param term character; the smooth term for which the derivative is required.
#' @param n numeric; the number of points to evaluate the derivative at.
#' @param eps numeric; the finite difference.
#' @param find the spot along the derivative to find; can be 'min', 'max', or
#'   'value', in which case the \code{value} parameter must be set.
#' @param value if \code{find} is set to 'value', this must be set to a specific
#'   value, like 0.
#' @param multimodal if the expected posterior might be multimodal, set this to
#'   \code{TRUE} so that the summary is a true, possibly discontinuous, highest
#'   posterior density interval.
#' @param adjust passed to \code{\link[stats]{density}} to adjust bandwidth.
#' @param prob the proportion of mass within the HPDI
#' @param summary_only set this to \code{FALSE} if you want all samples
#'   returned.
#'
#' @return An object of class \code{curvish.param} that contains the summary
#'   from \code{\link[HDInterval]{hdi}}, and the posterior samples. If
#'   \code{summary_only} is \code{TRUE}, the only the output from
#'   \code{\link[HDInterval]{hdi}} is returned.
#' @export
#'
#' @examples
posterior_x_at_fy <- function(object, term, n, eps, find = c('min', 'max', 'value'), value = NULL, multimodal = TRUE, adjust = 1, prob = .95, summary_only = FALSE){

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

  p <- curvish::derivatives(object, term = term, n = n, eps = eps, deriv_posterior = TRUE)
  interval <- p$deriv_posterior_summary[[data_term]]
  x <- apply(p$deriv_posterior, 1, function(a_sample) {
    the_spot <- find_spot(a_sample)
    closest_to_index <- abs(a_sample - the_spot) == min(abs(a_sample - the_spot))
    closest_x <- interval[which(closest_to_index)]
    if(length(closest_x) > 1){
      warning('More than one point matched. Sampling one at random.')
      closest_x <- sample(closest_x, size = 1)
    }
    return(closest_x)
  })
  adensity <- NULL
  x <- array(data = x, dim = list(length(x), 1), dimnames = list(NULL, data_term))
  if(!multimodal){
    summary_value <- HDInterval::hdi(x, credMass = prob)
  } else {
    adensity <- density(x, adjust = adjust)
    summary_value <- HDInterval::hdi(adensity, credMass = prob, allowSplit = TRUE)
  }
  if(summary_only){
    r <- summary_value
  } else {
    r <- list(param_posterior_sum = summary_value, param_posterior = x)
    attr(r,'class') <- c('curvish', 'curvish.param')
    attr(r, 'find') <- find
    attr(r, 'value') <- value
    attr(r, 'multimodal') <- multimodal
    attr(r, 'density') <- adensity
    attr(r, 'adjust') <- adjust
  }
  return(r)
}

#' @rdname posterior_x_at_fy
#' @export
posterior_x_at_miny <- function(...){
  return(posterior_x_at_fy(..., find = c('min'), value = NULL))
}

#' @rdname posterior_x_at_fy
#' @export
posterior_x_at_maxy <- function(...){
  return(posterior_x_at_fy(..., find = c('max'), value = NULL))
}

#' @rdname posterior_x_at_fy
#' @export
posterior_x_at_yequalto <- function(...){
  return(posterior_x_at_fy(..., find = c('value')))
}
