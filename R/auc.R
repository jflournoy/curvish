#' auc
#'
#' Find the area under the curve of the derivative. This function uses the
#' trapezoidal method.
#'
#' @param object An object of class \code{curvish.curve}
#' @param multimodal if the expected posterior might be multimodal, set this to
#'   \code{TRUE} so that the summary is a true, possibly discontinuous, highest
#'   posterior density interval.
#' @param adjust passed to \code{\link[stats]{density}} to adjust bandwidth.
#' @param prob the proportion of mass within the HPDI
#' @param summary_only set this to \code{FALSE} if you want all samples
#'   returned.
#'
#' @return
#' @export
#'
#' @examples
auc <- function(object, multimodal = FALSE, prob = .95, deriv = 1, abs = FALSE, adjust = 1, summary_only = FALSE){
  #note that if the derivative falls below zero, this will subtract area.
  if(!inherits(object, 'curvish.curve')){
    stop('Object must be of class "curvish.curve"')
  }
  if(is.null(object$deriv_posterior)){
    stop('Object must have derivative posterior. Rerun `curvish::derivatives()` with `deriv_posterior = TRUE`.')
  }
  data_term <- attr(object, which = 'term')$data_term

  if(deriv == 1){
    d <- object$deriv_posterior
    h <- diff(object$deriv_posterior_summary[[data_term]])
  } else if(deriv == 2){
    if(is.null(object$deriv2_posterior)) stop("No second derivative")
    d <- object$deriv2_posterior
    h <- diff(object$deriv2_posterior_summary[[data_term]])
  } else {
    d <- object$smooth_posterior
    h <- diff(object$smooth_posterior_summary[[data_term]])
  }

  posterior_auc <- apply(d, 1, function(iter){
    a <- iter[1:(length(iter) - 1)]
    b <- iter[-1]
    trapazoid_areas <- h * (a + b) / 2
    if(abs){
      s <- sum(abs(trapazoid_areas))
    } else {
      s <- sum(trapazoid_areas)
    }
    return(s)
  })

  adensity <- NULL
  x <- array(data = posterior_auc, dim = list(length(posterior_auc), 1), dimnames = list(NULL, paste0('auc_', data_term)))
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
    attr(r, 'find') <- 'AUC'
    attr(r, 'value') <- NULL
    attr(r, 'multimodal') <- multimodal
    attr(r, 'density') <- adensity
    attr(r, 'adjust') <- adjust
  }
  return(r)
}
