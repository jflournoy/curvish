#' derivatives
#'
#' Computes the derivative of a smooth term from the posterior of the smooth.
#'
#' @param object a \code{brmsfit} object
#' @param term character; the smooth term for which the derivative is required.
#' @param n numeric; the number of points to evaluate the derivative at.
#' @param eps numeric; the finite difference.
#'
#' @return
#' @export
#' @import brms
#'
#' @examples
derivatives <- function(object,
                        term,
                        # newdata,
                        n = 200,
                        eps = 1e-07, summary_only = TRUE){
  if(!grepl('s\\(.*\\)', term)){
    data_term <- term
    smooth_term <- sprintf('s(%s)', term)
  } else {
    data_term <- gsub('s\\((.*)\\)', '\\1', term)
    smooth_term <- term
  }
  srange <- range(object$data[data_term])
  interval <- list(seq(srange[[1]], srange[[2]], length.out = res))
  names(interval) <- data_term
  interval_eps <- interval
  interval_eps[[1]] <- interval_eps[[1]] + eps
  newd <- do.call(data.frame, c(interval, list(check.names = FALSE)))
  newd_eps <- do.call(data.frame, c(interval_eps, list(check.names = FALSE)))

  p <- brms::posterior_smooths(object = fit_b, smooth = smooth_term, newdata = newd)
  p_eps <- brms::posterior_smooths(object = fit_b, smooth = smooth_term, newdata = newd_eps)
  p_deriv <- (p_eps - p) / eps

  if(summary_only){
    p_deriv_sum <- do.call(rbind, apply(p_deriv, 2, function(x){
      q <- quantile(x, probs = c(.025, .975))
      m <- mean(x)
      data.frame(mean = m, l = q[[1]], u = q[[2]])
    }))

    p_deriv_sum[data_term] <- interval
    r <- p_deriv_sum
  } else {
    r <- list(posterior = p_deriv, interval = interval)
  }
  return(r)
}
