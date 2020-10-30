#' derivatives
#'
#' Computes the derivative of a smooth term from the posterior of the smooth.
#'
#' @param object a \code{brmsfit} object
#' @param term character; the smooth term for which the derivative is required.
#' @param n numeric; the number of points to evaluate the derivative at.
#' @param eps numeric; the finite difference.
#' @param prob numeric; the density of the posterior summary.
#' @param prob_outer numeric; the outer density of the posterior summary.
#' @param deriv_posterior return the posterior samples for the derivative as part of the object?
#'
#' @return a \code{curvish} object containing summaries of the posterior derivative, posterior smooths.
#' @export
#' @import brms
#'
#' @examples
derivatives <- function(object,
                        term,
                        # newdata,
                        n = 200,
                        eps = 1e-07,
                        prob = .95,
                        prob_outer = .99,
                        deriv_posterior = TRUE){
  if(!grepl('s\\(.*\\)', term)){
    data_term <- term
    smooth_term <- sprintf('s(%s)', term)
  } else {
    data_term <- gsub('s\\((.*)\\)', '\\1', term)
    smooth_term <- term
  }
  srange <- range(object$data[data_term])
  interval <- list(seq(srange[[1]], srange[[2]], length.out = n))
  names(interval) <- data_term
  interval_eps <- interval
  interval_eps[[1]] <- interval_eps[[1]] + eps
  newd <- do.call(data.frame, c(interval, list(check.names = FALSE)))
  newd_eps <- do.call(data.frame, c(interval_eps, list(check.names = FALSE)))

  p <- brms::posterior_smooths(object = object, smooth = smooth_term, newdata = newd)
  p_eps <- brms::posterior_smooths(object = object, smooth = smooth_term, newdata = newd_eps)
  p_deriv <- (p_eps - p) / eps

  probs <- sort(c(.5, unlist(lapply((1-c(prob, prob_outer))/2,
                                    function(x) c(0, 1) + c(1, -1)*x))))

  p_curve_sum <- do.call(rbind, apply(p, 2, function(x){
    q <- quantile(x, probs = probs)
    m <- mean(x)
    data.frame(mean = m, median = q[[3]], ll = q[[1]], l = q[[2]], u = q[[4]], uu = q[[5]])
  }))
  p_deriv_sum <- do.call(rbind, apply(p_deriv, 2, function(x){
    q <- quantile(x, probs = probs)
    m <- mean(x)
    data.frame(mean = m, median = q[[3]], ll = q[[1]], l = q[[2]], u = q[[4]], uu = q[[5]])
  }))

  p_deriv_sum[data_term] <- interval
  p_curve_sum[data_term] <- interval
  r <- list(deriv_posterior_summary = p_deriv_sum,  smooth_posterior_summary = p_curve_sum)
  if(deriv_posterior){
    r <- c(r, list(deriv_posterior = p_deriv))
  }
  attr(r,'class') <- 'curvish'
  attr(r,'class') <- 'curvish.curve'
  attr(r, 'term') <- list(data_term = data_term, smooth_term = smooth_term)
  attr(r, 'resolution') <- res
  attr(r, 'prob') <- c(prob, prob_outer)
  return(r)
}
