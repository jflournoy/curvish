#' a_pal
#'
#' @return
#'
#' @examples
a_pal <- function(){
  paste0('#', c('003E51', '007DBA', 'FFFFFF', 'D6D2C4', 'B78B20'))
}

#' @export
#' @import ggplot2
plot.curvish.curve <- function(x, robust = FALSE, deriv = TRUE){
  apal <- curvish::a_pal()
  if(deriv){
    d <- x$deriv_posterior_summary
  } else {
    d <- x$smooth_posterior_summary
  }
  if(robust){
    yvar <- 'median'
  } else {
    yvar <- 'mean'
  }
  xvar <- attr(x, 'term')$data_term

  aplot <- ggplot(data = d, aes_string(x = xvar, y = yvar, ymin = 'l', ymax = 'u')) +
    geom_ribbon(alpha = .5, fill = apal[[5]]) +
    geom_line(color = apal[[1]]) +
    theme_minimal()
  return(aplot)
}

#' @export
#' @import ggplot2
plot.curvish.param <- function(x, robust = FALSE, ...){
  d <- with(density(x$param_posterior, ...), data.frame(x, y))
  q_ <- lapply(1:dim(x$param_posterior_sum)[[1]], function(i){
    arow <- x$param_posterior_sum[i, ]
    r <- d$x < arow['end'] & d$x > arow['begin']
    r*i
  })
  q <- do.call(cbind, q_)
  d$CI <- apply(q, 1, any)
  d$CI_group <- apply(q, 1, sum)

  aplot <- ggplot(data = d, mapping = aes(x = x, y = y)) +
    geom_line() +
    geom_area(aes(group = CI_group, fill = CI)) +
    scale_fill_manual(aesthetics = c('fill', 'color'),
                      breaks = c(TRUE, FALSE),
                      values = apal[c(3,2)],
                      labels = c('Inside HDI', ''),
                      name = '') +
    theme_minimal()
  return(aplot)
}

#' is.curvish
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
is.curvish <- function(x) {
  inherits(x, 'curvish')
}

