#' contiguous_zeros
#'
#' @param x
#' @param colname
#' @param indexcol
#'
#' @return
#'
#' @examples
contiguous_zeros <- function(x, colname, indexcol = NULL){
  x$gt0__ <- as.numeric(x[colname] != 0)
  if(is.null(indexcol)){
    indexcol <- 'indexcol__'
    x[indexcol] <- 1:dim(x)[[1]]
  } else {
    x['indexcol__'] <- x[indexcol]
    indexcol <- 'indexcol__'
  }
  x$diff__ <- abs(x$gt0__ - c(x$gt0__[1], x$gt0__[-dim(x)[[1]]]))
  x$diff_cum__ <- cumsum(x$diff__)
  x$is_zero__ <- x[, colname] == 0
  block_lengths <- lapply(split(x, x$diff_cum__), function(d){
    n <- dim(d)[[1]]
    m <- mean(unlist(d[colname]))
    s <- min(d[indexcol])
    e <- max(d[indexcol])
    return(data.frame(n = n, mean = m, is_zero = unique(d[, 'is_zero__']), start = s, end = e))
  })
  block_lengths_df <- do.call(rbind, block_lengths)
  return(list(block_lengths = block_lengths_df, index = x$diff_cum__, is_zero = x$is_zero__))
}

#' a_pal
#'
#' @return
#'
#' @examples
a_pal <- function(){
  paste0('#', c('003E51', '007DBA', 'FFFFFF', 'D6D2C4', 'B78B20'))
}

#' plot.curvish.curve
#'
#' @param x An object of class \code{curvish.curve}.
#' @param robust This sets the parameter of central tendency to the median.
#' @param deriv Order of derivative to plot, if any (default 0 plots the smooth).
#' @param outer Use the outer credible interval bands?
#'
#' @export
#' @import ggplot2
plot.curvish.curve <- function(x, robust = FALSE, deriv = 0, outer = FALSE){
  apal <- curvish:::a_pal()
  if(deriv == 1){
    d <- x$deriv_posterior_summary
  } else if(deriv == 2){
    if(is.null(x$deriv2_posterior_summary)) stop("No second derivative")
    d <- x$deriv2_posterior_summary
  } else {
    d <- x$smooth_posterior_summary
  }
  if(robust){
    yvar <- 'median'
  } else {
    yvar <- 'mean'
  }
  xvar <- attr(x, 'term')$data_term

  aplot <- ggplot(data = d, aes_string(x = xvar, y = yvar, ymin = 'l', ymax = 'u'))
  if(outer){
    aplot <- aplot +
      geom_ribbon(aes_string(ymin = 'll', ymax = 'uu'),
                  alpha = .2, fill = apal[[5]])
  }
  aplot <- aplot +
    geom_ribbon(alpha = .5, fill = apal[[5]]) +
    geom_line(color = apal[[1]]) +
    theme_minimal()
  return(aplot)
}


#' plot.curvish.param
#'
#' @param x An object of class \code{curvish.param}.
#' @param robust If the type of summary is multimodal and this is \code{TRUE},
#'   this sets the parameter of central tendency to the median.
#' @param range Two-element numeric vector specifying limits on the X-axis range
#'   of the posterior to consider for the plot.
#' @param mode If the type of summary is multimodal, default is to use mode as
#'   the central tendency. If this is \code{FALSE}, will use the median of the
#'   posterior within the relevant portion of the HPDI.
#' @param adjust Will be passed to the \code{density} function. If not set, will
#'   attempt to take it from the object.
#' @param histogram Add a histogram overlay.
#'
#' @return
#' @export
#'
#' @examples
plot.curvish.param <- function(x, robust = FALSE, range = NULL, mode = TRUE, adjust = NULL, histogram = FALSE){
  apal <- curvish:::a_pal()
  if(!is.null(range)){
    top <- range[[2]]
    bottom <- range[[1]]
    x$param_posterior <- x$param_posterior[x$param_posterior < top & x$param_posterior > bottom]
  }

  p_cent <- NULL

  if(!is.null(attr(x, 'multimodal')) && attr(x, 'multimodal')){
    orig_adjust <- attr(x, 'adjust')
    if(is.null(adjust) & is.null(range)){
      adensity <- attr(x, 'density')
    } else {
      if(is.null(adjust))
        adjust <- 1
      if(orig_adjust != adjust){
        warning('Adjusting density bandwidth to be different from that used to compute the HPDI.')
      }
      adensity <- density(x$param_posterior, adjust = adjust)
    }

    d <- with(adensity, data.frame(x, y))

    q_ <- lapply(1:dim(x$param_posterior_sum)[[1]], function(i){
      arow <- x$param_posterior_sum[i, ]
      r <- d$x < arow['end'] & d$x > arow['begin']
      r*i
    })
    p_cent <- do.call(rbind, lapply(1:dim(x$param_posterior_sum)[[1]], function(i){
      arow <- x$param_posterior_sum[i, ]
      if(mode){
        cent_y <- max(d$y[d$x < arow['end'] & d$x > arow['begin']])
      } else {
        d_within <- x$param_posterior[x$param_posterior > arow['begin'] & x$param_posterior < arow['end']]
        cent <- median(d_within)
        cent_y <- d$y[which(abs(d$x - cent) == min(abs(d$x - cent)))]
      }
      which_y <- which(d$y == cent_y)
      which_ys <- (which_y - 1):(which_y + 1)
      which_ys <- which_ys[which_ys > 0]
      ys <- d$y[which_ys]
      xs <- d$x[which_ys]
      return(rbind(data.frame(x = xs, y = ys, type = 'areas', index = i),
                   data.frame(x = d$x[which_y], y = cent_y, type = 'points', index = i)))
    }))
    q <- do.call(cbind, q_)
    d$CI <- apply(q, 1, any)
  } else {
    if(is.null(adjust)){
      adjust = 1
    }
    adensity <- density(x$param_posterior, adjust = adjust)
    d <- with(adensity, data.frame(x, y))

    d$CI <- d$x < x$param_posterior_sum['upper',] & d$x > x$param_posterior_sum['lower',]
    if(robust){
      cent_x <- median(x$param_posterior)
    } else {
      cent_x <- mean(x$param_posterior)
    }
    cent_y <- d$y[which(abs(d$x - cent_x) == min(abs(d$x - cent_x)))]
    which_y <- which(d$y == cent_y)
    which_ys <- (which_y - 1):(which_y + 1)
    which_ys <- which_ys[which_ys > 0]
    ys <- d$y[which_ys]
    xs <- d$x[which_ys]
    p_cent <- rbind(data.frame(x = xs, y = ys, type = 'areas', index = 1),
                    data.frame(x = d$x[which_y], y = cent_y, type = 'points', index = 1))
  }

  d$CI_group <- curvish:::contiguous_zeros(d, colname = 'CI')$index

  prob <- attr(x$param_posterior_sum, 'credMass')

  aplot <- ggplot(data = d, mapping = aes(x = x, y = y)) +
    geom_line() +
    geom_area(aes(group = CI_group, fill = CI), alpha = .5) +
    geom_area(data = p_cent[p_cent$type == 'areas',],
              aes(x = x, y = y, group = index),
                 fill = apal[[1]], alpha = .8) +
    # geom_point(data = p_cent[p_cent$type == 'points',],
    #            aes(y = 0, x = x), color = apal[[5]]) +
    scale_fill_manual(aesthetics = c('fill', 'color'),
                      breaks = c(TRUE),
                      values = apal[c(3,2)],
                      labels = c(sprintf('%0.1f%% HDI', prob*100)),
                      name = '') +

    theme_minimal()
  if(histogram){
    aplot <- aplot + geom_histogram(data = data.frame(x = x$param_posterior[, 1]), aes(x = x, y = ..density..), alpha = .2)
  }
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

#' @rdname is.curvish
#' @export
is.curvish.param <- function(x) {
  inherits(x, 'curvish.param')
}

#' @rdname is.curvish
#' @export
is.curvish.curve <- function(x) {
  inherits(x, 'curvish.curve')
}
