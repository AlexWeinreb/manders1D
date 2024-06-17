#' Plot 1 individual's profile
#' Title
#'
#' @param data data.frame for a all individuals, with columns `individual`, `values.C2`, `values.C1`
#' @param ind_name one value in `data$individual` to filter on
#'
#' @return invisibly, the input data
plot_one_worm <- function(data, ind_name="#1"){

  X <- seq_along(data$values.C2)

  plot(X,
       data$values.C2 / max(data$values.C2),
       t="l", col = 'green4',
       xlab="x (pixels)",ylab="intensity of fluo (AU)",
       main=paste("Fluorescence of worm\n",ind_name))

  graphics::lines(X,
                  data$values.C1 / max(data$values.C1),
        col="magenta2")

  graphics::legend("bottomright",c("C2","C1"),lty=c(1,1),col=c("green4","magenta3"))

  invisible(data)
}


#' Plot 1 individual's profile with threshold
#'
#' @param C2 intensity channel
#' @param C2_peaks_thresholded_pos positions of peaks
#' @param threshold threshold
#'
#' @return NULL
plot_one_thresholded <- function(C2, C2_peaks_thresholded_pos, threshold){

  plot(seq_along(C2), C2,
       t="l", col = 'green4',
       xlab="pixels",ylab="Intensity of fluorescence (AU)")

  graphics::points(C2_peaks_thresholded_pos, C2[C2_peaks_thresholded_pos],
         t="p", col = 'darkred', lwd = 2, cex = 1.5)

  graphics::abline(h=threshold, col="grey30", lty = 'dashed', lwd = 3)

  invisible(NULL)
}

