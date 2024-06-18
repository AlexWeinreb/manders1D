#' Plot 1 individual's profile
#' Title
#'
#' @param data data.frame for a all individuals, with columns `individual`, `values.C2`, `values.C1`
#' @param ind_name one value in `data$individual` to filter on
#'
#' @return invisibly, the input data
plot_one_worm <- function(data, ind_name){

  X <- seq_along(data$values.C2)

  plot(X,
       data$values.C2 / max(data$values.C2),
       t="l", col = 'green4',
       xlab="x (pixels)",ylab="intensity of fluo (AU)",
       main=paste("Fluorescence of individual\n",ind_name))

  graphics::lines(X,
                  data$values.C1 / max(data$values.C1),
        col="magenta2")

  graphics::legend("bottomright",c("C2","C1"),lty=c(1,1),col=c("green4","magenta3"))

  invisible(data)
}




