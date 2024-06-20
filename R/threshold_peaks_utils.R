#' Find the peaks positions
#'
#' @param C2 channel intensity
#'
#' @return peak positions
find_peaks <- function(C2){

  which(splus2R::peaks(C2, span = 7))
}


#' Compute distance between peaks
#'
#' @param x_peaks positions of peaks
#' @param n number of pixels
#'
#' @return a vector of length `length(x_peaks)-1` with the distances between peaks
interpeak_dist <- function(x_peaks, n){

  x_peaks <- c(0, x_peaks, n)

  x_peaks[-1] - x_peaks[-(length(x_peaks))]

}


#' Determine threshold
#'
#' @param C2 channel intensity
#'
#' @return a threshold for `C2`
get_threshold <- function(C2){
  stats::quantile(C2, 0.5)
}


#' Define the peak domains
#'
#' @param np peak number
#' @param C2 channel intensity
#' @param C2_peaks_pos positions of peaks
#' @param C2_peaks_thresholded_pos positions of peaks above threshold
#' @param spanD.left left span of domain
#' @param spanD.right right span of domain
#' @param spread parameter
#'
#' @return a data.frame with the domains around peak
define_domain_per_peak <- function(np, C2, C2_peaks_pos, C2_peaks_thresholded_pos, spanD.left, spanD.right, spread = 3){

  pixel.size <- length(C2)
  C2_peaks_interp_dist <- interpeak_dist(C2_peaks_pos, length(C2))


  peak_nb <- which(C2_peaks_thresholded_pos == C2_peaks_thresholded_pos[np])

  bb <- ceiling(C2_peaks_thresholded_pos[np] + C2_peaks_interp_dist[peak_nb]/2)

  if(bb > pixel.size){
    x.coord <- c(floor(C2_peaks_thresholded_pos[np]-C2_peaks_interp_dist[peak_nb]/2):pixel.size)
  } else {
    x.coord <- c(floor(C2_peaks_thresholded_pos[np]-C2_peaks_interp_dist[peak_nb]/2):ceiling(C2_peaks_thresholded_pos[np]+C2_peaks_interp_dist[peak_nb]/2))  #the span to fit
  }




  halfmax <- ifelse((C2[x.coord] - min(C2[x.coord])) / (C2[C2_peaks_thresholded_pos[np]]- min(C2[x.coord])) > 0.5, 1, 0)
  decision.left <- spanD.left[x.coord] + halfmax
  decision.right <- spanD.right[x.coord] + halfmax

  # borders of C2 domain in the span
  bord.left <- utils::tail(which(decision.left[1:floor(C2_peaks_interp_dist[peak_nb]/2 - 1)] == 0), n=1)
  if(length(bord.left) == 0) bord.left <- 0

  bord.right <- which(decision.right[ceiling(C2_peaks_interp_dist[peak_nb]/2 + 1):length(decision.right)] == 0)[1] + floor(C2_peaks_interp_dist[peak_nb]/2)
  if(is.na(bord.right)) bord.right <- length(x.coord)


  #convert to worm coordinates

  border_init <- round(bord.left + x.coord[1] - spread)
  border_end <-  round(bord.right + x.coord[1] + spread)-1

  if(border_init < 1){
    border_init <- 1
  }

  if(border_end > pixel.size){
    border_end <- pixel.size
  }

  c(init = border_init, end = border_end)

}



#' Compute area under curve
#'
#' @param x x coordinates
#' @param y y coordinates
#' @param method integration method
#'
#' @details
#' If method is `trapz`, uses the Trapezoidal Integration from package `pracma`. If method is `sum`, simply
#' sums the heights of `y` (ignoring `x` altogether).
#'
#'
#' @return The area under the curve (single numeric value)
intensity_under_curve <- function(x, y, method = c("trapz", "sum")){
  method <- match.arg(method)

  if(method == "trapz"){

    return(pracma::trapz(x, y))
  } else if(method == "sum"){

    return(sum(y))
  } else{
    stop("Invalid method: ", method)
  }
}

