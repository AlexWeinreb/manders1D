find_peaks <- function(C2){

  which(splus2R::peaks(C2, span = 7))
}


interpeak_dist <- function(x_peaks, n){

  x_peaks <- c(0, x_peaks, n)

  x_peaks[-1] - x_peaks[-(length(x_peaks))]

}


get_threshold <- function(C2){
  quantile(C2, 0.5)
}


define_area_under_peak_C2 <- function(np, C2, C2_peaks_pos, C2_peaks_thresholded_pos, spanD.left, spanD.right, spread = 3){
  pixel.size <- length(C2)
  C2_peaks_interp_dist <- interpeak_dist(C2_peaks_pos, length(C2))


  n.peak <- which(C2_peaks_thresholded_pos == C2_peaks_thresholded_pos[np])   # np is the rank among the peaks above thresh, n.peak the rank in all local maxima

  bb <- ceiling(C2_peaks_thresholded_pos[np]+C2_peaks_interp_dist[n.peak]/2)
  if(bb > pixel.size){
    x.coord <- c(floor(C2_peaks_thresholded_pos[np]-C2_peaks_interp_dist[n.peak]/2):pixel.size)
  } else {
    x.coord <- c(floor(C2_peaks_thresholded_pos[np]-C2_peaks_interp_dist[n.peak]/2):ceiling(C2_peaks_thresholded_pos[np]+C2_peaks_interp_dist[n.peak]/2))  #the span to fit
  }




  halfmax <- ifelse((C2[x.coord] - min(C2[x.coord])) / (C2[C2_peaks_thresholded_pos[np]]- min(C2[x.coord])) > 0.5, 1, 0)
  decision.left <- spanD.left[x.coord] + halfmax
  decision.right <- spanD.right[x.coord] + halfmax

  # borders of C2 domain in the span
  bord.left <- tail(which(decision.left[1:floor(C2_peaks_interp_dist[n.peak]/2 - 1)] == 0), n=1)
  if(length(bord.left) == 0) bord.left <- 0

  bord.right <- which(decision.right[ceiling(C2_peaks_interp_dist[n.peak]/2 + 1):length(decision.right)] == 0)[1] + floor(C2_peaks_interp_dist[n.peak]/2)
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
