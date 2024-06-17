#A function to plot 1 worm's profile
plot_one_worm <- function(data, ind_name="#1"){
  data_filt <- filter(data,
                      individual == ind_name)


  plot(data_filt$X,
       data_filt$values.C2 / max(data_filt$values.C2),
       t="l", col = 'green3',
       xlab="x (pixels)",ylab="intensity of fluo (AU)",
       main=paste("Fluorescence of worm\n",ind_name))
  lines(data_filt$X,
        data_filt$values.C1 / max(data_filt$values.C1),
        col="magenta3")
  legend("bottomright",c("C2","C1"),lty=c(1,1),col=c("green3","magenta3"))
  invisible(data)
}


plot_one_thresholded <- function(C2, C2_peaks_thresholded_pos, threshold){
  plot(seq_along(C2), C2,
       t="l", col = 'green4',
       xlab="pixels",ylab="Intensity of fluorescence (AU)")
  points(C2_peaks_thresholded_pos, C2[C2_peaks_thresholded_pos],
         t="p", col = 'darkred', lwd = 2, cex = 1.5)
  abline(h=threshold, col="grey30", lty = 'dashed', lwd = 3)
}

