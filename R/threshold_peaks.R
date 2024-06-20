
#' Threshold the peaks for single individual
#'
#' @param data_filt data.frame subsetted for one individual
#' @param plot_path either NULL (do not save plots) or a path to a directory to save plots
#' @param dim.png dimensions of the plot to save, ignored if plot_path = FALSE
#'
#' @return a data.frame with 1 row, and columns for all computed variables
#'
threshold_peaks_one_indiv <- function(data_filt, plot_path = NULL, dim.png = c(1000,600), method_integration = "trapz"){


  ind_name <- unique(data_filt$individual_copy)

  stopifnot(length(ind_name) == 1L)


  C2 <- data_filt$values.C2

  C2_peaks_pos <- find_peaks(C2)

  if(length(C2_peaks_pos) < 1){
    stop("No peak found")
  }

  threshold = get_threshold(C2)

  C2_peaks_values <- C2[C2_peaks_pos]

  C2_peaks_thresholded_pos <- C2_peaks_pos[C2_peaks_values > threshold]

  if(length(C2_peaks_thresholded_pos) < 1){
    stop("No peak found above threshold")
  }


  #plot the worms' C2
  if(! is.null(plot_path)){
    grDevices::png(filename = paste0(plot_path, "/peaks_", ind_name, ".png"),
                   width = dim.png[1], height = dim.png[2])

    plot(C2,
         t="l", col = 'green4',
         xlab="pixels",ylab="Intensity of fluorescence (AU)",
         main = paste("Thresholding of individual\n", ind_name))

    graphics::points(C2_peaks_thresholded_pos, C2[C2_peaks_thresholded_pos],
                     t="p", col = 'darkred', lwd = 2, cex = 1.5)

    graphics::abline(h=threshold, col="grey30", lty = 'dashed', lwd = 3)

    grDevices::dev.off()
  }


  # Define span: left then right. We exclude the 1st and last peaks.


  dC2 <- diff(C2)
  C2_is_increasing <- ifelse(dC2 > 0.01, 1, 0) #thresholded derivative increase>0
  C2_is_decreasing <- ifelse(dC2 < -0.01, 1, 0) #thresholded derivative decrease>0

  l <- length(C2_is_increasing)
  spanD.left <- C2_is_increasing +
    c(0, C2_is_increasing[1:(l-1)]) +
    c(0,0, C2_is_increasing[1:(l-2)])

  spanD.right <- C2_is_decreasing +
    c(C2_is_decreasing[2:l], 0) +
    c(C2_is_decreasing[3:l], 0,0)



  peak_positions <- purrr::map_dfr(seq_along(C2_peaks_thresholded_pos),
                                   ~ define_domain_per_peak(.x,
                                                            C2,
                                                            C2_peaks_pos,
                                                            C2_peaks_thresholded_pos,
                                                            spanD.left,
                                                            spanD.right,
                                                            spread = 3))



  if(! is.null(plot_path)){
    grDevices::png(filename = paste0(plot_path, "/domains_worm_", ind_name, ".png"),
                   width=dim.png[1],height=dim.png[2])


    pal_alpha = grDevices::hcl.colors(nrow(peak_positions),
                                      palette = "Dark 2",
                                      alpha = .1)
    pal_noalpha = grDevices::hcl.colors(nrow(peak_positions),
                                        palette = "Dark 2")

    plot(C2,
         t='l', col = 'green4', lwd = 2,
         xlab = "X (pixel)", ylab = "Channel C2 intensity",
         main = paste("Domains of individual\n", ind_name))

    graphics::points(C2_peaks_thresholded_pos, C2_peaks_values*1.01,
                     pch = 25, col = pal_noalpha, bg = pal_noalpha,
                     cex = 1.5)

    graphics::rect(peak_positions$init, xright = peak_positions$end,
                   ybottom = 0, ytop = 1.1*max(C2_peaks_values),
                   col = pal_alpha, border = NA)

    grDevices::dev.off()
  }



  # Measure only the part under peaks
  x_under_peaks <- purrr::map2(peak_positions$init, peak_positions$end,
                               \(i,e) i:e) |>
    purrr::reduce(union)

  C1_worm <- data_filt$values.C1

  C1_under_peaks <- C1_worm
  C1_under_peaks[-x_under_peaks] <- 0
  # plot(C1_under_peaks)

  intensity_C1_colocalized <- intensity_under_curve(seq_along(C1_under_peaks), C1_under_peaks,
                                                    method = method_integration)
  intensity_C1_total <- intensity_under_curve(seq_along(C1_worm), C1_worm,
                                              method = method_integration)

  ratio <- intensity_C1_colocalized/intensity_C1_total

  size_C2_domains <- length(x_under_peaks)
  total_length_px <- length(C1_worm)
  size_C2_outside <- total_length_px - size_C2_domains

  intensity_C1_outside <- intensity_C1_total - intensity_C1_colocalized

  data.frame(
    intensity_C1_colocalized = intensity_C1_colocalized,
    intensity_C1_outside = intensity_C1_outside,
    intensity_C1_total = intensity_C1_total,

    size_C2_domains = size_C2_domains,
    size_C2_outside = size_C2_outside,
    size_C2_total = total_length_px,

    proportion_C1_intensity_in_domains = intensity_C1_colocalized / intensity_C1_total,
    proportion_C2_area_in_domains = size_C2_domains / total_length_px,

    # ratio of (C1coloc/C1noncoloc) / (C2domains/C2nondomains)
    ratio_domains_nondomains_intensity_area = intensity_C1_colocalized * size_C2_outside /
      (intensity_C1_outside * size_C2_domains),

    # intensity of area under curve of C1 signal within C2 domains divided by the size of the C2 domains
    mean_intensity_C1_in_domains = intensity_C1_colocalized / size_C2_domains,
    # intensity of area under curve of C1 signal outside C2 domains divided by the size outside C2 domains
    mean_intensity_C1_outside = intensity_C1_outside/size_C2_outside,

    nb_of_peaks = length(C2_peaks_thresholded_pos),
    C1_average_per_peak  =  intensity_C1_colocalized/length(C2_peaks_thresholded_pos)
  )

}
