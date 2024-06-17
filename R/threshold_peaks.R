
#' Threshold the peaks for single individual
#'
#' @param data_filt data.frame subsetted for one individual
#' @param plot_path either FALSE (do not save plots) or a path to a directory to save plots
#' @param dim.png dimensions of the plot to save, ignored if plot_path = FALSE
#'
#' @return a data.frame with 1 row, and columns for all computed variables
#'
threshold_peaks_one_indiv <- function(data_filt, plot_path = FALSE, dim.png = c(1000,600)){




  C2 <- data_filt$values.C2

  C2_peaks_pos <- find_peaks(C2)


  threshold = get_threshold(C2)

  C2_peaks_values <- C2[C2_peaks_pos]

  C2_peaks_thresholded_pos <- C2_peaks_pos[C2_peaks_values > threshold]




  # #plot the worms' C2
  # if(plot_path){
  #   grDevices::png(filename=paste(plot_path, "/peaks_",worm,".png",sep=""),
  #                  width=dim.png[1],height=dim.png[2])
  #   plot_one_thresholded(C2, C2_peaks_thresholded_pos, threshold)
  #
  #   grDevices::dev.off()
  # }


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
                                   ~ define_area_under_peak_C2(.x,
                                                               C2,
                                                               C2_peaks_pos,
                                                               C2_peaks_thresholded_pos,
                                                               spanD.left,
                                                               spanD.right,
                                                               spread = 3))





  ## Computing the ratio of C2 and non C2 areas
  colocalised.C1 <- 0
  size.C2.domains <- 0
  size.C2.outside.domains <- 0
  non.colocalised.C1 <- 0



  #right of last peak, excluding the last peak
  pos_last_peak <- C2_peaks_thresholded_pos[length(C2_peaks_thresholded_pos)]
  pos_penultimate_peak <- C2_peaks_thresholded_pos[length(C2_peaks_thresholded_pos)-1]

  last <- round((pos_last_peak - pos_penultimate_peak)/2 + pos_penultimate_peak)

  # if(save.plots){
  #   png(filename = paste(wd,"figures/C2_domains_worm_",worm,".png",sep=""),width=dim.png[1],height=dim.png[2])
  #   plot(all.C2[[worm]],t='l')
  #   points(xp,all.C2[[worm]][xp],col=vc)
  #   abline(v=inits,col=vc)
  #   abline(v=ends,col=vc)
  #   #abline(v=beg,col="red")
  #   #abline(v=last,col="red")
  #   dev.off()
  # }


  # look for overlaps
  inits <- peak_positions$init
  ends <- peak_positions$end
  if(ends[length(inits)-1]>inits[length(inits)]){ # should solve issues when the last peak overlaps, together with other loop below
    ends[length(inits)-1]=ends[length(inits)]
    inits[length(inits)]=inits[length(inits)-1]
  }
  if(ends[1]>inits[2]){ # should solve issues when the first peak overlaps
    ends[1]=ends[2]
    inits[2]=inits[1]
    print("overlaps in the first peak")
  }

  fluo_under_first_domain <- pracma::trapz(inits[1]:ends[1],
                                           data_filt$values.C1[inits[1]:ends[1]])
  fluo_under_last_domain <- pracma::trapz(inits[length(inits)]:ends[length(inits)],
                                          data_filt$values.C1[inits[length(inits)]:ends[length(inits)]])

  colocalised.C1 <- fluo_under_first_domain + fluo_under_last_domain

  size_first_domain <- length(inits[1]:ends[1])
  size_last_domain <- length(inits[length(inits)]:ends[length(inits)])

  size.C2.domains <- size_first_domain + size_last_domain



  for(domain in 2:(length(inits)-1)){
    if(ends[domain]>inits[domain+1]){ #should solve issues with overlapping peaks
      ends[domain]=ends[domain+1]
      inits[domain+1]=inits[domain]
      colocalised.C1 <- colocalised.C1+pracma::trapz(inits[domain]:ends[domain],
                                             data_filt$values.C1[inits[domain]:ends[domain]]) #synaptic receptors fluo # actually colocalized intensity (area under curve)
      size.C2.domains <- size.C2.domains+length(inits[domain]:ends[domain])  #length of the presyn domain. # size of C2-containing domains
      print("overlapping peaks")
    } else if(inits[domain]<ends[domain-1]){
      next
    } else {
      colocalised.C1 <- colocalised.C1+pracma::trapz(inits[domain]:ends[domain],
                                             data_filt$values.C1[inits[domain]:ends[domain]]) #synaptic receptors fluo # actually colocalized intensity (area under curve)
      size.C2.domains <- size.C2.domains+length(inits[domain]:ends[domain])  #length of the presyn domain. # size of C2-containing domains
    }
  }

  if(ends[length(inits)-1]>inits[length(inits)]){ # should solve issues when the last peak overlaps, together with other loop above
    colocalised.C1 <- colocalised.C1-pracma::trapz(inits[length(inits)]:ends[length(inits)-1],
                                           data_filt$values.C1[inits[length(inits)]:ends[length(inits)-1]])
    size.C2.domains <- size.C2.domains-length(inits[length(inits)]:ends[length(inits)-1])
    print("overlaps in the last peak")
  }

  pixel.size <- nrow(data_filt)

  intensity.C1.total <- pracma::trapz(1:pixel.size, data_filt$values.C1[1:pixel.size])
  non.colocalised.C1 <- intensity.C1.total-colocalised.C1
  size.C2.outside.domains <- pixel.size - size.C2.domains

  data.frame(
    ratio = colocalised.C1*size.C2.outside.domains/(non.colocalised.C1*size.C2.domains), # ratio of coloc / nonColoc
    ratio.area = colocalised.C1/non.colocalised.C1,
    ratio.total = colocalised.C1/(colocalised.C1+non.colocalised.C1), # ratio of coloc / total
    colocalised.C1.all = colocalised.C1,
    non.colocalised.C1.all = non.colocalised.C1,
    size.C2.domains.all = size.C2.domains,
    size.C2.outside.domains.all = size.C2.outside.domains,
    mean.intensity.C1.coloc.C2 = colocalised.C1/size.C2.domains, 		# intensity of area under curve of C1 signal within C2 domains divided by the size of the C2 domains
    mean.intensity.C1.non.coloc.C2 = non.colocalised.C1/size.C2.outside.domains,		# intensity of area under curve of C1 signal outside C2 domains divided by the size outside C2 domains
    synapse.all = length(C2_peaks_thresholded_pos),   # nb peaks

    C1.average.per.peak  =  colocalised.C1/length(C2_peaks_thresholded_pos), #average C1 intensity per peak
    C1.total  =  colocalised.C1+non.colocalised.C1, #sum of C1 intensities in each domains
    C1.total.true  =  pracma::trapz(inits[1]:ends[length(ends)],
                            data_filt$values.C1[inits[1]:ends[length(ends)]]) #total area under curve, should be similar to C1.total

  )

}
