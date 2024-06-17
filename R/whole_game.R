


whole_game <- function(wd, save.plots = TRUE, dim.png = c(1000,600), format="Berangere"){


  pixel.size <- 402 # enter number of pixel in x of images
  spread <- 3
  vc  <-  rainbow(7)










  #_________________________________________________________________________
  #                                                                         |
  # Main script -------------------------------------------------------------
  #_________________________________________________________________________|

  #initializations

  if(save.plots) if(length(list.dirs("figures"))!=1){
    cat("Warning: could not find the 'figures' subdirectory. Won't save figures.")
    save.plots <- FALSE
  }



  data <- read_condition("C:/Users/Alexis Weinreb/Projects/other/melissa/output/EN9174/")

  data |>
    tidyr::nest(.by = c(condition, individual)) |>
    dplyr::mutate(outs = purrr::map(data,
                                    threshold_peaks_one_indiv)) |>
    dplyr::select(-data) |>
    tidyr::unnest(outs) |>
    write.csv2(file.path(wd,"full.table.csv"))



}
