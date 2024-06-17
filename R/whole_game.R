


#' Run whole script
#'
#' @param input_dir where to find the csv files
#' @param output_dir where to save results
#' @param plot_path either FALSE (do not save plots) or a path to a directory to save plots
#' @param dim.png dimensions (in pixels) of the plots to save, ignored if plot_path = FALSE
#'
#' @return Saves in the output_dir, returns the same data
#' @export
whole_game <- function(input_dir, output_dir, plot_path = FALSE,
                       dim.png = c(1000,600)){



  if(plot_path){
    dir.create(plot_path)
  }



  genotypes <- list.files(input_dir)


  if(length(genotypes) < 1){
    stop("Did not find genotypes")
  }


  data <- purrr::map_dfr(genotypes,
                         ~ read_condition(file.path(input_dir, .x)))

  res <- data |>
    tidyr::nest(.by = c("condition", "individual")) |>
    dplyr::mutate(outs = purrr::map(data,
                                    threshold_peaks_one_indiv)) |>
    dplyr::select(-data) |>
    tidyr::unnest("outs")


  utils::write.csv2(res, file.path(output_dir, "full.table.csv"))


  invisible(res)
}
