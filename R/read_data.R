#' Read data for single condition
#'
#' @param dir_path Path to data
#' @param plot_path either NULL (do not save plots) or a path to a directory to save plots
#' @param dim.png dimensions of the plot to save, ignored if plot_path = FALSE
#'
#' @details
#' Path to the data files for a single condition (e.g. genotype).
#'
#' If it's the path to a directory, it is assumed to contain csv files
#'  for individuals, and the directory name is assumed to describe the condition.
#'
#'
#' @return data.frame with the input data in long format
#' @export
read_condition <- function(dir_path, plot_path = NULL, dim.png = c(1000,600)){

  if(! file.exists(dir_path)){
    stop("Path not found: ", dir_path)
  }

  if(dir.exists(dir_path)){
    format <- "directory"
  } else{
    stop("Only implemented format: one directory per condition")
  }


  # In this format, there is a directory for each condition, in this dir a csv file for each individual

  input_csvs <- list.files(path=dir_path,
                           pattern="*.csv",
                           recursive=TRUE,
                           full.names = TRUE)

  if(length(input_csvs) < 1){
    stop("Did not find csv files in ", dir_path)
  }


  condition <- basename(dir_path)

  individuals <- basename(input_csvs) |>
    tools::file_path_sans_ext()
  individuals <- gsub(" ", "_", individuals)


  data <- purrr::map2_dfr(input_csvs, individuals,
                          \(path, ind_name){
                            utils::read.csv(path) |>
                              tibble::add_column(individual = ind_name)
                          })
  data[["condition"]] <- condition
  data <- data[c("condition", "individual",
                 "values.C1", "values.C2")]





  if(! is.null(plot_path)){

    for(ind_name in unique(data$individual)){

      grDevices::png(filename=paste0(plot_path, "/whole_worm_", ind_name, ".png"),
                     width = dim.png[1], height = dim.png[2])


      plot_one_worm(data[data$individual == ind_name,], ind_name)

      grDevices::dev.off()
    }
  }



  data
}
