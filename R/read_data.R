#' Title
#'
#' @param dir_path Path to data
#' @param col_individuals If the path is to a csv file, the name of the column containing the individuals
#'
#' @details
#' Path to the data files for a single condition (e.g. genotype).
#'
#' If it's the path to a directory, it is assumed to contain csv files
#'  for individuals, and the directory name is assumed to describe the condition. If it's
#'  the path to a csv file, it is assumed to correspond to a single condition, and the
#'  individuals are given by a column `col_individuals`.
#'
#'
#' @return
#' @export
#'
#' @examples
read_condition <- function(dir_path, col_individuals, save.plots = FALSE){

  if(! file.exists(dir_path)){
    stop("Path not found: ", dir_path)
  }

  if(dir.exists(dir_path)){
    format <- "directory"
  } else{
    format <- "single_file"
  }






  if(format=="directory"){

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
      tools::file_path_sans_ext() |>
      stringr::str_replace_all(" ", "_")


    data <- purrr::map2_dfr(input_csvs, individuals,
                    \(path, ind_name){
                      read.csv(path) |>
                        tibble::add_column(individual = ind_name)
                    }) |>
      tibble::add_column(condition = condition) |>
      dplyr::select(condition, individual,
                    X, values.C1, values.C2)


    # plot_one_worm(data, ind_nametmp)



      if(save.plots){
        purrr::walk(plot_one_worm)
        png(filename=paste(wd,"figures/whole_worm_",gsub(" ","_",worm),".png",sep=""),width=dim.png[1],height=dim.png[2])
        plot_one_worm(data, ind_name)
        dev.off()
      }



  }


  if(format=="Haijun"){
    # In this format, there is 1 CSV file for each genotype. All the worms
    # are in this file, with the name of the worm in col "worm".
    genotypes <- c("wt","ok_259") # names of the genotypes
    in.files <- c("wt_28.csv","ok259_29.csv")  #the CSV input files (in the same order)

    if(length(genotypes)!=length(in.files)) cat("Error in input files: lengths must be the same.\n")


    # In Haijun's format, each file contains all the worms of 1 given phenotype
    # There is a col "worm" indicating which worm is each row
    data <- read.csv2(in.files[file.gen],header=TRUE)   #read data
    file.gen <- file.gen+1
    C2 <- data[,2]                      #read column 2
    postsyn <- data[,3]                     #read column3
    len <- length(C2)
    x <- c(1:len)
    xmax <- max(x)
    if(length(postsyn)!=len) cat("Warning: lengths do differ!\n")
    data$worms <- fill.missing(data$worms)
    if(length(data$worms)!=len) cat("Warning: lengths do differ!\n")

    worms <- levels(data$worms)             #names of the worms
    worms <- worms[-which(worms=="")]

    #Create list with 1 worm = 1 item
    all.C2 <- sapply(worms,function(wrm) split(data,data$worms)[[wrm]][,2])
    all.C1 <- sapply(worms,function(wrm) split(data,data$worms)[[wrm]][,3])


    if(save.plots){
      for(w in 1:length(worms)){
        png(filename=paste(wd,"figures/whole_worm_",w,".png",sep=""),width=dim.png[1],height=dim.png[2])
        plot_one_worm(w)
        dev.off()
      }
    } else{
      #display a random example of worm
      plot_one_worm(sample(data$individual, 1))
    }

  }

  data
}
