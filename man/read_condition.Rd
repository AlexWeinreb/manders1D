% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_data.R
\name{read_condition}
\alias{read_condition}
\title{Read data for single condition}
\usage{
read_condition(dir_path, plot_path = NULL, dim.png = c(1000, 600))
}
\arguments{
\item{dir_path}{Path to data}

\item{plot_path}{either NULL (do not save plots) or a path to a directory to save plots}

\item{dim.png}{dimensions of the plot to save, ignored if plot_path = FALSE}
}
\value{
data.frame with the input data in long format
}
\description{
Read data for single condition
}
\details{
Path to the data files for a single condition (e.g. genotype).

If it's the path to a directory, it is assumed to contain csv files
 for individuals, and the directory name is assumed to describe the condition.

Each csv file is expected to have columns named "values.C1" and "values.C2"
}
