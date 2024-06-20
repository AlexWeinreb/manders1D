
<!-- README.md is generated from README.Rmd. Please edit that file -->

# manders1D

Given two channels (from fluorescence microscopy), use one channel for
thresholding and define “domains”, then evaluate the intensity of the
second channel within the domains.

## Installation

You can install the development version of manders1D from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AlexWeinreb/manders1D")
```

## Example step by step

Given a directory (corresponding to a condition, e.g. a genotype), which
contains csv files (corresponding to individual biological replicates),
all the files in the directory can be read with:

``` r
library(manders1D)
genotype_dir <- "data/genotype1/"

data <- read_condition(genotype_dir)
```

The resulting file has columns `condition` (the name of the directory),
`individual` (the name of the csv file), `values.C1` and `values.C2`
(which had to be columns in the csv files).

``` r
head(data)
#>   condition                                                individual values.C1
#> 1    EN9174 EN9174_20231003_DNC_YA_01_T_w1SD561-single.TIF_DNC_AD.tif   295.255
#> 2    EN9174 EN9174_20231003_DNC_YA_01_T_w1SD561-single.TIF_DNC_AD.tif   409.266
#> 3    EN9174 EN9174_20231003_DNC_YA_01_T_w1SD561-single.TIF_DNC_AD.tif   413.830
#> 4    EN9174 EN9174_20231003_DNC_YA_01_T_w1SD561-single.TIF_DNC_AD.tif   548.406
#> 5    EN9174 EN9174_20231003_DNC_YA_01_T_w1SD561-single.TIF_DNC_AD.tif   554.100
#> 6    EN9174 EN9174_20231003_DNC_YA_01_T_w1SD561-single.TIF_DNC_AD.tif   641.145
#>   values.C2
#> 1         0
#> 2         0
#> 3         0
#> 4         0
#> 5         0
#> 6         0
```

Each individual can be processed by `threshold_peaks_one_individual()`:

``` r
first_individual <- unique(data$individual)[[1]]
data_single_ind <- data[ data$individual == first_individual, ]

quantif_single_individual <- threshold_peaks_one_individual(data_single_ind)
head(quantif_single_individual)
#>   intensity_C1_colocalized intensity_C1_outside intensity_C1_total
#> 1                 308905.4             193677.5           502582.9
#>   size_C2_domains size_C2_outside size_C2_total
#> 1              77             325           402
#>   proportion_C1_intensity_in_domains proportion_C2_area_in_domains
#> 1                          0.6146358                     0.1915423
#>   ratio_domains_nondomains_intensity_area mean_intensity_C1_in_domains
#> 1                                6.731922                     4011.759
#>   mean_intensity_C1_outside nb_of_peaks C1_average_per_peak
#> 1                  595.9307           7            44129.35
```

## Example in single command

Alternatively, the command `process_whole_directory()` runs all these
steps given an input an output directory:

``` r
library(manders1D)

process_whole_directory(input_dir = "data/input/genotype1",
                        output_dir = "data/output/genotype1")
```
