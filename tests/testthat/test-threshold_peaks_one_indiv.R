test_that("Uniform int under simple peak", {

  data_filt <- data.frame(values.C2 = dnorm(1:100, mean = 30),
                          values.C1 = dunif(1:100, min = 0, max = 100),
                          individual_copy = "my_test_ind")

  output <- threshold_peaks_one_individual(data_filt)

  size_dom <- 13L
  nb_pixels <- 100L

  correct_output <- tibble::tibble(
    intensity_C1_colocalized = size_dom * 1L/nb_pixels,
    intensity_C1_outside = (nb_pixels - size_dom) * 1L/nb_pixels,
    intensity_C1_total = nb_pixels * 1L/nb_pixels,
    size_C2_domains = size_dom,
    size_C2_outside = nb_pixels - size_dom,
    size_C2_total = nb_pixels,
    proportion_C1_intensity_in_domains = size_dom * 1L/nb_pixels,
    proportion_C2_area_in_domains = size_C2_domains / size_C2_total,
    ratio_domains_nondomains_intensity_area = (intensity_C1_colocalized/intensity_C1_outside) /
      (size_C2_domains/size_C2_outside),
    mean_intensity_C1_in_domains = intensity_C1_colocalized / size_C2_domains,
    mean_intensity_C1_outside = intensity_C1_outside / size_C2_outside,
    nb_of_peaks = 1L,
    C1_average_per_peak = intensity_C1_colocalized / nb_of_peaks
  ) |>
    as.data.frame()

  expect_equal(output,
               correct_output,
               tolerance = 0.05)

  # note: we need tolerance because trapz interpolates
})

test_that("Uniform int under simple peak, by sum", {

  data_filt <- data.frame(values.C2 = dnorm(1:100, mean = 30),
                          values.C1 = dunif(1:100, min = 0, max = 100),
                          individual_copy = "my_test_ind")

  output <- threshold_peaks_one_individual(data_filt, method_integration = "sum")

  size_dom <- 13L
  nb_pixels <- 100L

  correct_output <- tibble::tibble(
    intensity_C1_colocalized = size_dom * 1L/nb_pixels,
    intensity_C1_outside = (nb_pixels - size_dom) * 1L/nb_pixels,
    intensity_C1_total = nb_pixels * 1L/nb_pixels,
    size_C2_domains = size_dom,
    size_C2_outside = nb_pixels - size_dom,
    size_C2_total = nb_pixels,
    proportion_C1_intensity_in_domains = size_dom * 1L/nb_pixels,
    proportion_C2_area_in_domains = size_C2_domains / size_C2_total,
    ratio_domains_nondomains_intensity_area = (intensity_C1_colocalized/intensity_C1_outside) /
      (size_C2_domains/size_C2_outside),
    mean_intensity_C1_in_domains = intensity_C1_colocalized / size_C2_domains,
    mean_intensity_C1_outside = intensity_C1_outside / size_C2_outside,
    nb_of_peaks = 1L,
    C1_average_per_peak = intensity_C1_colocalized / nb_of_peaks
  ) |>
    as.data.frame()

  expect_equal(output,
               correct_output)
})




test_that("Uniform int under simple peak property", {


  for_all(mean.C2 = integer_bounded(5, 95,len = 1L),
          width.C2 = numeric_bounded(2, 4, len = 1L),
          property = function(mean.C2, width.C2){
            data_filt <- data.frame(values.C2 = dnorm(1:100, mean = mean.C2, sd = width.C2),
                                    values.C1 = dunif(1:100, min = 1, max = 100),
                                    individual_copy = "my_test_ind")

            output <- threshold_peaks_one_individual(data_filt)

            expect_equal(output$proportion_C1_intensity_in_domains, output$proportion_C2_area_in_domains,
                         tolerance = 0.05)
          })



})

test_that("Uniform int under simple peak property, sum method", {


  for_all(mean.C2 = integer_bounded(5, 95,len = 1L),
          width.C2 = numeric_bounded(2, 4, len = 1L),
          property = function(mean.C2, width.C2){
            data_filt <- data.frame(values.C2 = dnorm(1:100, mean = mean.C2, sd = width.C2),
                                    values.C1 = dunif(1:100, min = 1, max = 100),
                                    individual_copy = "my_test_ind")

            output <- threshold_peaks_one_individual(data_filt, method_integration = "sum")

            expect_equal(output$proportion_C1_intensity_in_domains, output$proportion_C2_area_in_domains)
          })

})


test_that("Uniform int under double peak property", {


  for_all(mean.C2 = integer_bounded(2, 94,len = 1L),
          width.C2 = numeric_bounded(2, 4, len = 1L),
          property = function(mean.C2, width.C2){
            data_filt <- data.frame(values.C2 = dnorm(1:101, mean = mean.C2, sd = width.C2) +
                                      dnorm(1:101, mean = mean.C2+2*width.C2, sd = width.C2),
                                    values.C1 = dunif(1:101, min = 1, max = 101),
                                    individual_copy = "my_test_ind")

            output <- threshold_peaks_one_individual(data_filt)

            expect_equal(output$proportion_C1_intensity_in_domains, output$proportion_C2_area_in_domains,
                         tolerance = 0.05)
          })



})



test_that("Uniform int under double peak property, sum methods", {


  for_all(mean.C2 = integer_bounded(3, 94,len = 1L),
          width.C2 = numeric_bounded(2, 4, len = 1L),
          property = function(mean.C2, width.C2){
            data_filt <- data.frame(values.C2 = dnorm(1:101, mean = mean.C2, sd = width.C2) +
                                      dnorm(1:101, mean = mean.C2+2*width.C2, sd = width.C2),
                                    values.C1 = dunif(1:101, min = 1, max = 101),
                                    individual_copy = "my_test_ind")

            output <- threshold_peaks_one_individual(data_filt, method_integration = "sum")

            expect_equal(output$proportion_C1_intensity_in_domains, output$proportion_C2_area_in_domains)
          })



})
