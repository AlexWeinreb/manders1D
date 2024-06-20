test_that("domain of simple peak", {
  C2 <- dnorm(1:100, mean = 30)
  C2_peaks_pos <- 30
  C2_peaks_thresholded_pos <- 30

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

  expect_identical(
    define_domain_per_peak(1, C2, C2_peaks_pos, C2_peaks_thresholded_pos, spanD.left, spanD.right, spread = 3),
    c(init = 24, end = 36)
  )


})


test_that("domain of double peak", {
  C2 <- dnorm(1:100, mean = 30) + dnorm(1:100, mean = 32, sd = .5)
  C2_peaks_pos <- 30
  C2_peaks_thresholded_pos <- 30

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

  expect_identical(
    define_domain_per_peak(1, C2, C2_peaks_pos, C2_peaks_thresholded_pos, spanD.left, spanD.right, spread = 3),
    c(init = 24, end = 37)
  )


})



