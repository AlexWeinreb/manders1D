test_that("Output matches older version", {
  ref_path <- "C:/Users/Alexis Weinreb/Projects/other/melissa/output/"

  expect_no_error(whole_game(ref_path, save.plots = TRUE))

  expect_identical(read.csv2(file.path(ref_path, "full.table.csv")) |>
                     dplyr::rename(
                       worm = individual,
                       genotype = condition,
                       intensity.non.colocalized.C1 = non.colocalised.C1.all,
                       intensity.colocalized.C1 = colocalised.C1.all,
                       size.syn.domains = size.C2.domains.all,
                       size.non.syn.domains = size.C2.outside.domains.all
                     ) |>
                     dplyr::select(-c(C1.average.per.peak, C1.total, C1.total.true)) |>
                     dplyr::relocate(worm, .before = genotype),
                   read.csv2(file.path(ref_path, "ref_full.table.csv")))
})
