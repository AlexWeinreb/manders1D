test_that("Output matches older version", {
  ref_path <- "C:/Users/Alexis Weinreb/Projects/other/melissa/"

  input_dir = file.path(ref_path, "input/csvs")
  output_dir = file.path(ref_path, "output")
  plot_path = NULL

  expect_no_error(process_whole_directory(input_dir, output_dir))

  expect_identical(read.csv2(file.path(ref_path, "output", "full.table.csv")),
                   read.csv2(file.path(ref_path, "output", "new_ref_full.table.csv")))
})

test_that("With plots saved", {
  ref_path <- "C:/Users/Alexis Weinreb/Projects/other/melissa/"

  save_figs <- "C:/Users/Alexis Weinreb/Projects/tests/figs"
  if(dir.exists(save_figs)) unlink(save_figs, recursive = TRUE)

  expect_no_error(process_whole_directory(input_dir = file.path(ref_path, "input/csvs"),
                                          output_dir = file.path(ref_path, "output"),
                                          plot_path = save_figs))

  expect_identical(read.csv2(file.path(ref_path, "output", "full.table.csv")),
                   read.csv2(file.path(ref_path, "output", "new_ref_full.table.csv")))
})
