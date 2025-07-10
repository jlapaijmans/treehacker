
# run tests with
# Rscript ./test_consensify.R

if(!require("testthat", quietly = TRUE)) install.packages("testthat")
if(!require("ape", quietly = TRUE)) install.packages("ape")


library(testthat)
# remove the old test fasta file if it still exists
if (file.exists("./*win*.fasta")){
  file.remove("./*win*.fasta")
}
if (dir.exists("./output_*")){
  dir.remove("./output_*")
}

# test with a simple set with full scaffolds
test_that("TREEhacker runs successfully",{
  res <- system2("../TREEhacker_1.0_raxml_DNA-BIN_NsFilter.sh"," fastafiles.txt test 20 20 0.5 DNA",
          stdout=TRUE)
  expect_true(file.exists("fasta1_20_20s_wins.fasta"))
  expect_true(dir.exists("output_TREEhackerFiles_test/"))
})