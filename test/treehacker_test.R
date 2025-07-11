
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

# test with a simple command that it runs at all
test_that("TREEhacker runs successfully",{
  res <- system2("../TREEhacker_1.0_raxml_DNA-BIN_NsFilter.sh"," fastafiles.txt test 20 20 0.5 BOTH",
          stdout=TRUE)
  # check if the output files are created
  expect_true(file.exists("fasta1_20_20s_wins.fasta"))
  # check if the output directory is created
  expect_true(dir.exists("output_TREEhackerFiles_test/"))
  # check if the tree files are created in the output directory
  expect_true(file.exists("*bestTree"))
  # remove output files
  file.remove("*_wins.fasta")
  file.remove("test*")
  # remove output directory
  unlink("output_TREEhackerFiles_test", recursive = TRUE)
})

# test with a simple set with full scaffolds
test_that("TREEhacker filters windows properly",{
  res <- system2("../TREEhacker_1.0_raxml_DNA-BIN_NsFilter.sh"," fastafiles.txt test 20 20 0.5 DNA",
          stdout=TRUE)
  # check if NO output files are created for scaffold_3_0-20 as it shouldnt pass 0.5 filter
  expect_false(file.exists("output_TREEhackerFiles_test/scaffold_3_0-20.concat.RN.fasta.raxml.bestTree"))
  file.remove("*_wins.fasta")
  file.remove("test*")
  # remove output directory
  unlink("output_TREEhackerFiles_test", recursive = TRUE)

  # but it should pass the 0.1 filter
  res <- system2("../TREEhacker_1.0_raxml_DNA-BIN_NsFilter.sh"," fastafiles.txt test 20 20 0.1 DNA",
          stdout=TRUE)
  expect_true(file.exists("output_TREEhackerFiles_test/scaffold_3_0-20.concat.RN.fasta.raxml.bestTree"))
  file.remove("*_wins.fasta")
  file.remove("test*")
  # remove output directory
  unlink("output_TREEhackerFiles_test", recursive = TRUE)
  })