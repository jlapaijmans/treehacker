
# run tests with
# Rscript ./treehacker_test.R

if(!require("testthat", quietly = TRUE)) install.packages("testthat")

library(testthat)
# remove the old test fasta file if it still exists
if (file.exists("./*win*.fasta")){
  unlink("./*win*.fasta*")
}
if (dir.exists("./output_*")){
  unlink("./output_*", recursive = TRUE)
}

# test with a simple command that it runs at all
test_that("TREEhacker runs successfully",{
  res <- system2("../TREEhacker.sh"," --threads 1 --parallel 1 fastafiles.txt test 20 20 0.5 BOTH",
          stdout=TRUE)
  # check if the output files are created
  expect_true(file.exists("fasta1_20_20s_wins.fasta"))
  # check if the master tree collection files are created
  expect_true(file.exists("test_raxml_trees.txt"))
  expect_true(file.exists("test_raxml_trees_BIN.txt"))
  # check if the topology summary files are created
  expect_true(file.exists("test_raxml_trees_TOPO.txt"))
  expect_true(file.exists("test_raxml_trees_TOPO_BIN.txt"))
  # check if master tree files have correct format (header + data)
  trees_dna <- readLines("test_raxml_trees.txt")
  trees_bin <- readLines("test_raxml_trees_BIN.txt")
  expect_true(length(trees_dna) > 1)  # Should have header + trees
  expect_true(length(trees_bin) > 1)  # Should have header + trees
  expect_true(grepl("^window\ttree$", trees_dna[1]))  # Check header format
  expect_true(grepl("^window\ttree$", trees_bin[1]))  # Check header format
  expect_true(grepl("^scaffold_.*\t.*\\);$", trees_dna[2]))  # Check data format
  expect_true(grepl("^scaffold_.*\t.*\\);$", trees_bin[2]))  # Check data format
  # check that output directory was cleaned up (should NOT exist after successful run)
  expect_false(dir.exists("output_TREEhackerFiles_test/"))
  # remove output files
  unlink("./*win*.fasta*")
  unlink("test*")
})

# test with a simple command that it runs at all with zipped input
test_that("TREEhacker runs successfully with zipped fasta.gz ",{
  res <- system2("../TREEhacker.sh"," --threads 1 --parallel 1 fastafiles_gz.txt test 20 20 0.5 BOTH",
          stdout=TRUE)
  # check if the output files are created
  expect_true(file.exists("fasta1_20_20s_wins.fasta"))
  # check if the master tree collection files are created
  expect_true(file.exists("test_raxml_trees.txt"))
  expect_true(file.exists("test_raxml_trees_BIN.txt"))
  # check if the topology summary files are created
  expect_true(file.exists("test_raxml_trees_TOPO.txt"))
  expect_true(file.exists("test_raxml_trees_TOPO_BIN.txt"))
  # check that output directory was cleaned up (should NOT exist after successful run)
  expect_false(dir.exists("output_TREEhackerFiles_test/"))
  # remove output files
  unlink("./*win*.fasta*")
  unlink("test*")
})

# test if the filtering happens correctly (example scaffold has 5 Ns out of 20 positions)
test_that("TREEhacker filters windows properly",{
  res <- system2("../TREEhacker.sh","  --threads 1 --parallel 1 fastafiles.txt test 20 20 0.5 DNA",
          stdout=TRUE)
  # check if master tree file is created and contains scaffold_3_0-20 (should pass 0.5 filter)
  expect_true(file.exists("test_raxml_trees.txt"))
  trees <- readLines("test_raxml_trees.txt")
  expect_true(any(grepl("scaffold_3_0-20", trees)))
  # check that output directory was cleaned up
  expect_false(dir.exists("output_TREEhackerFiles_test/"))
  unlink("*_wins.fasta*")
  unlink("test*")

  # but it should NOT pass the 0.1 filter
  res <- system2("../TREEhacker.sh"," --threads 1 --parallel 1 fastafiles.txt test 20 20 0.1 DNA",
          stdout=TRUE)
  # check if master tree file exists but should NOT contain scaffold_3_0-20
  if(file.exists("test_raxml_trees.txt")) {
    trees <- readLines("test_raxml_trees.txt")
    expect_false(any(grepl("scaffold_3_0-20", trees)))
  } else {
    # If no trees passed the filter, that's also acceptable
    expect_true(TRUE)
  }
  # check that output directory was cleaned up
  expect_false(dir.exists("output_TREEhackerFiles_test/"))
  unlink("*_wins.fasta*")
  unlink("test*")
  })
