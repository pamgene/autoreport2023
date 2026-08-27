# Runs all tests. From the repo root: Rscript tests/testthat.R
# Requires: testthat, tidyverse, flextable (install.packages(...) if missing).

library(testthat)
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(flextable)
})

source("R/00_GeneralFunctions.R")
source("R/01_BasicProcessing.R")

test_dir("tests/testthat", reporter = "summary")
