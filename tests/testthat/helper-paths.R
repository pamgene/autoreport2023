# Sourced automatically by testthat before tests run. test_dir() sets the working
# directory to tests/testthat/ for the duration of the run, so paths into the repo's
# real data/ folder need to go back up two levels.
repo_root <- normalizePath(file.path(getwd(), "..", ".."))

test_input_path <- function(...) file.path(repo_root, "data", "test inputs", ...)
qc_input_path <- test_input_path # kept as an alias; existing QC tests use this name

# Some tests need to run a file-discovery function (read_qc_dir, read_phosphosite_dir,
# read_kinase_dir) against a temp folder pre-populated with copies of real fixture files,
# so processing side effects (e.g. process_uka_allvsall() writing split files) don't touch
# the repo's real data/test inputs/ folder. `files` is a named vector: new filename (as it
# should appear in the temp input folder) -> source fixture path.
make_input_dir <- function(files, subfolder = NULL) {
  tmp <- tempfile("input_dir_")
  target <- if (is.null(subfolder)) tmp else file.path(tmp, subfolder)
  dir.create(target, recursive = TRUE)
  for (new_name in names(files)) {
    file.copy(files[[new_name]], file.path(target, new_name))
  }
  tmp
}

# Some functions (e.g. render_ref_qc_table()) read files via a path relative to the repo
# root (matching how knitr renders the real Rmd files, always from the repo root). Tests
# run with the working directory set to tests/testthat/, so wrap calls to such functions
# in with_repo_root() to temporarily switch there and back.
with_repo_root <- function(code) {
  old <- setwd(repo_root)
  on.exit(setwd(old))
  code()
}
