# Shared project root for all R scripts (works when cwd is New Thesis, R/, or R/cleaning/).
find_proj_root <- function() {
  cand <- normalizePath(getwd())
  for (i in 0:12) {
    if (dir.exists(file.path(cand, "data", "raw"))) {
      return(cand)
    }
    nc <- normalizePath(file.path(cand, ".."))
    if (identical(nc, cand)) {
      break
    }
    cand <- nc
  }
  normalizePath(getwd())
}
