# Load file from path (supports CSV and Excel)
load_file <- function(file_path) {
  if (! file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  if (grepl("\\.xlsx$", file_path)) {
    return(readxl::read_excel(file_path))
  } else if (grepl("\\.csv$", file_path)) {
    return(readr::read_csv(file_path, show_col_types = FALSE))
  } else {
    stop("File must be CSV or Excel (.xlsx) format")
  }
}

# Validate output path
validate_output_path <- function(path) {
  dir_path <- dirname(path)
  if (!dir.exists(dir_path)) {
    stop(paste("Output directory does not exist:", dir_path))
  }
  if (!file.access(dir_path, mode = 2) == 0) {
    stop(paste("Output directory is not writable:", dir_path))
  }
  return(TRUE)
}