# Load file from path (supports CSV and Excel)
load_file <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }
  
  # Accept .csv or .xlsx
  if (grepl("\\.xlsx?$", file_path, ignore.case = TRUE)) {
    return(readxl::read_excel(file_path))
  } else if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    return(readr::read_csv(file_path, show_col_types = FALSE))
  } else {
    stop("File must be CSV or Excel (.xlsx) format")
  }
}

# Validate output path (directory exists and writable)
validate_output_path <- function(path) {
  dir_path <- dirname(path)
  if (!dir.exists(dir_path)) {
    stop(paste("Output directory does not exist:", dir_path))
  }
  # file.access returns 0 on success; mode=2 checks write permission
  if (file.access(dir_path, mode = 2) != 0) {
    stop(paste("Output directory is not writable:", dir_path))
  }
  return(TRUE)
}
