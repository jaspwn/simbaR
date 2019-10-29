#' @export

fileRenameR <- function(my_vec) {
  ## set working directory by choosing any file in directory
  
  filename <- file.choose()
  
  ## create list of files in working directory
  
  file_list <- list.files(dirname(filename), full.names = TRUE)
  
  ## function that creates a new list of filenames w/o file extension
  
  rename.file <- function(x) {
    #dir_name <- dirname(filename)
    file_ext <- tools::file_ext(x)
    old_name <- tools::file_path_sans_ext(x)
    new_name <- paste0(old_name, my_vec, ".", file_ext)
  }
  
  ## renames files in file_list to remove file extension 
  
  file.rename(file_list, sapply(file_list, rename.file))
}





