# ************************************************************************
# * Name: list.dirs
# *
# * Input: - char (path)
# *        - char (pattern)
# *        - boolean (all.dirs)
# *        - boolean (full.names)
# *        - boolean (ignore.case)
# *
# * Return: char vector
# *
# * Features: This reads all the directory names for a given path
# *
# *
# * Ref: J. Ulrich (2011) How to obtain a list of directories within a
# *        directory, list.files(), but instead "list.dirs()", Stack
# *        Overflow, url: http://stackoverflow.com/questions/4749783/how-
# *        to-obtain-a-list-of-directories-within-a-directory-like-list-
# *        files-but-i
# ************************************************************************
list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
  all <- list.files(path, pattern, all.dirs,
                    full.names, recursive=FALSE, ignore.case)
  return(all)
}