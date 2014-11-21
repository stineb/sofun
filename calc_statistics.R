# ************************************************************************
# * Name: calc_statistics
# *
# * Input: numeric vector (data)
# *
# * Return: list of named statistics
# *
# * Features: Returns a list containing a variety of statistics 
# *           calculcated for a given vector of data 
# *
# ************************************************************************
calc_statistics <- function(data){
  # 
  stats.vals <- c()
  stats.vals$max <- max(data)
  stats.vals$min <- min(data)
  stats.vals$ave <- mean(data)
  stats.vals$std <- sd(data)
  stats.vals$skew <- sum((data - stats.vals$ave)^3)/((length(data)-1)*stats.vals$std^3)
  stats.vals$kurt <- sum((data - stats.vals$ave)^4)/((length(data)-1)*stats.vals$std^4) - 3
  stats.vals
}