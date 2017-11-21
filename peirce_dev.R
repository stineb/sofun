# R-Studio
#
# peirce_dev.R
#
# written by Tyler W. Davis
# Imperial College London
#
# created: 2013-07-15
# updated: 2014-11-27
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script removes outliers from a dataset based on observation pairs and
# a model fit using Peirce's criterion.
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# 01. updated documentation [14.11.27]
# 02. added remove_outliers function [14.11.27]
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### Define functions #########################################################
# /////////////////////////////////////////////////////////////////////////////

# ************************************************************************
# * Name: peirce_dev
# *
# * Input: - double, total number of observations (pN)
# *        - double, number of outliers to be removed (pn)
# *        - double, number of unknowns in model (pm)
# *
# * Return: double, squared error threshold (x2)
# *
# * Features: Returns the squared threshold error deviation for outlier 
# *           identification using Peirce's criterion based on Gould's
# *           methodology
# *
# * Ref: Gould, B. A. (1855), On Peirce's criterion for the rejection 
# *        of doubtful observations, with tables for facilitating its 
# *        application, Astronomical Journal, 4(11), 81-87.
# *      Peirce, B. (1852), Criterion for the rejection of doubtful
# *        observations, Astronomical Journal, 2(21), 161-163.
# *
# ************************************************************************
peirce_dev <- function(pN, pn, pm){
  # Check number of observations:
  if (pN > 1){
    # Calculate Q (Nth root of Gould's equation B):
    # Note: 1/N exponent is factored to each individual term to prevent
    # OverflowError with large N (e.g., >142)
    Q <- (pn^(pn/pN)*(pN - pn)^((pN - pn)/pN))/pN
    #
    # Initialize R values (as floats):
    Rnew <- 1.0
    Rold <- 0.0  # <- Necessary to prompt while loop
    #
    while(abs(Rnew - Rold) > (pN*2.0e-16)){
      # Calculate Lamda (1/(N-n)th root of Gould's equation A'):
      ldiv <- Rnew^pn
      if (ldiv == 0){
        ldiv <- 1.0e-6
      }
      Lamda <- ((Q^pN)/(ldiv))^(1.0/(pN - pn))
      #
      # Calculate x-squared (straight-forward Gould's equation C):
      x2 <- 1.0 + (pN - pm - pn)*(1.0 - Lamda^2)/pn
      #
      # If x2 goes negative, return as zero:
      if (x2 < 0){
        x2 <- 0
        Rold <- Rnew
      } else {
        #
        # Use x-squared to update R (Gould's equation D):
        # NOTE: error function (erfc) is replaced with pnorm:
        # source: 
        # http://stat.ethz.ch/R-manual/R-patched/library/stats/html/Normal.html
        Rold <- Rnew
        Rnew <- exp((x2 - 1)/2.0)*(2*pnorm(sqrt(x2)/sqrt(2)*sqrt(2), 
                                           lower=FALSE)
                                   )
      }
    }
  } else {
    x2 <- 0
  }
  return(x2)
}

# ************************************************************************
# * Name: remove_outliers
# *
# * Input: - numeric, observation abscissa data (x_data)
# *        - numeric, observation ordinate data (y_data)
# *        - numeric, modeled data (y_hat)
# *        - double, number of model fitting parameters (m)
# *
# * Return: data frame of outlier free x and y data
# *
# * Features: Returns outlier free data pairs
# *
# * Depends: peirce_dev()
# *
# * Ref: Peirce, B. (1852) Criterion for the rejection of doubtful
# *        observations, Astronomical Journal, 2(21), 161-163.
# *
# ************************************************************************
remove_outliers <- function(x_data, y_data, y_hat, m){
  # Define Peirce's parameters:
  peirce_N <- length(y_data)
  peirce_n <- 1
  peirce_m <- m
  #
  # Calculate the mean squared error of fit:
  se <- (y_data - y_hat)^2
  mse <- sum(se)/(peirce_N - peirce_m)
  #
  # Calculate Peirce's tolerance:
  peirce_x2 <- peirce_dev(peirce_N, peirce_n, peirce_m)
  peirce_delta2 <- mse*peirce_x2
  #
  # Count number of ouliers identified:
  outliers_index <- (se > peirce_delta2)
  outliers_found <- length(which(outliers_index == "TRUE"))
  #
  # Run a second time if no outliers were found:
  if (outliers_found == 0){
    peirce_n <- 2
    peirce_x2 <- peirce_dev(peirce_N, peirce_n, peirce_m)
    peirce_delta2 <- mse*peirce_x2
    outliers_index <- (se > peirce_delta2)
    outliers_found <- length(which(outliers_index == "TRUE"))
    #
    # Reset Peirce's n:
    peirce_n <- 1
  }
  # Loop Peirce's n until all outliers are accounted for:
  while (peirce_n <= outliers_found){
    peirce_n <- peirce_n + 1
    #
    # Check that n hasn't exceeded N:
    if (peirce_n >= peirce_N){
      peirce_n <- outliers_found + 1
    } else {
      peirce_x2 <- peirce_dev(peirce_N, peirce_n, peirce_m)
      peirce_delta2 <- mse*peirce_x2
      outlier_index <- (se > peirce_delta2)
      outliers_found <- length(which(outlier_index == "TRUE"))
    }
  }
  #
  # Find where exceedance occurs and remove outliers:
  non_outlier_index <- (se < peirce_delta2)
  x_data_ro <- x_data[non_outlier_index]
  y_data_ro <- y_data[non_outlier_index]
  outlier_free_data <- cbind(x_data_ro, y_data_ro)
  outlier_free_data <- as.data.frame(outlier_free_data)
  return(outlier_free_data)
}

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### Main program #############################################################
# /////////////////////////////////////////////////////////////////////////////
my_data = read.csv('peirce_example.csv')
my_data_ro = remove_outliers(my_data$x_obs, my_data$y_obs, my_data$y_mod, 3)

par(mar=c(4.5,4.5,1,1))
plot(my_data$x_obs, my_data$y_obs, pch=16, col='red', 
     xlab='x_data', ylab='y_data')
points(my_data_ro$x_data_ro, my_data_ro$y_data_ro, pch=16, col='gray50')
points(my_data$x_obs, my_data$y_mod, col='blue', pch='-')
legend('topright', c('outliers'), pch=c(16), col=c('red'), bty='n')
