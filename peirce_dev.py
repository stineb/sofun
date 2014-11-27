#!/usr/bin/python
#
# peirce_dev.py
#
# written by Tyler W. Davis
# Imperial College London
#
# created: 2013-07-15
# updated: 2013-11-27
#
# ------------
# description:
# ------------
# This script removes outliers from a dataset based on observation pairs and
# a model fit using Peirce's criterion.
#
# ----------
# changelog:
# ----------
# 01. Had to change the way Q was calculated (factor N exponent inside)
#     because it was throwing an OverflowError otherwise (N**N for N > 144)
# 02. Published on wikipedia.org [13.07.16]
# 03. Check for N > 1 [13.10.20]
# 04. Check for nan in Lamda calculation [13.10.20]
# 05. Return 0 for negative x2 values [13.10.20]
# 06. Updated peirce_dev variable names [14.11.27]
# 07. Added remove_outliers function [14.11.27]
#
################################################################################
## IMPORT MODULES:
################################################################################
import matplotlib.pyplot as plt
import numpy
import scipy.special

################################################################################
## FUNCTIONS:
################################################################################
def peirce_dev(pN, pn, pm):
    """
    Name:     peirce_dev
    Input:    - int, total number of observations (pN)
              - int, number of outliers to be removed (pn)
              - int, number of model unknowns (pm)
    Output:   float, squared error threshold (x2)
    Features: Returns the squared threshold error deviation for outlier 
              identification using Peirce's criterion based on Gould's
              methodology
    Ref:      Gould, B. A. (1855), On Peirce's criterion for the rejection of
                doubtful observations, with tables for facilitating its 
                application, Astronomical Journal, 4(11), 81-87.
              Peirce, B. (1852),  Criterion for the rejection of doubtful
                observations, Astronomical Journal, 2(21), 161-163.
    """ 
    # Assign floats to input variables:
    pN = float(pN)
    pn = float(pn)
    pm = float(pm)
    #
    # Check number of observations:
    if pN > 1:
        # Calculate Q (Nth root of Gould's equation B):
        # Note: 1/N exponent is factored to each individual term to prevent
        # OverflowError with large N (e.g., >142)
        Q = (pn**(pn/pN)*(pN - pn)**((pN - pn)/pN))/pN
        #
        # Initialize R values (as floats):
        Rnew = 1.0
        Rold = 0.0  # <- Necessary to prompt while loop
        #
        while ( abs(Rnew - Rold) > (pN*2.0e-16) ):
            # Calculate Lamda (1/(N-n)th root of Gould's equation A'):
            ldiv = Rnew**pn
            if ldiv == 0:
                ldiv = 1.0e-6
            Lamda = ((Q**pN)/(ldiv))**(1.0/(pN - pn))
            #
            # Calculate x-squared (straight-forward Gould's equation C):
            x2 = 1.0 + (pN - pm - pn)*(1.0 - Lamda**2)/pn
            #
            # If x2 goes negative, return as 0.0:
            if x2 < 0:
                x2 = 0.0
                Rold = Rnew
            else:
                # Use x-squared to update R (Gould's equation D):
                Rold = Rnew
                Rnew = numpy.exp((x2 - 1.0)/2.0)*scipy.special.erfc(
                    numpy.sqrt(x2)/numpy.sqrt(2.0)
                )
    else:
        x2 = 0.0
    return x2

def remove_outliers(x_data, y_data, y_hat, m):
    """
    Name:     remove_outliers
    Input:    - numpy.ndarray, observation data (x_data)
              - numpy.ndarray, observation data (y_data)
              - numpy.ndarray, modeled data (y_hat)
              - int, number of model parameters (m)
    Output:   - numpy.ndarray, outlier-free x data (x_data_ro)
              - numpy.ndarray, outlier-free y data (y_data_ro)
    Features: Returns x and y data with outliers removed based on Peirce's
              criterion.
    Depends:  peirce_dev
    Ref:      Peirce, B. (1852),  Criterion for the rejection of doubtful
                observations, Astronomical Journal, 2(21), 161-163.
    """
    # Define Peirce's variables:
    peirce_N = len(y_data)
    peirce_n = 1.0
    peirce_m = m
    #
    # Calculate the mean-squared error of fitted data:
    se = (y_data - y_hat)**2
    mse = se.sum()/(peirce_N - peirce_m)
    #
    # Calculate Peirce's tolerance:
    peirce_x2 = peirce_dev(peirce_N, peirce_n, peirce_m)
    peirce_delta2 = mse*peirce_x2
    #
    # Count number of outliers identified:
    outliers_index = numpy.where(se > peirce_delta2)[0]
    outliers_found = len(outliers_index)
    #
    # Run a second time if no outliers were found:
    if outliers_found == 0:
        peirce_n = 2
        peirce_x2 = peirce_dev(peirce_N, peirce_n, peirce_m)
        peirce_delta2 = mse*peirce_x2
        outliers_index = numpy.where(se > peirce_delta2)[0]
        outliers_found = len(outliers_index)
        #
        # Reset Peirce's n:
        peirce_n = 1
    #
    # Loop Peirce's n until all outliers are accounted for:
    while peirce_n <= outliers_found:
        peirce_n += 1
        #
        # Check that n hasn't exceeded N:
        if peirce_n >= peirce_N:
            peirce_n = outliers_found + 1
        else:
            peirce_x2 = peirce_dev(peirce_N, peirce_n, peirce_m)
            peirce_delta2 = mse*peirce_x2
            outliers_index = numpy.where(se > peirce_delta2)[0]
            outliers_found = len(outliers_index)
    #
    # Remove outliers:
    x_data_ro = numpy.delete(x_data, outliers_index)
    y_data_ro = numpy.delete(y_data, outliers_index)
    #
    return(x_data_ro, y_data_ro)

################################################################################
## MAIN PROGRAM:
################################################################################
my_file = 'peirce_example.csv'
my_data = numpy.loadtxt(my_file, delimiter=',', skiprows=1,
                        dtype={'names' : ('x_obs', 'y_obs', 'y_hat'),
                               'formats' : ('f4', 'f4', 'f4')})

x_ro, y_ro = remove_outliers(my_data['x_obs'], my_data['y_obs'], 
                             my_data['y_hat'], 3)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(my_data['x_obs'], my_data['y_obs'], 'ro', label='outliers')
ax1.legend(loc=2)
ax1.plot(x_ro, y_ro, 'bo')
ax1.plot(my_data['x_obs'], my_data['y_hat'], 'c_')
ax1.set_ylabel('y_data')
ax1.set_xlabel('x_data')
plt.show()
