# README
---------------
* LAST UPDATED: 2014-10-30
* TEAM: labprentice
* REPO: utilities (public)

## Description
This repository holds utility and general functionality scripts. 

## Files

### calc_statistics.R
* This R function calculates a variety of statistics for a given data vector.

### catfiles.pl
* This script concatenates multiple time series data files into one single file.

### diffdays.pl
* This script calculates the number of days between two given dates.

### list_dirs.R
* This function provides a list.dirs() function (similar to existing list.files() function) that returns directory names.

### rename.pl
* This script performs bulk file renames based on regular expression search and replace.

### rms.pl
* This script reads through a file (plain text or CSV) and removes all whitespace characters from lines.

### sort.pl
* This script performs lexicographical line sorting on a file.

### water_states.R
* Functions for calculating the temperature and pressure dependencies of:
    * the *density* of pure water (Tumlirz and Chen equations)
    * the *viscosity* of pure water (Huber and Vogel equations)

### zipall.pl
* This script performs bulk file compressions (e.g., gzip). It may be modified to run any bulk system commands (e.g., tar, gunzip, pdfcrop).

### area.R
* R function to calculate the surface area of a rectangular gridcell on a sphere.

### asc2cdf.jnl
* Ferret script to convert an ASCII file (gridded, geo-referenced data in text format) into NetCDF.

### cdf.write.R