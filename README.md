# README
---------------
* LAST UPDATED: 2016-10-23
* TEAM: labprentice
* REPO: utilities (public)

## Description
This repository holds utility and general functionality scripts.

## Files

### area.R
* R function to calculate the surface area of a rectangular gridcell on a sphere.

### asc2cdf.jnl
* Ferret script to convert an ASCII file (gridded, geo-referenced data in text format) into NetCDF.

### batch_audio_cvt.sh
* Performs batch audio conversions (e.g., lossless m4a to lossy mp3) using libav-tools and mediainfo (Debian Linux)

### biomisation.R
* R function returns the biome as a function of GDD, fractional plant cover of 9 PFTs (LPX) and vegetation height.

### calc_centroids.R
* Writes CSV file (ID, LAT, LON) of regular grid pixel centroids for given map extents and pixel resolution.

### calc_statistics.R
* This R function calculates a variety of statistics for a given data vector.

### catfiles.pl
* This script concatenates multiple time series data files into one single file.

### cdf.write.R
* R function to write a numeric array (maximum 4 dimensions: lon, lat, xxx, time) into a NetCDF file.

### diffdays.pl
* This script calculates the number of days between two given dates.

### etsrad
* Python (__etsrad.py__) and R (__etsrad.R__) scripts that calculate the half-hourly extraterrestrial solar radiation flux, daylight hours, sunset hour, daily solar irradation, half-hourly PAR, and daily PAR (as well as other parameters) for a given day and location.

### get.f.luc.R
* Get global CO2 emissions due to land use change from outputs of two LPX-Bern simulations (with and without landuse).

### get.nbp.R
* R function to get global NBP (total land-atmosphere C flux) from a LPX-Bern simulation.

### get_mean_seasonal_cycle.jnl
* Ferret script to calculate monthly climatology from a multi-annual time series.

### list_dirs.R
* This function provides a list.dirs() function (similar to existing list.files() function) that returns directory names.

### mean.bymonth.R
* R function to calculate monthly climatology from a multi-annual time series (analogue to 'get_mean_seasonal_cycle.jnl').

### mycolorbar.R
* R function to add a colorbar to a plot.

### ocr.py
* This script converts JPG images to text files (e.g., those pesky journal articles online where each page is an image file).
* Depends on [imagemagick](http://www.imagemagick.org/) (for image processing) and [tesseract](https://code.google.com/p/tesseract-ocr/) (for OCR) software packages.
* Includes __ocr_example.jpg__ for testing purposes.

### peirce_dev
* Python (__peirce_dev.py__) and R (__peirce_dev.R__) scripts that remove outliers from observation pairs based on a model fit using Peirce's criterion.
* Example data is available (__peirce_example.csv__)

### plot_biome.R
* R function creates a PDF with a nice map of biomes, given biome categorised using function 'biomisation.R'.

### regrid_landuse.R
* R function to regrid (remap) landuse data conserving total area in categories cropland, pasture, urban.

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
