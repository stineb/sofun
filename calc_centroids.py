#!/usr/bin/python
#
# reg_points.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-08-28 -- created
# 2014-12-01 -- last updated
#
# ------------
# description:
# ------------
# This script calculates the centroids for a regular grid and saves them in a
# vector-style CSV file (e.g., ID, LAT, LON) where the pixel numbering begins
# in the top (north) left (west) corner.
#
# ----------
# changelog:
# ----------
# 00. created [13.08.28]
# 01. updated comments [14.02.12]
# 02. updated function doc string [14.12.01]
#
###############################################################################
## DEFINE FUNCTIONS:
###############################################################################
def calc_centroids(gv, ev, out):
    """
    Name:     calc_centroids
    Input:    - tuple, grid values (i.e., the number of columns and rows, aka 
                pixels in x and y, and pixel resolution)
              - tuple, extent values (i.e., northing, southing, westing and 
                easting, in degrees)
              - str, output file name (with path) for CSV file
    Output:   None. 
    Features: Writes to CSV file the ID, latitude and longitude of a regular
              grid (centroid locations) given the extents and pixel resolution, 
              numbering from the top-left (north-west) corner.
    """
    # Create an output header:
    header = "id,lat,lon\n"
    #
    # Assign values:
    x_tot, y_tot, dx = gv
    y_max, y_min, x_min, x_max = ev
    #
    # Process centroids:
    # ct(0) = top left corner (x_min, y_max)
    for y in xrange(y_tot):
        # Latitude changes for each y (offset):
        lat = (y_max - 0.5*dx) - y*dx
        if y == 0:
            try:
                f = open(out, "w")
            except IOError:
                print "Cannot open", out, "for writing"
            else:
                f.write(header)
                f.close()
        for x in xrange(x_tot):
            # Row-major category ID:
            ct = x_tot*y + x
            # Longitude changes for each x (offset):
            lon = (x_min + 0.5*dx) + x*dx
            try:
                f = open(out, "a")
            except IOError:
                # If file writing not available, print to stdout:
                print ct, lat, lon
            else:
                f.write("%d,%0.3f,%0.3f\n" % (ct,lat,lon))
                f.close()

###############################################################################
## DEFINE VARIABLES:
###############################################################################
# Define number of cells (x and y):
grid_cells_x = 720
grid_cells_y = 360

# Define resolution:
grid_resolution = 0.05

# Define spatial extents:
ext_north = 90.0
ext_south = -90.0
ext_west = -180.0
ext_east = 180.0

# Define output file and directory:
output_dir = "/home/user/Desktop/"
output_file = "global_hd_grid.csv"

###############################################################################
## MAIN:
###############################################################################
# Calculate grid cell centroids and save to CSV file:
grid_vals = (grid_cells_x, grid_cells_y, grid_resolution)
ext_vals = (ext_north, ext_south, ext_west, ext_east)
output_location = output_dir + output_file
calc_centroids(grid_vals, ext_vals, output_location)

