#!/usr/bin/python
#
# ocr.py
#
# based on jpg2tiff.py and tif2txt.py
#
# written by Tyler W. Davis 
# created: 2011-02-17
# updated 2014-11-26
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script converts JPG files in the working directory to plain text 
# 1. Processes JPG images with a background noise filter (requires imagemagick) 
#    and saves the processed JPG files with a "_m" extension (i.e., preserving 
#    the original JPG files).
# 2. Converts processed JPG files to TIF image format (required for OCR).
# 3. Converts TIF images to plain text
#
# ~~~~~~
# notes:
# ~~~~~~
# These image processing functions depend on imagemagick 
#    http://www.imagemagick.org 
# which is a freeware for Windows, Mac and Linux. Other libraries that should 
# be installed include: libtiff4 (TIFF library), libjpeg62 (JPEG runtime 
# library), libpng12-0 (PNG runtime library), zliblg (runtime compression 
# library).
#
# The OCR processing requires the package tesseract, which was originally 
# developed by HP and is currently being maintained by Google. To download the 
# OCR software on Linux machines, use a package manager and search for 
# 'tesseract-ocr' and include any languages you require (e.g., 
# tesseract-ocr-eng for English). More information regarding the
# current state of tesseract OCR can be found at the Google code website:
#    http://code.google.com/p/tesseract-ocr/
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# 01. added threshold option to monochrome function [11.03.10]
# 02. added language option to totxt function [11.03.10]
# 03. updated function doc strings [14.11.24]
# 04. implemented glob for file searching [14.11.26]
#
###############################################################################
## IMPORT MODULES:
###############################################################################
import glob
import os

###############################################################################
# FUNCTIONS
###############################################################################
def findfiles(my_ext):
    """
    Name:     findfiles
    Input:    -str, directory name (my_dir)
              -str, file extension (my_ext)
    Output:   glob.glob list
    Features: Returns a list of file names the local directory based on the
              given search file extension
    """
    my_list = glob.glob("*" + my_ext)
    return (my_list)

def monochrome(myjpg, thresh=90):
    """
    Name:     monochrome
    Input:    -str, image file name (myjpg)
              -int, threshold value (thresh)
    Output:   None.
    Features: Processes JPG image with text with threshold filter to remove
              background noise.
    """
    myext = ".jpg"
    jpgbase = ""
    if myjpg.endswith(myext):
        #jpgbase holds the file name without the extension
        jpgbase = myjpg[:-len(myext)]  
    #
    mycmd = ("convert -threshold " + str(thresh) + "% " + 
             myjpg + " " + jpgbase + "_m.jpg")
    os.system(mycmd)

def totif(myjpg):
    """
    Name:     totif
    Input:    str, image file name (myjpg)
    Output:   None.
    Features: Converts JPG image to TIF format
    """
    myext = ".jpg"
    jpgbase = ""
    if myjpg.endswith(myext):
        jpgbase = myjpg[:-len(myext)]
    #
    mycmd = "convert " + myjpg + " " + jpgbase + ".tif"
    os.system(mycmd)

def totxt(mytif, lang = "eng"):
    """
    Name:     totxt
    Input:    -str, image file name (mytif)
              -str, OCR language (lang)
    Output:   None.
    Features: Converts TIF image file to TXT using tesseract OCR
    """
    myext = ".tif"
    mybase = ""
    if mytif.endswith(myext):
        mybase = mytif[:-len(myext)]
    #
    mycmd = "tesseract " + mytif + " " + mybase + " -l " + lang
    os.system(mycmd)

###############################################################################
# MAIN
###############################################################################
my_jpgs = findfiles('jpg')
if my_jpgs:
    # Process JPGs with noise filter:
    for name in my_jpgs:
        monochrome(name, 75)
    #
    # Convert JPGs to TIFs
    my_jpg_ms = findfiles('_m.jpg')
    for name in my_jpg_ms:
        totif(name)
    #
    # Convert TIFs to text:
    my_tifs = findfiles('_m.tif')
    for name in my_tifs:
        totxt(name)
else:
    print "Did not find any JPG files to process."
