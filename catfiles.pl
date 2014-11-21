#!/usr/bin/perl
#
# catfiles.pl
#
# written by Tyler W. Davis
# Imperial College London
#
# 2012-01-15 -- created
# 2014-10-27 -- last updated
#
# HOW TO RUN:
# > perl catfiles.pl
# > <ENTER OUTPUT FILE NAME>
# > <ENTER STRING FOR FILE SEARCHING>
#
# ------------
# description:
# ------------
# This script will concatenate all the files listed in the working directory 
# that match search string into a single file. Note that the lines in each 
# file must begin with a number (e.g., time series data files where the first
# field is the date), excluding the headerline.
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
## LOAD MODULES
# /////////////////////////////////////////////////////////////////////////////
use strict;
use warnings;

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
## VARIABLE INITIALIZATION
# /////////////////////////////////////////////////////////////////////////////
my $outname = "";     # Output file name prefix
my $output = "";      # Output file name (prefix + "_All_Data.txt"
my $ans = "";         # User response to if output file already exists
my $directory = ".";  # File directory for reading / writing
my @files;            # Array for holding files found in directory
my @filelist;         # Sorted files list
my $headerline = "";  # Headerline for output file
my $file = "";        # Individual file names (from filelist)
my @lines;            # Array of lines read from file
my $line = "";        # Individual line (from lines array)

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
## MAIN
# /////////////////////////////////////////////////////////////////////////////
# Prompt for name (used to name the outfile):
print "Enter output name: ";
$outname = <>;
chomp( $outname );

# Name the output file:
$output = $outname . "_All_Data.txt";

# Check to see if that file already exists:
$ans = "y";
if (-e $output) {
    print "That file already exists, continue (y/n)? ";
    $ans = <>;
    chomp( $ans );
}

# Continue if user says it's okay:
if ( lc( $ans ) eq "y" ) 
{
    # Prompt user for a file prefix:
    print "Enter search sting for files: ";
    my $prefix = <>;
    chomp( $prefix );
    
    # Open the directory and read all files matching regex
    # and sort the files you found by their filenames:
    opendir(DIR, $directory) or die $!;
    @files = grep  { /.*$prefix.*/ } readdir(DIR);
    @filelist = sort(@files);
    closedir(DIR);
    
    # Create the output file:
    open(OUT, ">$output") or die $!;
    
    # Open the first file in the working directory:
    open(ALTFILE,"<$filelist[1]") or die $!;
    
    # Print the headerline to the out file:
    $headerline = <ALTFILE>;
    print OUT $headerline;
    close(ALTFILE);
    close(OUT);
    
    # Check to see if the array has any files:
    if (@files)  
    {
        # Open outfile for writing/appending:
	    open(OUT, ">>$output") or die $!;
        
        # Iterate through each file found:
        foreach $file (@filelist) 
        {
            print "$file\n";
            
            open(FILE, "<$file") or die $!;
            @lines = <FILE>;
            foreach $line (@lines)
            {
                if ($line =~ /^\d+.*/)
                {
                    # Should avoid blank lines.
                    print OUT $line;
                }
            }
            close(FILE);
        }
        
        # Close out file:
        close(OUT);
    }
    else
    {
        print "Found no files.\n";
    }
}
else
{
    print "Quitting.\n";
}
