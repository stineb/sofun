#!/usr/bin/perl
#
# zipall.pl
#
# written by Tyler W. Davis
# Imperial College London
#
# 2012-03-01 -- created
# 2013-10-27 -- last updated
#
# HOW TO RUN:
# > perl zipall.pl
#
# ------------
# description:
# ------------
# This script reads all the files (based on given file extension)
# and performs a system call on each one (e.g., gzip).
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# LOAD MODULES:
# ////////////////////////////////////////////////////////////////////////////
use warnings;

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# MAIN
# ////////////////////////////////////////////////////////////////////////////
# Initialize script variables:
my ($command, $file, $fprefix, $tarfile);
my $dir = ".";
my $dh;
my @dirfiles;

# Read all files in directory:
# NOTE: change extension in grep based on needs
#       currently reads all csv files
opendir($dh, $dir) || die $!;
@dirfiles = grep { /.*\.{1}csv$/ } readdir($dh);
closedir $dh;

# Save number of files found:
my $numfiles = scalar( @dirfiles );

if ($numfiles > 0)
{
    # If files were found, process them:
    foreach $file (sort @dirfiles)
    {
        # Get file prefix, i.e. remove .csv:
        $fprefix = substr $file, 0, -4;

        # Define gzip file:
        $zfile = $fprefix . ".zip";

        # Define command:
        $command = "zip $zfile $file";
        #$command = "unzip $file";
        # $command = "gzip $file";
        # $command = "gunzip $file";

        # Debugging: check the command to see if it is right:
        print "$command\n";

        # Run the command:
        system( $command );
    }
}
else
{
    print "No files found.\n";
}
