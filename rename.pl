#!/usr/bin/perl
#
# rename.pl
#
# written by Tyler W. Davis
# Imperial College London
#
# 2011-06-06 -- created
# 2015-11-25 -- last updated
#
# HOW TO RUN:
# > perl rename.pl
#
# ------------
# description:
# ------------
# This script performs a file rename on all files in the working directory
# based on regular expression search and replace.
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# LOAD MODULES:
# ////////////////////////////////////////////////////////////////////////////
use warnings;

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# MAIN
# ////////////////////////////////////////////////////////////////////////////
# Initialize script variables:
my $dir = ".";

# Open directory and read all files:
# NOTE: change extension in grep based on needs
#       currently reads all txt files
opendir(MYDIR, $dir) or die $!;
my @files = grep {/.*\.{1}txt$/} readdir MYDIR;
closedir(MYDIR);

# Read through each file name you found and rename it appropriately:
for my $name (@files) {
    chomp($name);

    # Save original filename:
    my $oldname = $name;

    # Renaming scheme
    # e.g., s/SEARCH_FOR/REPLACE_WITH/;
    $name =~ s/_test/_master/;

    # Rename file:
    rename $oldname, $name;
    print "$oldname, $name\n";
}
