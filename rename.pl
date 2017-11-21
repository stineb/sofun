#!/usr/bin/perl
#
# rename.pl
#
# written by Tyler W. Davis
# Imperial College London
#
# 2011-06-06 -- created
# 2016-05-21 -- last updated
#
# HOW TO RUN:
# > perl rename.pl [-t]
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
use strict;
use warnings;

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# MAIN
# ////////////////////////////////////////////////////////////////////////////
# Initialize script variables:
my $dir = ".";

# Open directory and read all files:
# NOTE: change extension in grep based on needs
#       currently reads all jpg files
opendir(MYDIR, $dir) or die $!;
my @files = grep {/.*\.{1}jpg$/} readdir MYDIR;
closedir(MYDIR);

# Print warning if no files found
if ( @files < 1 ) {
    print "0 files found\n";
}

# Read through each file name you found and rename it appropriately:
for my $name (@files) {
    chomp($name);

    # Save original filename:
    my $oldname = $name;

    # Renaming scheme
    # e.g., s/SEARCH_FOR/REPLACE_WITH/;
    $name =~ s/[\s,]/_/g;
    $name =~ s/'//g;
    $name =~ s/&/_and_/g;
    $name =~ s/_-_/-/g;
    $name =~ s/_+/_/g;

    # Rename file:
    if ( @ARGV > 0 ) {
        if ( $ARGV[0] eq "-t" ) {
            print "$oldname, $name\n";
        } else {
            print "usage: perl rename.pl [-t]\n";
            print "       -t for testing or nothing to execute\n";
        }
    } else {
        rename $oldname, $name;
    }
}
