#!/usr/bin/perl
#
# sort.pl
#
# written by Tyler W. Davis
# Imperial College London
#
# 2012-06-07 -- created
# 2014-10-27 -- last updated
#
# HOW TO RUN:
# > perl sort.pl
# > <SELECT FILE BY NUMBER>
#
# ------------
# description:
# ------------
# This script reads a plain text file (or CSV) and returns a file with its 
# contents lexicographically sorted.  The headerline is preserved. The output
# file has the name '_sorted.txt' appended (original file preserved).
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# LOAD MODULES:
# ////////////////////////////////////////////////////////////////////////////
use warnings;

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# MAIN
# ////////////////////////////////////////////////////////////////////////////
# Initialize script variables:
my $counter = 0;
my $fileprefix = '';
my $writfile = "_sorted.txt";
my @lines = ();

# Ask user to select file for processing:
my $workingfile = &getFilename($fileprefix);

# Make certain file has contents and is readable:
# Check to see if the file exists and is readable:
if ( ( -s $workingfile ) && ( -r $workingfile ) )
{ 
	#print "$workingfile found!\n";
    
    # Update the writefile:
    $fileprefix = substr $workingfile, 0, -4;
    $writfile = $fileprefix . $writfile;
    
	open RFILE, "<$workingfile" or die $!;
	while (my $line = <RFILE>)
	{
		# Increment counter:
		$counter++;
		
		if ($counter == 1)
		{
			# Header line, print out:
			open WFILE, ">$writfile" or die $!;
			print WFILE $line;
			close WFILE;

		}
		else
		{
			chomp($line);
			# Save each line to an array:
			push(@lines,$line);
		}
	}
	# Read through array and print to file:
	open WFILE, ">>$writfile" or die $!;
	foreach my $item (sort(@lines))
	{
		print WFILE "$item\n";
	}
	close WFILE;
}

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# SUBROUTINES
# ////////////////////////////////////////////////////////////////////////////
sub getFilename
{
	# Reads the local directory and searches for text files based on the 
	# search prefix.
	local $prefix = $_[0];               # search prefix for text files
	local $dh;                           # file handle
	local $dir = ".";                    # directory path
	local @files = ();                   # array of found files
	local $numfiles;                     # number of files found
	local ($files, $selectedfile) = "";  # File names
	local $selection;                    # index of selected file
    
	# Read directory for files:
	opendir($dh, $dir) || die "cannot open directory $dir: $!";
	@files = grep{/^$prefix.*txt$/} readdir($dh);
	closedir($dh);
    
	# Sort and count the file names:
	@files = sort( @files );
	$numfiles = scalar(@files);
    
	if ($numfiles > 0)
	{
		# Initiate question and response variables:
		print "Select a file to process:\n";
        
		# Read through the files you found and print them:
		for (local $i = 1; $i < $numfiles+1; $i++)
		{
			$file = $files[$i-1];
			print "$i $file\n";
		}
		
		# User makes a selection
		$selection = <>;
		if ($selection >0 && $selection <= $numfiles)
		{
			$selectedfile = $files[$selection-1];
			print "You selected, $selectedfile\n";
			return $selectedfile;
		}
		else
		{
			die "You made a wrong selection.\n";
		}
	}
	else
	{
		die "No files found!\n";
	}
}
