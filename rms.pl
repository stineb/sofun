#!/usr/bin/perl
#
# rms.pl
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-03-28 -- created
# 2014-10-28 -- last updated
#
# HOW TO RUN:
# > perl rms.pl
# > <SELECT FILE TO PROCESS>
#
# ------------
# description:
# ------------
# This script reads through a file (e.g., CSV) and removes any whitespace 
# found on each line (preserves line returns).
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# LOAD MODULES:
# ////////////////////////////////////////////////////////////////////////////
use warnings;

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# MAIN
# ////////////////////////////////////////////////////////////////////////////
# Initialize variables
my $headerline = "";
my @content = ();

# Ask user to select a file:
my $readfile = &getFilename();

# Create the output filename (based on readfile):
my $writefile = $readfile;
$writefile =~ s/.txt/-out.txt/;

if ( (-s $readfile) && (-r $readfile) )
{
    # Open file for reading, save header and content and close file:
    open RFILE, "<$readfile" or die $!;
    $headerline = <RFILE>;
    @content = <RFILE>;
    close RFILE;
    print "Read " . scalar(@content) . " lines.\n";
    
    # Create writefile with header:
    &writeout($writefile, $headerline);
    
    # Open write file for appending and read through content:
    open WFILE, ">>$writefile" or die $!;
    foreach my $line (@content)
    {
        # Remove newline character:
        chomp($line);
        
        # Remove any whitespace characters globally:
        $line =~ s/\s+//g;
        
        # Append line to write file:
        print WFILE "$line\n";
    }
    close WFILE;
}

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# SUBROUTINES
# ////////////////////////////////////////////////////////////////////////////
sub getFilename{
    # Reads the local directory and searches for text files
    local $dh;                           # file handle
    local $dir = ".";                    # directory path
    local @files = ();                   # array of found files
    local $numfiles;                     # number of files found
    local ($files, $selectedfile) = "";  # File names
    local $selection;                    # index of selected file
    #
    # Read directory for files and close directory:
    opendir($dh, $dir) || die "cannot open directory $dir: $!";
    @files = grep{/^.*txt$/} readdir($dh);
    closedir($dh);
    #
    # Sort and count the file names:
    @files = sort( @files );
    $numfiles = scalar(@files);
    #
    if ($numfiles > 0)
    {
        print "Select a file to process:\n";
        #
        # Read through the files you found and print them:
        for (local $i = 1; $i < $numfiles+1; $i++)
        {
            $file = $files[$i-1];
            print "$i $file\n";
        }
        #
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
sub writeout{
    # Get local sent variables:
    local($filename, $outdata) = @_;
    #
    # Try to open output file:
    open WFILE, ">$filename" or die $!;
    print WFILE $outdata;
    close WFILE;
}
