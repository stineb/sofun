#!/usr/bin/perl
#
# diffdays.pl
#
# written by Tyler W. Davis
# Imperial College London
#
# 2012-09-03 -- created
# 2015-11-25 -- last modified
#
# HOW TO RUN:
# > perl diffdays.pl
# > <ENTER START DATE>
# > <ENTER END DATE>
# > <RUN AGAIN?>
#
# ------------
# description:
# ------------
# This script calculates the number of days between two given dates.

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# LOAD MODULES:
# ////////////////////////////////////////////////////////////////////////////
use warnings;
use Date::Calc qw(check_date Delta_Days);

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# MAIN
# ////////////////////////////////////////////////////////////////////////////
my @values = ();
my $ans = "";
my $date = "";

print "Please use the following format: \n";
print "  Today's date is: ";
my $td = (localtime)[3];
my $tm = (localtime)[4] + 1;
my $ty = (localtime)[5] + 1900;
print "$ty/$tm/$td \n";

# Allow user to check their input for correctness:
$ans = "Y";
while ( lc($ans) eq "y") {
    print "Enter start date: ";
    $date=<>;
    chomp $date;
    @values = split('/',$date);
    $check1 = check_date($values[0], $values[1], $values[2]);
    @date1 = ($values[0], $values[1], $values[2]);

    print "Enter end date: ";
    $date=<>;
    chomp $date;
    @values = split('/',$date);
    $check2 = check_date($values[0], $values[1], $values[2]);
    @date2 = ($values[0], $values[1], $values[2]);

    # Check the date (note this does not catch letters in the year):
    if ( $check1 & $check2 ) {
        $D = Delta_Days(@date1, @date2);
        print "Delta days: $D\n";

        print "Run again (Y/n)? ";
        $ans = <>;
        chomp $ans;
    }
    else {
        print " ... those dates were invalid, please try again.\n";
    }
}
