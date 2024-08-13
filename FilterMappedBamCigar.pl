#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_m );

# Usage
my $usage = "
FilterMappedBamCigar.pl - filters out reads with too many matching bp based on CIGAR string
Copyright (C) 2024 by Jacob A Tennessn 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: perl FilterMappedBamCigar.pl options
 required:
  -m  minimum number of matched bp
";

#############

# command line processing.
getopts('m:');
die $usage unless ($opt_m);

my $matchthresh = $opt_m if $opt_m;

while (<STDIN>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    my $mcount = 0;
    if (defined $data[10]) {
      if ((length($data[5]) > 1)&&($data[5] =~ /M/)) {
        my @matchdata = split "M", $data[5];
        foreach my $m (@matchdata) {
          if ($m =~ /\d$/) {
            my @matchdata2 = split /\D/, $m;
            $mcount += $matchdata2[-1];
          }
        }
      }
    }
    if ($mcount < $matchthresh) {
      print "$line\n";
    }
}
