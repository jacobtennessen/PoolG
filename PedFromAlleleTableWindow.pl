#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_f $opt_o $opt_l);

# Usage
my $usage = "
PedFromAlleleTableWindow.pl - converts tab-delimited R table to PED format
Copyright (C) 2024 by Jacob A Tennessen 

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

Usage: perl PedFromAlleleTableWindow.pl options
 required:
  -f  (path to) an allele table (numerical genotypes 0, 1, or 2)
  -o  (path to) output file baseline name
 optional:
  -l  comma-delimited list of samples to include (samples start at 0)
";

#############

# command line processing.
getopts('f:o:l:');
die $usage unless ($opt_f);
die $usage unless ($opt_o);

my ($genotypes,$outfile,@list);

$genotypes	= $opt_f if $opt_f;
$outfile = $opt_o if $opt_o;

if (defined $opt_l) {
  @list = split ",", $opt_l;
}

my @snps;

my %HoAgenos;

my @names;

my $snpcount = 0;

my $linecount = 0;

open(IN, "$genotypes") || die "can't open $genotypes\n";

while (<IN>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    my $pos = shift @data;
    my @realdata;
    if (defined $list[0]) {
      push @realdata, @data[@list];
    } else {
      push @realdata, @data;
    }
    if ($pos =~ /Site/) {
      foreach my $d (@realdata) {
        push @names, $d;
      }
      next;
    }
    $snpcount +=1;
    my @chromdata = split "_", $pos;
    push @snps, "$chromdata[0]\t$snpcount\t0\t$chromdata[1]";
    my $ind = 0;
    foreach my $g1 (@realdata) {
      if ($g1 =~ /0/) {
        push @{$HoAgenos{$names[$ind]}}, "1 1";
      } elsif ($g1 =~ /1/) {
        push @{$HoAgenos{$names[$ind]}}, "1 2";
      } elsif ($g1 =~ /2/) {
        push @{$HoAgenos{$names[$ind]}}, "2 2";
      } else {
        push @{$HoAgenos{$names[$ind]}}, "0 0";
      }
      $ind +=1;
    }
    $linecount +=1;
}

close (IN);

my $map = join "\n", @snps;

@snps = ();

unless ( open(MOUT, ">$outfile.map") ) {
    print "Cannot open file \"$outfile.map\" to write to!!\n\n";
    exit;
}
print MOUT $map;
close (MOUT);

$map = 1;

my @ped;

my $famcount = 1;

foreach my $n (@names) {
  my $genos = join " ", @{$HoAgenos{$n}};
  @{$HoAgenos{$n}} = ();
  push @ped, "$famcount $n 0 0 0 1 $genos";
  $famcount +=1;
}

my $ped = join "\n", @ped;

@ped = ();

unless ( open(POUT, ">$outfile.ped") ) {
    print "Cannot open file \"$outfile.ped\" to write to!!\n\n";
    exit;
}
print POUT "$ped\n";
close (POUT);

$ped = 1;

#########################
