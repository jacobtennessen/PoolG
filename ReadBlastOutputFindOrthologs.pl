#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw($opt_a $opt_b $opt_o );

# Usage
my $usage = "
ReadBlastOutputFindOrthologs.pl

Copyright (C) 2022 by Jacob A Tennessen

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

Usage: perl ReadBlastOutputFindOrthologs.pl options
  -a  (path to) blast output file 1
  -b  (path to) blast output file 2
  -o  (path to) outfile
";

#############

getopts('a:b:o:');

die $usage unless ($opt_a);
die $usage unless ($opt_b);
die $usage unless ($opt_o);

my $blastfileA = $opt_a;
my $blastfileB = $opt_b;
my $outfile = $opt_o;

my $query;

my %scoresA1;

my %scoresA2;

my %eA1;

my %eA2;

my %matchesA;

my $recording = 0;

open(FIRST, $blastfileA) || die "can't open $blastfileA\n";

while (<FIRST>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  unless (length($line) > 0) {
    next;
  }
  if ($line =~ /^Query=/) {
    my @qdata = split /\s+/, $line;
    $query = $qdata[1];
  } elsif ($line =~ /Sequences producing significant alignments/) {
    $recording = 1;
  } elsif ($line =~ /Lambda/) {
    $recording = 0;
  } elsif ($recording == 1) {
    my @sdata = split /\s+/, $line;
    if ((defined $sdata[-1])&&(defined $sdata[-2])&&($sdata[-1] =~ /\d/)&&($sdata[-2] =~ /\d/)) {
      if (defined $scoresA1{$query}) {
        $scoresA2{$query} = $sdata[-2];
        $eA2{$query} = $sdata[-1];
      } else {
        $scoresA1{$query} = $sdata[-2];
        $eA1{$query} = $sdata[-1];
        $matchesA{$query} = $sdata[0];
      }
    }
  }
}

close (FIRST);

my $qacount = scalar(keys %matchesA);

print "Found $qacount genes for file A.\n";

my %scoresB1;

my %scoresB2;

my %eB1;

my %eB2;

my %matchesB;

$recording = 0;

open(NEXT, $blastfileB) || die "can't open $blastfileB\n";

while (<NEXT>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  unless (length($line) > 0) {
    next;
  }
  if ($line =~ /^Query=/) {
    my @qdata = split /\s+/, $line;
    $query = $qdata[1];
  } elsif ($line =~ /Sequences producing significant alignments/) {
    $recording = 1;
  } elsif ($line =~ /Lambda/) {
    $recording = 0;
  } elsif ($recording == 1) {
    my @sdata = split /\s+/, $line;
    if ((defined $sdata[-1])&&(defined $sdata[-2])&&($sdata[-1] =~ /\d/)&&($sdata[-2] =~ /\d/)) {
      if (defined $scoresB1{$query}) {
        $scoresB2{$query} = $sdata[-2];
        $eB2{$query} = $sdata[-1];
      } else {
        $scoresB1{$query} = $sdata[-2];
        $eB1{$query} = $sdata[-1];
        $matchesB{$query} = $sdata[0];
      }
    }
  }
}

close (NEXT);

my $qbcount = scalar(keys %matchesB);

print "Found $qbcount genes for file B.\n";

my @out;

foreach my $genea (keys %scoresA1) {
  unless ((defined $scoresA2{$genea})&&(($scoresA2{$genea} >= $scoresA1{$genea})||($eA2{$genea} < $eA1{$genea}))) {
    my $geneb = $matchesA{$genea};
    #print "Looking for $geneb as a match to $genea\n";
    if (defined $scoresB1{$geneb}) {
      unless ((defined $scoresB2{$geneb})&&(($scoresB2{$geneb} >= $scoresB1{$geneb})||($eB2{$geneb} < $eB1{$geneb}))) {
        if ($genea =~ /^$matchesB{$geneb}$/) {
          push @out, "$genea\t$geneb";
        }
      }
    }
  }
}

my $result = join "\n", @out;

unless ( open(META, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print META "GeneA\tGeneB\n$result";
close (META);

################

sub by_number {
    if ($a < $b) {-1} elsif ($a > $b) {1} else {0}
}

