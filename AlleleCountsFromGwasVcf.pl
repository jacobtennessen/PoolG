#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_v $opt_o $opt_d $opt_l $opt_c $opt_z $opt_i $opt_q $opt_s);

# Usage
my $usage = "
AlleleCountsFromGwasVcf.pl
Converts a vcf into allele counts for R

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

Usage: perl AlleleCountsFromGwasVcf.pl options
 required:
  -v	(path to) vcf file
  -o  (path to) output file
 optional:
  -d  minimum allele count (across all samples) [default = 10]
  -l  comma-delimited list of sample positions to use (0-based, first sample column = 0)
  -c  comma-delimited list of chromosomes to use
  -z  file is zipped
  -i  exclude indels
  -q  minimum value for RPB, MQB, MQSB, BQB [default = none]
  -s  minimum size [default = none]
";

#############

getopts('v:o:d:l:c:q:s:zi');
die $usage unless ($opt_v);
die $usage unless ($opt_o);

my ($genotypes, $outfile, $mincount, @samplelist, @chromlist, $minmqb, $excludeindel, $minqual, $minsize);

$genotypes = $opt_v if $opt_v;

$outfile = $opt_o if $opt_o;

if (defined $opt_d) {
  $mincount = $opt_d;
} else {
  $mincount = 10;
}

if (defined $opt_q) {
  $minqual = $opt_q;
} else {
  $minqual = 0;
}

if (defined $opt_s) {
  $minsize = $opt_s;
} else {
  $minsize = 0;
}

if (defined $opt_l) {
  @samplelist = split ",", $opt_l;
}

if (defined $opt_c) {
  @chromlist = split ",", $opt_c;
}

if (defined $opt_i) {
  $excludeindel = 1;
}

if (defined $opt_z) {
  open(IN, "gunzip -c $genotypes |") || die "can’t open pipe to $genotypes";
} else {
  open(IN, $genotypes) || die "can't open $genotypes";
}

my %chroms;

foreach my $c (@chromlist) {
  $chroms{$c} = 1;
}

my $outputsize = 1000;

my @allnames;

my @allsamples;

my $linecount = 0;

my @out;

my $cleared;

while (<IN>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split "\t", $line;
  if ($data[0] =~ /CHROM/) {
    my $samplecount = 0;
    my @allsamplelist;
    for (my $fnx = 9; $fnx < (scalar(@data)); $fnx ++) {
        my @namedata = split /\//, $data[$fnx];
        my @namedata2 = split /\./, $namedata[-1];
        my $shortname = $namedata2[0];
        $shortname =~ s/-/_/g;
        push @allnames, $shortname;
        push @allsamples, $fnx;
        push @allsamplelist, $samplecount;
        $samplecount +=1;
    }
    unless (defined $samplelist[0]) {
      push @samplelist, @allsamplelist;
    }
    my @namelist;
    foreach my $n (@allnames[@samplelist]) {
      push @namelist, "Ref$n\tAlt$n";
    }
    my $namelist = join "\t", @namelist;
    push @out, "Site\t$namelist";
    next;
  } elsif ($line =~ /^#/) {
      next;
  }
  if (defined $chromlist[0]) {
    unless (defined $chroms{$data[0]}) {
      next;
    }
  }
  if ($data[4] =~ /\./) {
    next;
  }
  if ((defined $excludeindel)&&($data[7] =~ /INDEL/)) {
    next;
  }
  if ($minqual > 0) {
    my @rpb1 = split "RPB=", $data[7];
    unless (defined $rpb1[1]) {
      next;
    }
    my @rpb2 = split ";", $rpb1[1];
    my @mqb1 = split "MQB=", $data[7];
    unless (defined $mqb1[1]) {
      next;
    }
    my @mqb2 = split ";", $mqb1[1];
    my @mqsb1 = split "MQSB=", $data[7];
    unless (defined $mqsb1[1]) {
      next;
    }
    my @mqsb2 = split ";", $mqsb1[1];
    my @bqb1 = split "BQB=", $data[7];
    unless (defined $bqb1[1]) {
      next;
    }
    my @bqb2 = split ";", $bqb1[1];
    if (($rpb2[0] < $minqual)||($mqb2[0] < $minqual)||($mqsb2[0] < $minqual)||($bqb2[0] < $minqual)) {
      next;
    }
  }
  my @info = split ":", $data[8];
  my $adpos;
  for (my $in = 0; $in < (scalar(@info)); $in++) {
      if ($info[$in] =~ /AD/) {
          $adpos = $in;
      }
  }
  my @inds = @data[@allsamples[@samplelist]];
  my @ads;
  push @ads, "$data[0]_$data[1]";
  my $refcount = 0;
  my $altcount = 0;
  my $sizeflag = 0;
  foreach my $i (@inds) {
    my @idata = split ":", $i;
    my @ad = split ",", $idata[$adpos];
    my $adref = 0;
    my $size = 0;
    for (my $a = 0; $a < (scalar(@ad)); $a++) {
      unless ($a == 1) {
        if ($ad[$a] =~ /\d/) {
          $adref += $ad[$a];
        }
      }
      if ($ad[$a] =~ /\d/) {
        $size += $ad[$a];
      }
    }
    unless ((defined $ad[1])&&($ad[1] =~ /\d/)) {
      $ad[1] = 0;
    }
    push @ads, "$adref\t$ad[1]";
    $refcount += $adref;
    $altcount += $ad[1];
    if ($size < $minsize) {
      $sizeflag = 1;
    }
  }
  if (($refcount >= $mincount)&&($altcount >= $mincount)&&($sizeflag == 0)) {
    my $adlist = join "\t", @ads;
    push @out, $adlist;
  }
  if ((scalar (@out)) >= $outputsize) {
      my $result = join "\n", @out;
      if (defined $cleared) {
          unless ( open(META, ">>$outfile") ) {
              print "Cannot open file \"$outfile\" to write to!!\n\n";
              exit;
          }
          print META "$result\n";
          close (META);
      } else {
          unless ( open(META, ">$outfile") ) {
              print "Cannot open file \"$outfile\" to write to!!\n\n";
              exit;
          }
          print META "$result\n";
          close (META);
          $cleared = 1;
      }
      @out = ();
  }
  $linecount +=1;
}

close (IN);

if ((scalar(@out)) > 0) {
    my $result = join "\n", @out;
    if (defined $cleared) {
        unless ( open(META, ">>$outfile") ) {
            print "Cannot open file \"$outfile\" to write to!!\n\n";
            exit;
        }
        print META "$result\n";
        close (META);
    } else {
        unless ( open(META, ">$outfile") ) {
            print "Cannot open file \"$outfile\" to write to!!\n\n";
            exit;
        }
        print META "$result\n";
        close (META);
        $cleared = 1;
    }
    @out = ();
}
