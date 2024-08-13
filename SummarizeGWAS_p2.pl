#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_f $opt_o $opt_w $opt_m $opt_d );

# Usage
my $usage = "
SummarizeGWAS_p2.pl
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

Usage: perl SummarizeGWAS_p2.pl options
 required:
  -f	  comma-delimted list of (paths to) GWAS output files (product of R analysis)
  -o    (path to) output file
 optional:
  -w    window size [default = 10000]
  -m    min count per sample (ref + alt) [default = 1]
  -d    min dist [default =1000]

";

#############

# command line processing.
getopts('f:o:w:m:d:');
die $usage unless ($opt_f);
die $usage unless ($opt_o);

my (@freqfiles, $outfile, $window, $mincount, $mindist);

@freqfiles = split ",", $opt_f if $opt_f;
$outfile = $opt_o if $opt_o;

if (defined $opt_w) {
  $window = $opt_w;
} else {
  $window = 10000;
}

if (defined $opt_m) {
  $mincount = $opt_m;
} else {
  $mincount= 1;
}

if (defined $opt_d) {
  $mindist = $opt_d;
} else {
  $mindist= 1000;
}

my %HoHoApf;

my %HoHoApch;

my %HoHoApos;

foreach my $freqfile (@freqfiles) {

  open(IN, "$freqfile") || die "can't open $freqfile\n";

  while (<IN>) {
      my $line = $_;
      $line =~ s/\r|\n//g;
      my @data = split /\t/, $line;
      if ($data[0] =~ /\d/) {
        if (($data[2] >= $mincount)&&($data[3] >= $mincount)) {
          my @sitedata = split "_", $data[0];
          my $site = pop @sitedata;
          my $chrom = join "_", @sitedata;
          my $w = (int($site/$window))*$window;
          my $wstagup = (int(($site-($window*0.25))/$window))*$window + $window*0.25;
          my $wstag = (int(($site-($window*0.5))/$window))*$window + $window*0.5;
          my $wstagdown = (int(($site-($window*0.75))/$window))*$window + $window*0.75;
          push @{$HoHoApf{$chrom}{$w}}, $data[5];
          push @{$HoHoApch{$chrom}{$w}}, $data[6];
          push @{$HoHoApos{$chrom}{$w}}, $site;
          push @{$HoHoApf{$chrom}{$wstagup}}, $data[5];
          push @{$HoHoApch{$chrom}{$wstagup}}, $data[6];
          push @{$HoHoApos{$chrom}{$wstagup}}, $site;
          push @{$HoHoApf{$chrom}{$wstag}}, $data[5];
          push @{$HoHoApch{$chrom}{$wstag}}, $data[6];
          push @{$HoHoApos{$chrom}{$wstag}}, $site;
          push @{$HoHoApf{$chrom}{$wstagdown}}, $data[5];
          push @{$HoHoApch{$chrom}{$wstagdown}}, $data[6];
          push @{$HoHoApos{$chrom}{$wstagdown}}, $site;
        }
      }
  }

  close (IN);

}

my @out;

foreach my $chrom (sort keys %HoHoApos) {
  foreach my $win (sort by_number keys %{ $HoHoApos{$chrom} } ) {
    my $pcount = scalar(@{$HoHoApos{$chrom}{$win}});
    my $bestflower = 1;
    my $bestchlower = 1;
    my $bestfmean = 1;
    my $bestchmean = 1;
    my $pftotal = 0;
    my $pchtotal = 0;
    for (my $p1 = 0; $p1 < $pcount; $p1 ++) {
      $pftotal += log($HoHoApf{$chrom}{$win}[$p1]);
      $pchtotal += log($HoHoApch{$chrom}{$win}[$p1]);
      if (($HoHoApf{$chrom}{$win}[$p1] < $bestflower)||($HoHoApch{$chrom}{$win}[$p1] < $bestchlower)) {
        for (my $p2 = 0; $p2 < $pcount; $p2 ++) {
          if (($HoHoApf{$chrom}{$win}[$p2] < $bestflower)||($HoHoApch{$chrom}{$win}[$p2] < $bestchlower)) {
            if (abs($HoHoApos{$chrom}{$win}[$p1] - $HoHoApos{$chrom}{$win}[$p2]) >= $mindist) {
              my $tempflow = $HoHoApf{$chrom}{$win}[$p1];
              my $tempfhigh = $HoHoApf{$chrom}{$win}[$p2];
              if ($HoHoApf{$chrom}{$win}[$p2] > $tempflow) {
                $tempflow = $HoHoApf{$chrom}{$win}[$p2];
                $tempfhigh = $HoHoApf{$chrom}{$win}[$p1];
              }
              if ($tempflow < $bestflower) {
                $bestflower = $tempflow;
                $bestfmean = exp((log($tempflow)+log($tempfhigh))/2);
              }
              my $tempchlow = $HoHoApch{$chrom}{$win}[$p1];
              my $tempchhigh = $HoHoApch{$chrom}{$win}[$p2];
              if ($HoHoApch{$chrom}{$win}[$p2] > $tempchlow) {
                $tempchlow = $HoHoApch{$chrom}{$win}[$p2];
                $tempchhigh = $HoHoApch{$chrom}{$win}[$p1];
              }
              if ($tempchlow < $bestchlower) {
                $bestchlower = $tempchlow;
                $bestchmean = exp((log($tempchlow)+log($tempchhigh))/2);
              }
            }
          }
        }
      }
    }
    my $meanf = 0;
    my $meanch = 0;
    if ($pcount > 0) {
      $meanf = $pftotal/$pcount;
      $meanch = $pchtotal/$pcount;
    }
    push @out, "$chrom\t$win\t$bestflower\t$bestfmean\t$meanf\t$bestchlower\t$bestchmean\t$meanch\t$pcount";
  }
}

my $result = join "\n", @out;
unless ( open(OUT, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print OUT "Contig\tWindow\tSecondBestFisher\tFisherMean2\tLogFisherMeanAll\tSecondBestChi\tChiMean2\tLogChiMeanAll\tSNPs\n$result";
close (OUT);

###################

sub by_number {
    if ($a < $b) {-1} elsif ($a > $b) {1} else {0}
}
