#!/usr/bin/perl
use feature state;
use hanno_utility;
use dlpoly_utility;
use lammps_utility;
use Math::Trig;
use Cwd;
use strict;

open(my $fh,"<",$ARGV[0]);
my $ftot=0;
while(<$fh>) {
  $_=~s/(^\s+|\s+$)//;
  next if(/^#/);
  my @data=split(/\s+/,$_);
  $ftot+=20000*$data[2]+250000;
}
$ftot/=1e6;
print "$ftot\n";
close($fh);