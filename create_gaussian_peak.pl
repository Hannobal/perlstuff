#!/usr/bin/perl

use Math::Random qw(random_normal);
use hanno_utility;

$beauty = 1;

if($#ARGV<3) {
  print "input format:\n";
  print " 1. name of output file\n";
  print " 2. center of distribution\n";
  print " 3. standard deviation\n";
  print " 4. resolution\n";
  print "[5. beauty-factor (edits number of values generated)]";
  exit 1;
}

$filename = $ARGV[0];
$m0  = $ARGV[1];
$s0  = $ARGV[2];
$res = $ARGV[3];
if($#ARGV>3) { $beauty = $ARGV[4]; }

$offset  = $m0+10*$s0;
$numvals = int($beauty*5000*$s0/$res);
print "$numvals $offset\n";
foreach $f (random_normal($numvals,$m0,$s0)) {
  histogram_add_one(0,$f+$offset,$res);
}

open($fhhist,">",$filename) or die "**** error: can't open output file: $!\n";
printf $fhhist "# %15s %17s %17s\n","value","norm. cnt","count";
&histogram_normalize_integral(0,$res);
for($i=$histminindex[0];$i<=$histmaxindex[0];$i++) {
  printf $fhhist "%17.10g %17.10g %17.10g\n",$res*$i-$offset,$histogram[0][$i],$histnum[0][$i];
}
close($fhhist);