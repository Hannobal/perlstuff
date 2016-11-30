#!/usr/bin/perl

use hanno_utility;
use Math::Trig;
use Switch;

$countrow  = 2;
$lshowplot = 1;
$fitmin = -9e20;
$fitmax = 9e20;
$lfitmagnitude=0;

if($#ARGV<0 or $ARGV[0]=~/^-h/) {
  print "input format:\n";
  print " 1. file to fit\n";
  print "[-n   to not show plot after fit)]\n";
  print "[-r   row with counts (default:".($countrow+1).")]\n";
  print "[-m   to also fit magnitude]\n";
  print "[-min minimum boundary for fit]\n";
  print "[-max maximum boundary for fit]\n";
  exit 1;
}

$filename=$ARGV[0];

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-n/i) {
      $lshowplot = 0;
    } case (/^-r/i) {
      $i++;
      $countrow = $ARGV[$i]-1;
    } case (/^-min/i) {
      $i++;
      $fitmin = $ARGV[$i];
    } case (/^-max/i) {
      $i++;
      $fitmax = $ARGV[$i];
    } case (/^-m/i) {
      $lfitmagnitude=1;
#       $i++;
#       $magnitudeguess = $ARGV[$i];
    } else {
      die "**** error: unknown flag '$ARGV[$i]'!\n";
    }
  }
}

open($fh,"<",$filename) or die "**** error: can't open file $filename: $!\n";

$numlines = 0;
$maxval=-9e20;

#check if target row is integers or not
$isints=1;
while(<$fh>) {
  $_=~s/^\s+//;
  @linedata = split(/\s+/,$_);
  next if(not check_real($linedata[0]));
  if(not check_integer($linedata[$countrow])) {
    $isints=0; last;
  }
}
seek($fh,0,0);
$minx=9e20;
$maxx=-9e20;

while(<$fh>) {
  $_=~s/^\s+//; $_=~s/\s+$//;
  @linedata = split(/\s+/,$_);
  next if(not check_real($linedata[0]));
  next if($linedata[0]<$fitmin);
  next if($linedata[0]>$fitmax);
  if($isints) {
    for($i=0;$i<$linedata[$countrow];$i++) {
      push(@data,$linedata[0])
    }
  } else {
    $minx=min($minx,$linedata[0]);
    $maxx=max($maxx,$linedata[0]);
  }
  $maxval=max($maxval,$linedata[$countrow]);
  $numlines++;
}
if($numlines==0) {
  print "**** error: no valid lines found!\n"; exit 1;
}


if($isints) {
  ($av,$stdev) = calc_stdev(@data);
} else {
  $av=($maxx+$minx)/2.0;
  $stdev=($maxx-$minx)/8.0
}

# if($lfitmagnitude) {
  print "before fit (estimate): $av +/- $stdev\n";
# }

$test = `gnuplot -p -e 'm0=$av; s0=$stdev; filename="$filename"; lshowplot=$lshowplot; lfitmagnitude=$lfitmagnitude; fitmin=$fitmin; fitmax=$fitmax' /home/dietrich/scripts/perl/fit_gaussian.plt`;
print $test;