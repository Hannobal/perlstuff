#!/usr/bin/perl

use hanno_utility;
use Math::Trig;
use Switch;

$countrow  = 1;
$lshowplot = 1;
$fitmin = -9e20;
$fitmax = 9e20;
$fitblock = 0;
$lfitmagnitude=0;
$fitfunc = 0; # Gaussian

if($#ARGV<0 or $ARGV[0]=~/^-h/) {
  print "input format:\n";
  print " 1. file to fit\n";
  print "[-b   <int>   block to fit (default=0)]\n";
  print "[-r   <int>   row with counts (min=1, default=".($countrow+1).")]\n";
  print "[-f   <str>   function to use (Gaussian/Rayleigh/Weibull)]\n";
  print "[-min <real>  minimum boundary for fit]\n";
  print "[-max <real>  maximum boundary for fit]\n";
  print "[-n   do not show plot after fit)]\n";
  print "[-m   also fit magnitude]\n";
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
    } case (/^-b/i) {
      $i++;
      $fitblock = $ARGV[$i];
    } case (/^-max/i) {
      $i++;
      $fitmax = $ARGV[$i];
    } case (/^-m/i) {
      $lfitmagnitude=1;
    } case (/^-f/i) {
      $i++;
      if($ARGV[$i]=~/^g/i) {
        $fitfunc=0;
      }elsif($ARGV[$i]=~/^r/i) {
        $fitfunc=1;
      }elsif($ARGV[$i]=~/^w/i) {
        $fitfunc=2;
      }else{
        die "**** error: unknown function '$ARGV[$i]'!\n";
      }
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

$block=0;
$newblock=0;
$av=0;
$sum=0;
while(<$fh>) {
  $_=~s/^\s+//; $_=~s/\s+$//;
  @linedata = split(/\s+/,$_);
  if(not @linedata) {
    $block++ unless($newblock);
    last if($block>$fitblock);
    $newblock=1;
    next;
  } else {
    $newblock=0;
  }
  next if($_=~/^#/);
  next if($block!=$fitblock);
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
  $av  += $linedata[0]*$linedata[$countrow];
  $sum += $linedata[$countrow];
  $numlines++;
}
if($numlines==0) {
  print "**** error: no valid lines found!\n"; exit 1;
}

if($isints) {
  ($av,$stdev) = calc_stdev(@data);
} else {
  $av/=$sum;
  $stdev=($maxx-$minx)/8.0
}

# if($lfitmagnitude) {
  print "before fit (estimate): $av +/- $stdev\n";
# }

$countrow++;
$test = `gnuplot -p -e 'm0=$av; s0=$stdev; filename="$filename"; fitblock=$fitblock; fitrow=$countrow; lshowplot=$lshowplot; lfitmagnitude=$lfitmagnitude; fitmin=$fitmin; fitmax=$fitmax; fitfunc=$fitfunc' /home/dietrich/scripts/perl/fit_gaussian.plt`;
print $test;