#!/usr/bin/perl

use Math::Trig;
use dlpoly_utility;
use hanno_utility;
use POSIX qw(ceil floor);
# use strict;

if($#ARGV<1) {
  print "please specify input and ouptut file!\n";
  exit;
}

open($fhin,"<",$ARGV[0]) or die "**** error: can't open input file: $!\n";

$l=0;
while(<$fhin>) {
  $lines[$l]=$_;
  $_ =~ s/^\s+//;$_ =~ s/\s+$//;
  next if($_=~/^#/);
  @linedata = split(/\s+/,$_);
  if(defined($xold)) {
    $res=$linedata[0]-$xold;
    if(defined($resold)) {
#       print "$res vs $resold\n";
#       die "**** error: resolution chanes!\n" if($resold!=$res);
    }
    $resold=$res;
  }
  $xold=$linedata[0];
  for($i=1;$i<@linedata;$i++) {
    $factor[$i]+=$linedata[$i];
  }
  @{$data[$l]} = @linedata;
  $l++;
}

close($fhin);

for($i=1;$i<@factor;$i++) {
  $factor[$i]*=$res;
  print "$factor[$i]\n";
}

open($fhout,">",$ARGV[1]) or die "**** error: can't open output file: $!\n";
for($l=0;$l<@lines;$l++) {
  if($_=~/^\s*#/) {
    print $fhout $lines[$l];
  } else {
    print $fhout $data[$l][0];
    for($i=1;$i<@{$data[$l]};$i++) {
      print $fhout "  ",$data[$l][$i]/$factor[$i];
    }
    print $fhout "\n";
  }
}
close($fhout);