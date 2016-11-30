#!/usr/bin/perl

use Switch;
use dlpoly_utility;
use hanno_utility;

if($#ARGV<1) {
  print " 1. name of gaussian output file\n";
  print " 2. name of xyz file\n";
}

if(not open($fhg,"<",$ARGV[0])) {
  print "**** error: cannot open gaussian output file $ARGV[0]: $!\n";
  exit 1;
}

if(not open($fhout,">",$ARGV[1])) {
  print "**** error: cannot open output xyz file $ARGV[1]: $!\n";
  exit 1;
}

while($_=<$fhg>) {
  if(/Coordinates/ and /Angstrom/) {
    $_=<$fhg>;$_=<$fhg>;$_=<$fhg>;
    $natoms=0;
    @atoms=();
    until(/---------------/) {
      $_=~s/^\s+//;$_=~s/\s+$//;
      push(@atoms,[split(/\s+/,$_)]);
#       print join(" ",@{$atoms[$#atoms]}),"\n";
      $_=<$fhg>;
    }
    print $fhout $#atoms+1,"\n\n";
    for($i=0;$i<@atoms;$i++) {
      printf $fhout "%4u %15.6f %15.6f %15.6f\n",$atoms[$i][1],@{$atoms[$i]}[3..5];
    }
  }
}