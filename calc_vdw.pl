#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;

if($#ARGV<1) {
  print "the input format is:\n";
  print "1. gaff.dat\n";
  print "[ frcmod file (must contain \"frcmod\" in its name!)]\n";
  print "2. list of atom types in lower case letters\n";
  exit;
}

exit 1 if(read_gaff_file($ARGV[0]) != 0);

foreach $i ( 1..$#ARGV ) {
  if($ARGV[$i]=~/frcmod/) {
    exit 1 if(read_frcmod_file($ARGV[$i]) != 0);
  } else {
    $type[@type] = uc($ARGV[$i]);
  }
}

foreach $i ( 0..$#type ) {
  if(not exists $gaff_vdwparam{$type[$i]}[0]) {
    print "**** error: did not find entry for ".$type[$i]."\n";
    exit;
  }
}

$numvdw= @type * (@type + 1.0) / 2.0;

print "VDW $numvdw\n";
for($i=0;$i<@type;$i++) {
  for($j=$i;$j<@type;$j++) {
    ($epsilon, $sigma) = calc_vdw_params($type[$i], $type[$j]);
    printf "%-7s %-7s %-13s", uc($type[$i]), uc($type[$j]), "lj";
    printf "%11.6g %11.6g\n", $epsilon, $sigma;
  }
}
