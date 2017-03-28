#!/usr/bin/perl

use Switch;
use dlpoly_utility;
use hanno_utility;
use lammps_utility;

$startindex=1;
$separator=":";

if($#ARGV<0) {
  print "name of input FIELD-file\n";
  print "  -s <int> start with index 0 (default: $startindex)\n";
  print "  -p <str> separator (default: '$separator')\n";
  exit 1;
}

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-s$/) {
      $i++;
      $startindex=$ARGV[$i];
    } else {
      print "**** error: expected optional flag but got \"$ARGV[$i]\"!\n";
      exit 1;
    }
  }
}

exit 1 if(read_field_file($ARGV[0],0) != 0);

$curr=$startindex;
for($t=0;$t<$field_nummols[0];$t++) {
  printf "%-15s",$mol_name[0][$t];
  if($mol_numents[0][$t]==0) {
    print " no entities\n";
    next;
  }
  printf "%7u",$curr;
  $curr+=$mol_numents[0][$t]*$mol_numatoms[0][$t];
  print $separator,$curr-1,"\n";
}