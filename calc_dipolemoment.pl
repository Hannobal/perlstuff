#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;

if($#ARGV<1) {
  print "input format: <CONFIG filename> <FIELD filename>\n";
  exit;
}

$cfgname = $ARGV[0];
$fldname = $ARGV[1];
$test=&read_field_file($fldname,0);
if($test>0) { exit 1; }
if(read_config_file($cfgname,0,0)>0) { exit 1; }

@dipole = (0,0,0);
for($t=0;$t<$field_nummols[$fi];$t++) {
  for($m=0;$m<$mol_numents[$fi][$t];$m++) {
    for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
      $dipole[0] += $cdata[0][$t][$m][$a][0]*$mol_atomdata[0][$t][$a][2];
      $dipole[1] += $cdata[0][$t][$m][$a][1]*$mol_atomdata[0][$t][$a][2];
      $dipole[2] += $cdata[0][$t][$m][$a][2]*$mol_atomdata[0][$t][$a][2];
    }
  }
}
print join(" ",@dipole),"\n";
