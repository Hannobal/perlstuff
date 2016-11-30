#!/usr/bin/perl

use Switch;
use dlpoly_utility;
use hanno_utility;

if($#ARGV<1) {
  print "input format:\n";
  print "1. name of CONFIG file\n";
  print "2. name of corresponding FIELD file\n";
  exit;
}

$cfgfilename = $ARGV[0];
$fldfilename = $ARGV[1];

exit 1 if(read_field_file($fldfilename,0) != 0);
exit 1 if(read_config_file($cfgfilename,0,0) != 0);

$area=($cell[0][0][0]*$cell[0][1][1]+$cell[0][0][1]*$cell[0][1][0])/100;
print "area: $area square nanometers\n";
printf "%-15s %5s %15s\n","molecule", "#", "density";
$numpas = 0;
for($t=0;$t<$field_nummols[0];$t++) {
  next if($mol_name[0][$t] !~ /-PA/ and $mol_name[0][$t] !~ /-CA/);
  $numpas += $mol_numents[0][$t];
  $dens    = $mol_numents[0][$t]/$area;
  printf "%-15s %5u %15.8f\n", $mol_name[0][$t], $mol_numents[0][$t], $dens;
}
$dens = $numpas/$area;
print "-------------------------------------------------------------\n";
printf "%-15s %5u %15.8f\n", "TOTAL", $numpas, $dens;
