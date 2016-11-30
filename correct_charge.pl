#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;

if($#ARGV<2) {
  print "input format:\n";
  print "1. input FIELD filename\n";
  print "2. mol2 filename:\n";
  print "3. output FIELD filename:\n";
  exit;
}

$infld  = $ARGV[0];
$inmol2 = $ARGV[1];
$outfld = $ARGV[2];

exit 1 if(read_field_file($infld,0)!=0);
exit 1 if(read_mol2_file($inmol2,0) != 0);

$foundmolname=-1;
for($t=0;$t<@field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /$mol2_name[0]/) {
    $foundmolname = $t;
  }
}
if($foundmolname<0) {
  print "**** error: could not find molecule named $mol2_name[0]\n";
  exit 1;
}

if($mol2_numatoms[0] != $mol_numatoms[0][$foundmolname]) {
  print "**** error: number of atoms in mol2-file($mol2_numatoms[0]) does not fit ",
  "number of atoms in FIELD file ($mol_numatoms[0][$foundmolname]\n";
}

for($a=0;$a<$mol_numatoms[0][$foundmolname];$a++) {
  $mol_atomdata[0][$foundmolname][$a][2] = $mol2_atomdata[0][$a][8];
}
$mol_charge[0][$foundmolname] = $mol2_charge[0];

$totcharge=0;
for($t=0;$t<$field_nummols[0];$t++) {
  $totcharge += $mol_numents[0][$t]*$mol_charge[0][$t];
}

print "new total charge: $totcharge\n";

write_field_file($outfld,0);
