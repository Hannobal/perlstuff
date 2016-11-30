#!/usr/bin/perl
use strict vars;
use dlpoly_utility;
use hanno_utility;

my($i,$len,@d,$ia,$ib);

if($#ARGV<2) {
  print " 1. name of input FIELD-file\n";
  print " 2. name of input mol2-file (bonds marked with \"ar\" or \"asdf\"\n";
  print " 3. name of output FIELD-file\n";
}

exit 1 if(read_field_file($ARGV[0],0)!=0);
exit 1 if(read_mol2_file($ARGV[1],0)!=0);

for($i=0;$i<$mol2_numbonds[0];$i++) {
  if($mol2_bonddata[0][$i][3] eq "ar" or $mol2_bonddata[0][$i][3] eq "asdf") {
#     print "$mol_bonddata[0][0][$i][4]\n";
    $ia = $mol2_bonddata[0][$i][0];
    $ib = $mol2_bonddata[0][$i][1];
    $d[0] = $mol2_atomdata[0][$ia][0]-$mol2_atomdata[0][$ib][0];
    $d[1] = $mol2_atomdata[0][$ia][1]-$mol2_atomdata[0][$ib][1];
    $d[2] = $mol2_atomdata[0][$ia][2]-$mol2_atomdata[0][$ib][2];
    $mol_bonddata[0][0][$i][4] = vector_length(@d);
  }
}

write_field_file("$ARGV[3]",0);