#!/usr/bin/perl

use dlpoly_utility;

if($#ARGV<2) {
  print " 1. name of original CONFIG-file\n";
  print " 2. name of second CONFIG-file\n";
  print " 3. name of FIELD-file\n";
  print "[4. molecule name]\n";
  print "[5. index of molecule]\n";
  exit;
}

exit 1 if(read_field_file($ARGV[2],0) != 0);
exit 1 if(read_config_file($ARGV[0],0,0)!=0);
exit 1 if(read_config_file($ARGV[1],0,1)!=0);

$molname = $ARGV[3];
$mi      = $ARGV[4];

if(defined($molname)) {
  $ti=-1;
  for($t=0;$t<$field_nummols[0];$t++) {
    if($mol_name[0][$t] eq $molname) {
      $ti=$t;
      last;
    }
  }
  if($ti<0) {
    print "**** error: no molecule named $molname was found!\n";
    exit 1;
  }
}

print calc_rmsd(0,1,$ti,$mi),"\n";