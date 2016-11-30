#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;

if($#ARGV<1) {
  print " 1. name of input CONFIG-file\n";
  print " 2. name of input FIELD-file\n";
  print " [molecule name]\n";
  print " [molecule ids (from 0, e.g. '0 3-12 14')]\n";
  exit;
}

$incfg    = $ARGV[0];
$infld    = $ARGV[1];
$molname  = $ARGV[2];
$molidstr = join(" ",);
@molids   = parse_intlist(0,@ARGV[3..$#ARGV]);
$kB = 0.831451115;

exit 1 if(read_field_file($infld,0) != 0);
exit 1 if(read_config_file($incfg,0,0) != 0);

$ekintot     = 0;
$numatomstot = 0;
for($t=0;$t<$field_nummols[0];$t++) {
  $ekin[$t]     = 0;
  $numatoms[$t] = 0;
  if(defined($molname)) {
    next if(not $mol_name[0][$t]=~/$molname/i);
  }
  next if($mol_numents[0][$t]==0);
  if(@molids) {
    foreach $m (@molids) {
      if($m>=$mol_numents[0][$t]) {
	print "**** error: molecule $mol_name[0][$t] has only $mol_numents[0][$t] entities!\n";
	exit 1;
      }
      for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	next if($mol_atomdata[0][$t][$a][4]>0);
	$vxsq = $cdata[0][$t][$m][$a][3]*$cdata[0][$t][$m][$a][3];
	$vysq = $cdata[0][$t][$m][$a][4]*$cdata[0][$t][$m][$a][4];
	$vzsq = $cdata[0][$t][$m][$a][5]*$cdata[0][$t][$m][$a][5];
	$ekin[$t] += $mol_atomdata[0][$t][$a][1]*($vxsq+$vysq+$vzsq);
	$ekintot  += $mol_atomdata[0][$t][$a][1]*($vxsq+$vysq+$vzsq);
	$numatomstot++;
	$numatoms[$t]++;
      }
    }
  } else {
    for($m=0;$m<$mol_numents[0][$t];$m++) {
      for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	next if($mol_atomdata[0][$t][$a][4]>0);
	$vxsq = $cdata[0][$t][$m][$a][3]*$cdata[0][$t][$m][$a][3];
	$vysq = $cdata[0][$t][$m][$a][4]*$cdata[0][$t][$m][$a][4];
	$vzsq = $cdata[0][$t][$m][$a][5]*$cdata[0][$t][$m][$a][5];
	$ekin[$t] += $mol_atomdata[0][$t][$a][1]*($vxsq+$vysq+$vzsq);
	$ekintot  += $mol_atomdata[0][$t][$a][1]*($vxsq+$vysq+$vzsq);
	$numatomstot++;
	$numatoms[$t]++;
      }
    }
  }
  printf "%-20s %10.5f\n",$mol_name[0][$t],$ekin[$t]/(3*$numatoms[$t]-6)/$kB,"\n";
}
printf "%-20s %10.5f\n", "total", $ekintot/(3*$numatomstot-6)/$kB,"\n";