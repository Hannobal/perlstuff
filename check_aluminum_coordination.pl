#!/usr/bin/perl

use Switch;
use dlpoly_utility;
use hanno_utility;

if($#ARGV<1) {
  print "input format:\n";
  print " 1. name of CONFIG file\n";
  print " 2. name of corresponding FIELD file\n";
  print "[3. cutoff for bond length (default: 2.85)]\n";
  exit 99;
}

$cutoff = 2.85;
$cfgfilename = $ARGV[0];
$fldfilename = $ARGV[1];
$cutoff      = $ARGV[2] if($#ARGV>1);
if($cutoff <= 0) {
  print "**** error: cutoff must be greater than zero!\n";
  exit 99;
}
$cutoffsq=$cutoff*$cutoff;

exit 99 if(read_field_file($fldfilename,0) != 0);
exit 99 if(read_config_file($cfgfilename,0,0) != 0);

$numnonhex   = 0;
@nonhexa     = ();

#determine surface
$surfzmax=-9e20;
for($t=0;$t<$field_nummols[0];$t++) {
  next if(not $mol_name[0][$t] =~ /alumini?um/i);
  for($m=0;$m<$mol_numents[0][$t];$m++) {
    $surfzmax = $cdata[0][$t][$m][0][2] if($surfzmax<$cdata[0][$t][$m][0][2]);
  }
}
$surfzmax += 2;

for($t=0;$t<$field_nummols[0];$t++) {
  next if($mol_name[0][$t] !~ /alumini?um/i or $mol_name[0][$t] !~ /free/i);
  loopal:for($m=0;$m<$mol_numents[0][$t];$m++) {
    $coord=0;
    for($t2=0;$t2<$field_nummols[0];$t2++) {
      for($a2=0;$a2<$mol_numatoms[0][$t2];$a2++) {
	next if($mol_atomdata[0][$t2][$a2][0] ne 'O'
	    and $mol_atomdata[0][$t2][$a2][0] ne 'OX'
	    and $mol_atomdata[0][$t2][$a2][0] ne 'OA');
	for($m2=0;$m2<$mol_numents[0][$t2];$m2++) {
	  next if($cdata[0][$t2][$m2][$a2][2]>$surfzmax);
	  for($x=-1;$x<2;$x++) {
	    for($y=-1;$y<2;$y++) {
	      $dx = $cdata[0][$t2][$m2][$a2][0]-$cdata[0][$t][$m][0][0]+$x*$cell[0][0][0]+$y*$cell[0][1][0];
	      $dy = $cdata[0][$t2][$m2][$a2][1]-$cdata[0][$t][$m][0][1]+$x*$cell[0][0][1]+$y*$cell[0][1][1];
	      $dz = $cdata[0][$t2][$m2][$a2][2]-$cdata[0][$t][$m][0][2];
	      $distsq=$dx*$dx+$dy*$dy+$dz*$dz;
	      next if($distsq>$cutoffsq);
	      $coord++;
	      next loopal if($coord==6);
	    }
	  }
	} # end for a2
      } # end for m2
    } # end for t2
    if($coord != 6) {
      $numnonhex++;
      print "**** warning: aluminum t=$t m=$m i=$cdata[0][$t][$m][0][11] is $coord-fold coordinated!\n";
      push(@nonhexa,[$t,$m]);
    }
  } # end for m
} # end for t

exit $numnonhex;
