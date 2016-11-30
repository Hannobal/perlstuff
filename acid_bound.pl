#!/usr/bin/perl

use Switch;
use dlpoly_utility;
use hanno_utility;

# determines whether a PA is bound to the surface

if($#ARGV<1) {
  print "input format:\n";
  print " 1. name of CONFIG/HISTORY file\n";
  print " 2. name of corresponding FIELD file\n";
  print "[3. maximum Al-O bond length]\n";
  exit 2;
}

$cfgfilename = $ARGV[0];
$fldfilename = $ARGV[1];
if($#ARGV>1) {
  $dmaxsq = $ARGV[2]**2;
} else {
  $dmaxsq = 2.3**2;
}
# if($dmaxsq<2) {
#   print "**** error: $dmaxsq is not a realistic value for Al-O bond length...\n";
#   exit 1;
# }


exit 1 if(read_field_file($fldfilename,0)    != 0);
exit 1 if(read_config_file($cfgfilename,0,0) != 0);

$alfree=-1;
$alfrozen=-1;
@aluminium=();
for($t1=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /aluminum/i) {
    push(@aluminium,splice(@{$cdata[0][$t]},0,$mol_numents[0][$t]));
  }
}

$unbound=0;
for($t=0;$t<$field_nummols[0];$t++) {
  next if(not $mol_name[0][$t] =~ /-PA-/);
  loopmol:for($m=0;$m<$mol_numents[0][$t];$m++) {
    $bound=0;
    loopatom1:for($a=0;$a<$mol_numatoms[0][$t];$a++) {
      next if($mol_atomdata[0][$t][$a][0] ne 'O');
      if($config_key[0]>0) {
	for($x=-1;$x<2;$x++) {
	  for($y=-1;$y<2;$y++) {
	    $shiftx = $x*$cell[0][0][0]+$y*$cell[0][0][1];
	    $shifty = $x*$cell[0][1][0]+$y*$cell[0][1][1];
	    for($b=0;$b<@aluminium;$b++) {
	      $dx = $cdata[0][$t][$m][$a][0] - $aluminium[$b][0][0] + $shiftx;
	      $dy = $cdata[0][$t][$m][$a][1] - $aluminium[$b][0][1] + $shifty;
	      $dz = $cdata[0][$t][$m][$a][2] - $aluminium[$b][0][2];
	      $distsq = $dx*$dx + $dy *$dy + $dz*$dz;
	      if($distsq < $dmaxsq) {
		$bound = 1;
		next loopmol;
	      }
	    }
	  }
	}
      } # end if config_key
    } # end for $a
    if(not $bound) {
      print "**** $mol_name[0][$t] number $m is not bound!\n";
      $unbound++;
    }
  } #end for $m
} # end for $t

exit $unbound;
