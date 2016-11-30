#!/usr/bin/perl

use Storable qw(dclone);
use dlpoly_utility;
use hanno_utility;

if($#ARGV<2) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file\n";
  print "2. name of corresponding FIELD file\n";
  print "3. minimum distance\n";
  print "[-n to neglect periodic conditions]\n";
  exit;
}

$incfg = $ARGV[0];
$infld = $ARGV[1];
$dmin  = $ARGV[2];
if(not check_real($dmin)) {
  print "**** error: $dmin is not a number\n";
  exit 1;
}
$dminsq = $dmin**2;

$lperiodic=1;
if($#ARGV>2) {
  if($ARGV[3] =~ /^-[nN]/) {
    $lperiodic=0;
  } else {
    print "**** error: could not interpret input $ARGV[2]\n";
  }
}

exit 1 if(read_field_file($infld,0) != 0);
exit 1 if(read_config_file($incfg,0,0) != 0);

@cdata2=();
for($t=0;$t<$field_nummols[0];$t++) {
  push(@cdata2, @{dclone \@{$cdata[0][$t]}});
}

$mindistsq[0]=9999999;
$periodic_key[0]=0 if(not $lperiodic);
if($periodic_key[0] != 0) {
  if($periodic_key[0] == 6) {$zmax=1} else {$zmax=2};
  for($x=0;$x<2;$x++){
    for($y=0;$y<2;$y++){
      for($z=0;$z<$zmax;$z++){
# 	$i=0;
	print "$x $y $z\n";
	$pervec[0] = $x*$cell[0][0][0]+$y*$cell[0][1][0]+$z*$cell[0][2][0];
	$pervec[1] = $x*$cell[0][0][1]+$y*$cell[0][1][1]+$z*$cell[0][2][1];
	$pervec[2] = $x*$cell[0][0][2]+$y*$cell[0][1][2]+$z*$cell[0][2][2];
	for($m1=0;$m1<@cdata2;$m1++) {
	  for($a1=0;$a1<@{$cdata2[$m1]};$a1++) {
# 	    $i++;
# 	    $j=$i+@{$cdata2[$m1]};
	    for($m2=$m1+1;$m2<@cdata2;$m2++) {
	      for($a2=0;$a2<@{$cdata2[$m2]};$a2++) {
		$j++;
		$dx = $cdata2[$m1][$a1][0]-$cdata2[$m2][$a2][0]+$pervec[0];
		$dy = $cdata2[$m1][$a1][1]-$cdata2[$m2][$a2][1]+$pervec[1];
		$dz = $cdata2[$m1][$a1][2]-$cdata2[$m2][$a2][2]+$pervec[2];
		$distsq = $dx*$dx + $dy*$dy + $dz*$dz;
		if($distsq < $dminsq) {
		  printf "%8s %7u %8s %7u %10.5f\n",$cdata2[$m1][$a1][9],$cdata2[$m1][$a1][11],
			$cdata2[$m2][$a2][9],$cdata2[$m2][$a2][11],sqrt($distsq);
		}
		if($distsq<$mindistsq[0]) {
		  $mindistsq[0]=$distsq;
		  $mindistsq[1]=$m1;
		  $mindistsq[2]=$m2;
		  $mindistsq[3]=$a1;
		  $mindistsq[4]=$a2;
# 		  $mindistsq[6]=$i;
# 		  $mindistsq[7]=$j;
		}
	      }
	    }
	  }
	}
      }
    }
  }
} else {
  for($m1=0;$m1<@cdata2;$m1++) {
    for($a1=0;$a1<@{$cdata2[$m1]};$a1++) {
      for($m2=$m1+1;$m2<@cdata2;$m2++) {
	for($a2=0;$a2<@{$cdata2[$m2]};$a2++) {
	  $dx = $cdata2[$m1][$a1][0]-$cdata2[$m2][$a2][0]+$pervec[0];
	  $dy = $cdata2[$m1][$a1][1]-$cdata2[$m2][$a2][1]+$pervec[1];
	  $dz = $cdata2[$m1][$a1][2]-$cdata2[$m2][$a2][2]+$pervec[2];
	  $distsq = $dx*$dx + $dy*$dy + $dz*$dz;
	  if($distsq < $dminsq) {
	    printf "%8s %7u %8s %7u %10.5f\n",$cdata2[$m1][$a1][9],$cdata2[$m2][$a1][11],
		  $cdata2[$m2][$a2][9],$cdata2[$m2][$a2][11],sqrt($distsq);
	  }
	  if($distsq<$mindistsq[0]) {
	    $mindistsq[0]=$distsq;
	    $mindistsq[1]=$m1;
	    $mindistsq[2]=$m2;
	    $mindistsq[3]=$a1;
	    $mindistsq[4]=$a2;
	  }
	}
      }
    }
  }
}

$m1=$mindistsq[1];
$m2=$mindistsq[2];
$a1=$mindistsq[3];
$a2=$mindistsq[4];
print "minimum distance:";
printf "%8s %7u %8s %7u %10.5f\n",$cdata2[$m1][$a1][9],($mindistsq[6]),
       $cdata2[$m2][$a2][9],($mindistsq[7]),sqrt($mindistsq[0]);
