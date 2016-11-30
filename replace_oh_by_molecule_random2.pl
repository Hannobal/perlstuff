#!/usr/bin/perl

# replaces hydroxide ions by phosphonic acid molecules in a regular
# lattice and removes random molecules so that a certain density of
# molecules is achieved

use Math::Trig;
use Switch;

$doo = 2.755333333;
$doo3 = 3*$doo;
$tolerance = 0.2;
$minmoldist = 1.5*$doo;
$ohperiodicvec[0][0] =  $doo;
$ohperiodicvec[0][1] =  0;
$ohperiodicvec[1][0] =  $doo*cos(pi/3);
$ohperiodicvec[1][1] =  $doo*sin(pi/3);
$ohperiodicvec[2][0] = -$doo*cos(pi/3);
$ohperiodicvec[2][1] =  $doo*sin(pi/3);
$ohperiodicvec[3][0] = -$doo;
$ohperiodicvec[3][1] =  0;
$ohperiodicvec[4][0] = -$doo*cos(pi/3);
$ohperiodicvec[4][1] = -$doo*sin(pi/3);
$ohperiodicvec[5][0] =  $doo*cos(pi/3);
$ohperiodicvec[5][1] = -$doo*sin(pi/3);

#          x         
#         A          
#        /           
#     b/             
#     /              
#   /                
#  /        a        
# x---------------->x

if($#ARGV<9) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file with surface\n";
  print "2. name of corresponding FIELD-file\n";
  print "3. name of mol2-file with molecule\n";
  print "4. name of output CONFIG-file\n";
  print "5. name of output FIELD-file\n";
  print "6. reconst. vector 1 factor for a\n";
  print "7. reconst. vector 1 factor for b\n";
  print "8. reconst. vector 2 factor for a\n";
  print "9. reconst. vector 2 factor for b\n";
  print "10. requested density of molecules per squre nanometer\n";
  print "[11. mol shift vector in multiples of a]\n";
  print "[12. mol shift vector in multiples of b]\n";
  print "[13. (zero/mono/bi/tri)-dentat binding]\n";
  print "[14. additional distance to surface]\n";
  exit 1;
}

$factor[0][0]=$ARGV[5]; # first index: new vector, second index: old vector
$factor[0][1]=$ARGV[6];
$factor[1][0]=$ARGV[7];
$factor[1][1]=$ARGV[8];
$target_density = $ARGV[9];
if($#ARGV> 9){ $molshiftvec[0]=$ARGV[10]; } else { $molshiftvec[0]=0; }
if($#ARGV>10){ $molshiftvec[1]=$ARGV[11]; } else { $molshiftvec[1]=0; }
if(lc(substr($ARGV[12],0,4)) eq "zero") {
  $bindingmode=0;
} elsif(lc(substr($ARGV[12],0,4)) eq "mono") {
  $bindingmode=1;
} elsif(lc(substr($ARGV[12],0,2)) eq "bi") {
  $bindingmode=2;
} elsif(lc(substr($ARGV[12],0,3)) eq "tri") {
  $bindingmode=3;
} else {
  $bindingmode=undef;
}
if($#ARGV>12){ $surfdist=$ARGV[13]; } else { $surfdist=0; }

if($factor[0][0]==$factor[1][0] and $factor[0][1]==$factor[1][1]) {
  print "**** error: vectors must not be equal\n";
  exit 1;
}

if($target_density > 5) {
  print "**** error: target density must be smaller than 5\n";
  exit 1;
}

$molperiodicvec[0][0]=$factor[0][0]*$ohperiodicvec[0][0]+$factor[0][1]*$ohperiodicvec[1][0];
$molperiodicvec[0][1]=$factor[0][0]*$ohperiodicvec[0][1]+$factor[0][1]*$ohperiodicvec[1][1];
$molperiodicvec[1][0]=$factor[1][0]*$ohperiodicvec[0][0]+$factor[1][1]*$ohperiodicvec[1][0];
$molperiodicvec[1][1]=$factor[1][0]*$ohperiodicvec[0][1]+$factor[1][1]*$ohperiodicvec[1][1];

$filenameincfg  = $ARGV[0];
$filenameinfld  = $ARGV[1];
$filenamemol2   = $ARGV[2];
$filenameoutcfg = $ARGV[3];
$filenameoutfld = $ARGV[4];
$mol2name=$ARGV[2];
$mol2name =~ s/\.\S*$//; # remove file extension
$mol2name = substr($mol2name,rindex($mol2name,'/')+1); # remove path (everything before last "/")

# indices for cdata array:
# 0 name
# 1 x
# 2 y
# 3 z
# 4 velocity x
# 5 velocity y
# 6 velocity z
# 7 force x
# 8 force y
# 9 force z

#load molecule
open(MOL2, "<", $filenamemol2) or die "Can't open mol2-file $filenamemol2:\n$!\n";
while(<MOL2>) {
  if(/ATOM/) { last; }
}

$mol2charge = 0;
$mol2numox = 0;
$i=0;
while(<MOL2>) {
  $string = $_;
  if(/BOND/) {
    last;
  } else {
    ($mdata[$i][1],$mdata[$i][2],$mdata[$i][3],$mdata[$i][0],$mdata[$i][4]) =
    /^\s*\d+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)/; # position, type and charge
    $mol2charge += $mdata[$i][4];
    if(uc($mdata[$i][0]) eq "O" or uc($mdata[$i][0]) eq "OH") {
      $mol2oxpos[$mol2numox][0]=$mdata[$i][1];
      $mol2oxpos[$mol2numox][1]=$mdata[$i][2];
      $mol2oxpos[$mol2numox][2]=$mdata[$i][3];
      $mol2numox++;
    }
    $i++;
  }
}
close(MOL2);
@mol2oxpos = sort { $a->[2] <=> $b->[2] } @mol2oxpos;
$mol2charge = sprintf("%.0f",$mol2charge);
if(not defined($bindingmode)) {
  if($mol2charge==0) {
    $bindingmode=0;
  }elsif($mol2charge==-1) {
    $bindingmode=1;
  }elsif($mol2charge==-2) {
    $bindingmode=2;
  }
}

open(CONFIG, "<", $filenameincfg) or die "Can't open CONFIG-file $filenameincfg:\n$!\n";
open(FIELD, "<", $filenameinfld) or die "Can't open FIELD-file $filenameinfld:\n$!\n";
# open(TEST, ">", "test.xyz") or die "Can't open test.xyz\n$!\n";

#read header of CONFIG-file
$title=<CONFIG>;
$_=<CONFIG>;
$config_header[0]=$_;
($config_key,$periodic_key) = /^\s+(\S+)\s+(\S+)/;
if(not ($periodic_key==1 or $periodic_key==2 or $periodic_key==6)) {
  print "sorry, this script is only implemented for rectangular cells\n";
  exit 1;
}
switch($periodic_key) {
  case [1,2,6] {
    for($j=0;$j<=2;$j++) {
      $_=<CONFIG>;
      $config_header[$j+1]=$_;
      ($vec[$j][0], $vec[$j][1], $vec[$j][2]) = /\s*(\S*)\s*(\S*)\s*(\S*)/;
      $size[$j]=$vec[$j][$j]/2;
      $margin[$j] = $size[$j]-4;
    }
  } case 0 {
    print "**** error: periodic conditions needed for calculation\n";
    exit 1;
  } else {
    print "sorry, this script is only implemented for rectangular cells\n";
    exit 1;
  }
}

# read FIELD file and rest of CONFIG file
while(<FIELD>) {
  if(/molecules/ or /MOLECULES/) { last; }
}
($field_nummols) = /^\s*\S+\s+(\S+)/;
#   print "found $field_nummols molecules\n";
if($field_nummols==0) {
  print "no molecules found in FIELD file $filenameinfld\n";
  exit 1;
}

$zmin=7; #cutoff for surface atoms
$xmin=9999; # to find hydroxide in lower left corner
$ymin=9999; # to find hydroxide in lower left corner
$foundox=-1;
$foundhydrox=-1;
$foundal=-1;
$foundmol2name=-1;
$molindex=1;
for($t=0;$t<$field_nummols;$t++) {
  $j = 0; #number of atoms counted for the molecule
  $_=<FIELD>;
  $mol_name[$t]=$_;
  $mol_name[$t] =~ s/^\s+//;
  $mol_name[$t] =~ s/\s+$//;
  if($mol_name[$t] eq $mol2name) {
    $foundmol2name=$t
  }
  if(/^\s*oxygen/ or /^\s*OXYGEN/ or /^\s*Oxygen/ ) {
    if(/free/ or /Free/ or /FREE/) {
      $foundox=$t;
    }
  } elsif(/^\s*hydroxide/ or /^\s*HYDROXIDE/ or /^\s*Hydroxide/ ) {
    $foundhydrox=$t;
  } elsif(/^\s*aluminum/ or /^\s*ALUMINUM/ or /^\s*Aluminum/ ) {
    $foundal=$t;
  }
  $_=<FIELD>;
  ($mol_numents[$t])=/^\s*\S+\s+(\S+)/;
  $_=<FIELD>;
  ($mol_numatoms[$t])=/^\s*\S+\s+(\S+)/;
#   print "found ".$mol_numatoms[$t]." in $t. molecule named ".$mol_name[$t]." with ".$mol_numents[$t]." entities\n";
  while(<FIELD>) {
    if(/FINISH/ or /finish/ or /CLOSE/ or /close/ or /BONDS/
       or /bonds/ or /ANGLES/ or /angles/ or /DIHEDRALS/ or /dihedrals/) {
      if(not (/FINISH/ or /finish/)) {
	while(<FIELD>) { 
	  if(/FINISH/ or /finish/) { last; }
	}
      }
      last;
    }
    ($mol_atomname[$t][$j],$mol_atommass[$t][$j],$mol_atomcharge[$t][$j])=/^\s*(\S+)\s+(\S+)\s+(\S+)/;
#     printf "%8s %3.5f %3.5f\n", $mol_atomname[$t][$j],$mol_atommass[$t][$j],$mol_atomcharge[$t][$j];
    $j++;
  }
  if($j!=$mol_numatoms[$t]) {
    print "**** error: ".$mol_numatoms[$t]." expected for molecule ".$mol_name[$t]." but ".$j." were found\n";
    exit 1;
  }
  # read atom coordinates for the molecule from the CONFIG file
  for($m=0;$m<$mol_numents[$t];$m++) {
    for($a=0;$a<$mol_numatoms[$t];$a++) {
      $_=<CONFIG>;
      ($cdata[$t][$m][$a][0]) = /^\s*(\S*)/;
      if($cdata[$t][$m][$a][0] ne $mol_atomname[$t][$a]) {
	print "**** error: atom name ".$cdata[$t][$m][$a][0]." in CONFIG-file $filenameincfg does not match\n";
	print "name ".$mol_atomname[$t][$a]." in FIELD file $filenameinfld for molecule ".$mol_name[$t]."\n";
	exit 1;
      }
      $_=<CONFIG>;
      ($cdata[$t][$m][$a][1],$cdata[$t][$m][$a][2],$cdata[$t][$m][$a][3]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
      if($config_key>0) { #read velocity
	$_=<CONFIG>;
	($cdata[$t][$m][$a][4],$cdata[$t][$m][$a][5],$cdata[$t][$m][$a][6]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
	if($config_key>1) { #read force
	  $_=<CONFIG>;
	  ($cdata[$t][$m][$a][7],$cdata[$t][$m][$a][8],$cdata[$t][$m][$a][9]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
	}
      }
      $cdata[$t][$m][$a][10]=$molindex;
    }
    $molindex++;
  }
}
close(CONFIG, FIELD);

# find maximum z value of surface
$zsurf=-9999;
$candindex=0;
for($t=0;$t<$field_nummols;$t++) {
  if(lc($mol_name[$t]) =~ /hydroxide/) {
    for($m=0;$m<$mol_numents[$t];$m++) {
      for($a=0;$a<$mol_numatoms[$t];$a++) {
	if(uc($cdata[$t][$m][$a][0]) eq "OX") {
	  if($cdata[$t][$m][$a][3]>$zsurf) {
	    $zsurf=$cdata[$t][$m][$a][3];
	  }
	}
      }
    }
    $zmintol=$zsurf-$tolerance;
    $zmaxtol=$zsurf+$tolerance;
    for($m=0;$m<$mol_numents[$t];$m++) {
      for($a=0;$a<$mol_numatoms[$t];$a++) {
	if((uc($cdata[$t][$m][$a][0]) eq "OX")
	   and $cdata[$t][$m][$a][3]>$zmintol
	   and $cdata[$t][$m][$a][3]<$zmaxtol) {
	  if($cdata[$t][$m][$a][1]<$xmin) {
	    $xmin=$cdata[$t][$m][$a][1];
	  }
	  if($cdata[$t][$m][$a][2]<$ymin) {
	    $ymin=$cdata[$t][$m][$a][2];
	  }
	  $candidate[$candindex][0]=$cdata[$t][$m][$a][10];
	  $candidate[$candindex][1]=$cdata[$t][$m][$a][1];
	  $candidate[$candindex][2]=$cdata[$t][$m][$a][2];
	  $candidate[$candindex][3]=$cdata[$t][$m][$a][3];
	  $candidate[$candindex][4]=$t;
	  $candidate[$candindex][5]=$m;
	  $candindex++;
	}
      }
    }
    last;
  }
}

print "found ".($#candidate+1)."candidates\n";
$surface_area=($vec[0][0]*$vec[1][1]-$vec[1][0]*$vec[0][1])/100; # in squarenm
$nummols_needed = sprintf("%.0f",$surface_area*$target_density);

if($foundox==-1) {
  print "**** error: could not find entry for oxygen in FIELD file $filenameinfld!\n";
  exit 1;
}
if($foundhydrox==-1) {
  print "**** error: could not find entry for hydroxide ions in FIELD file $filenameinfld!\n";
  exit 1;
}
if($foundmol2name==-1) {
  print "**** error: could not find entry for molecule $mol2name in FIELD file $filenameinfld!\n";
  exit 1;
}

$imaxx = int(2*$size[0]/$doo+2)+2;
$imaxy = int(2*$size[1]/$ohperiodicvec[1][1])+2;

$i=0;
for($ix=-$imaxx;$ix<=$imaxx;$ix++) {
  for($iy=-$imaxy;$iy<=$imaxy;$iy++) {
    $molpos[$i][0]=$xmin+$ix*$molperiodicvec[0][0]+$iy*$molperiodicvec[1][0]
		    +$molshiftvec[0]*$ohperiodicvec[0][0]+$molshiftvec[1]*$ohperiodicvec[1][0];
    $molpos[$i][1]=$ymin+$ix*$molperiodicvec[0][1]+$iy*$molperiodicvec[1][1]
		    +$molshiftvec[0]*$ohperiodicvec[0][1]+$molshiftvec[1]*$ohperiodicvec[1][1];
    $dmin=9999;
    for($j=0;$j<=$#candidate;$j++) {
      for($x=-1;$x<2;$x++) {
	for($y=-1;$y<2;$y++) {
	  $dx=$molpos[$i][0]+$mol2oxpos[0][0]-$candidate[$j][1]-$x*$vec[0][0]-$y*$vec[1][0];
	  $dy=$molpos[$i][1]+$mol2oxpos[0][1]-$candidate[$j][2]-$x*$vec[0][1]-$y*$vec[1][1];
	  $dist=sqrt($dx**2+$dy**2);
	  if($dist<$dmin) {
	    $newpos[0] = $candidate[$j][1]+$x*$vec[0][0]+$y*$vec[1][0];
	    $newpos[1] = $candidate[$j][2]+$x*$vec[0][1]+$y*$vec[1][1];
	    $dmin=$dist;
	  }
	}
      }
    }
    $molpos[$i][0] = $newpos[0];
    $molpos[$i][1] = $newpos[1];
    if($molpos[$i][0] > -$size[0] and $molpos[$i][0] < $size[0] and 
      $molpos[$i][1] > -$size[1] and $molpos[$i][1] < $size[1]) {
      $possible=1;
      for($j=0;$j<$#molpos;$j++) {
	for($x=-1;$x<2;$x++) {
	  for($y=-1;$y<2;$y++) {
	    $dx=$molpos[$i][0]-$molpos[$j][0]-$x*$vec[0][0]-$y*$vec[1][0];
	    $dy=$molpos[$i][1]-$molpos[$j][1]-$x*$vec[0][1]-$y*$vec[1][1];
	    $dist=sqrt($dx**2+$dy**2);
	    if($dist<5) {
	    }
	    if($dist<$minmoldist) {
	      $possible=0;
	      break;
	    }
	  }
	  break if not $possible;
	}
	break if not $possible;
      }
      if($possible) {
	$mol_deltanuments[$foundmol2name]++;
	$i++;
      }
    }
  }
}
$#molpos=$i-1;

$ohremoved=0;
$nummols_remove = $i - $nummols_needed;
if($nummols_remove<0) {
  print "**** error: density cannot be achieved with the specified lattice!\n";
  exit 1;
}

print "randomly removing $nummols_remove molecules\n";
for($i=0;$i<$nummols_remove;$i++) {
  $j = rand(@molpos);
  splice(@molpos,$j,1);
  $mol_deltanuments[$foundmol2name]--;
}

$end_density = @molpos/$surface_area;
print "final density: $end_density molecule per square nanometer\n";

# for($i=0;$i<100;$i++) {
#   @possible=();
#   for($m=0;$m<@molpos;$m++) {
#     $force[$m][0] = 0;
#     $force[$m][1] = 0;
#     for($n=0;$n<@molpos;$n++) {
#       if($m==$n) { next; }
#       for($x=-1;$x<2;$x++) {
# 	for($y=-1;$y<2;$y++) {
# 	  $dx=$molpos[$m][0]-$molpos[$n][0]-$x*$vec[0][0]-$y*$vec[1][0];
# 	  $dy=$molpos[$m][1]-$molpos[$n][1]-$x*$vec[0][1]-$y*$vec[1][1];
# 	  $dist[$m][$n]=sqrt($dx**2+$dy**2);
# 	  if($dist[$m][$n]<$doo3) {
# 	    $force[$m][0] += $dx;
# 	    $force[$m][1] += $dy;
# 	    $force[$m][3] += $dist[$m][$n];
# 	  }
# 	}
#       }
#     }
#     $direction = sprintf("%.0f", atan2($force[$m][1],$force[$m][0])/pi*3);
#     $direction += 6 if $direction<0;
#     $newmolpos[0] = $molpos[$m][0] + $ohperiodicvec[$direction][0];
#     $newmolpos[1] = $molpos[$m][1] + $ohperiodicvec[$direction][1];
#   }
# }
# undef @dist;

open(TEST, ">", "test.xyz");
print "shuffling molecule positions\n";
for($i=0;$i<100;$i++) {
  print TEST $#molpos+1, "\n\n";
  if($i%10==0) {print "$i percent done\n";}
  for($m=0;$m<@molpos;$m++) {
    $k = int(rand(@ohperiodicvec));
    $newmolpos[0] = $molpos[$m][0] + $ohperiodicvec[$k][0];
    $newmolpos[1] = $molpos[$m][1] + $ohperiodicvec[$k][1];
    $possible =1;
    for($n=0;$n<@molpos;$n++) {
      if($m==$n) { next; }
      for($x=-1;$x<2;$x++) {
	for($y=-1;$y<2;$y++) {
	  $dx=$newmolpos[0]-$molpos[$n][0]-$x*$vec[0][0]-$y*$vec[1][0];
	  $dy=$newmolpos[1]-$molpos[$n][1]-$x*$vec[0][1]-$y*$vec[1][1];
	  $dist=sqrt($dx**2+$dy**2);
	  if($dist<$minmoldist) {
	    $possible = 0;
	    break;
	  }
	}
	break if not $possible;
      }
      break if not $possible;
    }
    if($possible) {
#       print "moved molecule $m\n";
      $molpos[$m][0] = $newmolpos[0];
      $molpos[$m][1] = $newmolpos[1];
      if($molpos[$m][0] > $size[0]) {
	$molpos[$m][0] -= 2*$size[0]
      } elsif($molpos[$m][0] < -$size[0]) {
	$molpos[$m][0] += 2*$size[0]
      }
      if($molpos[$m][1] > $size[1]) {
	$molpos[$m][1] -= 2*$size[1]
      } elsif($molpos[$m][1] < -$size[1]) {
	$molpos[$m][1] += 2*$size[1]
      }
    }
    printf TEST "%3s %10.5f %10.5f %10.5f\n", "C", @{$molpos[$m]}
  }
}
close(TEST);

$ohremoved=0;

for($i=0;$i<=$#molpos;$i++) {
  $xpos=0;
  $ypos=0;
  $zpos=0;
  for($b=0;$b<$bindingmode;$b++) {
    $dmin=99999;
    for($j=0;$j<=$#candidate;$j++) {
      for($x=-1;$x<2;$x++) {
	for($y=-1;$y<2;$y++) {
	  $dx=$molpos[$i][0]-$candidate[$j][1]-$x*$vec[0][0]-$y*$vec[1][0]+$mol2oxpos[$b][0];
	  $dy=$molpos[$i][1]-$candidate[$j][2]-$x*$vec[0][1]-$y*$vec[1][1]+$mol2oxpos[$b][1];
	  $dist=sqrt($dx**2+$dy**2);
	  if($dist<$dmin) {
	    $remove[0] = $j;
	    $remove[1] = $candidate[$j][1]+$x*$vec[0][0]+$y*$vec[1][0];
	    $remove[2] = $candidate[$j][2]+$x*$vec[0][1]+$y*$vec[1][1];
	    $remove[3] = $candidate[$j][3];
	    $dmin=$dist;
	  }
	}
      }
    }
    $xpos += $remove[1];
    $ypos += $remove[2];
    $zpos += $remove[3];
    $mol_deltanuments[$candidate[$remove[0]][4]]--;
    $ohremoved++;
    @removedmols[$#removedmols+1] = splice(@candidate,$remove[0],1);
  }
  if($bindingmode>0) {
    $molpos[$i][0]=$xpos/$bindingmode;
    $molpos[$i][1]=$ypos/$bindingmode;
    $molpos[$i][2]=$zpos/$bindingmode+$surfdist;
  } else {
    $molpos[$i][2]=$zsurf+$surfdist;
  }
}

@removedmols = sort { $a->[0] <=> $b->[0] } @removedmols;

print "replaced $ohremoved hydroxides by $mol_deltanuments[$foundmol2name] molecules in $bindingmode"."dentate coordination\n";
# print "".($mol_deltanuments[$foundmol2name]/$area)." molecules per square nanometer\n";

open(OUTCFG,">",$ARGV[3]);
$index=1;
printf OUTCFG $title;
printf OUTCFG "%10u", $config_key; # OUTCFG file key
printf OUTCFG "%10u", $periodic_key; # periodic boundary key
printf OUTCFG "\n%19.12e %19.12e %19.12e", $vec[0][0], $vec[0][1], $vec[0][2];
printf OUTCFG "\n%19.12e %19.12e %19.12e", $vec[1][0], $vec[1][1], $vec[1][2];
printf OUTCFG "\n%19.12e %19.12e %19.12e", $vec[2][0], $vec[2][1], $vec[2][2];

$i=0; # index for @removedmols
$index=0; # index for molecules in CONFIG
for($t=0;$t<$field_nummols;$t++) {
  for($m=0;$m<$mol_numents[$t];$m++) {
    if($cdata[$t][$m][0][10] != $removedmols[$i][0]) {
      for($a=0;$a<$mol_numatoms[$t];$a++) {
	printf OUTCFG "\n%-8s %9u", $cdata[$t][$m][$a][0], $index;
	printf OUTCFG "\n%19.12e",$cdata[$t][$m][$a][1];
	printf OUTCFG  " %19.12e",$cdata[$t][$m][$a][2];
	printf OUTCFG  " %19.12e",$cdata[$t][$m][$a][3];
	for($k=0;$k<$config_key;$k++) {
	  printf OUTCFG " \n%19.12e %19.12e %19.12e", 0, 0, 0;
	}
	$index++;
      }
    } else {
      $i++;
    }
  }
  if($mol_name[$t] eq $mol2name) {
    for($m=0;$m<=$#molpos;$m++) {
      for($a=0;$a<=$#mdata;$a++) {
	printf OUTCFG "\n%-8s %9u", uc($mdata[$a][0]), $index;
	printf OUTCFG "\n%19.12e",$mdata[$a][1]+$molpos[$m][0];
	printf OUTCFG  " %19.12e",$mdata[$a][2]+$molpos[$m][1];
	printf OUTCFG  " %19.12e",$mdata[$a][3]+$molpos[$m][2];
	for($k=0;$k<$config_key;$k++) {
	  printf OUTCFG " \n%19.12e %19.12e %19.12e", 0, 0, 0;
	}
	$index++;
      }
    }
  }
}
close(OUTCFG);

open(FIELD, "<", $ARGV[1]) or die "Can't open FIELD-file:\n$!";
open(OUTFLD, ">", $ARGV[4]) or die "Can't open FIELD-file:\n$!";

while(<FIELD>) {
  print OUTFLD $_;
  if(/molecules/ or /MOLECULES/) { last; }
}
for($t=0;$t<$field_nummols;$t++) {
  $_=<FIELD>; # molecule name
  print OUTFLD $_;
  $_=<FIELD>; # nummols
  print OUTFLD "NUMMOLS ".($mol_numents[$t]+$mol_deltanuments[$t])."\n";
  while(<FIELD>) {
    print OUTFLD $_;
    if(/finish/ or /FINISH/) { last; }
  }
}
while(<FIELD>) {
  print OUTFLD $_;
}
close(FIELD, OUTFLD);
exit 0;
