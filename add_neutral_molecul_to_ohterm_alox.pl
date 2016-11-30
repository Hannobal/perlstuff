#!/usr/bin/perl

# adds a molecule (not necessarily neutral) from a mol2-file to an aluminum
# oxide surface in a periodic lattice and with a certain distance to the surface.
# the FIELD file must already contain the molecule to be added.
# (apply combine_fields.pl to add the molecule before using this script)

use Math::Trig;
use Switch;

$distmolsurf=1.5;
$doo=2.75533;
$ohperiodicvec[0][0]=$doo;
$ohperiodicvec[0][1]=0;
$ohperiodicvec[1][0]=$doo*cos(pi/3);
$ohperiodicvec[1][1]=$doo*sin(pi/3);

#          x         
#         A          
#        /           
#     b/             
#     /              
#   /                
#  /        a        
# x---------------->x

if($#ARGV<8) {
  print "input format:\n";
  print "  1. name of CONFIG/REVCON-file with surface\n";
  print "  2. name of corresponding FIELD-file\n";
  print "  3. name of mol2-file with molecule\n";
  print "  4. name of output CONFIG-file\n";
  print "  5. name of output FIELD-file\n";
  print "  6. reconst. vector 1 factor for a\n";
  print "  7. reconst. vector 1 factor for b\n";
  print "  8. reconst. vector 2 factor for a\n";
  print "  9. reconst. vector 2 factor for b\n";
  print "[10. distance to surface]\n";
  print "[11. shift in x- and y-direction]\n";
  exit;
}

$factor[0][0]=$ARGV[5]; # first index: new vector, second index: old vector
$factor[0][1]=$ARGV[6];
$factor[1][0]=$ARGV[7];
$factor[1][1]=$ARGV[8];

if($factor[0][0]==$factor[1][0] and $factor[0][1]==$factor[1][1]) {
  print "error: vectors must not be equal\n";
  exit;
}

if($#ARGV>8) {
  $distmolsurf=$ARGV[9];
}

if($#ARGV>9) {
  $molshiftvec[0]=$ARGV[10];
}
if($#ARGV>10) {
  $molshiftvec[1]=$ARGV[11];
}

$molperiodicvec[0][0]=$factor[0][0]*$ohperiodicvec[0][0]+$factor[0][1]*$ohperiodicvec[1][0];
$molperiodicvec[0][1]=$factor[0][0]*$ohperiodicvec[0][1]+$factor[0][1]*$ohperiodicvec[1][1];
$molperiodicvec[1][0]=$factor[1][0]*$ohperiodicvec[0][0]+$factor[1][1]*$ohperiodicvec[1][0];
$molperiodicvec[1][1]=$factor[1][0]*$ohperiodicvec[0][1]+$factor[1][1]*$ohperiodicvec[1][1];

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
open(MOL2, "<", $ARGV[2]) or die "Can't open mol2-file:\n$!\n";
while(<MOL2>) {
  if(/ATOM/) { last; }
}

$i=0;
while(<MOL2>) {
  $string = $_;
  if(/BOND/) {
    last;
  } else {
    ($mdata[$i][1],$mdata[$i][2],$mdata[$i][3],$mdata[$i][0]) =
    /^\s*\d+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    $i++;
  }
}
close(MOL2);

open(CONFIG, "<", $ARGV[0]) or die "Can't open CONFIG-file:\n\n$!";
open(FIELD, "<", $ARGV[1]) or die "Can't open FIELD-file:\n$!\n";

#read header of CONFIG-file
$title=<CONFIG>;
$_=<CONFIG>;
$config_header[0]=$_;
($config_key,$periodic_key) = /^\s+(\S+)\s+(\S+)/;
if(not ($periodic_key==1 or $periodic_key==2 or $periodic_key==6)) {
  print "sorry, this script is only implemented for rectangular cells\n";
  exit;
}
switch($periodic_key) {
  case [1,2,6] {
    for($j=0;$j<=2;$j++) {
      $_=<CONFIG>;
      $config_header[$j+1]=$_;
      ($vec[$j][0], $vec[$j][1], $vec[$j][2]) = /\s*(\S*)\s*(\S*)\s*(\S*)/;
      $size[$j]=$vec[$j][$j]/2;
    }
  } case 0 {
    print "error: periodic conditions needed for calculation\n";
    exit;
  } else {
    print "sorry, this script is only implemented for rectangular cells\n";
    exit;
  }
}

# read FIELD file and rest of CONFIG file
while(<FIELD>) {
  if(/molecules/ or /MOLECULES/) { last; }
}
($field_nummols) = /^\s*\S+\s+(\S+)/;
#   print "found $field_nummols molecules\n";
if($field_nummols==0) {
  print "no molecules found in FIELD file\n";
  exit;
}

$mol2name=substr($ARGV[2],0,-5);
$i=0; #atom index
$zmax=-999; # to find topmost hydroxide
$xmin=9999; # to find hydroxide in lower left corner
$ymin=9999; # to find hydroxide in lower left corner
$foundhydrox=-1;
$foundal=-1;
$foundmol2name=-1;
for($t=0;$t<$field_nummols;$t++) {
  $j = 0; #number of atoms counted for the molecule
  $_=<FIELD>;
  $mol_name[$t]=$_;
  $mol_name[$t] =~ s/^\s+//;
  $mol_name[$t] =~ s/\s+$//;
  if($mol_name[$t] eq $mol2name) {
    $foundmol2name=$t
  }
  $startatomindex[$t]=$i;
  if(/^\s*hydroxide/ or /^\s*HYDROXIDE/ or /^\s*Hydroxide/ ) {
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
    print "error: ".$mol_numatoms[$t]." expected for molecule ".$mol_name[$t]." but ".$j." were found\n";
    exit;
  }
  # read atom coordinates for the molecule from the CONFIG file
  for($m=0;$m<$mol_numents[$t];$m++) {
    for($a=0;$a<$mol_numatoms[$t];$a++) {
      $_=<CONFIG>;
      ($cdata[$i][0]) = /^\s*(\S*)/;
      if($cdata[$i][0] ne $mol_atomname[$t][$a]) {
	print "error: atom name".$cdata[$i][0]." in CONFIG-file does not match\n";
	print "name ".$mol_atomname[$t][$a]." in FIELD file for molecule ".$mol_name[$t]."\n";
	exit;
      }
      $_=<CONFIG>;
      ($cdata[$i][1],$cdata[$i][2],$cdata[$i][3]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
      if($config_key>0) { #read velocity
	$_=<CONFIG>;
	($cdata[$i][4],$cdata[$i][5],$cdata[$i][6]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
	if($config_key>1) { #read force
	  $_=<CONFIG>;
	  ($cdata[$i][7],$cdata[$i][8],$cdata[$i][9]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
	}
      }
      if(uc($cdata[$i][0]) eq 'OX' and $t==$foundhydrox) {
	if($cdata[$i][1]<$xmin) {
	  $xmin=$cdata[$i][1];
	}
	if($cdata[$i][2]<$ymin) {
	  $ymin=$cdata[$i][2];
	}
      } elsif(uc($cdata[$i][0]) eq 'HG' and $t==$foundhydrox) {
	if($cdata[$i][3]>$zmax) {
	  $zmax=$cdata[$i][3];
	}
      }
      $i++;
    }
  $endatomindex[$t]=$i;
  }
}
close(CONFIG, FIELD);

if($foundhydrox==-1) {
  print "error: could not find entry for hydroxide ions in FIELD file!\n";
  exit;
}
if($foundmol2name==-1) {
  print "error: could not find entry for molecule name $mol2name in FIELD file!\n";
  exit;
}

$imax = 100;

$i=0;
for($x=-$imax;$x<=$imax;$x++) {
  for($y=-$imax;$y<=$imax;$y++) {
    $molpos[$i][0]=$xmin+$x*$molperiodicvec[0][0]+$y*$molperiodicvec[1][0]
                  +$molshiftvec[0]*$ohperiodicvec[0][0]+$molshiftvec[1]*$ohperiodicvec[1][0];
    $molpos[$i][1]=$ymin+$x*$molperiodicvec[0][1]+$y*$molperiodicvec[1][1]
                  +$molshiftvec[0]*$ohperiodicvec[0][1]+$molshiftvec[1]*$ohperiodicvec[1][1];
    $molpos[$i][2]=$zmax+$distmolsurf;
    if($molpos[$i][0] > -$size[0] and $molpos[$i][0] <= $size[0] and 
       $molpos[$i][1] > -$size[1] and $molpos[$i][1] <= $size[1]) {
      $i++;
    }
  }
}
$#molpos=$i-1;

open(OUTCFG,">",$ARGV[3]);
$index=1;
printf OUTCFG $title;
printf OUTCFG "%10u", $config_key; # OUTCFG file key
printf OUTCFG "%10u", $periodic_key; # periodic boundary key
printf OUTCFG "\n%19.12e %19.12e %19.12e", $vec[0][0], $vec[0][1], $vec[0][2];
printf OUTCFG "\n%19.12e %19.12e %19.12e", $vec[1][0], $vec[1][1], $vec[1][2];
printf OUTCFG "\n%19.12e %19.12e %19.12e", $vec[2][0], $vec[2][1], $vec[2][2];
for($i=0;$i<=$#molpos;$i++) {
  for($j=0;$j<=$#mdata;$j++) {
    printf OUTCFG "\n%-8s %9u", uc($mdata[$j][0]), $index;
    printf OUTCFG "\n%19.12e",$mdata[$j][1]+$molpos[$i][0];
    printf OUTCFG  " %19.12e",$mdata[$j][2]+$molpos[$i][1];
    printf OUTCFG  " %19.12e",$mdata[$j][3]+$molpos[$i][2];
    for($k=0;$k<$config_key;$k++) {
      printf OUTCFG " \n%19.12e %19.12e %19.12e", 0, 0, 0;
    }
    $index++;
  }
}
for($i=0;$i<=$#cdata;$i++) {
  printf OUTCFG "\n%-8s %9u", $cdata[$i][0], $index;
  printf OUTCFG "\n%19.12e %19.12e %19.12e",$cdata[$i][1],$cdata[$i][2],$cdata[$i][3];
  for($k=0;$k<$config_key;$k++) {
    printf OUTCFG " \n%19.12e %19.12e %19.12e", 0, 0, 0;
  }
  $index++;
}
close(OUTCFG);

open(FIELD, "<", $ARGV[1]) or die "Can't open FIELD-file:\n$!";
open(OUTFLD, ">", $ARGV[4]) or die "Can't open FIELD-file:\n$!";

while(<FIELD>) {
  if(/$mol2name/) {
    print OUTFLD $_;
    $_=<FIELD>;
    $n = $#molpos+1;
    print OUTFLD "NUMMOLS $n\n";
  } else {
    print OUTFLD $_;
  }
}

close(FIELD, OUTFLD);
print "added ".($#molpos+1)." molecules\n";
