#!/usr/bin/perl
use Math::Trig;
use Switch;

$doo=2.74375;
$tolerance=0.2; # tolerance for surface oxygens
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

if($#ARGV<9) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file with surface\n";
  print "2. name of corresponding FIELD-file (must contain molecule)\n";
  print "3. name of mol2-file with molecule\n";
  print "4. name of output CONFIG-file\n";
  print "5. name of output FIELD-file\n";
  print "6. reconst. vector 1 factor for a\n";
  print "7. reconst. vector 1 factor for b\n";
  print "8. reconst. vector 2 factor for a\n";
  print "9. reconst. vector 2 factor for b\n";
  print "10 surface z\n";
  print "[11. mol shift vector in multiples of a]\n";
  print "[12. mol shift vector in multiples of b]\n";
  exit 1;
}

$factor[0][0]=$ARGV[5]; # first index: new vector, second index: old vector
$factor[0][1]=$ARGV[6];
$factor[1][0]=$ARGV[7];
$factor[1][1]=$ARGV[8];
$molshiftvec[0]=$ARGV[10];
$molshiftvec[1]=$ARGV[11];

if($factor[0][0]==$factor[1][0] and $factor[0][1]==$factor[1][1]) {
  print "error: vectors must not be equal\n";
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
$zsurf=$ARGV[9];
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

$bindingmode=0;
$i=0;
while(<MOL2>) {
  $string = $_;
  if(/BOND/) {
    last;
  } else {
    ($mdata[$i][1],$mdata[$i][2],$mdata[$i][3],$mdata[$i][0]) =
    /^\s*\d+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    if(uc($mdata[$i][0]) eq "O") {
      if($mdata[$i][3]<$tolerance and $mdata[$i][3]>-$tolerance) {
	$bindingmode++;
      }
    }
    $i++;
  }
}
close(MOL2);

open(CONFIG, "<", $filenameincfg) or die "Can't open CONFIG-file $filenameincfg:\n$!\n";
open(FIELD, "<", $filenameinfld) or die "Can't open FIELD-file $filenameinfld:\n$!\n";

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
    }
  } case 0 {
    print "error: periodic conditions required for calculation\n";
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
  print "error: no molecules found in FIELD file $filenameinfld\n";
  exit 1;
}

$xmin=9999; # to find hydroxide in lower left corner
$ymin=9999; # to find hydroxide in lower left corner
$foundmol2name=-1;
$molindex=1;
$candindex=0;
for($t=0;$t<$field_nummols;$t++) {
  $j = 0; #number of atoms counted for the molecule
  $_=<FIELD>;
  $mol_name[$t]=$_;
  $mol_name[$t] =~ s/^\s+//;
  $mol_name[$t] =~ s/\s+$//;
  if($mol_name[$t] eq $mol2name) {
    $foundmol2name=$t;
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
    print "error: ".$mol_numatoms[$t]." expected for molecule ".$mol_name[$t]." but ".$j." were found\n";
    exit 1;
  }
  # read atom coordinates for the molecule from the CONFIG file
  for($m=0;$m<$mol_numents[$t];$m++) {
    for($a=0;$a<$mol_numatoms[$t];$a++) {
      $_=<CONFIG>;
      ($cdata[$t][$m][$a][0]) = /^\s*(\S*)/;
      if($cdata[$t][$m][$a][0] ne $mol_atomname[$t][$a]) {
	print "error: atom name ".$cdata[$t][$m][$a][0]." in CONFIG-file $filenameincfg does not match\n";
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
      # select suitable PA molecules 
      if(($cdata[$t][$m][$a][0] eq "O" )#or $cdata[$t][$m][$a][0] eq "OX")
      and $cdata[$t][$m][$a][3]<$zsurf+$tolerance and $cdata[$t][$m][$a][3]>$zsurf-$tolerance) {
	if($cdata[$t][$m][$a][1]<$xmin) {
	  $xmin=$cdata[$t][$m][$a][1];
	}
	if($cdata[$t][$m][$a][2]<$ymin) {
	  $ymin=$cdata[$t][$m][$a][2];
	}
	if($candidate[$candindex][0] != $molindex) {
	  $candidate[$candindex][0]=$molindex;
	  $candidate[$candindex][1]=$cdata[$t][$m][$a][1];
	  $candidate[$candindex][2]=$cdata[$t][$m][$a][2];
	  $candidate[$candindex][3]=$cdata[$t][$m][$a][3];
 	  $candidate[$candindex][4]=$t;
 	  $candidate[$candindex][5]=$m;
 	  $candidate[$candindex][6]=1; #number of oxygen atoms on surface
	  $candindex++;
	} else {
	  $candidate[$candindex][6]++;
	  if($cdata[$t][$m][$a][1]<$candidate[$candindex][1]) {
	    $candidate[$candindex][1]=$cdata[$t][$m][$a][1];
	    $candidate[$candindex][2]=$cdata[$t][$m][$a][2];
	    $candidate[$candindex][3]=$cdata[$t][$m][$a][3];
	  }
	}
      }
    }
    $molindex++;
  }
}
close(CONFIG, FIELD);

print "".($#candidate+1)." candidates for replacement found\n";

if($foundmol2name<0) {
  print "error: molecule named $mol2name was not found in FIELD file $filenameinfld!\n";
  exit 1;
}

$imax = 100;
$molsreplaced=0;
# find desired positions
$i=0;
for($x=-$imax;$x<=$imax;$x++) {
  for($y=-$imax;$y<=$imax;$y++) {
    $molpos[$i][0]=$xmin+$x*$molperiodicvec[0][0]+$y*$molperiodicvec[1][0]
		    +$molshiftvec[0]*$ohperiodicvec[0][0]+$molshiftvec[1]*$ohperiodicvec[1][0];
    $molpos[$i][1]=$ymin+$x*$molperiodicvec[0][1]+$y*$molperiodicvec[1][1]
		    +$molshiftvec[0]*$ohperiodicvec[0][1]+$molshiftvec[1]*$ohperiodicvec[1][1];
    if($molpos[$i][0] > -$size[0] and $molpos[$i][0] < $size[0] and 
      $molpos[$i][1] > -$size[1] and $molpos[$i][1] < $size[1]) {
      $dmin=99999;
      for($j=0;$j<=$#candidate;$j++) {
	for($xx=-1;$xx<2;$xx++) {
	  for($yy=-1;$yy<2;$yy++) {
	    $dx=$molpos[$i][0]-$candidate[$j][1]+$xx*$vec[0][0]+$yy*$vec[1][0];
	    $dy=$molpos[$i][1]-$candidate[$j][2]+$xx*$vec[0][1]+$yy*$vec[1][1];
	    $dist=sqrt($dx**2+$dy**2);
	    if($dist<$dmin) {
	      $remove[0] = $j;
	      $dmin=$dist;
	    }
	  }
	}
      }
      if($dmin>1) {
	print "sorry, no molecule found at position ".$molpos[$i][0].", ".$molpos[$i][1].".\n";
	print "please check your periodicity and shift vectors.\n";
	print $remove[1]." ".$remove[2]." $dmin\n";
      } elsif($bindingmode != $candidate[$remove[0]][6]) {
	print "molecule at position ".$candidate[$remove[0]][1].", ".$candidate[$remove[0]][2]." is bound with\n";
	print $candidate[$remove[0]][6]." oxygen atoms instead of $bindingmode.\n";
      } else {
	$mol_deltanuments[$foundmol2name]++;
	$mol_deltanuments[$candidate[$remove[0]][4]]--;
	$molpos[$i][0]=$candidate[$remove[0]][1];
	$molpos[$i][1]=$candidate[$remove[0]][2];
	$molpos[$i][2]=$candidate[$remove[0]][3];
	@removedmols[$#removedmols+1] = splice(@candidate,$remove[0],1);
	$molsreplaced++;
	$i++;
      }
    }
  }
}

$#molpos=$i-1;

@removedmols = sort { $a->[0] <=> $b->[0] } @removedmols;

open(OUTCFG,">",$ARGV[3]) or die "Can't open output CONFIG-file $filenameoutcfg:\n$!";
$index=1;
printf OUTCFG $title;
printf OUTCFG "%10u", $config_key; # OUTCFG file key
printf OUTCFG "%10u", $periodic_key; # periodic boundary key
printf OUTCFG "\n%19.12e %19.12e %19.12e", $vec[0][0], $vec[0][1], $vec[0][2];
printf OUTCFG "\n%19.12e %19.12e %19.12e", $vec[1][0], $vec[1][1], $vec[1][2];
printf OUTCFG "\n%19.12e %19.12e %19.12e", $vec[2][0], $vec[2][1], $vec[2][2];

$i=0; # index for @removedmols
$index=1; # index for molecules in CONFIG
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

open(FIELD, "<", $ARGV[1]) or die "Can't open FIELD-file $filenameinfld:\n$!";
open(OUTFLD, ">", $ARGV[4]) or die "Can't open FIELD-file:$filenameoutfld\n$!";

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
print "replaced $molsreplaced molecules\n";
exit 0;
