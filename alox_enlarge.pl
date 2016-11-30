#!/usr/bin/perl
use Switch;

# repeats the content of a DL_POLY CONFIG file in x- and y- direction
# and cuts a surface with a certain slab thickness

if($#ARGV<6) {
  print "1. dl_poly REVCON-file\n";
  print "2. corresponding dl_poly FIELD-file\n";
  print "3. repetitions in x-direction (integer)\n";
  print "4. repetitions in y-direction (integer)\n";
  print "5. thickness of the slab in Angstrom\n";
  print "6. output-CONFIG with ending\n";
  print "7. output-FIELD with ending\n";
  exit;
}

$filenameincfg=$ARGV[0];
$filenameinfld=$ARGV[1];
$repx=$ARGV[2];
$repy=$ARGV[3];
$thickness=$ARGV[4];
$filenameoutcfg=$ARGV[5];
$filenameoutfld=$ARGV[6];
if($repx<1 or $repy<1) {
  print "repetitions must be >= 1\n";
  exit;
}
if($thickness<2) {
  print "You don't really want a slab of only $ARGV[3] angstroms thick...\n";
  exit;
}

open(CONFIG, "<", $filenameincfg) or die "Can't open input CONFIG-file: $!";
open(FIELD,  "<", $filenameinfld) or die "Can't open input FIELD-file: $!";
$title=<CONFIG>;
$_=<CONFIG>;
($config_key,$periodic_key) = /^\s+(\S+)\s+(\S+)/;
if(not ($periodic_key==1 or $periodic_key==2 or $periodic_key==6)) {
  print "sorry, this script is only implemented for rectangular cells\n";
  exit;
}
if($periodic_key==0) {
  print "error: periodic conditions required for calculation\n";
  exit;
} else {
  for($j=0;$j<=2;$j++) {
    $_=<CONFIG>;
    ($vec[$j][0], $vec[$j][1], $vec[$j][2]) = /\s*(\S*)\s*(\S*)\s*(\S*)/;
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
$molindex=1;
@maxpos=(-99999,-99999,-99999);
@minpos=( 99999, 99999, 99999);
for($t=0;$t<$field_nummols;$t++) {
  $j = 0; #number of atoms counted for the molecule
  $_=<FIELD>;
  $mol_name[$t]=$_;
  $mol_name[$t] =~ s/^\s+//;
  $mol_name[$t] =~ s/\s+$//;
  $_=<FIELD>;
  ($mol_numents[$t])=/^\s*\S+\s+(\S+)/;
  $_=<FIELD>;
  ($mol_numatoms[$t])=/^\s*\S+\s+(\S+)/;
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
      ($cdata[$t][$m][$a][0]) = /^\s*(\S*)/;
      if($cdata[$t][$m][$a][0] ne $mol_atomname[$t][$a]) {
	print "error: atom name ".$cdata[$t][$m][$a][0]." in CONFIG-file does not match\n";
	print "name ".$mol_atomname[$t][$a]." in FIELD file for molecule ".$mol_name[$t]."\n";
	exit;
      }
      $_=<CONFIG>;
      ($cdata[$t][$m][$a][1],$cdata[$t][$m][$a][2],$cdata[$t][$m][$a][3]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
      for($i=1;$i<4;$i++) {
	if($cdata[$t][$m][$a][$i]>$maxpos[$i]) {
	  $maxpos[$i]=$cdata[$t][$m][$a][$i];
	}
	if($cdata[$t][$m][$a][$i]<$minpos[$i]) {
	  $minpos[$i]=$cdata[$t][$m][$a][$i];
	}
      }
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

$size[1] = ($vec[0][0]+$vec[1][0]+$vec[2][0])/2;
$size[2] = ($vec[0][1]+$vec[1][1]+$vec[2][1])/2;
$size[3] = ($vec[0][2]+$vec[1][2]+$vec[2][2])/2;
$cmax=3;
# image molecules (make whole if extended over periodic boundaries)
for($t=0;$t<$field_nummols;$t++) {
  if($mol_numatoms[$t]>1) {
    for($m=0;$m<$mol_numents[$t];$m++) {
      for($a=1;$a<$mol_numatoms[$t];$a++) {
	$a_old = $a-1;
	for($c=1;$c<=$cmax;$c++) {
	  $d[$c] = $cdata[$t][$m][$a][$c]-$cdata[$t][$m][$a_old][$c];
	  if(abs($d[$c])>$size[$c]) {
	    if($cdata[$t][$m][$a][$c]>0) {
	      $cdata[$t][$m][$a][$c] -= 2*$size[$c];
	    } else {
	      $cdata[$t][$m][$a][$c] += 2*$size[$c];
	    }
	  }
	}
      }
    }
  }
}

$zmin=$maxpos[3]-$thickness;
$dx = -(($rep1-1)*$vec[0][0]+($rep2-1)*$vec[1][0]+($rep3-1)*$vec[2][0])/2;
$dy = -(($rep1-1)*$vec[0][1]+($rep2-1)*$vec[1][1]+($rep3-1)*$vec[2][1])/2;
$dz = -(($rep1-1)*$vec[0][2]+($rep2-1)*$vec[1][2]+($rep3-1)*$vec[2][2])/2;

open(OUTCFG,">",$filenameoutcfg);
$index=1;
printf OUTCFG $title;
printf OUTCFG "%10u", $config_key; # OUTCFG file key
printf OUTCFG "%10u", $periodic_key; # periodic boundary key
printf OUTCFG "\n%19.12e %19.12e %19.12e", $repx*$vec[0][0], $repy*$vec[0][1], $vec[0][2];
printf OUTCFG "\n%19.12e %19.12e %19.12e", $repx*$vec[1][0], $repy*$vec[1][1], $vec[1][2];
printf OUTCFG "\n%19.12e %19.12e %19.12e",       $vec[2][0],       $vec[2][1], $vec[2][2];

for($t=0;$t<$field_nummols;$t++) {
  for($m=0;$m<$mol_numents[$t];$m++) {
    for($x=0;$x<$repx;$x++) {
      for($y=0;$y<$repy;$y++) {
	for($a=0;$a<$mol_numatoms[$t];$a++) {
	  if($cdata[$t][$m][$a][3]>$zmin) {
	    printf OUTCFG "\n%-8s %9u", $cdata[$t][$m][$a][0], $index;
	    printf OUTCFG "\n%19.12e",$cdata[$t][$m][$a][1]+$x*$vec[0][0]+$y*$vec[1][0]+$dx;
	    printf OUTCFG  " %19.12e",$cdata[$t][$m][$a][2]+$x*$vec[0][1]+$y*$vec[1][1]+$dy;
	    printf OUTCFG  " %19.12e",$cdata[$t][$m][$a][3];
	    if($config_key>0) { #print velocity
	      printf OUTCFG "\n%19.12e",$cdata[$t][$m][$a][4];
	      printf OUTCFG  " %19.12e",$cdata[$t][$m][$a][5];
	      printf OUTCFG  " %19.12e",$cdata[$t][$m][$a][6];
	      if($config_key>1) { #print force
		$_=<CONFIG>;
		printf OUTCFG "\n%19.12e",$cdata[$t][$m][$a][7];
		printf OUTCFG  " %19.12e",$cdata[$t][$m][$a][8];
		printf OUTCFG  " %19.12e",$cdata[$t][$m][$a][9];
	      }
	    }
	    $index++;
	  }
	}
      }
    }
  }
}

close(OUTCFG);
open(FIELD, "<",$filenameinfld) or die "Can't open input FIELD-file:\n$!";
open(OUTFLD, ">", $filenameoutfld) or die "Can't open output FIELD-file:\n$!";
while(<FIELD>) {
  print OUTFLD $_;
  if(/molecules/ or /MOLECULES/) { last; }
}

for($t=0;$t<$field_nummols;$t++) {
  $_=<FIELD>; # molecule name
  print OUTFLD $_;
  $_=<FIELD>; # nummols
  print OUTFLD "NUMMOLS ".($mol_numents[$t]*$repx*$repy)."\n";
  while(<FIELD>) {
    print OUTFLD $_;
    if(/finish/ or /FINISH/) { last; }
  }
}
while(<FIELD>) {
  print OUTFLD $_;
}
close(FIELD, OUTFLD);