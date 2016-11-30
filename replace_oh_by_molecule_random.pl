#!/usr/bin/perl
use Math::Trig;
use Switch;

$doo=2.755333333;
$tolerance=0.2;
$nmaxtrys = 5;

#          x         
#         A          
#        /           
#     b/             
#     /              
#   /                
#  /        a        
# x---------------->x

if($#ARGV<6) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file with surface\n";
  print "2. name of corresponding FIELD-file\n";
  print "3. name of mol2-file with molecule\n";
  print "4. name of output CONFIG-file\n";
  print "5. name of output FIELD-file\n";
  print "6. desired density of molecules per square nanometer (max 5)\n";
  print "7. minimum distance between molecules\n";
  print "[8. (zero/mono/bi/tri)-dentat binding]\n";
  print "[9. additional distance to surface]\n";
  exit 1;
}

$targetdensity = $ARGV[5];
$molmindist    = $ARGV[6];

if(lc(substr($ARGV[7],0,4)) eq "zero") {
  $bindingmode=0;
} elsif(lc(substr($ARGV[7],0,4)) eq "mono") {
  $bindingmode=1;
} elsif(lc(substr($ARGV[7],0,2)) eq "bi") {
  $bindingmode=2;
} elsif(lc(substr($ARGV[7],0,3)) eq "tri") {
  $bindingmode=3;
} else {
  $bindingmode=undef;
}
if($#ARGV>7){ $surfdist=$ARGV[8]; } else { $surfdist=0; }

$filenameincfg  = $ARGV[0];
$filenameinfld  = $ARGV[1];
$filenamemol2   = $ARGV[2];
$filenameoutcfg = $ARGV[3];
$filenameoutfld = $ARGV[4];
# $mol2name=$ARGV[2];
# $mol2name =~ s/\.\S*$//; # remove file extension
# $mol2name = substr($mol2name,rindex($mol2name,'/')+1); # remove path (everything before last "/")

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
$_=<MOL2>;
$mol2name=<MOL2>;
$mol2name=~s/^\s+//;$mol2name=~s/\s+$//;
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

$area=($vec[0][0]*$vec[1][1]-$vec[1][0]*$vec[0][1])/100;
$ohremoved=0;

$imol    = 0;
$density = 0;
$ntrys   = 0;
$best_density = 0;
@candidatecopy = @candidate;

if($bindingmode==1) {
  while($density<$targetdensity) {
    if($density > $best_density) {
      $best_density     = $density;
      @best_molpos      = @molpos;
      @best_removedmols = @removedmols;
      $best_ohremoved   = $ohremoved;
    }
    if($#candidate<0) {
      if($ntrys>$nmaxtrys) {
	print "**** error: desired density could not be reached within $nmaxtrys trys\n";
	last;
      }
      print "**** faild in iteration $ntrys out of $nmaxtrys (density=$density, best=$best_density)\n";
      @molpos      = ();
      @removedmols = ();
      $density     = 0;
      $imol        = 0;
      $ohremoved   = 0;
      @candidate   = @candidatecopy;
      $ntrys++;
    }
    $ioh = int(rand($#candidate+1));
    $molpos[$imol][0] = $candidate[$ioh][1];
    $molpos[$imol][1] = $candidate[$ioh][2];
    $molpos[$imol][2] = $candidate[$ioh][3];
    $ohremoved++;
    @removedmols[$#removedmols+1] = splice(@candidate,$ioh,1);
    for($j=0;$j<=$#candidate;$j++) {
      if($j<0 or $imol<0) {print "$imol $j\n";}
      $deleted = 0;
      for($x=-1;$x<2;$x++) {
	for($y=-1;$y<2;$y++) {
	  $dx=$molpos[$imol][0]-$candidate[$j][1]+$x*$vec[0][0]+$y*$vec[1][0]+$mol2oxpos[0][0];
	  $dy=$molpos[$imol][1]-$candidate[$j][2]+$x*$vec[0][1]+$y*$vec[1][1]+$mol2oxpos[0][1];
	  $dist=sqrt($dx**2+$dy**2);
	  if($dist<$molmindist) {
	    splice(@candidate,$j,1);
	    $deleted=1;
	    $j--;
	    last;
	  }
	}
	last if($deleted);
      }
    }
    $density=$imol/$area;
    $imol++;
  }
  $density = $best_density;
  $mol_deltanuments[$foundmol2name] =  $best_ohremoved;
  $mol_deltanuments[$foundhydrox]   = -$best_ohremoved;
  $ohremoved = $best_ohremoved;
} else {
  print "**** error: algorithm for bindingmode $bindingmode is not implemented yet\n";
  exit 1;
}

@removedmols = sort { $a->[0] <=> $b->[0] } @best_removedmols;
@molpos = @best_molpos;

print "replaced $ohremoved hydroxides by $mol_deltanuments[$foundmol2name] molecules in $bindingmode"."dentate coordination\n";
# print "surface area: $area squarenaometers\n";
print "$density molecules per square nanometer\n";

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
