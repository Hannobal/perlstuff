#!/usr/bin/perl

# adds a new molecule in a random position above the topmost atom in a DL_POLY
# CONFIG file and assigns a speed vector to the molecule so that the randomly
# assigned binding atom approaches the ranomly selected binding site

use Math::Trig;
use Switch;
use strict vars;

my $zmin = 6;
my $vscale = 10;
our $kB  = 1.3806488e-23; #Boltzmann constant in J/K
our $u   = 1.660538921e-27; #unified atomic mass unit in kg
my ($i, $j, $k, $t, $m, $at, $c);
my ($mol2name, $mol2charge, $mol2numox, @mdata, @mol2ox);
my ($title, @config_header, $periodic_key, $config_key, @cdata, @cellvec, @size);
my (@mol_atomcharge, @mol_atommass, @mol_atomname, @mol_deltanuments, $temperature);
my ($field_nummols, @mol_mass, @mol_name, @mol_charge, @mol_numatoms, @mol_numents);
my ($xmin, $ymin, $molindex, $ztop, @pos, @rot_matrix, @molpos, $mzmin, $moldist, $distance, @vector);
my ($thermal_energy, @direction, $bindingatom, $bindingsite, $angle, @unityvector, $veclen, $index);
my ($candindex, @candidate, $foundmol2name, $avvelocity, $lrotation);

if($#ARGV<4) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file with surface\n";
  print "2. name of corresponding FIELD-file\n";
  print "3. name of mol2-file with molecule\n";
  print "4. name of output CONFIG-file\n";
  print "5. name of output FIELD-file\n";
  print "[6. minimum distance between topmost atom and molecule in Angstroms]\n";
  print "[7. temperature of added molecule in K]\n";
  print "[8. toggle random rotation on/off (default:on)]\n";
  exit 1;
}

if($#ARGV>4){ $moldist=$ARGV[5]; } else { $moldist=10; }
if($#ARGV>5){ $temperature=$ARGV[6]; } else { $temperature=300; }
if($#ARGV>6){
  if($ARGV[7]=~/^y/i or $ARGV[7]=~/^on/i) {
    $lrotation=1;
  } elsif($ARGV[7]=~/^n/i or $ARGV[7]=~/^off/i) {
    $lrotation=0;
  } else {
    print "**** error: input for \"toggle random rotation on/off\" could not be understood\n";
    exit 1;
  }
} else { $lrotation=1; }
my $filenameincfg  = $ARGV[0];
my $filenameinfld  = $ARGV[1];
my $filenamemol2   = $ARGV[2];
my $filenameoutcfg = $ARGV[3];
my $filenameoutfld = $ARGV[4];
my $mol2name       = $ARGV[2];
$mol2name = substr($mol2name,0,rindex($mol2name,'.')); # remove file type
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
  my $string = $_;
  if(/BOND/) {
    last;
  } else {
    ($mdata[$i][1],$mdata[$i][2],$mdata[$i][3],$mdata[$i][0],$mdata[$i][4]) =
    /^\s*\d+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)/; # position, type and charge
    $mol2charge += $mdata[$i][4];
    if(uc($mdata[$i][0]) eq "O" or uc($mdata[$i][0]) eq "OH") {
      $mol2ox[$mol2numox] = $i;
      $mol2numox++;
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
      ($cellvec[$j][0], $cellvec[$j][1], $cellvec[$j][2]) = /\s*(\S*)\s*(\S*)\s*(\S*)/;
      $size[$j]=$cellvec[$j][$j]/2;
    }
  } case 0 {
    print "error: periodic conditions needed for calculation\n";
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

$xmin=9999; # to find hydroxide in lower left corner
$ymin=9999; # to find hydroxide in lower left corner
$foundmol2name=-1;
$molindex=1;
for($t=0;$t<$field_nummols;$t++) {
  $j = 0; #number of atoms counted for the molecule
  $_=<FIELD>;
  $mol_name[$t]=$_;
  $mol_name[$t] =~ s/^\s+//;
  $mol_name[$t] =~ s/\s+$//;
  $mol_charge[$t] = 0;
  $mol_mass[$t]   = 0;
  if($mol_name[$t] eq $mol2name) {
    $foundmol2name = $t;
  }
  $_=<FIELD>;
  ($mol_numents[$t]) = /^\s*\S+\s+(\S+)/;
  $_=<FIELD>;
  ($mol_numatoms[$t]) = /^\s*\S+\s+(\S+)/;
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
    $mol_mass[$t]   += $mol_atommass[$t][$j];
    $mol_charge[$t] += $mol_atomcharge[$t][$j];
#     printf "%8s %3.5f %3.5f\n", $mol_atomname[$t][$j],$mol_atommass[$t][$j],$mol_atomcharge[$t][$j];
    $j++;
  }
  if($j!=$mol_numatoms[$t]) {
    print "error: ".$mol_numatoms[$t]." expected for molecule ".$mol_name[$t]." but ".$j." were found\n";
    exit 1;
  }
  # read atom coordinates for the molecule from the CONFIG file
  for($m=0;$m<$mol_numents[$t];$m++) {
    for($at=0;$at<$mol_numatoms[$t];$at++) {
      $_=<CONFIG>;
      ($cdata[$t][$m][$at][0]) = /^\s*(\S*)/;
      if($cdata[$t][$m][$at][0] ne $mol_atomname[$t][$at]) {
	print "error: atom name ".$cdata[$t][$m][$at][0]." in CONFIG-file $filenameincfg does not match\n";
	print "name ".$mol_atomname[$t][$at]." in FIELD file $filenameinfld for molecule ".$mol_name[$t]."\n";
	exit 1;
      }
      $_=<CONFIG>;
      ($cdata[$t][$m][$at][1],$cdata[$t][$m][$at][2],$cdata[$t][$m][$at][3]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
      if($config_key>0) { #read velocity
	$_=<CONFIG>;
	($cdata[$t][$m][$at][4],$cdata[$t][$m][$at][5],$cdata[$t][$m][$at][6]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
	if($config_key>1) { #read force
	  $_=<CONFIG>;
	  ($cdata[$t][$m][$at][7],$cdata[$t][$m][$at][8],$cdata[$t][$m][$at][9]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
	}
      }
      $cdata[$t][$m][$at][10]=$molindex;
    }
    $molindex++;
  }
}
close(CONFIG, FIELD);

if($foundmol2name==-1) {
  print "error: could not find entry for molecule $mol2name in FIELD file $filenameinfld!\n";
  exit 1;
}

# find maximum z value and candidates for docking
$candindex=0;
$ztop = -9999;
for($t=0;$t<$field_nummols;$t++) {
  $mol_deltanuments[$t]=0;
  for($m=0;$m<$mol_numents[$t];$m++) {
    for($at=0;$at<$mol_numatoms[$t];$at++) {
      if($cdata[$t][$m][$at][3]>$ztop) {
	$ztop=$cdata[$t][$m][$at][3];
      }
      if((uc($cdata[$t][$m][$at][0]) eq "O")
	  and $cdata[$t][$m][$at][3]>$zmin) {
# 	if($cdata[$t][$m][$at][1]<$xmin) {
# 	  $xmin=$cdata[$t][$m][$at][1];
# 	}
# 	if($cdata[$t][$m][$at][2]<$ymin) {
# 	  $ymin=$cdata[$t][$m][$at][2];
# 	}
	$candidate[$candindex][0]=$cdata[$t][$m][$at][10];
	$candidate[$candindex][1]=$cdata[$t][$m][$at][1];
	$candidate[$candindex][2]=$cdata[$t][$m][$at][2];
	$candidate[$candindex][3]=$cdata[$t][$m][$at][3];
	$candidate[$candindex][4]=$t;
	$candidate[$candindex][5]=$m;
	$candindex++;
      }
    }
  }
}

# rotate molecule by random angle around x- y and z-axis
if($lrotation) {
  @unityvector = ([1,0,0],[0,1,0],[0,0,1]);
  for($c=0;$c<3;$c++) {
    $angle = rand(2*pi);
    @rot_matrix = &gen_rot_matrix($unityvector[$c], $angle);
    for($at=0;$at<$mol_numatoms[$foundmol2name];$at++) {
      @pos = ($mdata[$at][1],$mdata[$at][2],$mdata[$at][3]);
      ($mdata[$at][1],$mdata[$at][2],$mdata[$at][3]) = &rotate(\@rot_matrix, \@pos);
    }
  }
}
# determine minimum z-position of molelcule
$mzmin=9999;
for($at=0;$at<$mol_numatoms[$foundmol2name];$at++) {
  if($mdata[$at][3]<$mzmin) {
    $mzmin=$mdata[$at][3];
  }
}
$bindingsite   = int(rand($#candidate+1));
$bindingatom   = int(rand($mol2numox));
$molpos[2]     = $ztop+$moldist-$mzmin;#-$candidate[$bindingsite][3];
$molpos[1]     = 2*rand($size[1])-$size[1];#-$candidate[$bindingsite][2];
$molpos[0]     = 2*rand($size[0])-$size[0];#-$candidate[$bindingsite][1];
for($at=0;$at<$mol_numatoms[$foundmol2name];$at++) {
  $mdata[$at][1] += $molpos[0];
  $mdata[$at][2] += $molpos[1];
  $mdata[$at][3] += $molpos[2];
  for($c=0;$c<2;$c++) {
    if($mdata[$at][$c+1] < -$size[$c]) {
      $mdata[$at][$c+1] += 2*$size[$c];
    } elsif($mdata[$at][$c+1] > $size[$c]) {
      $mdata[$at][$c+1] -= 2*$size[$c];
    }
  }
}
$direction[0]  = $candidate[$bindingsite][1]-$mdata[$mol2ox[$bindingatom]][1];
$direction[1]  = $candidate[$bindingsite][2]-$mdata[$mol2ox[$bindingatom]][2];
$direction[2]  = $candidate[$bindingsite][3]-$mdata[$mol2ox[$bindingatom]][3];
$distance      = sqrt($direction[0]**2 + $direction[1]**2 + $direction[2]**2);
$direction[0] /= $distance;
$direction[1] /= $distance;
$direction[2] /= $distance;

$thermal_energy = 3/2*$kB*$temperature;
$avvelocity=sqrt(2*$thermal_energy/($mol_mass[$foundmol2name]*$u))*1e-2*$vscale; # in Angstroms/ps
for($at=0;$at<$mol_numatoms[$foundmol2name];$at++) {
  $mdata[$at][5]=$avvelocity*$direction[0];
  $mdata[$at][6]=$avvelocity*$direction[1];
  $mdata[$at][7]=$avvelocity*$direction[2];
}
# @vector = (rand(2)-1,rand(2)-1,rand(2)-1);
# $veclen = sqrt($vector[0]**2+$vector[1]**2+$vector[2]**2);
# $vector[0] /= $veclen;
# $vector[1] /= $veclen;
# $vector[2] /= $veclen;

$mol_deltanuments[$foundmol2name]++;

# open(TEST,">","test.xyz");
# print TEST "".(2+$mol_numatoms[$foundmol2name])."\n";
# for($at=0;$at<$mol_numatoms[$foundmol2name];$at++) {
#   printf  TEST "\n%4s %10.5f %10.5f %10.5f",$mdata[$at][4],$mdata[$at][1],$mdata[$at][2],$mdata[$at][3];
#   printf  TEST " %10.5f %10.5f %10.5f",$mdata[$at][5],$mdata[$at][6],$mdata[$at][7];
# }
#   printf  TEST "\n%4s %10.5f %10.5f %10.5f",$candidate[$bindingsite][0],$candidate[$bindingsite][1],$candidate[$bindingsite][2],$candidate[$bindingsite][3];
#   printf  TEST " %10.5f %10.5f %10.5f",0,0,0;
#   printf  TEST "\n%4s %10.5f %10.5f %10.5f","C",$mdata[$mol2ox[$bindingatom]][1],$mdata[$mol2ox[$bindingatom]][2],$mdata[$mol2ox[$bindingatom]][3];
#   printf  TEST " %10.5f %10.5f %10.5f",0,0,0;


open(OUTCFG,">",$filenameoutcfg);
$index=1;
printf OUTCFG $title;
printf OUTCFG "%10u", 2; # OUTCFG file key
printf OUTCFG "%10u", $periodic_key; # periodic boundary key
printf OUTCFG "\n%19.12e %19.12e %19.12e", $cellvec[0][0], $cellvec[0][1], $cellvec[0][2];
printf OUTCFG "\n%19.12e %19.12e %19.12e", $cellvec[1][0], $cellvec[1][1], $cellvec[1][2];
printf OUTCFG "\n%19.12e %19.12e %19.12e", $cellvec[2][0], $cellvec[2][1], $cellvec[2][2];

$i=0; # index for @removedmols
$index=1; # index for molecules in CONFIG
for($t=0;$t<$field_nummols;$t++) {
  if($t==$foundmol2name) {
    for($at=0;$at<=$#mdata;$at++) {
      printf OUTCFG "\n%-8s %9u", uc($mdata[$at][0]), $index;
      printf OUTCFG "\n%19.12e",$mdata[$at][1];
      printf OUTCFG  " %19.12e",$mdata[$at][2];
      printf OUTCFG  " %19.12e",$mdata[$at][3];
      printf OUTCFG "\n%19.12e",$mdata[$at][5];
      printf OUTCFG  " %19.12e",$mdata[$at][6];
      printf OUTCFG  " %19.12e",$mdata[$at][7];
      printf OUTCFG " \n%19.12e %19.12e %19.12e", 0,0,0;
      $index++;
    }
  }
  for($m=0;$m<$mol_numents[$t];$m++) {
    for($at=0;$at<$mol_numatoms[$t];$at++) {
      printf OUTCFG "\n%-8s %9u", $cdata[$t][$m][$at][0], $index;
      printf OUTCFG "\n%19.12e",$cdata[$t][$m][$at][1];
      printf OUTCFG  " %19.12e",$cdata[$t][$m][$at][2];
      printf OUTCFG  " %19.12e",$cdata[$t][$m][$at][3];
      printf OUTCFG "\n%19.12e",$cdata[$t][$m][$at][4];
      printf OUTCFG  " %19.12e",$cdata[$t][$m][$at][5];
      printf OUTCFG  " %19.12e",$cdata[$t][$m][$at][6];
      printf OUTCFG "\n%19.12e",$cdata[$t][$m][$at][7];
      printf OUTCFG  " %19.12e",$cdata[$t][$m][$at][8];
      printf OUTCFG  " %19.12e",$cdata[$t][$m][$at][9];
      $index++;
    }
  }
}
close(OUTCFG);

open(FIELD, "<", $filenameinfld) or die "Can't open FIELD-file:\n$!";
open(OUTFLD, ">", $filenameoutfld) or die "Can't open FIELD-file:\n$!";

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

sub gen_rot_matrix {
  # args: axis (3 dim vector), angle (radian)
  my @axis = @{$_[0]};
  my $angle = $_[1];
  my @rot_matrix;
  if($#axis == 2) {
    # 3 dimensional rotation matrix
    $rot_matrix[0][0] = cos($angle) + ($axis[0]**2)*(1 - cos($angle));
    $rot_matrix[1][0] = $axis[1]*$axis[0]*(1 - cos($angle)) + $axis[2]*sin($angle);
    $rot_matrix[2][0] = $axis[2]*$axis[0]*(1 - cos($angle)) - $axis[1]*sin($angle);
    $rot_matrix[0][1] = $axis[0]*$axis[1]*(1 - cos($angle)) - $axis[2]*sin($angle);
    $rot_matrix[1][1] = cos($angle) + ($axis[1]**2)*(1 - cos($angle));
    $rot_matrix[2][1] = $axis[2]*$axis[1]*(1 - cos($angle)) + $axis[0]*sin($angle);
    $rot_matrix[0][2] = $axis[0]*$axis[2]*(1 - cos($angle)) + $axis[1]*sin($angle);
    $rot_matrix[1][2] = $axis[1]*$axis[2]*(1 - cos($angle)) - $axis[0]*sin($angle);
    $rot_matrix[2][2] = cos($angle) + ($axis[2]**2)*(1 - cos($angle));
    my $rows1 = @rot_matrix;
    my $cols1 = @{$rot_matrix[0]};
    return @rot_matrix;
  } elsif ($#axis == 1) {
    # 2 dimensional rotation matrix
    $rot_matrix[0][0] =   cos($angle);
    $rot_matrix[0][1] = - sin($angle);
    $rot_matrix[1][0] =   sin($angle);
    $rot_matrix[1][1] =   cos($angle);
  } else {
    print "error in subroutine gen_rot_matrix: axis vector has @axis dimensions instead of 2 or 3!\n";
    exit 1;
  }
}

sub rotate {
  my @matrix = @{$_[0]};
  my @vector = @{$_[1]};
  my @result;
  my $rows1 = @matrix;
  my $cols1 = @{$matrix[0]};
  my $rows2 = @vector;
  my $cols2 = @{$vector[0]};
  if($cols1 != $rows2) {
    print "error in subroutine rotate: dimensions of matrix don't match vector!\n";
    print "r1 $rows1 c1 $cols1 r2 $rows2 c2 $cols2\n";
    exit 1;
  }
  if($cols2 != 0) {
    print "error in subroutine rotate: vector must not be a matrix!\n";
    exit 1;
  }
  $result[0] = $vector[0]*$matrix[0][0] + $vector[1]*$matrix[0][1] + $vector[2]*$matrix[0][2];
  $result[1] = $vector[0]*$matrix[1][0] + $vector[1]*$matrix[1][1] + $vector[2]*$matrix[1][2];
  $result[2] = $vector[0]*$matrix[2][0] + $vector[1]*$matrix[2][1] + $vector[2]*$matrix[2][2];
  return @result;
}

sub isnumber {
  for(my $i=0; $i<@_; $i++) {
    if($_[$i]*1 ne $_[$i]) {
      return 0;
    }
  }
  return 1;
}
