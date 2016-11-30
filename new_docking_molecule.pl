#!/usr/bin/perl
use Math::Trig;
use Switch;

use dlpoly_utility;
use Storable qw(dclone);

if($#ARGV<6) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file with surface\n";
  print "2. name of corresponding FIELD-file\n";
  print "3. name of mol2 file with molecule to dock\n";
  print "4. distance to surface of new molecule\n";
  print "5. name of solvent molecule\n";
  print "6. name of output CONFIG-file\n";
  print "7. name of output FIELD-file\n";
  print "[8. velocity vector in A/ps]\n";
  print "[9. maximum angular deviation from normal for velocity in deg]\n";
  exit 1;
}

$filenameincfg  = $ARGV[0];
$filenameinfld  = $ARGV[1];
$filenamemol2   = $ARGV[2];
$surfdist       = $ARGV[3];
$solvname       = $ARGV[4];
$filenameoutcfg = $ARGV[5];
$filenameoutfld = $ARGV[6];
if($#ARGV>6) {$velocity = $ARGV[7]} else {$velocity=1};
if($#ARGV>7) {$maxangle = $ARGV[8]} else {$maxangle=20};

exit 1 if(read_mol2_file($filenamemol2,0)!=0);

# read FIELD file
exit 1 if(read_field_file($filenameinfld,0)!=0);
$foundmol = -1;
$foundsolv = -1;
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /$mol2_name[0]/ and not ($mol_name[0][$t] =~ /frozen/i or $mol_name[0][$t] =~ /fix/i)) {
    $foundmol=$t;
  } elsif($mol_name[0][$t] =~ /$solvname/) {
    $foundsolv=$t;
  }
}

if($foundmol<0) {
  print "**** error: molecule named '$mol2_name[0]' was not found in FIELD-file!\n";
  exit 1;
}

if($foundsolv<0) {
  print "**** warning: no molecule named '$solvname' was found in FIELD-file!\n";
}

exit 1 if(read_config_file($filenameincfg,0,0)!=0);


$zsurf=-9e20;
# freeze everything in the FIELD file
for($t=0;$t<$field_nummols[0];$t++) {
  next if($t==$foundmol);
  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
    $mol_atomdata[0][$t][$a][4] = 1;
  }
  # remove solvent
  if($t==$foundsolv) {
    @{$cdata[0][$t]} = ();
    $mol_numents[0][$t] = 0;
  }
  if($mol_name[0][$t] =~ /hydroxide/i) {
    for($m=0;$m<$mol_numents[0][$t];$m++) {
      $zsurf=$cdata[0][$t][$m][0][2] if($zsurf<$cdata[0][$t][$m][0][2]);
    }
  }
}

if($zsurf==-9e20) {
  print "**** error: could not determine surface!\n";
  exit 1;
}

####### place new molecule ########
$m = 0;
unshift(@{$cdata[0][$foundmol]},undef);
@{$cdata[0][$foundmol][$m]} = mol2_to_cdata(0);
$mol_numents[0][$foundmol]++;
# rotate around z-axis
@rotmatrix = &gen_rot_matrix([0,0,1], rand(2*pi));
rotate_molecule(\@{$cdata[0][$foundmol][$m]}, \@rotmatrix);
# determine molecule boundaries
@mol2_minpos = ( 9e20,  9e20,  9e20);
@mol2_maxpos = (-9e20, -9e20, -9e20);
for($a=0;$a<$mol2_numatoms[0];$a++) {
  for($c=0;$c<3;$c++) {
    $mol2_maxpos[$c]=$mol2_atomdata[0][$a][$c] if($mol2_atomdata[0][$a][$c]>$mol2_maxpos[$c]);
    $mol2_minpos[$c]=$mol2_atomdata[0][$a][$c] if($mol2_atomdata[0][$a][$c]<$mol2_minpos[$c]);
  }
}
# shift to random position
$shiftz = $zsurf + $surfdist - $mol2_minpos[2];
# @shiftvec = (rand(2*$size[0][0])-$size[0][0], rand(2*$size[0][1])-$size[0][1], $shiftz);
@shiftvec = (0,0,$shiftz);
move_mol(\@{$cdata[0][$foundmol][$m]},\@shiftvec);
# add velocity to the molecule
@velvec = (0,0,-$velocity);
@rotmatrix = gen_rot_matrix([1,0,0], rand($maxangle/180*pi)); # deviation from normal
@velvec = rotate_vector(\@velvec, \@rotmatrix);
@rotmatrix = gen_rot_matrix([0,0,1], rand(2*pi)); # rotate around z-axis
@velvec = rotate_vector(\@velvec, \@rotmatrix);
add_mol_velocity(\@{$cdata[0][$foundmol][$m]},\@velvec);

write_config_file($filenameoutcfg,0,$config_title[0]);
write_field_file($filenameoutfld,0);
exit 0;