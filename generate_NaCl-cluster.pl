#!/usr/bin/perl

use dlpoly_utility;

@nm = ('Na','Cl');
$imax = 6;
$dist = 5.6402;
$m=0;
for($z=0;$z<$imax;$z++) {
  for($y=0;$y<$imax;$y++) {
    for($x=0;$x<$imax;$x++) {
      $cdata[0][0][$m][0][0] = $x*$dist;
      $cdata[0][0][$m][0][1] = $y*$dist;
      $cdata[0][0][$m][0][2] = $z*$dist;
      $cdata[0][0][$m][0][9] = 'Na';
      $cdata[0][1][$m][0][0] = ($x+0.5)*$dist;
      $cdata[0][1][$m][0][1] = $y*$dist;
      $cdata[0][1][$m][0][2] = $z*$dist;
      $cdata[0][1][$m][0][9] = 'Cl';
      $m++;
      $cdata[0][0][$m][0][0] = ($x+0.5)*$dist;
      $cdata[0][0][$m][0][1] = ($y+0.5)*$dist;
      $cdata[0][0][$m][0][2] = $z*$dist;
      $cdata[0][0][$m][0][9] = 'Na';
      $cdata[0][1][$m][0][0] = $x*$dist;
      $cdata[0][1][$m][0][1] = ($y+0.5)*$dist;
      $cdata[0][1][$m][0][2] = $z*$dist;
      $cdata[0][1][$m][0][9] = 'Cl';
      $m++;
      $cdata[0][0][$m][0][0] = ($x+0.5)*$dist;
      $cdata[0][0][$m][0][1] = $y*$dist;
      $cdata[0][0][$m][0][2] = ($z+0.5)*$dist;
      $cdata[0][0][$m][0][9] = 'Na';
      $cdata[0][1][$m][0][0] = $x*$dist;
      $cdata[0][1][$m][0][1] = $y*$dist;
      $cdata[0][1][$m][0][2] = ($z+0.5)*$dist;
      $cdata[0][1][$m][0][9] = 'Cl';
      $m++;
      $cdata[0][0][$m][0][0] = $x*$dist;
      $cdata[0][0][$m][0][1] = ($y+0.5)*$dist;
      $cdata[0][0][$m][0][2] = ($z+0.5)*$dist;
      $cdata[0][0][$m][0][9] = 'Na';
      $cdata[0][1][$m][0][0] = ($x+0.5)*$dist;
      $cdata[0][1][$m][0][1] = ($y+0.5)*$dist;
      $cdata[0][1][$m][0][2] = ($z+0.5)*$dist;
      $cdata[0][1][$m][0][9] = 'Cl';
      $m++;
    }
  }
}
$field_title[0] = 'sodium chloride particle';
$field_units[0] = 'kj';
$field_nummols[0] = 2;
$mol_name[0][0] = 'sodium ion';
$mol_name[0][1] = 'chloride ion';
$mol_numents[0][0] = ($imax**3)*4;
$mol_numents[0][1] = ($imax**3)*4;
$mol_numatoms[0][0] = 1;
$mol_numatoms[0][1] = 1;
@{$mol_atomdata[0][0][0]}=('Na',22.98976928,2,1,0,0);
@{$mol_atomdata[0][1][0]}=('Cl',35.45,-2,1,0,0);
$field_numvdw[0]=1;
@{$field_vdwdata[0][0]}=('Na','Na','lj',0,1);

$periodic_key[0] = 0;
$config_key[0] = 2;
@cell = [[100,0,0],[0,100,0],[0,0,100]];

write_config_file("CONFIG",0,$field_title[0]);
write_field_file("FIELD",0);