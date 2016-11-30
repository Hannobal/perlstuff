#!/usr/bin/perl

use Math::Trig;
use dlpoly_utility;
use POSIX qw(ceil floor);
# use strict;

$nmax = 1500;
$xmax = 100;
$ymax = 100;
$dminsq = 4;

$xmax2=$xmax/2;
$ymax2=$ymax/2;

read_field_file("FIELD",0);

$mol_numents[0][0]=0;

@{$cdata[0][0]}=();
for($i=0;$i<$nmax;$i++) {
  print "$i\n";
  newpos:while(1) {
    $pos[0]=rand($xmax);
    $pos[1]=rand($ymax);
    for($j=0;$j<@{$cdata[0][0]};$j++) {
      $dx = abs($cdata[0][0][$j][0][0] - $pos[0]);
      $dy = abs($cdata[0][0][$j][0][1] - $pos[1]);
      $dx-= $xmax if($dx>$xmax2);
      $dy-= $ymax if($dy>$ymax2);
      next newpos if($dx*$dx+$dy*$dy)<$dminsq;
    }
    last;
  };
  push(@{$cdata[0][0]},[[$pos[0],$pos[1],0,0,0,0,0,0,0,"X"]]);
  $mol_numents[0][0]++;
}

write_config_file("CONFIG",0,"title",">",0);

