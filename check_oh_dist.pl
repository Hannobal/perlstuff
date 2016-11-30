#!/usr/bin/perl

if($#ARGV<1) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file\n";
  print "2. minimum distance\n";
  exit;
}

$dmax=$ARGV[1];
if(not ($dmax*1) eq $dmax) {
  print "error: $dmax is not a number\n";
  exit;
}

if($#ARGV>1) {
  if($ARGV[2] =~ /^-[nN]/) {
    $periodic=0;
  } else {
    print "error: could not interpret input $ARGV[2]\n";
  }
} else {
  $periodic=1;
}

$anzges = 0; #number of atoms counted in CONFIG
open(CONFIG, "<", $ARGV[0]) or die "Can't open CONFIG-File: $!";

$title=<CONFIG>;
$_=<CONFIG>;
($CONFIGkey,$periodic_key) = /^\s*(\S+)\s+(\S+)/;
if($periodic_key!=0) { # skip lines for box size if necessary
  for($i=0;$i<=2;$i++) {
    $_=<CONFIG>;
    ($vec[$i][0], $vec[$i][1], $vec[$i][2]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
  }
}
while(<CONFIG>) {
  ($name[$anzges],$index[$anzges]) = /^\s*(\S+)\d*\s+(\S+)/;
#   print $name[$anzges]."\n";
  $_=<CONFIG>;
  ($posx[$anzges],$posy[$anzges],$posz[$anzges]) = /\s*(\S*)\s*(\S*)\s*(\S*)/;
  for($i=0;$i<$CONFIGkey;$i++) { #skip lines for speed and force if necessary
    $_ = <CONFIG>;
  }
  $anzges++;
}

$maxdist[0]=-9999999;
for($i=0;$i<$anzges;$i++){
  if($name[$i] eq "OX") {
    $dmin=99999;
    for($x=-1;$x<2;$x++){
      for($y=-1;$y<2;$y++){
	$j=$i+1;
	if($name[$j] ne "HG") {
	  print "error for atom ".($i+1).": HG expected!";
	  exit;
	}
	$dx=$posx[$i]-$posx[$j]+$x*$vec[0][0]+$y*$vec[1][0]+$z*$vec[2][0];
	$dy=$posy[$i]-$posy[$j]+$x*$vec[0][1]+$y*$vec[1][1]+$z*$vec[2][1];
	$dz=$posz[$i]-$posz[$j]+$x*$vec[0][2]+$y*$vec[1][2]+$z*$vec[2][2];
	$dist=sqrt($dx**2+$dy**2+$dz**2);
	if($dist<$dmin) {
	  $dmin=$dist;
	}
      }
    }
    if($dmin>$dmax) {
      printf "%8s %7u %8s %7u %10.5f\n",$name[$i],($i+1),$name[$j],($j+1),$dmin;
    }
    if($dmin>$maxdist[0]) {
      $maxdist[0]=$dmin;
      $maxdist[1]=$i;
      $maxdist[2]=$j;
    }
  }
}

print "minimum distance:";
printf "%8s %7u %8s %7u %10.5f\n",$name[$maxdist[1]],($maxdist[1]+1),$name[$maxdist[2]],($maxdist[2]+1),$maxdist[0];
