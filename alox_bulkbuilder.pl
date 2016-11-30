#!/usr/bin/perl

# repeats the aluminum oxide surface from a DL_POLY CONFIG file by
# integer numbers of the periodicity vectors to yield a larger slab

if($#ARGV<4) {
  print "1. dl_poly REVCON-file\n";
  print "2. output-file with ending\n";
  print "3. repetitions in x-, y- and z-direction (integer)\n";
  exit;
}
if($ARGV[2]<1 or $ARGV[3]<1 or $ARGV[4]<1) {
  print "repetitions must be >= 1\n";
  exit;
}
$rep[0]=$ARGV[2];
$rep[1]=$ARGV[3];
$rep[2]=$ARGV[4];
open(REVCON, "<", $ARGV[0]) or die "Can't open REVCON-file: $!";
open(OUTPUT, ">", $ARGV[1]) or die "Can't open output file: $!";

# indices for data array:
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

$title=<REVCON>;
$_=<REVCON>;
($REVCONkey,$periodic_key) = /^\s*(\S+)\s+(\S+)/;
if($periodic_key!=0) { # skip lines for box size if necessary
  for($j=0;$j<=2;$j++) {
    $_=<REVCON>;
    ($vec[$j][0], $vec[$j][1], $vec[$j][2]) = /\s*(\S+)\s+(\S+)\s+(\S+)/;
  }
} else {
  print "error: periodic information is missing in CONFIG file!";
  exit;
}

$numold=0;
while(<REVCON>) {
  ($data[$numold][0]) = /^\s*(\S*)/;
  $_=<REVCON>;
  ($data[$numold][1],$data[$numold][2],$data[$numold][3]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
  if($REVCONkey>0) {
    $_ = <REVCON>;
    ($data[$numold][4],$data[$numold][5],$data[$numold][6]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
    if($REVCONkey>1) {
      $_ = <REVCON>;
      ($data[$numold][7],$data[$numold][8],$data[$numold][9]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
    }
  }
  $numold++;
}
$numnew=0;
@maxpos=(0,-99999,-99999,-99999);
@minpos=(0, 99999, 99999, 99999);
for($i=0;$i<$numold;$i++){
  for($z=0;$z<$rep[2];$z++) {
    for($y=0;$y<$rep[1];$y++) {
      for($x=0;$x<$rep[0];$x++) {
	$newdata[$numnew][0]=$data[$i][0];
	$newdata[$numnew][1]=$data[$i][1]+$x*$vec[0][0]+$y*$vec[1][0]+$z*$vec[2][0];
	$newdata[$numnew][2]=$data[$i][2]+$x*$vec[0][1]+$y*$vec[1][1]+$z*$vec[2][1];
	$newdata[$numnew][3]=$data[$i][3]+$x*$vec[0][2]+$y*$vec[1][2]+$z*$vec[2][2];
	$newdata[$numnew][4]=$data[$i][4];
	$newdata[$numnew][5]=$data[$i][5];
	$newdata[$numnew][6]=$data[$i][6];
	$newdata[$numnew][7]=$data[$i][7];
	$newdata[$numnew][8]=$data[$i][8];
	$newdata[$numnew][9]=$data[$i][9];
	for($j=1;$j<3;$j++){
	  if($newdata[$newnum][$j]>$maxpos[$j]){
	    $maxpos[$j]=$newdata[$newnum][$j];
	  }
	  if($newdata[$newnum][$j]<$minpos[$j]){
	    $minpos[$j]=$newdata[$newnum][$j];
	  }
	}
	$numnew++;
      }
    }  
  }
}
@newdata = sort { $a->[3] <=> $b->[3] } @newdata; #sort by z-coordiante
@newdata = sort { $a->[0] cmp $b->[0] } @newdata; #sort by name

printf OUTPUT "%-80s\n", "title";
printf OUTPUT "%10u", "2"; # OUTPUT file key
printf OUTPUT "%10u\n", "6"; # periodic boundary key
printf OUTPUT "%20.12e%20.12e%20.12e\n", $rep[0]*$vec[0][0], $rep[0]*$vec[0][1], $rep[0]*$vec[0][2];
printf OUTPUT "%20.12e%20.12e%20.12e\n", $rep[1]*$vec[1][0], $rep[1]*$vec[1][1], $rep[1]*$vec[1][2];
printf OUTPUT "%20.12e%20.12e%20.12e\n", $rep[2]*$vec[2][0], $rep[2]*$vec[2][1], $rep[2]*$vec[2][2];
$j=131;
for($i=0;$i<$numnew;$i++){
  printf OUTPUT "%-8s%10u\n" ,$newdata[$i][0],$i+1;
  printf OUTPUT "%20.12f%20.12f%20.12f\n",$newdata[$i][1],$newdata[$i][2],$newdata[$i][3]; # xyz coordiantes
  printf OUTPUT "%20.12f%20.12f%20.12f\n",$newdata[$i][4],$newdata[$i][5],$newdata[$i][6]; # xyz of velocity
  printf OUTPUT "%20.12f%20.12f%20.12f\n",$newdata[$i][7],$newdata[$i][8],$newdata[$i][9]; # xyz of force
}
close(REVCON, OUTPUT, INFO);