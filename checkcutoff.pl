#!/usr/bin/perl

if($#ARGV<3) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file\n";
  print "2. cutoff in Angstroms\n";
  print "3. index of atom\n";
  print "4. name of output xyz-file\n";
  exit;
}

$cutoff=$ARGV[1];
if(not ($cutoff*1) eq $cutoff) {
  print "error: $cutoff is not a number\n";
  exit;
}

$center=$ARGV[2]-1;

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
  ($data[$anzges][0],$data[$anzges][4]) = /^\s*(\S+)\d*\s+(\S+)/;
#   print $name[$anzges]."\n";
  $_=<CONFIG>;
  ($data[$anzges][1],$data[$anzges][2],$data[$anzges][3]) = /\s*(\S*)\s*(\S*)\s*(\S*)/;
  for($i=0;$i<$CONFIGkey;$i++) { #skip lines for speed and force if necessary
    $_ = <CONFIG>;
  }
  $anzges++;
}

if($anzges-1<$center) {
  print "error: index $center is out of range of atoms!";
  exit;
}

if($periodic_key==2) {
  for($x=-1;$x<2;$x++) {
    for($y=-1;$y<2;$y++) {
      for($z=-1;$z<2;$z++) {
	for($i=0;$i<$anzges;$i++){
	  if($i != $center) {
	    $dx=$data[$center][1]-$data[$i][1]+$x*$vec[0][0]+$y*$vec[1][0]+$z*$vec[2][0];
	    $dy=$data[$center][2]-$data[$i][2]+$x*$vec[0][1]+$y*$vec[1][1]+$z*$vec[2][1];
	    $dz=$data[$center][3]-$data[$i][3]+$x*$vec[0][2]+$y*$vec[1][2]+$z*$vec[2][2];
	    $dist=sqrt($dx**2+$dy**2+$dz**2);
	    if($dist<=$cutoff) {
	      push(@inside,$data[$i]);
	      $inside[$#inside][1]-=$x*$vec[0][0]+$y*$vec[1][0]+$z*$vec[2][0];
	      $inside[$#inside][2]-=$x*$vec[0][1]+$y*$vec[1][1]+$z*$vec[2][1];
	      $inside[$#inside][3]-=$x*$vec[0][2]+$y*$vec[1][2]+$z*$vec[2][2];
	    }
	  }
	}
      }
    }
  }
} elsif($periodic_key==6) {
  for($x=-1;$x<2;$x++) {
    for($y=-1;$y<2;$y++) {
      for($i=0;$i<$anzges;$i++){
	if($i != $center) {
	  $dx=$data[$center][1]-$data[$i][1]+$x*$vec[0][0]+$y*$vec[1][0]+$z*$vec[2][0];
	  $dy=$data[$center][2]-$data[$i][2]+$x*$vec[0][1]+$y*$vec[1][1]+$z*$vec[2][1];
	  $dz=$data[$center][3]-$data[$i][3]+$x*$vec[0][2]+$y*$vec[1][2]+$z*$vec[2][2];
	  $dist=sqrt($dx**2+$dy**2+$dz**2);
	  if($dist<=$cutoff) {
	    push(@inside,$data[$i]);
	    $inside[$#inside][1]-=$x*$vec[0][0]+$y*$vec[1][0];
	    $inside[$#inside][2]-=$x*$vec[0][1]+$y*$vec[1][1];
	    $inside[$#inside][3]-=$x*$vec[0][2]+$y*$vec[1][2];
	  }
	}
      }
    }
  }
} elsif($periodic_key==0) {
  for($i=0;$i<$anzges;$i++){
    if($i != $center) {
      $dx=$data[$center][1]-$data[$i][1]+$x*$vec[0][0]+$y*$vec[1][0]+$z*$vec[2][0];
      $dy=$data[$center][2]-$data[$i][2]+$x*$vec[0][1]+$y*$vec[1][1]+$z*$vec[2][1];
      $dz=$data[$center][3]-$data[$i][3]+$x*$vec[0][2]+$y*$vec[1][2]+$z*$vec[2][2];
      $dist=sqrt($dx**2+$dy**2+$dz**2);
      if($dist<=$cutoff) {
	push(@inside,$data[$i]);
      }
    }
  }
} else {
  print "sorry, the periodicity for periodic key $periodic is not implemented.";
  exit;
}

open(XYZ, ">", $ARGV[3]) or die "Can't open output XYZ-File: $!";
print XYZ ($#inside+2)."\n\n";
printf XYZ "%-8s %10.5f %10.5f %10.5f", $data[$center][0], $data[$center][1], $data[$center][2], $data[$center][3];
for($i=0;$i<@inside;$i++) {
  printf XYZ "\n%-8s %10.5f %10.5f %10.5f", $inside[$i][0], $inside[$i][1], $inside[$i][2], $inside[$i][3];
  $anz{$inside[$i][0]}++
}
foreach $key (keys %anz) {
  print "$key $anz{$key}\n";
}