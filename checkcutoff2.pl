#!/usr/bin/perl

if($#ARGV<4) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file\n";
  print "2. index of atom\n";
  print "3. min cutoff\n";
  print "4. max cutoff\n";
  print "5. step\n";
  exit;
}

$center=$ARGV[1]-1;
$mincutoff=$ARGV[2];
$maxcutoff=$ARGV[3];
$step=$ARGV[4];
if($step<=0) {
  print "error: step must be greater than zero\n";
  exit 1;
}
if($maxcutoff<$mincutoff) {
  print "error: maximum cutoff must be greater than minimum cutoff\n";
  exit 1;
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

print "centeratom: $data[$center][0]\n";

for($cutoff=$maxcutoff;$cutoff>$mincutoff;$cutoff-=$step) {
  %anz=();
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
		$anz{$data[$i][0]}++;
# 		push(@inside,$data[$i]);
# 		$inside[$#inside][1]-=$x*$vec[0][0]+$y*$vec[1][0]+$z*$vec[2][0];
# 		$inside[$#inside][2]-=$x*$vec[0][1]+$y*$vec[1][1]+$z*$vec[2][1];
# 		$inside[$#inside][3]-=$x*$vec[0][2]+$y*$vec[1][2]+$z*$vec[2][2];
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
	      $anz{$data[$i][0]}++;
# 	      push(@inside,$data[$i]);
# 	      $inside[$#inside][1]-=$x*$vec[0][0]+$y*$vec[1][0];
# 	      $inside[$#inside][2]-=$x*$vec[0][1]+$y*$vec[1][1];
# 	      $inside[$#inside][3]-=$x*$vec[0][2]+$y*$vec[1][2];
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
	  $anz{$data[$i][0]}++;
	}
      }
    }
  } else {
    print "sorry, the periodicity for periodic key $periodic is not implemented.";
    exit;
  }
  $ratio=$anz{"AL"}/$anz{"OA"};
  printf "%5.2f %2s %4d %2s %4d %5s %10.4f\n", $cutoff, "AL", $anz{"AL"}, "OA", $anz{"OA"}, "ratio", $ratio;
}

# open(XYZ, ">", $ARGV[3]) or die "Can't open output XYZ-File: $!";
# print XYZ ($#inside+2)."\n\n";
# printf XYZ "%-8s %10.5f %10.5f %10.5f", $data[$center][0], $data[$center][1], $data[$center][2], $data[$center][3];
# for($i=0;$i<@inside;$i++) {
#   printf XYZ "\n%-8s %10.5f %10.5f %10.5f", $inside[$i][0], $inside[$i][1], $inside[$i][2], $inside[$i][3];
#   $anz{$inside[$i][0]}++
# }
# foreach $key (keys %anz) {
#   print "$key $anz{$key}\n";
# }