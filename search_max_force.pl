#!/usr/bin/perl

if($#ARGV==-1) {
  print "input format:\n";
  print " 1. name of CONFIG/REVCON-file\n";
  print "[2. -o only organic molecules (no Al2O3)]\n";
  exit;
}

if($#ARGV>0) {
  if($ARGV[1] =~ /^-[oO]/) {
    $noalox=1;
  }else {
    print "error: could not interpret flag $ARGV[1]\n";
    exit;
  }
} else {
  $noalox=0;
}

open(CONFIG, "<", $ARGV[0]) or die "Can't open CONFIG-File: $!";

$anzges = 0; #number of atoms counted in CONFIG
$title=<CONFIG>;
$_=<CONFIG>;
($config_key,$periodic_key) = /^\s*(\S+)\s+(\S+)/;
if($config_key<2) {
  print "no force data in CONFIG file (config_key<2)!\n";
  exit;
}
if($periodic_key!=0) { # skip lines for box size if necessary
  for($j=0;$j<=2;$j++) {
    $_=<CONFIG>;
    ($vec[$j][0], $vec[$j][1], $vec[$j][2]) = /^\s+(\S+)\s+(\S+)\s+(\S+)/;
  }
}
$maxforce=-1;
while(<CONFIG>) {
  ($name[$anzges]) = /^\s*([a-zA-Z]+)/;
  $_=<CONFIG>;
  ($pos[$anzges][0],$pos[$anzges][1],$pos[$anzges][2]) = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
  $_=<CONFIG>;
  ($vel[$anzges][0],$vel[$anzges][1],$vel[$anzges][2]) = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
  $_=<CONFIG>;
  ($frc[$anzges][0],$frc[$anzges][1],$frc[$anzges][2]) = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
  if((uc($name[$anzges]) ne "AL" and uc($name[$anzges]) ne "OA") or $noalox==0) {
    $absforce=sqrt($frc[$anzges][0]**2+$frc[$anzges][1]**2+$frc[$anzges][2]**2);
    if($absforce>$maxforce) {
      $i = $anzges;
      $maxforce = $absforce;
    }
  }
  $anzges++;
}

print "maximum force with absolute length ".$maxforce." for atom ".($i+1).":\n";
printf "%-8s", $name[$i];
printf "%10u\n",$i+1; #atom index
printf " %19.12e %19.12e %19.12e\n",$pos[$i][0],$pos[$i][1],$pos[$i][2];
printf " %19.12e %19.12e %19.12e\n",$vel[$i][0],$vel[$i][1],$vel[$i][2];
printf " %19.12e %19.12e %19.12e\n",$frc[$i][0],$frc[$i][1],$frc[$i][2];