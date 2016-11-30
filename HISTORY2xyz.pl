#!/usr/bin/perl

use dlpoly_utility;

if($#ARGV<2) {
  print "1. name of HISTORY-file\n";
  print "2. name of XYZ-file\n";
  print "3. (c/f/v) for coordinate (c) + force (f) or velocity (v) information\n";
  print "[4. timestep (default: last frame)]\n";
  print "[5. divide force by charge? (y/n)]\n";
  print "[6. scale factor for force/velocity vector]\n";
  print "[7. list of atom types to include]\n";
  exit;
}

if($ARGV[2] =~ /^f/i) {
  print "print forces\n";
  $print_key=2;
} elsif($ARGV[2] =~ /^v/i) {
  print "print velocities\n";
  $print_key=1;
} elsif($ARGV[2] =~ /^c/i) {
  $print_key=0;
} else {
  print "error: input for additional force/velocity information could not be interpreted.\n";
  print "       possible values are \"c\", \"f\" and \"v\".\n";
  exit;
}

if($#ARGV>2) {
  $wantedframe = $ARGV[3];
  $foundframe  = 0;
} else {
  $wantedframe = -1;
}
if($ARGV[4] =~ /^y/i) {
  $divide = 1;
} else {
  $divide = 0;
}

if($#ARGV>4) {
  $scalefactor = $ARGV[5];
} else {
  $scalefactor = 1;
}

$numtypes = 0;
for($i=6;$i<@ARGV;$i++) {
  $includetype[$numtypes] = $ARGV[$i];
  $numtypes++;
}

%type = (
  "al", "Al",
  "c", "C",
  "c1", "C",
  "c2", "C",
  "c3", "C",
  "ca", "C",
  "cs", "C",
  "cx", "C",
  "h", "H",
  "ha", "H",
  "h1", "H",
  "hc", "H",
  "hg", "H",
  "ho", "H",
  "hr", "H",
  "o",  "O",
  "hn", "H",
  "oa", "O",
  "oh", "O",
  "os", "O",
  "ox", "O",
  "p5", "P",
  "p", "P",
  "n", "N",
  "n1", "N",
  "n2", "N",
  "n3", "N",
  "n4", "N",
  "nb", "N",
  "nc", "N",
  "nd", "N",
  "ne", "N",
  "nf", "N",
  "nh", "N",
  "no", "N",
  "ss", "C",
  "pt", "Pt",
  "na", "Na",
);

open($fhhist, "<", $ARGV[0]) or die "Can't open HISTORY-File: $!";

$err=0;
while($err==0) {
  $err=read_history_timestep($fhhist,-1,0);
  exit 1 if($err>0);
  last if($frame_number[0]==$wantedframe);
}
close($fhhist);

open($fhxyz,">",$ARGV[1]);
write_xyz_timestep($fhxyz,0,$config_title[0]);
close($fhxyz);
