#!/usr/bin/perl

# count the number of surface hydroxides availablele for docking

use hanno_utility;
use dlpoly_utility;
use aloxsam_utility;
use Switch;

if($#ARGV<1) {
  print "input format:\n",
  " 1. CONFIG file\n",
  " 2. FIELD file\n",
  " optional arguments:\n",
  " -m <real>     minimum distance between binding sites (default: 3.5 A)\n",
  " -l <4*int>    lattice vectors\n",
  " -c <2*real>   center of lattice\n",
  " -t <real>     tolerance for hydroxide z-positions (default: 2 A)\n",
  " -b [-1/0/1]   find OH on top (1), bottom (-1) or both (0) (default: top)\n",
  " -xyz <string> print xyz file\n";
  " -p            print a list of all sites\n";
  " -r <int>      print a list of n randomly selected sites\n";
  exit;
}

$cfgname=$ARGV[0];
$fldname=$ARGV[1];
$tolerance=2;
$top=1;
$mindist=3.5;
@gridvec=();
@gridcenter=();
$lxyz=0;
$lprint=0;
$rand=-1;
for($i=2;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-m/) { # minimum PA distance
      $i++;
      $mindist = $ARGV[$i];
      if(not check_real($mindist)) {
	print "**** error: expected real number after flag -m but got $mindist!\n"; exit 1;
      }
    } case (/^-l/) {
      @gridvec = @ARGV[$i+1..$i+4];
      if(not check_integer(@gridvec)) {
	print "**** error: lattice vectors must be four integer numbers!\n"; exit 1;
      }
      $i+=4;
    } case (/^-c/) {
      @gridcenter = @ARGV[$i+1..$i+2];
      if(not check_real(@gridcenter)) {
	print "**** error: lattice center must be two real numbers!\n"; exit 1;
      }
      $i+=2;
    } case (/^-b/) {
      $i++;
      $top = $ARGV[$i];
      if(not $top=~/^[123]$/) {
	print "**** error: expected -1, 0 or 1 for top/bottom!\n"; exit 1;
      }
    } case (/^-xyz/) {
      $i++;
      $xyzname = $ARGV[$i];
      $lxyz    = 1;
    } case (/^-xyz/) {
      $i++;
      $xyzname = $ARGV[$i];
      $lxyz    = 1;
    } case (/^-p/) {
      $lprint=1;
    } case (/^-r/) {
      $i++;
      $rand = $ARGV[$i];
      if(not check_integer($rand)) {
	print "**** error: argument of flag -r must be an integer!\n"; exit 1;
      }
    } else {
      print "**** error: expected optional flag but got \"$ARGV[$i]\"!\n"; exit 1;
    }
  }
}

exit 0 if(read_field_file($fldname,0)!=0);
exit 0 if(read_config_file($cfgname,0,0)!=0);

if(@gridvec) {
  @candidates=find_hydroxide_candidates(0,0,$top,$mindist,2,@gridvec,@gridcenter);
} else {
  @candidates=find_hydroxide_candidates(0,0,$top,$mindist);
}
$numcands=$#candidates+1;
if($lxyz) {
  for($t=0;$t<$field_nummols[0];$t++) {
    last if ($mol_name[0][$t]=~/hydroxide/i);
  }
  open($fh,">",$xyzname) or die "**** error: can't open xyz file $xyzname: $!\n";
  print $fh "$numcands\n\n";
  for($i=0;$i<@candidates;$i++) {
    printf $fh "H %10.4f %10.4f %10.4f\n", @{$cdata[0][$t][$candidates[$i]][0]}[0..2];
  }
}
if($rand>0 and $rand<@candidates) {
  @ids=();
  for($i=0;$i<$rand;$i++) {
    $j = int(rand(@candidates));
    push(@ids,splice(@candidates,$j,1));
  }
  print join(" ",sort {$a<=>$b} @ids);
} elsif($lprint or $rand>=@candidates) {
  print join(" ",@candidates);
} else {
  print $numcands;
}
print "\n";