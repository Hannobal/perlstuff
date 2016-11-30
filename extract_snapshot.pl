#!/usr/bin/perl

# calculates the tilt angles, dihedral angles and chain lengths of carbon backbone of
# certain phosponic acids from DL_POLY HISTORY files and their develeopment in time

use Math::Trig;
use Cwd;
use Switch;
use hanno_utility;
use dlpoly_utility;
use aloxsam_utility;
use Storable qw(dclone);
@includetypes = ();

use Math::Trig;

$histrestiltdeg = 2;
$histresdihdeg  = 5;
$histreslen     = 2;
$startframe = 0;
$endframe   = 9e20;
$histname   = "HISTORY";
$fldname    = "FIELD";
$outdir     = $ARGV[0];

if($#ARGV<0) {
  print "input format:\n";
  print " 1. output filename\n";
  print " 2. frame number\n";
  print "-xyz   output as xyz\n";
  exit;
}

$lxyz=0;
$outfilename=$ARGV[0];
$targetframe=$ARGV[1];

for($i=2;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case(/^-xyz/i) {
      $lxyz = 1;
    } else {
      &error(2,"could not interpret flag $ARGV[$i]");
    }
  }
}

@directories=get_numeric_directories(".");
if(@directories) {
  @directories=get_numeric_directories(".",$targetframe,$targetframe);
  if(not @directories) {
    push(@directories,".") if(-e "./HISTORY");
  }
}

&error(3,"no suitable directories found") if(not @directories);

open($fhhist,"<","$directories[0]/HISTORY") or die "**** error: can't open file $directories[0]/HISTORY!\n";
if(not find_history_timestep($fhhist,$targetframe)) {
  print "**** error: can't find timestep $targetframe in file $directories[0]/HISTORY\n"; exit 1;
}
if(read_history_timestep($fhhist,-1,0)!=0) {
  print "**** error: can't read timestep $targetframe in file $directories[0]/HISTORY\n"; exit 1;
}
close($fhhist);

open($fhout,">",$outfilename) or die "**** error: can't open file $outfilename: $!\n";
if($lxyz) {
  write_xyz_timestep($fhout,0,$config_title[0]);
} else {
  write_history_timestep($fhout	,-1,0);
}
close($fhout);