#!/usr/bin/perl

# extracts the box dimensions from the DL_POLY HISTORY file
# and calculates average values 

use dlpoly_utility;
use hanno_utility;
use Switch;

if($#ARGV<0) {
  print "the input format is:\n";
  print " 1. output directory for the output files\n";
  print "-s <int>     start frame\n";
  print "-e <int>     end frame\n";
  exit 1;
}

$startframe = 0;
$endframe   = 9e20;
$outdir     = $ARGV[0];

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case(/^-s/i) {
      $i++;
      $startframe = $ARGV[$i];
      &error(1,"start frame must be an integer") if not check_integer($startframe);
    } case(/^-e/i) {
      $i++;
      $endframe = $ARGV[$i];
      &error(1,"end frame must be an integer") if not check_integer($endframe);
    } else {
      print "**** error: unkown flag \"$ARGV[$i]\"!\n"; exit 1;
    }
  }
}

&error(2,"end frame must be greater than start frame") if($endframe<$startframe);

@directories=get_numeric_directories(".",$startframe,$endframe);
push(@directories,".") if(-e "./HISTORY");
&error(3,"no suitable directories found") if(not @directories);


$f=0;
@data=();
foreach $dir (@directories) {
  exit 1 if(not open($fhhist,"<","$dir/HISTORY"));
  if($startframe>0 and $firsthist) {
    &error(5,$startframe) if(not find_history_timestep($fhhist,$startframe));
  }
  $firsthist=0;
  $err=0;
  while($err==0) {
    $err=read_history_timestep($fhhist,-1,0);
    last   if($err<0);
    exit 1 if($err>0);
    last if($frame_number[0]>$endframe);
    push(@data,[$frame_number[0],
      $cell[0][0][0],$cell[0][0][1],$cell[0][0][2],
      $cell[0][1][0],$cell[0][1][1],$cell[0][1][2],
      $cell[0][2][0],$cell[0][2][1],$cell[0][2][2]]);
    $f++;
  }
  
  close($fhhist);
}
print "asdf\n";
mkdir($outdir) if(not -d $outdir);

open($fhout,">","$outdir/BOX_SIZE");
printf $fhout "#%7s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n","frame",
  "cell 1x","cell 1y","cell 1z","cell 2x","cell 2y","cell 2z","cell 3x","cell 3y","cell 3z";
for($f=0;$f<@data;$f++) {
  printf $fhout "%8u %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",@{$data[$f]};
}
