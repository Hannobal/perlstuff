#!/usr/bin/perl

use hanno_utility;
use dlpoly_utility;
use Switch;
use Storable qw(dclone);
$startframe      = 0;
$endframe        = 9e20;
use Math::Trig;

$lessoutput = 0;
if($#ARGV<0) {
  print "input format:\n";
  print " 1. output directory for the output files\n";
  print "-s <int>     start frame\n";
  print "-e <int>     end frame\n";
  print "-lessoutput  do not print current frame\n";
  exit;
}

$outdir=$ARGV[0];

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case(/^-s/i) {
      $i++;
      $startframe = $ARGV[$i];
      &error(1) if not check_integer($startframe);
    } case(/^-e/i) {
      $i++;
      $endframe = $ARGV[$i];
      &error(2) if not check_integer($endframe);
    } case(/^-lessoutput/i) {
      $lessoutput=1;
    } else {
      print "**** error: could not interpret flag \"$ARGV[$i]\"!\n"; exit 1;
    }
  }
}

@directories=get_numeric_directories(".");
if(@directories) {
  @directories=get_numeric_directories(".",$startframe,$endframe);
}
push(@directories,".") if(-e "./$histname");
if(not @directories) {
  print "**** error: did not find HISTORY files!\n"; exit 1;
}

exit 1 if(read_field_file("$directories[0]/FIELD",0)!=0);
mkdir($outdir) if(not -d $outdir);
if(not open($fhout,">","$outdir/CELLVECTORS")) {
  print "**** error: could not open output file $fhout!\n";
}
printf $fhout "#%9s %17s %17s %17s %17s %17s %17s\n","frame","length_vec_1","length_vec_2","length_vec_3","alpha","beta","gamma";


$firsthist=1;
$f=0;
if($lessoutput) {
  print "begin analysis...\n";
} else {
  print "\n" 
}
foreach $dir (@directories) {
  exit 1 if(not open($fhhist,"<","$dir/HISTORY"));
  if($startframe>0 and $firsthist) {
    if(not find_history_timestep($fhhist,$startframe)) {
      print "**** error: did not find frame $startframe in file $dir/HISTORY!\n"; exit 1;
    }
  }
  $firsthist=0;
  $err=0;
  while($err==0) {
    $err=read_history_timestep($fhhist,0,0);
    last   if($err<0);
    exit 1 if($err>0);
    last if($frame_number[0]>$endframe);
    print "\ranalyzing timestep $frame_number[0]" unless($lessoutput);
    @celldata = calc_cell_abc(0);
    printf $fhout "%10u %17.10f %17.10f %17.10f %17.10f %17.10f %17.10f\n", $frame_number[0],@celldata;
  }
  close($fhhist);
}