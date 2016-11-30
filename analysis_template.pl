#!/usr/bin/perl

use hanno_utility;
use dlpoly_utility;
use Switch;
use Storable qw(dclone);
use Math::Trig;

# $|=1;

$lessoutput      = 0;
$startframe      = 0;
$endframe        = 9e20;
$outdir          = "analysis";

if($#ARGV<0) {
  print "required input:\n";
  
  print "optional flags:\n";
  print "  -o <str>     output folder (default: $outdir)\n";
  print "  -s <int>     start frame\n";
  print "  -e <int>     end frame\n";
  print "  -lessoutput  do not print current frame\n";
}

for($i=0;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    # required input
    case (/^-XXX/i) {
    # optional flags
    } case(/^-o/i) {
      $i++;
      $outdir=$ARGV[$i];
    } case(/^-s/i) {
      $i++;
      $startframe = $ARGV[$i];
      &error(1,"Start frame must be an integer number!") if not check_integer($startframe);
    } case(/^-e/i) {
      $i++;
      $endframe = $ARGV[$i];
      &error(2,"End frame must be an integer number!") if not check_integer($endframe);
    } case(/^-lessoutput/i) {
      $lessoutput=1;
    } else {
      &error(10,"Could not interpret flag $ARGV[$i]!");
    }
  }
}

@directories=get_numeric_directories(".");
if(@directories) {
  @directories=get_numeric_directories(".",$startframe,$endframe);
}
push(@directories,".") if(-e "./HISTORY");
&error(7,"No suitable directories found for analysis!") if(not @directories);

exit 1 if(read_field_file("$directories[0]/FIELD",0)!=0);

mkdir($outdir) if(not -d $outdir);

$firsthist=1;
print "begin analysis...\n" if($lessoutput);

foreach $dir (@directories) {
  &error(11,"Could not open file $dir/HISTORY!") if(not open($fhhist,"<","$dir/HISTORY"));
  if($startframe>0 and $firsthist) {
    &error(9,"Did not find frame $startframe in $dir/HISTORY!") if(not find_history_timestep($fhhist,$startframe));
  }
  $firsthist=0;
  $err=0;
  while($err==0) {
    $err=read_history_timestep($fhhist,0,0);
    last if($err!=0);
    last if($frame_number[0]>$endframe);
    print "\ranalyzing timestep $frame_number[0]" unless($lessoutput);
    # do stuff here with the frame
  }
  close($fhhist);
}

print "\n" unless($lessoutput);

# write histograms and other output



sub error {
  if(@_>1) {
    print "**** error: $_[1]\n";
  } else {
    print "**** unknown error!\n";
  }
  exit $_[0] if(@_);
  exit 99;
}