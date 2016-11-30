#!/usr/bin/perl

use hanno_utility;
use dlpoly_utility;
use Switch;
use Storable qw(dclone);
use Math::Trig;
use Math::Matrix;

$lessoutput = 0;
$outdir     = $ARGV[0];
$doo        = 2.755333333;
$tolz       = 1.5;
$angle      = 0;

if($#ARGV<0) {
  print "input format:\n";
  print " 1. output directory for the output files\n";
  print "-t <real>    tolerance for z-position (default: $tolz)\n";
  print "-r <real>    grid rotation in degree (default: $angle = one O-O vector in x-direction)\n";
  print "-s <int>     start frame\n";
  print "-e <int>     end frame\n";
  print "-lessoutput  do not print current frame\n";
  exit;
}

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    # general flags
    case (/^-t/i) {
      $i++;
      $tolz = $ARGV[$i];
      &error(1,"Tolerance must be a real number > 0") if($tolz<=0);
    } case (/^-r/i) {
      $i++;
      &error(1,"Grid rotation must be a real number > 0") if($ARGV[$i]<=0);
      $angle = $ARGV[$i]/180.0*pi;
    } case(/^-s/i) {
      $i++;
      $startframe = $ARGV[$i];
      &error(1,"Start frame must be an integer number") if not check_integer($startframe);
    } case(/^-e/i) {
      $i++;
      $endframe = $ARGV[$i];
      &error(1,"End frame must be an integer number") if not check_integer($endframe);
    } case(/^-lessoutput/i) {
      $lessoutput=1;
    } else {
      &error(1,$ARGV[$i]);
    }
  }
}

@directories=get_numeric_directories(".");
if(@directories) {
  @directories=get_numeric_directories(".",$startframe,$endframe);
}
push(@directories,".") if(-e "./HISTORY");
&error(2,"No suitable directories found for analysis!") if(not @directories);

exit 1 if(read_field_file("$directories[0]/FIELD",0)!=0);
$toh = -1;
$to  = -1;
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t]=~/hydrox/i) {
    $toh = $t;
  } elsif($mol_name[0][$t]=~/oxygen/i and not $mol_name[0][$t]=~/(frozen|fix)/i) {
    $to  = $t;
  }
}
if($toh<0) {
  &error(3,"Did not find hydroxide molecules or oxide atoms in FIELD file!") if($to<0);
  print "**** warning: did not find hydroxide ions in FIELD file\n";
  print "              using $mol_name[0][$to] (",$to+1,") instead\n";
  $toh = $to;
}

exit 1 if(read_config_file("$directories[0]/CONFIG",0,0)!=0);

# refine the o-o distance
$ix  = int($cell[0][0][0]/$doo);
$doo = $cell[0][0][0]/$ix;
$aoh = sqrt(3.0)/2.0*$doo*$doo;
@{$dvec[0]} = ($doo,0);
@{$dvec[1]} = ($doo*cos(pi/3.0),$doo*sin(pi/3.0));

# find the top layer of hydroxide molecules
$zmax=-9.0e20;
for($m=0;$m<$mol_numents[0][$toh];$m++) {
  if($zmax<$cdata[0][$toh][$m][0][2]) {
    $zmax=$cdata[0][$toh][$m][0][2];
  }
}

$zmax-=$tolz;
@moh=();
for($m=0;$m<$mol_numents[0][$toh];$m++) {
  next if($zmax>$cdata[0][$toh][$m][0][2]);
  push(@moh,$m);
}

$numvac = sprintf("%.0f",$cell[0][0][0]*$cell[0][1][1]/$aoh-@moh);
if($numvac==1) {
  print "found $numvac vacancy\n";
} else {
  print "found $numvac vacancies\n";
}

# find the grid center
$mat = new Math::Matrix([cos($angle),sin($angle)],[cos($angle+pi/3.0),sin($angle+pi/3.0)]);
$mat = $mat->multiply_scalar($doo)->transpose();
$invmat = $mat->invert();
$center = new Math::Matrix([$cdata[0][$toh][$moh[0]][0][0]],[$cdata[0][$toh][$moh[0]][0][1]]);
$shift = new Math::Matrix([0],[0]);
$deviaton = 0;
foreach $m (@moh) {
  $vec = new Math::Matrix([$cdata[0][$toh][$m][0][0]],[$cdata[0][$toh][$m][0][1]]);
  $vec = $invmat*($vec-$center);
  $vec2 = new Math::Matrix([sprintf("%.0f",$vec->[0][0])],[sprintf("%.0f",$vec->[1][0])]);
  $shift = $shift+$vec-$vec2;
  $deviaton+=($mat*($vec-$vec2))->absolute();
}
print "deviaton: $deviaton\n";
$shift = $shift->multiply_scalar(1.0/@moh);
$center = $center+$mat*$shift;

# mkdir($outdir) if(not -d $outdir);
# &error(4) if(not open($fhavdihedral,">","$outdir/AV_CCCC_DIHEDRAL"));
# print $fhavdihedral "# average dihedral angle in degrees within $radiuscutoff A from center\n";
# print $fhavdihedral "# molecules used for analysis: ",join(" ",@includetypes),"\n";
# printf $fhavdihedral "#%9s %17s %17s\n","frame","dihedral","stdev";

$firsthist=1;
$f=0;
if($lessoutput) {
  print "begin analysis...\n";
} else {
  print "\n" 
}
foreach $dir (@directories) {
  exit 1 if(not open($fhhist,"<","$dir/$histname"));
  if($startframe>0 and $firsthist) {
    &error(5,$startframe) if(not find_history_timestep($fhhist,$startframe));
  }
  $firsthist=0;
  $err=0;
  while($err==0) {
    $err=read_history_timestep($fhhist,0,0);
    last   if($err<0);
    exit 1 if($err>0);
    last if($frame_number[0]>$endframe);
    print "\ranalyzing timestep $frame_number[0]" unless($lessoutput);
    
    $f++;
  }
  close($fhhist);
}


sub calc_deviation {
  my @pos  = @{$_[0]};
  my @vecs = @{$_[1]};
  my($i,@dev,$int);
  if($#_>1) {
    my @center = @{$_[2]};
    for($i=0;$i<@pos;$i++) {
#       $int=sprintf("%.0f",$pos[$i]/
    }
  } else {
    
  }
}

sub error {
  if(@_>1) {
    print "**** error: $_[1]\n";
  } else {
    print "**** unknown error\n";
  }
  exit $_[0];
}
