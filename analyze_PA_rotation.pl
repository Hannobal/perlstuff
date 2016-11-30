#!/usr/bin/perl

# calculates the rotational orientation angles of phosponic acids from DL_POLY
# HISTORY files

use Math::Trig;
use Cwd;
use Switch;
use hanno_utility;
use dlpoly_utility;
use aloxsam_utility;
use strict;

my($t,$m,$a,$b,$i,$j,@includetypes,@tpas,$tpamax,@optimangle,$histresrad,$talfree,
@directories,$dir,@p5,@c3,@numrot,$numpas,$fhhist,@fhhisto,@fhhistdiff,@fhnumrot,
$angle,$diffangle,@d,$iall,$change,$lessoutput,$firsthist,$f,$err,$tolrad,@oldangle);

my $histresdeg = 1;
my $toldeg     = 20;
my $startframe = 0;
my $endframe   = 9e20;
my $histname   = "HISTORY";
my $fldname    = "FIELD";
my $outdir     = $ARGV[0];
our $twopi      = 2*pi;
our $pisixth = pi/6.0;

if($#ARGV<0) {
  print "input format:\n";
  print " 1. <output-directory>\n";
  print "-r <real>    histogram resolution for PA orientation angles (default: $histresdeg deg)]\n";
  print "-n <real>    tolerance for rotation (default: $toldeg deg)]\n";
  print "-t n*<str>   list of molecules to analyze\n";
  print "-s <2*int>   start frame\n";
  print "-e <2*int>   end frame\n";
  print "-lessoutput  do not print current frame\n";
  exit 1;
} 

$outdir = $ARGV[0];

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    # general flags
    case (/^-r/i) {
      $i++;
      $histresdeg = $ARGV[$i];
      &error(1,"histogram resolution for orientation angle must be >=0") if($histresdeg<=0);
    } case (/^-n/i) {
      $i++;
      $toldeg = $ARGV[$i];
      &error(1,"tolerance for orientation angle must be >=0") if($toldeg<=0);
    } case(/^-s/i) {
      $i++;
      $startframe = $ARGV[$i];
      &error(1,"start frame must be an integer") if not check_integer($startframe);
    } case(/^-e/i) {
      $i++;
      $endframe = $ARGV[$i];
      &error(1,"end frame must be an integer") if not check_integer($endframe);
    } case(/^-t/i) {
      $i++;
      until($ARGV[$i]=~/^-/ or $i>$#ARGV) {
	push(@includetypes,$ARGV[$i]);
	$i++;
      }
      $i--;
    } case(/^-lessoutput/i) {
      $lessoutput=1;
    } else {
      &error(2,"could not interpret flag $ARGV[$i]");
    }
  }
}

&error(1,"end frame must not be smaller than start frame") if($endframe<$startframe);
$histresrad = $histresdeg/180*pi;
$tolrad     = $toldeg/180*pi;

############### find the relevant directories to read #########################
@directories=get_numeric_directories(".",$startframe,$endframe);
push(@directories,".") if(-e "./$histname");
&error(3,"no suitable directories found") if(not @directories);

############### analyze the FIELD file ########################################

exit 1 if(read_field_file("$directories[0]/$fldname",0)!=0);
@tpas=();
$talfree=-1;
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t]=~/alumin.*free/) {
    $talfree=$t;
    next;
  }
  next if(not $mol_name[0][$t]=~/-(PA|CA)-\d?d/);
  next if($mol_numents[0][$t]==0);
  if(@includetypes) {
    next if(not contains(@includetypes,$mol_name[0][$t]));
  }
  push(@tpas,$t);
  findp5:for($a=0;$a<$mol_numatoms[0][$t];$a++) {
    if($mol_atomdata[0][$t][$a][0]=~/^P5$/i) {
      $p5[$#tpas] = $a;
      foreach $b (@{$mol_bondatoms[0][$t][$a]}) {
	if($mol_atomdata[0][$t][$b][0]=~/^C3$/i) {
	  $c3[$#tpas] = $b;
	  last findp5;
	}
      }
    }
  }
  &error(4,"did not find P5 atom for molecule $mol_name[0][$t]") if(not defined($p5[$#tpas]));
  &error(4,"did not find C3 atom for molecule $mol_name[0][$t]") if(not defined($c3[$#tpas]));
}

&error(4,"did not find free aluminum") if($talfree<0);
&error(4,"no free aluminum atoms defined") if($mol_numents[0][$talfree]==0);
&error(4,"no molecules to analyze") if(not @tpas);
$tpamax=$#tpas+1;
if(@tpas>1) {
  push(@{$mol_name[0]},"all");
  push(@tpas,$#{$mol_name[0]});
  $iall=$#tpas;
}

$numpas = 0;
for($i=0;$i<$tpamax;$i++) {
  $numpas += $mol_numents[0][$tpas[$i]];
  $numrot[$i] = 0;
}

############### read the HISTORY file from each directory #####################

mkdir($outdir) if(not -d $outdir);
for($i=0;$i<@tpas;$i++) {
  $t=$tpas[$i];
  &error(4) if(not open($fhnumrot[$i],">","$outdir/NUM_PA_ROTATIONS_".$mol_name[0][$t]));
  print {$fhnumrot[$i]} "#number of rotations of PAs\n";
  &error(4) if(not open($fhhisto[$i],">","$outdir/HIST_PA_ROTATIONS_".$mol_name[0][$t]));
  print {$fhhisto[$i]} "#histogram of rotation angles of PAs\n";
  &error(4) if(not open($fhhistdiff[$i],">","$outdir/HIST_PA_ROTATIONS_DIFF_".$mol_name[0][$t]));
  print {$fhhistdiff[$i]} "#histogram of rotation angle difference of PAs with respect to their optimum\n";
}

if($lessoutput) {
  print "begin analysis...\n";
}

$firsthist = 1;
$f         = 0;
foreach $dir (@directories) {
  print "analyzing directory $dir\n" if($lessoutput);
  &error(2,"cannot open $dir/$histname") if(not open($fhhist,"<","$dir/$histname"));
  if($startframe>0 and $firsthist) {
    &error(5,$startframe) if(not find_history_timestep($fhhist,$startframe));
  }
  $firsthist=0;
  $err=0;
  while($err==0) {
    $err=read_history_timestep($fhhist,0,0);
    last if($err!=0);
    last if($frame_number[0]>$endframe);
    print "\ranalyzing timestep $frame_number[0]" unless($lessoutput);
    if(not @optimangle) {
      for($i=0;$i<$tpamax;$i++) {
	@{$optimangle[$i]} = calc_optim($tpas[$i],$talfree);
      }
    }
    for($i=0;$i<$tpamax;$i++) {
      $t=$tpas[$i];
      for($m=0;$m<$mol_numents[0][$t];$m++) {
	# vector from P5 to C3 versus x-axis (in xy-plane)
	@d = calc_dvec_orthocell(0, $t, $m, $c3[$i], 0, $t, $m, $p5[$i]);
	$angle = atan2($d[1],$d[0]);
	$angle += $twopi if($angle<0);
	$diffangle = $angle-$optimangle[$i][$m];
	$diffangle += $twopi if($diffangle<0);
	$diffangle -= $twopi if($diffangle>$twopi);
	histogram_add_one($i,$angle,$histresrad);
	histogram_add_one(@tpas+$i,$diffangle,$histresrad);
	histogram_add_one($iall,$angle,$histresrad) if(@tpas>1);
	histogram_add_one(@tpas+$iall,$diffangle,$histresrad) if(@tpas>1);
	unless($firsthist) {
	  $change=$angle-$oldangle[$i][$m];
	  $change+=$twopi if($angle < -(pi));
	  if(abs($change)>$tolrad) {
	    $numrot[$i]++;
	    $numrot[$iall]++ if(@tpas>1);
	  }
	}
	$oldangle[$i][$m]=$angle;
      }
    }
    $f++;
  }
  close($fhhist);
}
  
############### print histograms ##############################################
  
for($i=0;$i<@tpas;$i++) {
  printf {$fhnumrot[$i]} "%4u %10u\n",$dir,$numrot[$i];
  histogram_normalize_integral($i,$histresdeg);
  for($j=$histminindex[$i];$j<=$histmaxindex[$i];$j++) {
    printf {$fhhisto[$i]} "%17.10g %17.10g %10u\n",$histresdeg*$j,$histogram[$i][$j],$histnum[$i][$j];
  }
  print {$fhhisto[$i]} "\n\n";
  histogram_clear($i);
  my $i2=$i+@tpas;
  histogram_normalize_integral($i2,$histresdeg);
  for($j=$histminindex[$i2];$j<=$histmaxindex[$i2];$j++) {
    printf {$fhhistdiff[$i]} "%17.10g %17.10g %10u\n",$histresdeg*$j,$histogram[$i2][$j],$histnum[$i2][$j];
  }
  print {$fhhistdiff[$i]} "\n\n";
  histogram_clear($i2);
}

if($lessoutput) {
  print "last frame analyzed: $frame_number[0]\n";
} else {
  print "\n";
}

close(@fhhisto, @fhnumrot, @fhhistdiff);

sub calc_optim {
  my $tpa = $_[0];
  my $tal = $_[1];
  my @optimangle=();
  for(my $m=0;$m<$mol_numents[0][$tpa];$m++) {
    my $mindsq=999;
    my $al1 = -1;
    my $al2 = -1;
    my $o = -1;
    for(my $a=0;$a<$mol_numatoms[0][$tpa];$a++) {
      next if($mol_atomdata[0][$tpa][$a][0] !~ /^O$/);
      for(my $n=0;$n<$mol_numents[0][$tal];$n++) {
	  my $dsq = calc_dsq_orthocell(0,$tal,$n,0, $tpa,$m,$a);
	  if($dsq<$mindsq) {
	    $al2=$al1;
	    $al1=$n;
	    $o=$a;
	    $mindsq=$dsq;
	  }
	}
    }
    # set al1 to be the upper atom
    if($cdata[0][$tal][$al2][0][2]>$cdata[0][$tal][$al1][0][2]) {
      my $tmp=$al1;
      $al1=$al2;
      $al2=$tmp;
    }
    my @vec = calc_dvec_orthocell(0,$tal,$al1,0,  0,$tpa,$m,$o);
    my $angle=(round(atan2($vec[1],$vec[0])/$pisixth)+2)*$pisixth;
    push(@optimangle,$angle);
  }
  return @optimangle;
}

sub error {
  if(defined $_[1]) {
    print "**** error $_[1]\n";
  } else {
    print "**** unknown error!\n";
  }
  exit $_[0]
}