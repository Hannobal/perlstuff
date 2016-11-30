#!/usr/bin/perl

# calculates the rotational orientation angles of phosponic acids from DL_POLY
# HISTORY files and their develeopment in time for force_docking_z-type systems

use Math::Trig;
use Cwd;
use Switch;
use hanno_utility;
use dlpoly_utility;
use aloxsam_utility;

$histresdeg = 1;
$toldeg     = 10;
$startframe = 0;
$endframe   = 9e20;
$histname   = "HISTORY";
$fldname    = "FIELD";
$outdir     = $ARGV[0];
$subdir     = "run";
$lnumsubdir = 0;
$twopi      = 2*pi;

if($#ARGV<0) {
  print "input format:\n";
  print " 1. <output-directory>\n";
  print "-r <real>    histogram resolution for PA orientation angles (default: $histresdeg deg)]\n";
  print "-n <real>    tolerance for rotation (default: $toldeg deg)]\n";
  print "-d  <str>    subdirectory to analyze (default: $subdir, \"num\" for last numeric directory)\n";
  print "-t n*<str>   list of molecules to analyze\n";
  print "-f <2*int>   start and end frame\n";
  print "-s <2*int>   start and end step\n";
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
    } case(/^-f/i) {
      $i++;
      $startframe = $ARGV[$i];
      &error(1,"start frame must be an integer") if not check_integer($startframe);
      $i++;
      $endframe = $ARGV[$i];
      &error(1,"end frame must be an integer") if not check_integer($endframe);
    } case(/^-s/i) {
      $i++;
      $startstep = $ARGV[$i];
      &error(1,"start step must be an integer") if not check_integer($startstep);
      $i++;
      $endstep = $ARGV[$i];
      &error(1,"end step must be an integer") if not check_integer($endstep);
    } case(/^-d/i) {
      $i++;
      $subdir=$ARGV[$i];
      if($subdir eq "num") {
	$lnumsubdir=1;
      } else {
	$lnumsubdir=0;
      }
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
print "reading directories\n";
@directories=get_numeric_directories(".",$startstep,$endstep);
&error(3,"no suitable directories found") if(not @directories);
if($lnumsubdir) {
  @numdirs = get_numeric_directories($directories[0]);
  &error(3,"no numeric subdirectories found") if(@numdirs);
  if(-d "$directories[0]/min") {
    exit 1 if(read_field_file("$directories[0]/$numdirs[$#numdirs]/min/FIELD",0)!=0);
  } else {
    exit 1 if(read_field_file("$directories[0]/$numdirs[$#numdirs]/FIELD",0)!=0);
  }
} else {
  &error(3,"directory $directories[0]/$subdir does not exist") if(not -d "$directories[0]/$subdir");
  @numdirs = get_numeric_directories("$directories[0]/$subdir");
  if(@numdirs) {
    exit 1 if(read_field_file("$directories[0]/$subdir/$numdirs[0]/FIELD",0)!=0);
  } else {
    exit 1 if(read_field_file("$directories[0]/$subdir/FIELD",0)!=0);
  }
}

############### analyze the FIELD file ########################################

@tpas=();
for($t=0;$t<$field_nummols[0];$t++) {
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

&error(4,"no molecules to analyze") if(not @tpas);
$tpamax=$#tpas+1;
if(@tpas>1) {
  push(@{$mol_name[0]},"all");
  push(@tpas,$#{$mol_name[0]});
  $iall=$#tpas;
}

############### read the HISTORY file from each directory #####################

mkdir($outdir) if(not -d $outdir);
for($i=0;$i<@tpas;$i++) {
  $t=$tpas[$i];
  &error(4) if(not open($fhnumrot[$i],">","$outdir/NUM_PA_ROTATIONS_".$mol_name[0][$t]));
  print {$fhnumrot[$i]} "#number of rotations of PAs\n";
  &error(4) if(not open($fhhisto[$i],">","$outdir/HIST_PA_ROTATIONS_".$mol_name[0][$t]));
  print {$fhhisto[$i]} "#histogram of rotation angles of PAs\n";
}

if($lessoutput) {
  print "begin analysis...\n";
}

foreach $dir (@directories) {
  $nummols = $dir;
  if($lnumsubdir) {
    @numdirs = get_numeric_directories("$dir");
    $dir .= "/$numdirs[$#numdirs]" if(@numdirs>0);
    $dir .= "/min" if(-d "$dir/min");
  } else {
    $dir .= "/$subdir";
  }
  $firsthist = 1;
  $f         = 0;
  $numpas    = 0;
  for($i=0;$i<$tpamax;$i++) {
    $numpas += $mol_numents[0][$tpas[$i]];
    $numrot[$i] = 0;
  }
  print "analyzing directory $dir\n";
  @directories2=get_numeric_directories($dir,$startframe,$endframe);
  push(@directories2,".") if(-e "$dir/HISTORY");
  if(not @directories2) {
    print "**** warning: did not find HISTORY files in $dir\n";
    next;
  }
  &error(6) if(0!=update_numents("$dir/$directories2[0]/FIELD",0));
  foreach $dir2 (@directories2) {
    &error(2,"$dir/$histname") if(not open($fhhist,"<","$dir/$histname"));
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
      for($i=0;$i<$tpamax;$i++) {
	$t=$tpas[$i];
	for($m=0;$m<$mol_numents[0][$t];$m++) {
	  @d = calc_dvec_orthocell(0, $t, $m, $p5[$i], 0, $t, $m, $c3[$i]);
	  $angle = atan2($d[1],$d[0]);
	  $angle += $twopi if($angle<0);
	  histogram_add_one($i,$angle,$histresrad);
	  histogram_add_one($iall,$angle,$histresrad) if(@tpas>1);
	  unless($firsthist) {
	    $change=$angle-$oldangle[$i][$m];
	    $change+=$twopi if($angle<-pi);
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
  $firsthist=0;
  
############### print histograms ##############################################
  
  for($i=0;$i<@tpas;$i++) {
    printf {$fhnumrot[$i]} "%4u %10u\n",$dir,$numrot[$i];
    print {$fhhisto[$i]} "# directory $dir\n";
    histogram_normalize_integral($i,$histresdeg);
    for($j=$histminindex[$i];$j<=$histmaxindex[$i];$j++) {
      printf {$fhhisto[$i]} "%17.10g %17.10g %10u\n",$histresdeg*$j,$histogram[$i][$j],$histnum[$i][$j];
    }
    print {$fhhisto[$i]} "\n\n";
    histogram_clear($i);
  }
  
  if($lessoutput) {
    print "last frame analyzed: $frame_number[0]\n";
  } else {
    print "\n";
  }
}

close(@fhhisto, @fhnumrot);

sub update_numents {
  # arguments: filename, fi
  my($fh,$t,@new);
  return 1 if(not open($fh,"<",$_[0]));
  $t=0;
  while($_=<$fh>) {
    if(/NUMMOLS\s+(\d+)/i) {
      $new[$t]=$1;
      $t++;
    }
  }
  return 1 if($field_nummols[$_[1]] != $t);
  @{$mol_numents[$_[1]]}=@new;
  return 0;
}

sub error {
  if(defined $_[1]) {
    print "**** error $_[1]\n";
  } else {
    print "**** unknown error!\n";
  }
  exit $_[0]
}