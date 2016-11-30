#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;
use Switch;

# extracts the information from DL_POLY STATIS files distributed
# over several folders with numbered names

$fldname = undef;
$outdir = $ARGV[0];
$startframe = 0;
$endframe   = 9e20;
$resolution = 1;
$lheatup    = 0;
$growth     = 0;
$fldname    = undef;
$ctrlname   = undef;
$ldefctrl   = 0;

use Cwd;
if($#ARGV<0) {
  print "the input format is:\n";
  print " 1. output directory for the output files\n";
  print "-f           FIELD file name (default: \"FIELD\" in first directory\n";
  print "-c           CONTROL file name (default: \"CONTROL\" in first directory\n";
  print "-h           toggle heatup mode\n";
  print "-ga          growth mode: last numeric/main + min + run\n";
  print "-gr          growth mode: run only\n";
  print "-s <int>     start frame\n";
  print "-e <int>     end frame\n";
  print "-r <int>     only print every n-th frame\n";
  print "-p <real>    ignore timestamp in STATIS and use this step size instead for 1st col\n";
  exit 1;
}

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case(/^-h/i) {
      $lheatup = 1;
    } case(/^-f/i) {
      $i++;
      $fldname = $ARGV[$i];
    } case(/^-c/i) {
      $i++;
      $ctrlname = $ARGV[$i];
      $ldefctrl = 1;
    } case(/^-s/i) {
      $i++;
      $startframe = $ARGV[$i];
      &error(1,"start frame must be an integer") if not check_integer($startframe);
    } case(/^-e/i) {
      $i++;
      $endframe = $ARGV[$i];
      &error(1,"end frame must be an integer") if not check_integer($endframe);
    } case (/^-r/i) {
      $i++;
      $resolution = $ARGV[$i];
      &error(1,"resolution must be an integer") if not check_integer($resolution);
      &error(1,"resolution must be greater than zero") if($resolution<=0);
    } case (/^-r/i) {
      $i++;
      $resolution = $ARGV[$i];
      &error(1,"resolution must be an integer") if not check_integer($resolution);
      &error(1,"resolution must be greater than zero") if($resolution<=0);
    } case (/^-ga/i) {
      $growth = 1;
    } case (/^-gr/i) {
      $growth = 2;
    } case (/^-p/i) {
      $i++;
      $step = $ARGV[$i];
      &error(1,"step size must be a real number") if not check_real($step);
    } else {
      &error(1,"cannot interpret flag $ARGV[$i]");
    }
  }
}

&error(2,"end frame must be greater than start frame") if($endframe<$startframe);

if(not $growth) {
  @directories=get_numeric_directories(".",$startframe,$endframe);
  push(@directories,".") if(-e "./STATIS");
} else {
  @directories=();
  @dirs=get_numeric_directories(".",$startframe,$endframe);
  foreach $dir (@dirs) {
    if($growth==1) {
      if(-e "$dir/HISTORY") {
	push(@directories,$dir);
      } else {
	@subdirs=get_numeric_directories($dir);
	next if(not @subdirs);
	push(@directories,"$dir/$subdirs[$#subdirs]");
      }
      push(@directories,"$dir/min") if(-d "$dir/min");
    }
    next if(not -d "$dir/run");
    push(@directories,"$dir/run") if(-d "$dir/run");
  }
}
&error(3,"no suitable directories found") if(not @directories);

############### read the FIELD file from the first directory ##################
if(not defined($fldname)) {
  $fldname = "$directories[0]/FIELD";
}
exit 1 if(read_field_file($fldname,0)!=0);

############### generate the output files #####################################
if(not -d $outdir) {
  mkdir $outdir or &error(4,"cannot create output directory:\n$!");
}
open($fhout, ">", "$outdir/STATIS_DATA") or &error(4,"cannot open file STATIS_DATA:\n$!");
print_statis_data_header($fhout,0);

############### read the STATIS file from each directory ######################
$time=0;
$time_old=0;
$starttime=0;
$f=-1;
foreach $dir (@directories) {
  print "reading $dir/STATIS\n";
  open($fhstatis, "<", "$dir/STATIS") or print "file STATIS was not found in $dir\n";
  if($ldefctrl) {
    open(CONTROL, "<", $ctrlname) or print "file CONTROL was not found in $dir\n";
  } else {
    if($growth) {
      if($dir=~/min/i) {
	open(CONTROL, "<", "CONTROL_min") or print "file CONTROL_min was not found\n";
      } elsif($dir=~/run/i) {
	open(CONTROL, "<", "CONTROL_run") or print "file CONTROL_run was not found\n";
      } else {
	open(CONTROL, "<", "CONTROL_growth") or print "file CONTROL_growth was not found\n";
      }
    } else {
      open(CONTROL, "<", "$dir/CONTROL") or print "file CONTROL was not found in $dir\n";
    }
  }
  $i=0;
  # read the time step interval out of the CONTROL file
  while(<CONTROL>) {
    if(/\s*timestep/) {
      ($timestep) = /\s*timestep\s+(\S*)/;
    } elsif(/traj/) {
      ($trajinterval) = /\s*traj\S*\s+\S+\s+(\S*)/;
    }
  }
  close(CONTROL);
  while(1) {
    last if(read_statis_timestep($fhstatis,0) != 0);
    last if($sdata[0][0]>$endframe and not $lheatup and not $growth);
    $f++;
    next if($f % $resolution != 0);
    if(defined($step)) {
      $time = $f*$step;
    } else {
      $time  = $sdata[0][1];
      $time += $starttime;
      if($time<$time_old) {
	$time += $time_old-$starttime;
	$starttime = $time_old;
      }
    }
    if($lheatup or $sdata[0][0]>=$startframe) {
      print_statis_data($fhout,0,$time);
    }
    $time_old = $time;
  } # end while
  close($fhstatis);
} # end for directories
close($fhout);


sub error {
  print $_[1],"\n";
  exit $_[0];
}
