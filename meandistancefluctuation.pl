#!/usr/bin/perl
if($#ARGV<1) {
  print "1. HISTORY file to analyze\n";
  print "2. name of output file\n";
  print "3. optional: -n list of atom names (e.g. \"-n Al C3 H1)\"\n";
  print "             -i list of atom indices (e.g. \"-i 1-5 78 120-154\")\n";
  print "4. optional: -t timerange in framenumbers (e.g. \"-t 5000-200000\")\n";
  exit;
}

open(HISTORY, "<", $ARGV[0]) or die "cannot open HISTORY file: $!";
open(OUTPUT, ">", $ARGV[1]) or die "cannot open output file: $!";

$checktime=0;
$checknames=0;
$checkindices=0;

for($i=2;$i<=$#ARGV;$i++) {
  # read list of atom names
  if ($ARGV[$i] eq "-n") {
    $checknames=0;
    $i++;
    while(not substr($ARGV[$i],0,1) eq "-" and $i<$#ARGV) {
      $namelist[@namelist]=$ARGV[$i];
      $i++;
    }
  # read list of indices
  } elsif ($ARGV[$i] eq "-i") {
    $checkindices=1;
    $i++;
    while(not substr($ARGV[$i],0,1) eq "-" and $i<=$#ARGV) {
      if($ARGV[$i] =~ /^\d+$/) {
	$indexlist[@indexlist]=$ARGV[$i];
	$i++;
      } elsif($ARGV[$i] =~ /^\d+-\d+$/) {
	@argument=split('-',$ARGV[$i]);
	@argument=sort(@argument);
	for($j=$argument[0];$j<=$argument[1];$j++) {
	  $indexlist[@indexlist]=$j;
	}
	$i++
      } else {
	print "argument \"".$ARGV[$i]."\" is not an integer!\n";
	exit;
      }
    }
  # read list of timeframes
  } elsif ($ARGV[$i] eq "-t") {
    $checktime=1;
    $i++;
    while(not substr($ARGV[$i],0,1) eq "-" and $i<=$#ARGV) {
      if($ARGV[$i] =~ /^\d+$/) {
	$timelist[@timelist]=$ARGV[$i];
	$i++;
      } elsif($ARGV[$i] =~ /^\d+-\d+$/) {
	@argument=split('-',$ARGV[$i]);
	@argument=sort(@argument);
	for($j=$argument[0];$j<=$argument[1];$j++) {
	  $timelist[@timelist]=$j;
	}
	$i++
      } else {
	print "argument \"".$ARGV[$i]."\" is not an integer!\n";
	exit;
      }
    }
  } else {
    print "argument \"".$ARGV[$i]."\" could not be interpreted!\n";
    exit;
  }
}

if(not $checktime) {
  $evalframe=1
}

$timesteptot=0;
$frame=-1;
$timetot=0;
while(<HISTORY>) {
# reading frame out of history file
  if(/timestep/) {
    $frame++;
    #print "evaluating frame $frame\n";
    ($timestep[$frame], $numatoms[$frame], $trajkey[$frame], $periodkey[$frame], $integration[$i]) =
    /^timestep\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    if($numatoms[$frame] /= $field_numatoms) {
      print "error in HISOTRY file: $field_numatoms expected, $numatoms[$frame] given in HISTORY!\n";
      exit;
    }
    $_=<HISTORY>;
    ($cell[$frame][0][0],$cell[$frame][0][1],$cell[$frame][0][2])=/^\s*(\S+)\s+(\S+)\s+(\S+)/;
    ($cell[$frame][1][0],$cell[$frame][1][1],$cell[$frame][1][2])=/^\s*(\S+)\s+(\S+)\s+(\S+)/;
    ($cell[$frame][2][0],$cell[$frame][2][1],$cell[$frame][2][2])=/^\s*(\S+)\s+(\S+)\s+(\S+)/;
    for($i=0;$i<=$numatoms[$frame];$i++) {
      $_=<HISTORY>;
      ($name[$frame][$i], $index[$frame][$i], $mass[$frame][$i],$charge[$frame][$i]) = 
      /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
      $_=<HISTORY>;
      ($coords[$frame][$i][0],$coords[$frame][$i][1],$coords[$frame][$i][2]) =
      /^\s*(\S+)\s+(\S+)\s+(\S+)/;
      if($historykey > 0) {
        $_=<HISTORY>;
        ($vel[$frame][$i][0],$vel[$frame][$i][1],$vel[$frame][$i][2]) =
        /^\s*(\S+)\s+(\S+)\s+(\S+)/;
         if($historykey > 0) {
          $_=<HISTORY>;
          ($force[$frame][$i][0],$force[$frame][$i][1],$force[$frame][$i][2]) =
          /^\s*(\S+)\s+(\S+)\s+(\S+)/;
        }
      }
    }

#  evaluating frame
    #check whether or not to evaluate the frame
    if($checktime) {
      $evalframe=0;
      for($i=0$i<@timelist;$i++) {
	if($timestep[$frame]==$timelist[$i]) {
	  $evalframe=1
	}
      }
    }
    if($evalframe) {
      for($i=0;$i<=$numatoms[$frame];$i++) {
      
      }
    }
  }
}

