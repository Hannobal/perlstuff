#!/usr/bin/perl

# use strict vars;
use hanno_utility;
use dlpoly_utility;
use aloxsam_utility;

if($#ARGV<1) {
  print "input format:\n";
  print " 1. CONFIG file\n";
  print " 2. FIELD file\n";
  print "[3. minimum distance, default: 3.5]\n";
  print "[4. top/bottom (1/-1) default: 1]\n";
  exit 1;
}

$incfg = $ARGV[0];
$infld = $ARGV[1];

if($#ARGV>1) {
  $minpadist=$ARGV[2];
  if(not check_real($minpadist)) {
    print "****error: minimum distance must be a real number!\n";
    exit 1;
  }
} else {
  $minpadist = 3.5;
}

if($#ARGV>2) {
  $top=$ARGV[3];
  if(not ($top eq "-1" or $top eq "0" or $top eq "1")) {
    print "****error: minimum distance must be a real number!\n";
    exit 1;
  }
} else {
  $top=1;
}

exit 1 if(read_field_file($infld,0)!=0);
exit 1 if(read_config_file($incfg,0,0)!=0);

@candidates=find_hydroxide_candidates(0,0,$top,$minpadist);
print "found ",$#candidates+1," candidates";