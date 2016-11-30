#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;
use Switch;

if($#ARGV==-1) {
  print "the input format is:\n";
  print "<output-filename>\n";
  print "-s <int>   start frame\n";
  print "-e <int>   end frame\n";
  print "-r reverse (cool down)\n";
  exit 1;
}

$start = 1;
$end   = 9e20;
for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-s/) {
      $i++;
      $start = $ARGV[$i];
      if(not check_integer($start)) {
	print "**** error: start frame must be integer number!\n";
	exit 1;
      }
    }
    case (/^-e/) {
      $i++;
      $end = $ARGV[$i];
      if(not check_integer($end)) {
	print "**** error: end frame must be integer number!\n";
	exit 1;
      }
    }
    case (/^-r/) {
      $lreverse=1;
    }
    else {
      print "**** error: expected optional flag but got \"$ARGV[$i]\"!\n";
      exit 1;
    }
  }
}

if($lreverse) {
  @directories = reverse(get_numeric_directories(".",$start,$end));
} else {
  @directories = get_numeric_directories(".",$start,$end);
}

push(@directories,'.') if(-e './HISTORY');

open(HISTORY_OUT, ">", $ARGV[0]) or die "Cant open $ARGV[0]: $!";
$error = 0;
foreach $dir (@directories) {
  if(-e "./$dir/HISTORY") {
    open(HISTORY_IN, "<", "./$dir/HISTORY");
  } elsif(-e "./$i/HISTORY.dlpolyhist") {
    open(HISTORY_IN, "<", "./$dir/HISTORY.dlpolyhist");
  } else {
    print "did not find HISTORY file for $dir K!\n";
    next;
  }
  print "writing out $dir K\n";
  while(<HISTORY_IN>) {
    print HISTORY_OUT $_;
  }
}
