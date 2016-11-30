#!/usr/bin/perl

use hanno_utility;
use dlpoly_utility;
use lammps_utility;
use Switch;

use strict vars;

if($#ARGV<0) {
  print "input format:\n";
  print " 1. output directory\n";
  print "-f <str>     input filename (default: log.lammps)\n";
  print "-s <int>     start frame\n";
  print "-e <int>     end_frame\n";
  print "-h           toggle heatup mode\n";
  print "-r <int>     only print every n-th frame\n";
  print "-p <real>    ignore first row and use this step size\n";
  exit;
}

my $step       = 0;
my $outdir     = $ARGV[0];
my $infile     = "log.lammps";
my $resolution = 1;
my $lheatup    = 0;
my $startframe = 0;
my $endframe   = 9e20;
my($i,$j,@directories,$fhout,$dir);
for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case(/^-h/i) {
      $lheatup = 1;
    } case(/^-f/i) {
      $i++;
      $infile = $ARGV[$i];
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
    } case (/^-p/i) {
      $i++;
      $step = $ARGV[$i];
      &error(1,"step size must be a real number") if not check_real($step);
    } else {
      &error(1,"cannot interpret flag $ARGV[$i]");
    }
  }
}

@directories=get_numeric_directories(".",$startframe,$endframe);
push(@directories,".") if(-e "./$infile");

if(not @directories) {
  print "found no directories containing $infile.\n"; exit 1;
}

if(not -d $outdir) {
  mkdir $outdir or &error(4,"cannot create output directory:\n$!");
}
open($fhout, ">", "$outdir/STATIS_DATA") or &error(4,"cannot open file STATIS_DATA:\n$!");

my $firstfile=1;
my $f=-1;
foreach $dir (@directories) {
  print "reading $dir/$infile\n";
  if(not -e "$dir/$infile") {
    print "**** warning: file $dir/$infile was not found\n";
    next;
  }
  exit 1 if(read_logfile("$dir/$infile",0)!=0);
  if(not @{$ldata[0]}) {
    print "**** warning: no thermo data found in $dir/$infile\n";
    next;
  }
  if($firstfile) {
    for($i=0;$i<@{$ldata_names[0]};$i++) {
      printf $fhout "%14s",$ldata_names[0][$i];
    }
    print $fhout "\n";
    $firstfile=0;
  }
  for($i=0;$i<@{$ldata[0]};$i++) {
    $f++;
    next if($f % $resolution != 0);
    if($step==0) {
      for($j=0;$j<@{$ldata[0][$i]};$j++) {
        printf $fhout "%14s",$ldata[0][$i][$j];
      }
    } else {
      printf $fhout "%14.6g",$step*$f;
      for($j=1;$j<@{$ldata[0][$i]};$j++) {
        printf $fhout "%14s",$ldata[0][$i][$j];
      }
    }
    print $fhout "\n";
  }
}

sub error {
  print $_[1],"\n";
  exit $_[0];
}