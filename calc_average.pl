#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;

if($#ARGV<1) {
  print "format:\n";
  print " 1. input file\n";
  print " 2. list of rows to average\n";
  exit 1;
}

@targetrows = @ARGV[1..$#ARGV];

for($i=0;$i<@targetrows;$i++) {
  if($targetrows[$i]<1 or not check_integer($targetrows[$i])) {
    print "**** error: rows must be integers starting from 1\n";
    exit 1;
  }
  $targetrows[$i]--;
  $average[$i]=0;
}

open(IN, '<',$ARGV[0]) or die "**** error: can't open input file $ARGV[0]:\n$!\n";

$numdata=0;
$linenum=0;
while(<IN>) {
  $linenum++;
  $line=$_;
  $_ =~ s/^\s+//;  $_ =~ s/\s+$//;
  @data = split(/\s+/);
  if(check_real($data[0])) {
    $numdata++;
    for($i=0;$i<@targetrows;$i++) {
      if(not check_real($data[$targetrows[$i]])) {
	print "**** error on line $linenum: expected data but found '$data[$targetrows[$i]]'\n";
	exit 1;
      }
      $average[$i]+=$data[$targetrows[$i]];
    }
  }
}

for($i=0;$i<@targetrows;$i++) {
  $targetrows[$i]++;
  $average[$i] /= $numdata;
  print "average($targetrows[$i]) = $average[$i]\n";
}

close(IN,OUT)
