#!/usr/bin/perl

use hanno_utility;

if($#ARGV<2) {
  print "format:\n";
  print " 1.  input file\n";
  print " 2.  output file\n";
  print " 3.  list of rows to logarithmize in the following format:\n";
  print "     [+/-]<rownumber>[+/-shift] (e.g. -2+1.643e5)\n";
  exit 1;
}

@targetrows = @ARGV[2..$#ARGV];

for($i=0;$i<@targetrows;$i++) {
  $_ = $targetrows[$i];
  if(/^[+-]?(\d+)$/) {
    $targetrows[$i] = $1-1;
    $shift[$i]      = 0;
    if(/^-/) { $negative[$i]=1; } else { $negative[$i]=0; }
  } elsif (/^[+-]?(\d+)([+-]?\d+\.?\d*([eEdD][+-]?\d+)?)$/) {
    $targetrows[$i] = $1-1;
    $shift[$i]      = $2;
    if(/^-/) { $negative[$i]=1; } else { $negative[$i]=0; }
  } else {
    print "**** error: input format '$_' for row #",($i+1)," was not recognized!";
    exit 1;
  }
}

open(IN, '<',$ARGV[0]) or die "**** error: can't open input file $ARGV[0]:\n$!\n";
open(OUT,'>',$ARGV[1]) or die "**** error: can't open output file $ARGV[1]:\n$!\n";

while(<IN>) {
  $line=$_;
  $_ =~ s/^\s+//;  $_ =~ s/\s+$//;
  @data = split(/\s+/);
  if($line =~ /^\s+#/ or not check_real($data[0])) {
    print OUT $line;
  } else {
    for($i=0;$i<@targetrows;$i++) {
      if($negative[$i]) {
        $data[$targetrows[$i]] = log(-$data[$targetrows[$i]]+$shift[$i]);
      } else {
        $data[$targetrows[$i]] = log($data[$targetrows[$i]]+$shift[$i]);
      }
    }
    for($i=0;$i<@data;$i++) {
      printf OUT " %20.13e",$data[$i];
    }
    print OUT "\n";
  }
}

close(IN,OUT)
