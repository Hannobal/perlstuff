#!/usr/bin/perl

use hanno_utility;

if($#ARGV<1) {
  print "input format: <input filename> <output filename> [write mode ('>'/'>>')]\n";
  exit;
}

$infilename  = $ARGV[0];
$outfilename = $ARGV[1];
$writemode   = '>';
if($#ARGV>1) {
  $writemode = $ARGV[2]
}

if($infilename eq $outfilename) {
  print "**** error: input and output filenames must not be equal!";
}

open(INFILE, "<", $infilename) or die "**** error: can't open input file '$infilename':\n$!";
open(OUTFILE, ">", $outfilename) or die "**** error: can't open input file '$outfilename':\n$!";

$firstline=1;
while(<INFILE>) {
    $line = $_;
    $_ =~ s/^\s*(.*)\s*$/$1/;
    @indata = split(/\s+/, $_);
  if(/^\s*#/ or not check_real($indata[0])) { # comment line
    print OUTFILE $line;
  } else {
    if($#indata==-1) { #empty line
      print OUTFILE $_;
      $firstline=1;
    } else {
      if($firstline) {
	@indata_old = @indata;
	$firstline  = 0;
      } else {
	printf OUTFILE "%20.13e", ($indata[0]+$indata_old[0])/2.0;
	for($i=1;$i<@indata;$i++) {
	  printf OUTFILE " %20.13e", ($indata[$i]-$indata_old[$i])/($indata[0]-$indata_old[0]);
	}
	print OUTFILE "\n";
      }
    }
  }
}
