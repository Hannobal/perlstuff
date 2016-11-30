#!/usr/bin/perl

use hanno_utility;

if($#ARGV<2) {
  print "input format:\n";
  print "1. input file name\n";
  print "2. output file name\n";
  print "3. histogram step size\n";
  exit 1;
}

$infilename=$ARGV[0];
$outfilename=$ARGV[1];
$histstep=$ARGV[2];
if(not check_real($histstep)) {
  print "**** error: input for histogram step size ($histstep) is not a real number!\n";
  exit 1;
}

open(INPUT,"<",$infilename) or die "**** error: Can't open input file!";

$numblocks = -1;
$newblock  = 1;
while(<INPUT>) {
  $_ =~ s/^\s+//;
  $_ =~ s/\s+$//;
  if(length($_)==0) {
    $newblock = 1;
  } elsif(not /^#/) {
    if($newblock) {
      $numblocks++;
      $newblock = 0;
      $numcols[$numblocks]     = 0;
      $numdatasets[$numblocks] = 0;
    }
    @{$data[$numblocks][$numdatasets[$numblocks]]} = split(/\s+/,$_);
    if(@{$data[$numblocks][$numdatasets[$numblocks]]} > $numcols[$numblocks]) {
      $numcols[$numblocks] = @{$data[$numblocks][$numdatasets[$numblocks]]};
    }
    $numdatasets[$numblocks]++;
  } else {
    push(@{$comment[$numblocks]},$_);
  }
}
close(INPUT);

$maxnumhistentries=0;
for($b=0;$b<=$numblocks;$b++) {
  for($c=1;$c<$numcols[$b];$c++) {
  $histmaxindex[$b][$c] = -99999999999999;
  $histminindex[$b][$c] =  99999999999999;
    for($s=0;$s<$numdatasets[$b];$s++) {
      $i = int($data[$b][$s][$c]/$histstep);
      if($i>$histmaxindex[$b][$c]) {
	$histmaxindex[$b][$c] = $i;
      }
      if($i<$histminindex[$b][$c]) {
	$histminindex[$b][$c] = $i;
      }
      $histogram[$b][$c]{$i}++;
    }
  }
}

open(OUTPUT,">",$outfilename) or die "**** error: Can't open input file!";
for($i=$histminindex[0][1];$i<$histmaxindex[0][1];$i++) {
  print "asdf\n";
# for($c=1;$c<5;$c++) {
#   for($b=0;$b<=$numblocks;$b++) {
    printf OUTPUT "%20.10g %8d %8d\n", $i*$histstep,$histogram[0][1]{$i};
#   }
}
close(OUTPUT);
