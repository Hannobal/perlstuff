#!/usr/bin/perl

use Math::Trig;
use Math::MatrixReal;
use dlpoly_utility;
use hanno_utility;
use POSIX qw(ceil floor);

open($fhin,"<",$ARGV[0]);
open($fhout,">",$ARGV[1]);

while(<$fhin>) {
  if($_=~/^\s*#/) {
    print $fhout $_;
  } else {
    $_=~s/^\s+//;$_=~s/\s+$//;
    @linedata=split(/\s+/,$_);
    if(not @old) {
      @old=@linedata; next;
    }
    printf $fhout "%17.10g", $linedata[0];
    for($i=1;$i<@linedata;$i++) {
      printf $fhout " %17.10g", $linedata[$i]-$old[$i];
      $old[$i]=$linedata[$i];
    }
    print $fhout "\n";
  }
}

close($fhout,$fhin);