#!/usr/bin/perl
use PDL;
use strict(vars);

# use hanno_utility;
# use dlpoly_utility;
# use Math::Trig;
# use aloxsam_utility;
# use List::Util 'shuffle';
# use Storable qw(dclone);

&loop;

sub loop {
  # how to loop over a piddle
  # !!! note: indices are reversed !!!
  my($test,@dims,$i,$j);
  $test = pdl [[1,2,3],[4,5,6],[7,8,9]];
  @dims=$test->dims;
  for($i=0;$i<$dims[0];$i++) {
    for($j=0;$j<$dims[1];$j++) {
      print $test->index($j)->index($i),"\n";
    }
  }
}