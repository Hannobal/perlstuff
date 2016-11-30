#!/usr/bin/perl
use feature state;
use hanno_utility;
use dlpoly_utility;
use lammps_utility;
use Math::Trig;
use Cwd;

# use strict;

# my $str="10001000    291.19786    833.21939    432.27389    +1655.0325E2   -1410.3391\n";
# print "asdf\n" if($str=~/^[0-9\s\-\+eE.]+$/);
# 
@arr=(0,1,2,3,4,5,6,7,8,9);
splice(@arr,4,0,splice(@arr,5,1));
print join(" ",@arr),"\n";