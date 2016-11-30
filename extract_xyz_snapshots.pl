#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;
use Switch;

# join restarted DLPOLY HISTORY files

if($#ARGV<2) {
  print "the input format is:\n";
  print " 1. input filename\n";
  print " 2. output filename\n";
  print " 3. frame (0=first -1=last)\n";
  exit 1;
}

$infilename  = $ARGV[0];
$outfilename = $ARGV[1];
$targetframe = $ARGV[2];
&error(1) if(not check_integer($targetframe));

&error(2,$infilename) if not(open($fhin,"<",$infilename));

$f=0;
$stat=0;
while($stat==0) {
  exit 1 if(read_xyz_timestep($fhin,-1,0)>0);
  &error(4,$infilename,$targetframe) if($stat<0 and $targetframe>=0);
  last if($f==$targetframe and $targetframe>=0);
  $f++;
}
close($fhin);

&error(2,$outfilename) if not(open($fhout,">",$outfilename));
write_xyz_timestep($fhout,0,$config_title[0]);
close($fhout);

sub error {
  if($_[0]==1) {
    print "**** error: target frame must be an integer number!\n";
  } elsif($_[0]==2) {
    print "**** error: could not open file $_[1]: $!\n";
  } elsif($_[0]==4) {
    print "**** error: did not find timestep $_[2] in file $_[1]!\n";
  } else {
    print "**** unknown error!\n";
    exit 99;
  }
  exit $_[0]
}