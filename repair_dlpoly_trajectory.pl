#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;
use POSIX;

# extracts the information from DL_POLY HISTORY files distributed
# over several folders with numbered names

use Cwd;
use Switch;

if($#ARGV<1) {
  print "the input format is:\n";
  print " 1. name of input HISTORY file\n";
  print " 2. name of corresponding FIELD file\n";
  print "[3. name of output HISTORY FILE]\n";
  exit 1;
}

$inhist  = $ARGV[0];
$infld   = $ARGV[1];
if($#ARGV>1){
  $outhist = $ARGV[2];
  if($inhist eq $outhist) {
    print "**** error: input and output HISTORY file must not be the same!\n";
    exit 1;
  }
}

&read_field_file($infld,0);
die "*** error: cannot open HISTORY file $inhist: $!\n" if(not open(HISTORY,"<",$inhist));

if(defined($outhist)) {
  # die "*** error: cannot open HISTORY file $outhist: $!\n" if(not open(OUT,">",$outhist));
}

$linenumber=1;
while(<HISTORY>) {
  last if(/timestep/);
  $linenumber++;
  print OUT $_ if(defined($outhist));
}

timestep : while(1) {
#   print $_;
  if(not /timestep\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
    print "**** error in line $linenumber:\n$_";
    exit 1;
  }
  $framenumber  = $1;
  $numexpect    = $2;
  $config_key   = $3;
  $periodic_key = $4;
  $numexpect*=(2+$config_key);
  $numexpect+=3 if($periodic_key>0);
  if(defined($outhist)) {
    @newlines=();
    ($newlines[0]) = /(timestep.*)/;
    $newlines[0] .= "\n";
  }
  for($i=0;$i<$numexpect;$i++) {
    last if(not $_=<HISTORY>);
    $linenumber++;
    if(/timestep/) {
      print "**** warning: timestep $framenumber ended too early ($i/$numexpect lines)!\n";
      print $_;
      next timestep;
    }
    push(@newlines,$_) if(defined($outhist));
  }
  # write output
  if(defined($outhist)) {
    for($i=0;$i<@newlines;$i++) {
      print OUT $newlines[$i];
    }
  }
  last if(not $_=<HISTORY>);
  $linenumber++;
}

close(HISTORY);
close(OUT) if(defined($outhist));
# 
