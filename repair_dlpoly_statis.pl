#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;
use POSIX;
use Cwd;
use Switch;

if($#ARGV<0) {
  print "the input format is:\n";
  print " 1. name of input STATIS file\n";
  print "[2. name of output STATIS FILE]\n";
  exit 1;
}

$instat = $ARGV[0];
if($#ARGV>0){
  $outstat = $ARGV[1];
  if($instat eq $outstat) {
    print "**** error: input and output STATIS file must not be the same!\n";
    exit 1;
  }
}

die "*** error: cannot open STATIS file $inhist: $!\n" if(not open(STATIS,"<",$instat));

$newlines[0]=<STATIS>;
$linenumber=1;
$line=$newlines[0];
$line =~ s/^\s+//; $line =~ s/\s+$//;
@linedata  = split(/\s+/,$line);
if(not check_real(@linedata)) {
  $newlines[1]=<STATIS>;
  $newlines[2]=<STATIS>;
  $line=$newlines[2];
  $line =~ s/^\s+//; $line =~ s/\s+$//;
  @linedata  = split(/\s+/,$line);
  $framestep = $linedata[0];
  $numenttot = $linedata[2];
  $timestep  = $linedata[1]/$linedata[0];
  $linenumber=3;
}

if(defined($outstat)) {
  die "*** error: cannot open STATIS file $outstat: $!\n" if(not open(OUT,">",$outstat));
  print OUT @newlines;
}

timestep : while(1) {
  if($#linedata != 2 or not check_real(@linedata)) {
    print "**** error on line $linenumber:\n$_";
    exit 1;
  }
  $framenumber = $linedata[0];
  $numexpect   = ceil($linedata[2]/5.0)-1;
  $numentlast  = $linedata[2]%5;
  if(defined($outstat)) {
    @newlines=();
    ($newlines[0]) = $_;
  }
  for($i=0;$i<=$numexpect;$i++) {
    last if(not $_=<STATIS>);
    $line=$_;
    $line =~ s/^\s+//; $line =~ s/\s+$//;
    @linedata   = split(/\s+/,$line);
    $linenumber++;
    $numints = contains_int(@linedata);
    if($numints>0) {
      #this line definitely contains the header of a new frame!
      print "**** warning: frame $framenumber ended too early ($i/$numexpect lines)!\n";
      if($numints>=2) {
	# then we can reconstruct the next frame
	for($i=$#linedata;$i>=0;$i++) {
	  if(check_integer($linedata[$i])) {
	    $last = $i;
	    last;
	  }
	}
	$_ = sprintf("%10u%14.6e%10u\n",@linedata[$last-2..$last]);
	@linedata=@linedata[$last-2..$last];
      } else {
	print "**** error: cannot reconstruct frame on line $linenumber:\n$_";
	exit 1;
      }
      next timestep;
    }
#     print join(" ",@linedata),"\n";
    push(@newlines,$_) if(defined($outstat));
  }
  # write output
  if(defined($outstat)) {
    for($i=0;$i<@newlines;$i++) {
      print OUT $newlines[$i];
    }
  }
  last if(not $_=<STATIS>);
  $line=$_;
  $line =~ s/^\s+//; $line =~ s/\s+$//;
  @linedata   = split(/\s+/,$line);
  $linenumber++;
}
print "done\n";

close(STATIS);
close(OUT) if(defined($outstat));

sub contains_int {
  $numints=0;
  for(my $i=0;$i<@_;$i++) {
    if ($_[$i] =~ /^[+-]?\d+$/ ) {
      $numints++;
    }
  }
  print "numints $numints\n";
  return $numints;
}
