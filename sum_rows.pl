#!/usr/bin/perl

# calculate average values and standard deviation of datasets from a data file

use hanno_utility;
use Switch;

if($#ARGV<0) {
  print "format:\n";
  print " 1. input file\n";
  print "-r <row number(s)> (standard: all, first=1)\n";
  print "   or expression in quotes '' with row as e.g. \'\$1+5*\$2\'\n";
  print "-s <start value for x>\n";
  print "-e <end   value for x>\n";
  exit 1;
}

@targetrows = ();
for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-r/) {
      while($i+1<@ARGV and not $ARGV[$i+1]=~/^-/) {
	$i++;
	push(@targetrowsstr,$ARGV[$i]);
	if(not check_integer($ARGV[$i])) {
	  $expr = $ARGV[$i];
	} else {
	  $expr = '$'.$ARGV[$i];
	}
	$expr =~ s/\$(\d+)/\$linedata[\$$1]/g;
	$expr =~ s/\$(\d+)/($1-1)/eg;
	push(@targetrows,$expr);
      }
    } case (/^-s/) {
      $i++;
      die "**** error: start value must be real!\n" if(not check_real($ARGV[$i]));
      $startval = $ARGV[$i];
    } case (/^-e/) {
      $i++;
      die "**** error: end value must be real!\n" if(not check_real($ARGV[$i]));
      $endval = $ARGV[$i];
    } else {
      print "**** error: input $ARGV[$i] could not be interpreted!\n";
    }
  }
}

open(IN, '<',$ARGV[0]) or die "**** error: can't open input file $ARGV[0]:\n$!\n";

$numdata=0;
while(<IN>) {
  $line=$_;
  $_ =~ s/^\s+//;  $_ =~ s/\s+$//;
  @linedata = split(/\s+/);
  next if(not @linedata or $linedata[0] =~ /^#/);
  next if(not check_real($linedata[0]));
  if(defined($startval)) {
    next if($linedata[0]<$startval);
  }
  if(defined($endval)) {
    last if($linedata[0]>$endval);
  }
  if(@targetrows) {
    for($i=0;$i<@targetrows;$i++) {
      $sum[$i] += eval($targetrows[$i]);
    }
  } else {
    for($i=0;$i<@linedata;$i++) {
      $sum[$i] += $linedata[$i];
    }
  }
  $numdata++;
}

print "results from $ARGV[0]:\n";
printf "%4s  %20s  %20s  %20s\n","row","sum","# data points";

if(@targetrows) {
  for($i=0;$i<@targetrows;$i++) {
    printf "%10s  %20.13g  %20u\n",$targetrowsstr[$i],$sum[$i],$numdata;
  }
} else {
  for($i=0;$i<@sum;$i++) {
    printf "%4u  %20.13g  %20u\n",$i+1,$sum[$i],$numdata;
  }
}

close(IN);

sub max {
    my ($max, $next, @vars) = @_;
    return $max if not $next;
    return max( $max > $next ? $max : $next, @vars );
}