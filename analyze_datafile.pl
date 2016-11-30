#!/usr/bin/perl

use hanno_utility;
use Switch;

$startframe = -9e20;
$endframe   = 9e20;

if($#ARGV<1) {
  print "the input format is:\n",
  " 1. input file\n",
  "-s <int>  start value in first column\n",
  "-e <int>  end value in first column\n",
  "-<expr>   <resolution> <output filename>\n\n",
  "expr can be a row number starting from 1 (e.g. 2) or an expression with row e.g. '\$2+\$1/2'\n",
  "you may have to put the argument in quotes (e.g. '-(\$2+\$1)/2)'\n";

  exit 1;
}

@cols=();

$infilename = $ARGV[0];
for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-s$/) { # start frame
      $i++;
      $startframe = $ARGV[$i];
    } case (/^-e$/) { # end frame
      $i++;
      $endframe = $ARGV[$i];
    } else {
      push(@colstr,substr($ARGV[$i],1));
      if($ARGV[$i]=~/^-(\d+)$/) {
        push(@cols,"\$linedata[".($1-1)."]");
      } else {
        $expr = substr($ARGV[$i],1); # cut the "-" from the flag
        $expr =~ s/\$(\d+)/\$linedata[\$$1]/g;
        $expr =~ s/\$(\d+)/($1-1)/eg;
        push(@cols,$expr);
      }
    }
  }
}
# print "min ",join(" ",@mincols),"\n";
# print "max ",join(" ",@maxcols),"\n";
# print "Dev ",join(" ",@stdevcols),"\n";
foreach $i (@mincols,@maxcols,@stdevcols) {
  if(not check_integer($i) or $i<1) {
    print "**** errror: columns must be integer numbers greater than zero!\n";
    exit 1;
  }
}

open(INPUT,"<",$infilename) or die "**** error: Can't open input file!";
$line=0;
@coldata=();
while(<INPUT>) {
  $_ =~ s/^\s+//;
  $_ =~ s/\s+$//;
  next if(length($_)==0);
  next if(/^#/);
  @linedata = split(/\s+/,$_);
  next if(not check_real($linedata[0]));
  next if($linedata[0]<$startframe);
  next if($linedata[0]>$endframe);
  for($i=0;$i<@cols;$i++) {
    push(@{$coldata[$i]},eval($cols[$i]));
  }
}
close(INPUT);

$maxlen=4;
$maxlen=max($maxlen,$_) foreach @colstr;

printf "%-$maxlen\s %17s %17s %17s %17s\n","expr","min","max","average","stdev";
for($i=0;$i<@cols;$i++) {
  $minval=$coldata[$i][0];
  $maxval=$coldata[$i][0];
  for($j=1;$j<@{$coldata[$i]};$j++) {
    $minval=$coldata[$i][$j] if($coldata[$i][$j]<$minval);
    $maxval=$coldata[$i][$j] if($coldata[$i][$j]>$maxval);
  }
  printf "%-$maxlen\s %17.10g %17.10g %17.10g %17.10g\n",
  $colstr[$i],$minval,$maxval,calc_stdev(@{$coldata[$i]});
}
