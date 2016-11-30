#!/usr/bin/perl

# calculate average values and standard deviation of datasets from a data file

use hanno_utility;
use Symbol;
use Switch;

if($#ARGV<0) {
  print "format:\n";
  print " 1. name of output file\n";
  print " 2. list of input files\n";
  print "-r <row number(s)> (standard: all, first=1, format: file.row (e.g. 1.2)\n";
  print "   or expression in quotes '' with file.row e.g. '\$2.5+\$1.5'\n";
  print "-s <start value for x>\n";
  print "-e <end   value for x>\n";
  print "-a append to file\n";
  exit 1;
}

$outfilename = $ARGV[0];
$lappend     = 0;

@targetrows = ();
for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-r/) {
      while($i+1<@ARGV and not $ARGV[$i+1]=~/^-/) {
	$i++;
	if(check_integer($ARGV[$i])) {
	  $ARGV[$i]--;
	  for($j=0;$j<@infilenames;$j++) {
	    push(@targetrows,"\$linedata[$j][$ARGV[$i]]");
	  }
	} elsif($ARGV[$i] =~ /^\d+.\d+$/) {
	  $ARGV[$i] =~ s/(\d+).(\d+)/\$linedata[\$$1][\$$2]/g;
	  $ARGV[$i] =~ s/\$(\d+)/($1-1)/eg;
	  push(@targetrows,$ARGV[$i]);
	} else {
	  $expr = $ARGV[$i];
	  $expr =~ s/\$(\d+).(\d+)/\$linedata[\$$1][\$$2]/g;
	  $expr =~ s/\$(\d+)/($1-1)/eg;
	  push(@targetrows,$expr);
	}
      }
    } case (/^-a/) {
      $lappend=1;
    } case (/^-s/) {
      $i++;
      die "**** error: start value must be real!\n" if(not check_real($ARGV[$i]));
      $startval = $ARGV[$i];
    } case (/^-e/) {
      $i++;
      die "**** error: end value must be real!\n" if(not check_real($ARGV[$i]));
      $endval = $ARGV[$i];
    } case(/^-/) {
      print "**** error: unknown flag $ARGV[$i]!\n";
    } else {
      push(@infilenames,$ARGV[$i]);
    }
  }
}

# print "input files:\n",join("\n",@infilenames),"\n";
# print "output rows:\n",join("\n",@targetrows),"\n";

# open files
if(not @infilenames) {
  print "**** error: no input files specified!\n";
  exit 1;
}
for($i=0;$i<@infilenames;$i++) {
  local *FILE;
  open(FILE, '<',$infilenames[$i]) or die "**** error: can't open input file $infilenames[$i]:\n$!\n";
  push(@fh,*FILE);
}
if($lappend) {
  open($fhout,">>",$outfilename) or die "**** error: can't open output file $outfilename:\n$!\n";
} else {
  open($fhout,">",$outfilename) or die "**** error: can't open output file $outfilename:\n$!\n";
}

#read data
$stop = 0;
$linenr = 0;
readlines:while($stop==0) {
  $linenr++;
  @linedata = ();
  $comment  =  0;
  for($i=0;$i<@fh;$i++) {
    $file = $fh[$i];
    if(not $_=<$file>) {
      $stop=1;
      last readlines;
    }
    $line=$_;
    $_ =~ s/^\s+//;  $_ =~ s/\s+$//;
    @{$linedata[$i]} = split(/\s+/);
    if(not @{$linedata[$i]} or $linedata[$i][0] =~ /^#/ or not check_real($linedata[$i][0])) {
      print $fhout $line;
      $comment=1;
    }
  }
  next if($comment);
  for($i=0;$i<@fh;$i++) {
    if($linedata[$i][0] != $linedata[0][0]) {
      print "**** warning: first rows don't match on line $linenr!\n";
      last;
    }
  }
  if(defined($startval)) {
    next if($linedata[0][0]<$startval);
  }
  if(defined($endval)) {
    last if($linedata[0][0]>$endval);
  }
  if(@targetrows) {
    for($i=0;$i<@targetrows;$i++) {
      printf $fhout "%20.13g",eval($targetrows[$i]);
    }
  } else {
    printf $fhout "%20.13g",$linedata[0][0];
    for($i=0;$i<@fh;$i++) {
      for($j=1;$j<@{$linedata[$i]};$j++) {
	printf $fhout "  %20.13g",$linedata[$i][$j];
      }
    }
  }
  print $fhout "\n";
}

close(@fh,$fhout);