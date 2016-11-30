#!/usr/bin/perl

use hanno_utility;
use strict vars;
use Switch;

our($outfile,$lsortnum,$row,$fhin,$fhout,@lines,$cmtsign,@comment);
my($i);

$outfile  = "";
$lsortnum = 0;
$row      = -1;
$cmtsign  = '#';

if($#ARGV<0) {
  print "input format:\n";
  print "1. name of input file\n";
  print " -o <str>   specify output file\n";
  print " -r <int>   specify row (starting from 1)\n";
  print " -c <str>   comment sign (default: $cmtsign)\n";
  print " -n         sort numerically instead of lexically\n";
  exit 1;
}

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-o/) {
      $i++;
      $outfile=$ARGV[$i];
      if($outfile eq $ARGV[0]) {
	print "**** error: input and output file must not be equal!\n";
	exit 1;
      }
    } case (/^-c/) {
      $i++;
      $cmtsign=$ARGV[$i];
    } case (/^-r/) {
      $i++;
      $row=$ARGV[$i];
      if(not check_integer($row) or $row<=0) {
	print "**** error: row must be an integer number > 0!";
	exit 1;
      }
    } case (/^-n/) {
      $lsortnum=1;
    }
  }
}

open($fhin,"<",$ARGV[0]) or die "**** error: can't open input file $ARGV[0]: $!";
if($outfile) {
  open($fhout,">",$outfile) or die "**** error: can't open output file $outfile: $!";
}

$i=0;
while($_=<$fhin>) {
  if(/\s*$cmtsign/) {
    push(@comment,$_); next;
  }
  if(length(trim($_))==0) {
    &sort_and_print(1);
    $i=0;
    next;
  }
  $lines[$i][0]=$_;
  if($row>0 or $lsortnum) {
    $_=~s/^\s+//; $_=~s/\s+$//;
    push(@{$lines[$i]},split(/\s+/));
  }
  $i++;
}
&sort_and_print();
close($fhin,$fhout);


sub sort_and_print {
  my($i);
  if($row>0) {
    if($lsortnum) {
      @lines = sort { $a->[$row] <=> $b->[$row] } @lines;
    } else {
      @lines = sort { $a->[$row] cmp $b->[$row] } @lines;
    }
  } else {
    if($lsortnum) {
      @lines = sort { $a->[1] <=> $b->[1] } @lines;
    } else {
      @lines = sort { $a->[0] cmp $b->[0] } @lines;
    }
  }
  for($i=0;$i<@comment;$i++) {
    if($outfile) {
      print $fhout $comment[$i];
    } else {
      print $comment[$i];
    }
  }
  for($i=0;$i<@lines;$i++) {
    if($outfile) {
      print $fhout $lines[$i][0];
    } else {
      print $lines[$i][0];
    }
  }
  if(@_) {
    if($outfile) {
      print $fhout "\n";
    } else {
      print "\n";
    }
  }
  undef @lines;
  undef @comment;
}


sub ltrim { my $s = shift; $s =~ s/^\s+//;       return $s };
sub rtrim { my $s = shift; $s =~ s/\s+$//;       return $s };
sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };