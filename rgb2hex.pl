#!/usr/bin/perl

if($#ARGV < 2) {
  print "input format:\n";
  print "<red> <green> <blue> [alpha]\n";
  exit 0;
}

$str="";
foreach $v (@ARGV) {
  $str .= sprintf "%02x", $v;
}
print "$str\n";