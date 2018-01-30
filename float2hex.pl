#!/usr/bin/perl

# use Math::Round;

if($#ARGV < 2) {
  print "input format:\n";
  print "<red> <green> <blue> [alpha]\n";
  exit 0;
}

$str="";
foreach $v (@ARGV) {
  die "**** error: values must be in range [0:1]!" if($v>1);
  $str .= sprintf "%02x", sprintf("%.0f",$v*255);
}
print "$str\n";