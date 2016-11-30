#!/usr/bin/perl
if($#ARGV<2) {
  print "1. name of input mol2-file\n";
  print "1. name of output mol2-file\n";
  print "2. file with charges from gaussian\n";
  exit;
}
open(IN_MOL2, "+<", $ARGV[0]) or die "Cant open $ARGV[0]: $!";
open(OUT_MOL2, ">", $ARGV[1]) or die "Cant open $ARGV[1]: $!";
open(CHARGES, "<", $ARGV[2]) or die "Cant open $ARGV[2]: $!";
$i=0;
$_="asdf";
until(/ATOM/) {
  $_ = <IN_MOL2>;
  $string = $_;
  print OUT_MOL2 $_;
}
print "start reading atoms\n";
$atoms=0;
$bonds=0;
while(<IN_MOL2>) {
  $string = $_;
  if(/BOND/) {
    last;
  } else {
    ($at_index[$atoms], $name[$atoms], $x[$atoms], $y[$atoms], $z[$atoms],
     $at_type[$atoms], $resnum[$atoms], $resname[$atoms]) =
    /\s*(\d*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)/;
    $atoms++;
  }
}
print "start reading bonds\n";
while(<IN_MOL2>) {
  $string = $_;
  if(/SUBSTRUCTURE/) {
    last;
  } else {
    ($bd_index[$bonds], $bd_id1[$bonds], $bd_id2[$bonds], $bd_type[$bonds]) =
    /\s*(\d*)\s*(\S*)\s*(\S*)\s*(\S*)\s*/;
    $bonds++;
  }
}
$i=0;
print "start reading charges\n";
while(<CHARGES>){
  $string = $_;
  ($charge[$i]) = /\s*\S*\s*\S*\s*(\S*)/;
  $i++;
}
print "writing output\n";
for($i=0;$i<@at_index;$i++) {
  printf OUT_MOL2 "%7u ", $at_index[$i];
  printf OUT_MOL2 "%-8s ",$name[$i];
  printf OUT_MOL2 "%9.4f ",$x[$i];
  printf OUT_MOL2 "%9.4f ",$y[$i];
  printf OUT_MOL2 "%9.4f ",$z[$i];
  printf OUT_MOL2 "%-4s ",$at_type[$i];
  printf OUT_MOL2 "%6u ",$resnum[$i];
  printf OUT_MOL2 "%-7s ",$resname[$i];
  printf OUT_MOL2 "%9.6f\n",$charge[$i];
}
print OUT_MOL2 "@<TRIPOS>BOND\n";
for($i=0;$i<@bd_index;$i++) {
  printf OUT_MOL2 "%6u ", $bd_index[$i];
  printf OUT_MOL2 "%4u ", $bd_id1[$i];
  printf OUT_MOL2 "%4u ", $bd_id2[$i];
  printf OUT_MOL2 "%-4s\n", $bd_type[$i];
}
print OUT_MOL2 "@<TRIPOS>SUBSTRUCTURE\n";
while(<IN_MOL2>) {
  $string = $_;
  print OUT_MOL2 $string;
}
close IN_MOL2;
close OUT_MOL2;
close CHARGES;