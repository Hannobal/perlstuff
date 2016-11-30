#!/usr/bin/perl
if($#ARGV<1) {
  print "1. name of input mol2-file\n";
  print "2. name of output CONFIG-file\n";
  print "[3. CONFIG-key]\n";
  print "[4. start index for atom numbers]\n";
  exit;
}
open(IN_MOL2, "+<", $ARGV[0]) or die "Cant open $ARGV[0]: $!";
open(CONFIG, ">", $ARGV[1]) or die "Cant open $ARGV[1]: $!";
open(POS, "<", "pos") or die "Cant open pos: $!";
$configkey = $ARGV[2];
print $configkey."\n";
$startindex = $ARGV[3];
$i=0;
$_="asdf";
until(/ATOM/) {
  $_ = <IN_MOL2>;
  $string = $_;
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
print "start reading pos\n";
while(<POS>){
  $string = $_;
  ($vecx[$i],$vecy[$i],$vecz[$i]) = /\s*(\S+)\s+(\S+)\s+(\S+)/;
  $i++;
}
print "writing output\n";
$j=0;
for($a=0;$a<@vecx;$a++) {
  for($i=0;$i<@at_index;$i++) {
    printf CONFIG "%8s", uc($at_type[$i]);
    printf CONFIG " %9u\n", $i+$startindex;
    printf CONFIG " %19.12e",   $x[$i]+$vecx[$a];
    printf CONFIG " %19.12e",   $y[$i]+$vecy[$a];
    printf CONFIG " %19.12e\n", $z[$i]+$vecz[$a];
    for($j=0;$j<$configkey;$j++) {
      printf CONFIG " %19.12e %19.12e %19.12e\n", 0, 0, 0;
    }
  }
}
close IN_MOL2;
close CONFIG;
close POS;