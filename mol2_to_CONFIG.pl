#!/usr/bin/perl
if($#ARGV<1) {
  print "1. name of input mol2-file\n";
  print "2. name of output CONFIG-file\n";
  print "[3. CONFIG-key]\n";
  print "[4. start index for atom numbers]\n";
  print "[5. shift vector]\n";
  exit;
}
open(IN_MOL2, "+<", $ARGV[0]) or die "Cant open $ARGV[0]: $!";
open(CONFIG, ">", $ARGV[1]) or die "Cant open $ARGV[1]: $!";
$configkey = $ARGV[2];
$startindex = $ARGV[3];
$vec[0] = $ARGV[4];
$vec[1] = $ARGV[5];
$vec[2] = $ARGV[6];
$i=0;

until(/ATOM/) {
  $_ = <IN_MOL2>;
  $string = $_;
}

$atoms=0;
$bonds=0;
while(<IN_MOL2>) {
  $string = $_;
  if(/BOND/) {
    last;
  } else {
    ($at_index[$atoms], $name[$atoms], $x[$atoms], $y[$atoms], $z[$atoms],
     $at_type[$atoms], $resnum[$atoms], $resname[$atoms]) =
    /^\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    $atoms++;
  }
}

while(<IN_MOL2>) {
  $string = $_;
  if(/SUBSTRUCTURE/) {
    last;
  } else {
    ($bd_index[$bonds], $bd_id1[$bonds], $bd_id2[$bonds], $bd_type[$bonds]) =
    /^\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/;
    $bonds++;
  }
}
$i=0;
while(<CHARGES>){
  $string = $_;
  ($charge[$i]) = /^\s*\S+\s+\S+\s+(\S+)/;
  $i++;
}

printf CONFIG $ARGV[0]."\n";
printf CONFIG "%10u", $configkey; # CONFIG file key
printf CONFIG "%10u\n", 2; # periodic boundary key
printf CONFIG " %19.12e %19.12e %19.12e\n", 100,  0,  0;
printf CONFIG " %19.12e %19.12e %19.12e\n",   0,100,  0;
printf CONFIG " %19.12e %19.12e %19.12e\n",   0,  0,100;
$j=0;
for($i=0;$i<@at_index;$i++) {
  printf CONFIG "%-8s", uc($at_type[$i]);
  printf CONFIG " %9u\n", $i+$startindex;
  printf CONFIG " %19.12e",   $x[$i]+$vec[0];
  printf CONFIG " %19.12e",   $y[$i]+$vec[1];
  printf CONFIG " %19.12e\n", $z[$i]+$vec[2];
  for($j=0;$j<$configkey;$j++) {
    printf CONFIG " %19.12e %19.12e %19.12e\n", 0, 0, 0;
  }
}
close IN_MOL2;
close CONFIG;