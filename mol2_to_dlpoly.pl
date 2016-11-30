#!/usr/bin/perl
$mol2_file = $ARGV[0]
$gaff_file = $ARGV[1]
#while($i<=$#list) {
#  print "$list[$i] ";
#  $i++;
#}
#print "\n";
open(MOL2, "<", $ARGV[0]) or die "Cant open $ARGV[0]: $!";
until(/ATOM/) {
  $string = $_;
  $_ = <MOL2>;
}
while(<MOL2>) {
  $string = $_;
  if(/BOND/) {
    last;
  } else {
    ($id, $name, $pos_x, $pos_y, $pos_z, $type) =
    /\s*(\d*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s/;
    $x{$id} = $pos_x;
    $y{$id} = $pos_y;
    $z{$id} = $pos_z;
    $atom_type{$id} = $type;
  }
}
while(<MOL2>) {
  $string = $_;
  if(/SUBSTRUCTURE/) {
    last;
  } else {
    ($bond_id, $id1, $id2, $type) =
    /\s*(\d*)\s*(\S*)\s*(\S*)\s*(\S*)\s*/;
    $bond_id1{$bond_id} = $id1;
    $bond_id2{$bond_id} = $id2;
    $bond_type{$bond_id} = $type;
    $j=0;
    for($i=0; $i<=$#list; $i++) {
      if($list[$i]==$id1 or $list[$i]==$id2) {
        $j++;
      }
    }
    if($j==2) {
      $dx = $x{$id1}-$x{$id2};
      $dy = $y{$id1}-$y{$id2};
      $dz = $z{$id1}-$z{$id2};
      $bond_length{$bond_id} = sqrt($dx**2+$dy**2+$dz**2);
      printf "%6d", $bond_id;
      printf "%5s", $id1;
      printf "%5s", $id2;
      printf " %-4s  ", $type;
      printf "%7.5f\n", $bond_length{$bond_id};
      #print " $bond_id $id1 $id2 $type $bond_length{$bond_id}\n";
    }
  }
}