#!/usr/bin/perl
use Math::Trig;

sub calc_angle {
  # $_[0] is absolute of vector
  # $_[1] is absolute of projected vector into plane
  # $_[2] is vector component out of projection plane
  if($_[1] == 0) {
    $angle = pi/2;
  } elsif($_[2] == 0) {
    $angle = 0;
  } else {
    $angle = acos($_[1]/($_[0]));
  }
  if($angle != 0) {
    $angle = $angle*180.0/pi*$_[2]/abs($_[2]);
  }
  return $angle;
}

if($#ARGV<4) {
  print "the input format is:\n";
  print "1. name of CONFIG/REVCON-file containing Al2O3-surface\n";
  print "2. name of mol2-file containing phosphonic acid\n";
  print "3. name of output-file\n";
  print "4. radius of island\n";
  print "5. distance from surface\n";
  print "6. leave center free (0=false, 1=true) (optional, default=0)\n";
  exit;
}
if($ARGV[3]<7.0) {
  print"That's pretty small for an island. Maybe you shoud make it bigger...\n";
  exit;
}
open(REVCON, "<", $ARGV[0])  or die "Failed to open HISTORY-file containing surface: $!\n";
open(MOLECULE, "<", $ARGV[1])  or die "Failed to open mol2-file: $!\n";
#open(CONFIG, ">", $ARGV[3]) or die "Failed to create outupt-file: $!\n";
$radius = $ARGV[3];
$distance = $ARGV[4];
$center_free = $ARGV[5];

# read out surface from REVCON-file

$title=<REVCON>;
$_=<REVCON>;
($REVCONkey,$periodic_key) = /\s*(\S*)\s*(\S*)\s*\S*\s*\S*/;
for($j=0;$j<=2;$j++) {
  $_=<REVCON>;
  ($vec[$j][0], $vec[$j][1], $vec[$j][2]) = /\s*(\S*)\s*(\S*)\s*(\S*)/;
}
$anzges = 0; #number of atoms counted in REVCON
$xmin = 0;
$ymin = 0;
$zmax = -9999;
while(<REVCON>) {
  $anzges++;
  ($name) = /(\S*)\s*\S*/;
  $_=<REVCON>;
  ($x{$name}[$anz{$name}],$y{$name}[$anz{$name}],$z{$name}[$anz{$name}]) = /\s*(\S*)\s*(\S*)\s*(\S*)/;
  #print "$z{$name}[$anz{$name}] ";
  if($z{$name}[$anz{$name}] > $zmax) {
    $zmax = $z{$name}[$anz{$name}];
  }
  for($i=0;$i<$REVCONkey;$i++) { #überspringe ggf. Zeilen für Geschwindigkeit und Kraft
    $_ = <REVCON>;
  }
  $anz{$name}++;
}
# read mol2 file
$_ = <MOLECULE>;
until(/ATOM/) {
  $string = $_;
  $_ = <MOLECULE>;
  print $_;
}
while(<MOLECULE>) {
  $string = $_;
  if(/BOND/) {
    last;
  } else {
    ($id, $name, $pos_x, $pos_y, $pos_z, $type) =
    /\s*(\d*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s/;
    $mol_pos[$id][0] = $pos_x;
    $mol_pos[$id][1] = $pos_y;
    $mol_pos[$id][2] = $pos_z;
    $mol_type[$id] = $type;
  }
}
while(<MOLECULE>) {
  $string = $_;
  if(/SUBSTRUCTURE/) {
    last;
  } else {
    ($bond_id, $id1, $id2, $type) =
    /\s*(\d*)\s*(\S*)\s*(\S*)\s*(\S*)\s*/;
    $mol_bond_id1[$bond_id] = $id1;
    $mol_bond_id2[$bond_id] = $id2;
    $mol_bond_type[$bond_id] = $type;
    #search for phosphonic acid group and remember oxygens
    if($mol_type[$id1] =~ m/p5/ or $mol_type[$id2] =~  m/p5/) {
      if($mol_type[$id1] =~ m/oh/) {
        $p5_oh=$id1;
      } elsif($mol_type[$id2] =~ m/oh/) {
        $p5_oh=$id2;
      } elsif($mol_type[$id1] =~ m/o/) {
        $p5_o[@p5_o]=$id1;
      } elsif($mol_type[$id2] =~  m/o/) {
        $p5_o[@p5_o]=$id2;
      }
    }
  }
}
for($i=0;$i<@p5_o;$i++) {
}

# rotate and align template phosphonic acid molecule
$dx = $mol_x[$p5_o[1]] - $mol_x[$p5_o[0]];
$dy = $mol_y[$p5_o[1]] - $mol_y[$p5_o[0]];
$dz = $mol_z[$p5_o[1]] - $mol_z[$p5_o[0]];
$angle_xy = calc_angle(sqrt($dx*$dx+$dy*$dy+$dz*$dz), sqrt($dx*$dx+$dy*$dy), $dz);
$angle_xz = calc_angle(sqrt($dx*$dx+$dy*$dy+$dz*$dz), sqrt($dx*$dx+$dz*$dz), $dy);
$angle_yz = calc_angle(sqrt($dx*$dx+$dy*$dy+$dz*$dz), sqrt($dy*$dy+$dz*$dz), $dx);

close(REVCON,MOLECULE);