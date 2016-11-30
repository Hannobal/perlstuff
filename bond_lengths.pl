#!/usr/bin/perl
$kcal = 4.184; #kJ
if($#ARGV==-1) {
  print "the input format is:\n";
  print "1. name of gaff-file\n";
  print "2. name of mol2-file\n";
  print "3. name of output file (optional)\n";
  print "4. list of atom numbers whose bond lengths shall be evaluated\n\n";
  print "example with output file:\n";
  print "perl bond_length.pl gaff.dat example.mol2 length_example 1-5 8 12 17-20 23\n";
  print "example without output file (output is written to console):\n";
  print "perl bond_length.pl gaff.dat example.mol2 1-5 8 12 17-20 23\n\n";
  exit;
}
$i=2;
$outf = 0;
if ($ARGV[$i] == 0 && $$ARGV[$i] ne "0")  {
  open(OUTFILE, ">", $ARGV[$i]) or die "Cant open $ARGV[$i]: $!";
  $outf = 1;
  $i++;
}
while($i<=$#ARGV ) {
  if($ARGV[$i] =~ /-/) {
    ($j, $last) = split(/-/, $ARGV[$i]);
    if($last<$j) {
      die "$j is larger than $last";
    }
    while($j <= $last) {
      $list[@list] = $j;
      $j++;
    }
  } else {
    $list[@list] = $ARGV[$i];
  }
  $i++;
}
$i = 0;
open(GAFF, "<", $ARGV[0]) or die "Cant open $ARGV[0]: $!";
open(MOL2, "<", $ARGV[1]) or die "Cant open $ARGV[0]: $!";
#---------------------read bond parameters from gaff.dat------------------------
$_ = <GAFF>;
until(length($_)<3) {
  $_ = <GAFF>;
}
$_ = <GAFF>;
$_ = <GAFF>;
until(length($_)<3) {
  ($t1, $t2, $frc, $bl) = /(\S*)\s*-(\S*)\s*(\d*\.\d*)\s*(\d*\.\d*)/;
  $bondlength{$t1}{$t2} = $bl;
  $force{$t1}{$t2} = $frc;
  $_ = <GAFF>;
}
#---------------------read selected atoms from mol2 file------------------------
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
    $x[$id] = $pos_x;
    $y[$id] = $pos_y;
    $z[$id] = $pos_z;
    $atom_type[$id] = $type;
    $atom_name[$id] = $name;
  }
}
#---read bonds between selected atoms from mol2 file calculate bond length and print to output----------
while(<MOL2>) {
  $string = $_;
  if(/SUBSTRUCTURE/) {
    last;
  } else {
    ($bond_id, $id1, $id2, $type) =
    /\s*(\d*)\s*(\S*)\s*(\S*)\s*(\S*)\s*/;
    $bond_id1[$bond_id] = $id1;
    $bond_id2[$bond_id] = $id2;
    $bond_type[$bond_id] = $type;
    #check whether bond is between two selected atoms (then j=2 -> see below)
    $j=0;
    for($i=0; $i<=$#list; $i++) {
      if($list[$i]==$id1 or $list[$i]==$id2) {
        $j++;
      }
    }
    if($j==2) {
      #calculate bond length
      $dx = $x[$id1]-$x[$id2];
      $dy = $y[$id1]-$y[$id2];
      $dz = $z[$id1]-$z[$id2];
      $bond_length[$bond_id] = sqrt($dx**2+$dy**2+$dz**2);
      #get force constant from gaff.dat
      if($force{$atom_type[$id1]}{$atom_type[$id2]} == 0) {
        if($force{$atom_type[$id2]}{$atom_type[$id1]} == 0) {
          $frc = -1.0;
        } else {
          $frc = $force{$atom_type[$id2]}{$atom_type[$id1]}*$kcal*2;
        }
      } else {
        $frc = $force{$atom_type[$id1]}{$atom_type[$id2]}*$kcal*2;
      }
      if($outf) {
        print OUTFILE "harm";
        printf OUTFILE " %6s", $id1;
        printf OUTFILE " %6s", $id2;
        printf OUTFILE " %12.5f", $frc;
        printf OUTFILE " %12.5f", $bond_length[$bond_id];
        printf OUTFILE "       ; %-4s-%-4s", $atom_name[$id1],$atom_name[$id2];
        printf OUTFILE " %-4s-%-4s", $atom_type[$id1],$atom_type[$id2];
        printf OUTFILE " %-4s", $type;
        print OUTFILE "\n";
      } else {
        print "harm";
        printf " %6s", $id1;
        printf " %6s", $id2;
        printf " %12.5f", $frc;
        printf " %12.5f", $bond_length[$bond_id];
        printf "       ; %-4s-%-4s", $atom_name[$id1],$atom_name[$id2];
        printf " %-4s-%-4s", $atom_type[$id1],$atom_type[$id2];
        printf " %-4s", $type;
        print "\n";
      }
    }
  }
}
close GAFF or die "error closing $ARGV[0]: $!";
close MOL2 or die "error closing $ARGV[1]: $!";