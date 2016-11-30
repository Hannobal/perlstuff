#!/usr/bin/perl

# moves an atom in the atom list of a FIELD or mol2 file
# and adjusts bond, angle and dihedral parameters

use hanno_utility;
use dlpoly_utility;

if($#ARGV<3) {
  print "input format:\n",
  " 1. name of input FIELD/mol2 file\n",
  " 2. name of output FIELD/mol2 file\n",
  " 3. name of molecule\n",
  " 4. index of atom (starting from 1)\n",
  "[5. new index for atom (starting from 1)]\n";
  exit 1;
}

$infile=$ARGV[0];
$outfile=$ARGV[1];
$name=$ARGV[2];
$index=$ARGV[3];
if(not check_integer($index)) {
  print "**** error: atom index must be integer number!\n"; exit 1;
}

$new=0;
if($#ARGV>3) {
  $new=$ARGV[4];
  if(not check_integer($new)) {
    print "**** error: new atom index must be integer number!\n"; exit 1;
  }
}

$index--;
$new--;
if($infile=~/.mol2/) {
  exit 1 if (read_mol2_file($infile,0)!=0);
  if(not $mol2_name[0] eq $name) {
    print "**** error: molecule name does not fit!\n"; exit 1;
  }
  if($new<0 or $new>=$mol2_numatoms[0]) {
    $new = $mol2_numatoms[0]-1;
  }
  if($new>$index) {
    $shift=-1;
  } else {
    $shift=1;
  }
  ($lower,$higher) = sort {$a <=> $b} ($index,$new);
  splice(@{$mol2_atomdata[0]},$new,0,splice(@{$mol2_atomdata[0]},$index,1));
  for($b=0;$b<$mol2_numbonds[0];$b++) {
    if($mol2_bonddata[0][$b][0]==$index) {
      $mol2_bonddata[0][$b][0]=$new;
    } elsif(($mol2_bonddata[0][$b][0]>$lower and $mol2_bonddata[0][$b][0]<$higher) or $mol2_bonddata[0][$b][0]==$new) {
      $mol2_bonddata[0][$b][0] += $shift;
    }
    if($mol2_bonddata[0][$b][1]==$index) {
      $mol2_bonddata[0][$b][1]=$new;
    } elsif(($mol2_bonddata[0][$b][1]>$lower and $mol2_bonddata[0][$b][1]<$higher) or $mol2_bonddata[0][$b][1]==$new) {
      $mol2_bonddata[0][$b][1] += $shift;
    }
#     @{$mol2_bonddata[0]} = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$mol2_bonddata[0]} !!! wrong !!!
  }
  write_mol2_file($outfile,0);
  
  
} else {
  exit 1 if (read_field_file($infile,0)!=0);
  read_field_file($infile,0);
  for($t=0;$t<$field_nummols[0];$t++) {
    last if ($mol_name[0][$t] eq $name);
  }
  if($t==$field_nummols[0]){
    print "**** error: no molecule named $name found!\n"; exit 1;
  }
  if($new<0 or $new>=$mol_numatoms[0][$t]) {
    $new = $mol_numatoms[0][$t]-1;
  }
  if($new>$index) {
    $shift=-1;
  } else {
    $shift=1;
  }
  print "atom name: $mol_atomdata[0][$t][$index][0]\n";
  ($lower,$higher) = sort {$a <=> $b} ($index,$new);
  splice(@{$mol_atomdata[0][$t]},$new,0,splice(@{$mol_atomdata[0][$t]},$index,1));
  for($b=0;$b<$mol_numbonds[0][$t];$b++) {
    for($a=0;$a<=2;$a++) {
      if($mol_bonddata[0][$t][$b][$a]==$index) {
	$mol_bonddata[0][$t][$b][$a] = $new;
      } elsif(($mol_bonddata[0][$t][$b][$a]>$lower and $mol_bonddata[0][$t][$b][$a]<$higher) or $mol_bonddata[0][$t][$b][$a]==$new) {
	$mol_bonddata[0][$t][$b][$a] += $shift;
      }
    }
#     @{$mol_bonddata[0][$t]} = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$mol_bonddata[0][$t]} !!! wrong !!!
  }
  for($b=0;$b<$mol_numangles[0][$t];$b++) {
    for($a=0;$a<=3;$a++) {
      if($mol_angledata[0][$t][$b][$a]==$index) {
	$mol_angledata[0][$t][$b][$a] = $new;
      } elsif(($mol_angledata[0][$t][$b][$a]>$lower and $mol_angledata[0][$t][$b][$a]<$higher) or $mol_angledata[0][$t][$b][$a]==$new) {
	$mol_angledata[0][$t][$b][$a] += $shift;
      }
    }
#     @{$mol_angledata[0][$t]} = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3] } @{$mol_angledata[0][$t]} !!! wrong !!!
  }
  for($b=0;$b<$mol_numdihedrals[0][$t];$b++) {
    for($a=0;$a<=4;$a++) {
      if($mol_dihedraldata[0][$t][$b][$a]==$index) {
	$mol_dihedraldata[0][$t][$b][$a] = $new;
      } elsif(($mol_dihedraldata[0][$t][$b][$a]>$lower and $mol_dihedraldata[0][$t][$b][$a]<$higher) or $mol_dihedraldata[0][$t][$b][$a]==$new) {
	$mol_dihedraldata[0][$t][$b][$a] += $shift;
      }
    }
#     @{$mol_dihedraldata[0][$t]} = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3]  || $a->[4] <=> $b->[4] } @{$mol_dihedraldata[0][$t]} !!! wrong !!!
  }
  write_field_file($outfile,0);
}