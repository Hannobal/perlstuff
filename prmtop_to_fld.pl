#!/usr/bin/perl

# for prmtop file format confer http://ambermd.org/formats.html

use dlpoly_utility;
use hanno_utility;
use Math::Trig;
use strict;

my(@data,$section,$fhin,$name,$i,$j,$k);
our(@pointer,@bond_param,@bond_length,@bond_index,@angle_param,@angle_value,@angle_index,@dih_param,
@dih_period,@dih_phase,@dih_index,@attype,@charge,@mass,@attype_index,@vdw_index,@vdw,@scee_scale,
@scbn_scale,@vdw_acoeff,@vdw_bcoeff);

our $engconv=4.184; # kcal to kJ

if($#ARGV<0) {
  print "please specify common name for .prmtop and .fld file\n";
  exit 1;
}

$name=$ARGV[0];

open($fhin,"<","$name.prmtop") or die "**** error: can't open input file:\n$!";

# read the prmtop file
while(<$fhin>) {
  if(/\%FLAG (\S+)/) {
    assign_data($section,\@data);
    $section=$1;
    @data=();
    next;
  }
  next if(/\%FORMAT/);
  $_=~s/^\s+//;$_=~s/\s+$//;
  push(@data,split(" ",$_));
}
close($fhin);

# now generate the FIELD data
$field_units[0]="kJ/mol";
$field_title[0]=$name;
$field_nummols[0]=1;
$mol_numents[0][0]=1;
$mol_numatoms[0][0]=$pointer[0];
$mol_numbonds[0][0]=$pointer[2]+$pointer[3];
$mol_numangles[0][0]=$pointer[4]+$pointer[5];
$mol_numdihedrals[0][0]=$pointer[6]+$pointer[7];
$field_numvdw[0]=0;

@{$mol_atomdata[0][0]}=();
for($i=0;$i<@attype_index;$i++) {
  $mol_atomdata[0][0][$i][0] = uc($attype[$i]);
  $mol_atomdata[0][0][$i][1] = $mass[$i];
  $mol_atomdata[0][0][$i][2] = $charge[$i]/18.2223;
  $mol_atomdata[0][0][$i][3] = 1;
  $mol_atomdata[0][0][$i][4] = 0;
  $mol_atomdata[0][0][$i][5] = 0;
}

$j=0;
@{$mol_bonddata[0][0]}=();
for($i=0;$i<@bond_index;$i+=3) {
  $mol_bonddata[0][0][$j][0]="harm";
  $mol_bonddata[0][0][$j][1]=$bond_index[$i]/3;
  $mol_bonddata[0][0][$j][2]=$bond_index[$i+1]/3;
  $mol_bonddata[0][0][$j][3]=$bond_param[$bond_index[$i+2]-1]*$engconv*2;
  $mol_bonddata[0][0][$j][4]=$bond_length[$bond_index[$i+2]-1];
  $j++;
}
if(@{$mol_bonddata[0][0]}) {
  @{$mol_bonddata[0][0]} = sort {$a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]} @{$mol_bonddata[0][0]};
}

$j=0;
@{$mol_angledata[0][0]}=();
for($i=0;$i<@angle_index;$i+=4) {
  $mol_angledata[0][0][$j][0]="harm";
  $mol_angledata[0][0][$j][1]=$angle_index[$i]/3;
  $mol_angledata[0][0][$j][2]=$angle_index[$i+1]/3;
  $mol_angledata[0][0][$j][3]=$angle_index[$i+2]/3;
  $mol_angledata[0][0][$j][4]=$angle_param[$angle_index[$i+3]-1]*$engconv*2;
  $mol_angledata[0][0][$j][5]=$angle_value[$angle_index[$i+3]-1]*180/pi;
  if($mol_angledata[0][0][$j][1]>$mol_angledata[0][0][$j][3]) {
    my $tmp=$mol_angledata[0][0][$j][1];
    $mol_angledata[0][0][$j][1]=$mol_angledata[0][0][$j][3];
    $mol_angledata[0][0][$j][3]=$tmp;
  }
  $j++;
}
if(@{$mol_angledata[0][0]}) {
  @{$mol_angledata[0][0]} = sort {$a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3]} @{$mol_angledata[0][0]};
}

$j=0;
@{$mol_dihedraldata[0][0]}=();
for($i=0;$i<@dih_index;$i+=5) {
  $mol_dihedraldata[0][0][$j][0]="cos";
  $mol_dihedraldata[0][0][$j][1]=$dih_index[$i]/3;
  $mol_dihedraldata[0][0][$j][2]=$dih_index[$i+1]/3;
  $mol_dihedraldata[0][0][$j][3]=$dih_index[$i+2]/3;
  $mol_dihedraldata[0][0][$j][4]=$dih_index[$i+3]/3;
  $mol_dihedraldata[0][0][$j][5]=$dih_param[$dih_index[$i+4]-1]*$engconv;
  $mol_dihedraldata[0][0][$j][6]=$dih_phase[$dih_index[$i+4]-1]*180/pi;
  $mol_dihedraldata[0][0][$j][7]=$dih_period[$dih_index[$i+4]-1];
  $mol_dihedraldata[0][0][$j][8]=1.0/$scee_scale[$dih_index[$i+4]-1];
  $mol_dihedraldata[0][0][$j][9]=1.0/$scbn_scale[$dih_index[$i+4]-1];
  print join("|",@{$mol_dihedraldata[0][0][$j]}),"\n";
  if($mol_dihedraldata[0][0][$j][1]>$mol_dihedraldata[0][0][$j][4]) {
    my $tmp=$mol_dihedraldata[0][0][$j][1];
    $mol_dihedraldata[0][0][$j][1]=$mol_dihedraldata[0][0][$j][4];
    $mol_dihedraldata[0][0][$j][4]=$tmp;
    $tmp=$mol_dihedraldata[0][0][$j][2];
    $mol_dihedraldata[0][0][$j][2]=$mol_dihedraldata[0][0][$j][3];
    $mol_dihedraldata[0][0][$j][3]=$tmp;
  }
  $j++;
}
if(@{$mol_dihedraldata[0][0]}) {
  @{$mol_dihedraldata[0][0]} = sort {$a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3] || $a->[3] <=> $b->[4]} @{$mol_dihedraldata[0][0]};
}

# find the indices for the atom types
my @asdf=();
for($i=0;$i<@attype;$i++) {
  next if(contains(@{$field_atomtypes[0]},$attype[$i]));
  push(@{$field_atomtypes[0]},$attype[$i]);
  push(@asdf,[$attype[$i],$attype_index[$i]]);
}
@asdf= sort {$a->[0] cmp $b->[0]} @asdf;

$k=0;
for($i=0;$i<@asdf;$i++) {
#   print "$asdf[$i][0] $asdf[$i][1]\n";
  for($j=$i;$j<@asdf;$j++) {
#     print @asdf*($asdf[$i][1]-1)+$asdf[$j][1]-1," ";
    my $n=$vdw_index[@asdf*($asdf[$i][1]-1)+$asdf[$j][1]-1]-1;
    $field_vdwdata[0][$k][0]=uc($asdf[$i][0]);
    $field_vdwdata[0][$k][1]=uc($asdf[$j][0]);
    $field_vdwdata[0][$k][2]="lj";
    @{$field_vdwdata[0][$k]}[3..4]=lj_ab_to_sigeps($vdw_acoeff[$n],$vdw_bcoeff[$n]);
    $field_vdwdata[0][$k][3]*=$engconv;
    $k++;
  }
}
$field_numvdw[0]=$k;
write_field_file("$name.fld",0);

sub assign_data {
  my @data=@{$_[1]};
  # omitted flags:
  # ATOM_NAME 
  # ATOMIC_NUMBER 
  # NUMBER_EXCLUDED_ATOMS
  # RESIDUE_LABEL
  # RESIDUE_POINTER
  # all others
  my($i);
  if($_[0] eq "TITLE") {
    $mol_name[0][0]=join(" ",@data);
  } elsif($_[0] eq "POINTERS") {
    @pointer=@data;
  } elsif($_[0] eq "ATOM_TYPE_INDEX") {
    @attype_index=@data;
  } elsif($_[0] eq "CHARGE") {
    @charge=@data;
  } elsif($_[0] eq "MASS") {
    @mass=@data;
  } elsif($_[0] eq "NONBONDED_PARM_INDEX") {
    @vdw_index=@data;
  } elsif($_[0] eq "BOND_FORCE_CONSTANT") {
    @bond_param=@data;
  } elsif($_[0] eq "BOND_EQUIL_VALUE") {
    @bond_length=@data;
  } elsif($_[0] eq "ANGLE_FORCE_CONSTANT") {
    @angle_param=@data;
  } elsif($_[0] eq "ANGLE_EQUIL_VALUE") {
    @angle_value=@data;
  } elsif($_[0] eq "DIHEDRAL_FORCE_CONSTANT") {
    @dih_param=@data;
  } elsif($_[0] eq "DIHEDRAL_PERIODICITY") {
    @dih_period=@data;
  } elsif($_[0] eq "DIHEDRAL_PHASE") {
    @dih_phase=@data;
  } elsif($_[0] eq "SCEE_SCALE_FACTOR") {
    @scee_scale=@data;
  } elsif($_[0] eq "SCNB_SCALE_FACTOR") {
    @scbn_scale=@data;
  } elsif($_[0] eq "LENNARD_JONES_ACOEF") {
    @vdw_acoeff=@data;
  } elsif($_[0] eq "LENNARD_JONES_BCOEF") {
    @vdw_bcoeff=@data;
  } elsif($_[0] eq "SCEE_SCALE_FACTOR") {
    @scee_scale=@data;
  } elsif($_[0] eq "AMBER_ATOM_TYPE") {
    @attype=@data;
  } elsif($_[0] =~/^BONDS_/) {
    push(@bond_index,@data);
  } elsif($_[0] =~/^ANGLES_/) {
    push(@angle_index,@data);
  } elsif($_[0] =~/^DIHEDRALS_/) {
    push(@dih_index,@data);
  }
}

sub error {
  if(@_>0) {
    print "**** error: $_[1]\n";
  }
  exit $_[0];
}