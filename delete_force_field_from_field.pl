#!/usr/bin/perl

# deletes all bond, angle, dihedral, vdw, etc. information from a FIELD file
# so that only the electrostatic interactions are left

use dlpoly_utility;
use hanno_utility;

if($#ARGV<1) {
  print "input format:\n";
  print "1. name of input FIELD file\n";
  print "2. name of output FIELD file\n";
  exit 1;
}

$infilename=$ARGV[0];
$outfilename=$ARGV[1];

if($infilename eq $outfilename) {
  print "error: target FIELD and output FIELD file must not be the same";
  exit 1;
}

&read_field_file($infilename,0);

for($t=0;$t<$field_nummols[0];$t++) {
  $mol_numbonds[0][$t] = 0;
  @{$mol_bonddata[0][$t]} = ();
  @{$mol_bondatoms[0][$t]} = ();
  $mol_numconstraints[0][$t] = 0;
  @{$mol_constraintdata[0][$t]} = ();
  $mol_numangles[0][$t] = 0;
  @{$mol_angledata[0][$t]} = ();
  $mol_numdihedrals[0][$t] = 0;
  @{$mol_dihedraldata[0][$t]} = ();
  $mol_numinversions[0][$t] = 0;
  @{$mol_inversiondata[0][$t]} = ();
}

$field_numtbp[0]=0;
$field_numfbp[0]=0;
$field_numextern[0]=0;
$field_numvdw[0]=1;
@{$field_vdwdata[0][0]} = ($field_atomtypes[0][0],$field_atomtypes[0][0],'lj',0,1);

write_field_file($outfilename,0);
