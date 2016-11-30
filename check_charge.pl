#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;

if($#ARGV<0) {
  print "input format:\n",
  " 1. filename (*.mol2/FIELD/*.fld)>\n",
  "[2. list of atom numbers (mol2) only]\n";
  exit;
}

$filename = $ARGV[0];
@atomids = parse_intlist(-1,@ARGV[1..$#ARGV]);

if($filename =~ /mol2/) {
  exit 1 if(read_mol2_file($ARGV[0],0) != 0);
  $molcharge = 0;
  for($a=0;$a<@{$mol2_atomdata[0]};$a++) {
    next if(@atomids and not contains(@atomids,$a));
    $molcharge += $mol2_atomdata[0][$a][8];
  }
  print "$molcharge\n";
} elsif($filename =~ /.fld/ or $filename =~ /FIELD/) {
  exit 1 if(read_field_file($ARGV[0],0) != 0);
  $totcharge = 0;
  printf "%15s %14s %14s\n", "molecule","charge","# molecules";
  for($t=0;$t<$field_nummols[0];$t++) {
    $molcharge = 0;
    for($a=0;$a<$mol_numatoms[0][$t];$a++) {
      $molcharge += $mol_atomdata[0][$t][$a][2];
    }
    printf "%15s %14.10f %14u\n", $mol_name[0][$t],$molcharge,$mol_numents[0][$t];
    $totcharge += $mol_numents[0][$t]*$molcharge;
  }
  print "---------------------------------------------\n";
  printf "%15s %29.10f\n", "total charge",$totcharge;
} else {
  print "**** error: file format not recognized!\n";
}
