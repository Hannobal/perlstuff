#!/usr/bin/perl

use hanno_utility;
use dlpoly_utility;
use strict vars;

if($#ARGV<3) {
  print " 1. name of input FIELD/mol2-file\n";
  print " 2. name of output FIELD/mol2-file\n";
  print " 3. molecule name\n";
  print " 4. new charge\n";
  print "[list of atom ids (starting from 1) or names for scaling]\n";
  exit 1;
}
my $infile    = $ARGV[0];
my $outfile   = $ARGV[1];
my $molname   = $ARGV[2];
my $newcharge = $ARGV[3];

my($i,$t,$a,$tm,$diff,$rescharge);
our(@atomids,@atomnames);

if(not check_real($newcharge)) {
  print "**** error: expected real number as new charge but got $newcharge\n";
  exit 1;
}

for($i=4;$i<@ARGV;$i++) {
  if(check_integer($ARGV[$i])) {
    push(@atomids,$ARGV[$i]-1);
  } else {
    push(@atomnames,$ARGV[$i]);
  }
}

if($infile=~/.mol2/) {
  exit 1 if(read_mol2_file($infile,0) != 0);
  $rescharge=0;
  for($a=0;$a<$mol2_numatoms[0];$a++) {
    if(@atomids or @atomnames) {
      next if(not checkfit2($a));
    }
    $rescharge+=abs($mol2_atomdata[0][$a][8]);
  }
  if($rescharge==0) {
    print "**** error: no charges found to scale!\n";
    exit 1;
  }
  $diff=$newcharge-$mol2_charge[0];
  for($a=0;$a<$mol2_numatoms[0];$a++) {
    if(@atomids or @atomnames) {
      next if(not checkfit2($a));
    }
    print "scaling atom ",$a+1," named $mol2_atomdata[0][$a][3]\n";
    $mol2_atomdata[0][$a][8] += $diff*abs($mol2_atomdata[0][$a][8]/$rescharge);
  }
  write_mol2_file($outfile,0);
  
} else {

  exit 1 if(read_field_file($infile,0) != 0);
  $tm=-1;
  $rescharge=0;
  for($t=0;$t<$field_nummols[0];$t++) {
    next if($mol_name[0][$t]!~/$molname/i);
    $tm=$t;
    for($a=0;$a<$mol_numatoms[0][$t];$a++) {
      if(@atomids or @atomnames) {
	next if(not checkfit($t,$a));
      }
      print "$mol_atomdata[0][$t][$a][0] $mol_atomdata[0][$t][$a][2]\n";
      $rescharge+=abs($mol_atomdata[0][$t][$a][2]);
    }
    last;
  }
  if($tm<0) {
    print "**** error: molecule named '$ARGV[$i]' was not found in FIELD file!\n";
    exit 1;
  }
  if($rescharge==0) {
    print "**** error: no charges found to scale!\n";
    exit 1;
  }
  $diff=$newcharge-$mol_charge[0][$tm];
  print "$rescharge\n";
  for($a=0;$a<$mol_numatoms[0][$tm];$a++) {
    if(@atomids or @atomnames) {
      next if(not checkfit($tm,$a));
    }
    print "scaling atom ",$a+1," named $mol_atomdata[0][$tm][$a][0]\n";
    $mol_atomdata[0][$tm][$a][2] += $diff*abs($mol_atomdata[0][$tm][$a][2]/$rescharge);
  }

  write_field_file($outfile,0);
}

sub checkfit {
  my $t=$_[0];
  my $a=$_[1];
  my($i,$name);
  foreach $i (@atomids) {
    return 1 if($i == $a);
  }
  foreach $name (@atomnames) {
    return 1 if(uc($mol_atomdata[0][$t][$a][0]) eq uc($name));
  }
  return 0;
}
sub checkfit2 {
  my $a=$_[0];
  my($i,$name);
  foreach $i (@atomids) {
    return 1 if($i == $a);
  }
  foreach $name (@atomnames) {
    return 1 if(uc($mol2_atomdata[0][$a][4]) eq uc($name));
  }
  return 0;
}