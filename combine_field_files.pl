#!/usr/bin/perl

use Storable qw(dclone);
use dlpoly_utility;
use hanno_utility;

if($#ARGV<4) {
  print "input format:\n";
  print "1. name of target FIELD file with surface\n";
  print "2. name of source FIELD file with molecule\n";
  print "3. name of molecule in source FIELD file\n";
  print "4. name of output FIELD file\n";
  print "5. name of GAFF parameter file\n";
  print "[6. additional frcmod file]\n";
  print "[7. (y/n) add additional frozen version of molecule?]\n";
  exit 1;
}

if($#ARGV>4) {
  $frcmodfilename=$ARGV[5];
}
if($#ARGV>5) {
  if($ARGV[6] =~ /^[yn]/i) {
    if($ARGV[6] =~ /^[y]/i) {
      $addfrozenversion = 1;
    } else {
      $addfrozenversion = 0;
    }
  } else {
    print "error: could not interpret input \"$ARGV[5]\"!";
    exit 1;
  }
}

$targetfilename=$ARGV[0];
$sourcefilename=$ARGV[1];
$molname=$ARGV[2];
$outputfilename=$ARGV[3];
$gafffilename=$ARGV[4];

$number1 = 1;

if($targetfilename eq $outputfilename) {
  print "error: target FIELD and output FIELD file must not be the same";
  exit 1;
}
if($sourcefilename eq $outputfilename) {
  print "error: source FIELD and output FIELD file must not be the same";
  exit 1;
}


# first check whether molecule is already added
exit 1 if(read_field_file($targetfilename,0) != 0);
$foundmolname=-1;
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] eq $molname) {
    $foundmolname = $t;
    last;
  }
}
if($foundmolname>=0) {
  print "FIELD file $targetfilename already contains molecule $molname\n";
  exit 0;
}


exit 1 if(read_field_file($sourcefilename,1) != 0);
$foundmolname=-1;
for($t=0;$t<$field_nummols[1];$t++) {
  if($mol_name[1][$t] eq $molname) {
    $foundmolname = $t;
    last;
  }
}
if($foundmolname < 0) {
  print "**** error: molecule name '$molname' was not found in FIELD file $sourcefilename!";
  exit 1;
}

exit 1 if(read_gaff_file($gafffilename) != 0);
if(defined($frcmodfilename)) {
  exit 1 if(read_frcmod_file($frcmodfilename) != 0);
}
unshift(@{$mol_atomdata[0]},\@{dclone \@{$mol_atomdata[1][$foundmolname]}});
unshift(@{$mol_bonddata[0]},\@{dclone \@{$mol_bonddata[1][$foundmolname]}});
unshift(@{$mol_angledata[0]},\@{dclone \@{$mol_angledata[1][$foundmolname]}});
unshift(@{$mol_dihedraldata[0]},\@{dclone \@{$mol_dihedraldata[1][$foundmolname]}});
unshift(@{$mol_inversiondata[0]},\@{dclone \@{$mol_inversiondata[1][$foundmolname]}});
unshift(@{$mol_constraintdata[0]},\@{dclone \@{$mol_constraintdata[1][$foundmolname]}});
unshift(@{$mol_tetherdata[0]},\@{dclone \@{$mol_tetherdata[1][$foundmolname]}});
unshift(@{$mol_rigiddata[0]},\@{dclone \@{$mol_rigiddata[1][$foundmolname]}});
unshift(@{$mol_bondatoms[0]},\@{dclone \@{$mol_bondatoms[1][$foundmolname]}});
unshift(@{$mol_numents[0]},$mol_numents[1][$foundmolname]);
unshift(@{$mol_numatoms[0]},$mol_numatoms[1][$foundmolname]);
unshift(@{$mol_numbonds[0]},$mol_numbonds[1][$foundmolname]);
unshift(@{$mol_numangles[0]},$mol_numangles[1][$foundmolname]);
unshift(@{$mol_numdihedrals[0]},$mol_numdihedrals[1][$foundmolname]);
unshift(@{$mol_numinversions[0]},$mol_numinversions[1][$foundmolname]);
unshift(@{$mol_numconstraints[0]},$mol_numconstraints[1][$foundmolname]);
unshift(@{$mol_numtether[0]},$mol_numtether[1][$foundmolname]);
unshift(@{$mol_name[0]},$mol_name[1][$foundmolname]);
unshift(@{$mol_numrigid[0]},$mol_numrigid[1][$foundmolname]);
$field_nummols[0]++;
$mol_numents[0][0]=0;

forname0: for($i=0;$i<$mol_numatoms[1][$foundmolname];$i++) {
  $name1 = $mol_atomdata[1][$foundmolname][$i][0];
  foreach $name2 (values @{$field_atomtypes[0]}) {
    next forname0 if($name1 eq $name2);
  }
  unshift(@{$field_atomtypes[0]},$name1);
  $field_numatomtypes[0]++;
  foreach $name2 (values @{$field_atomtypes[0]}) {
    $exists = 0;
    for($j=0;$j<$field_numvdw[1];$j++) {
      if(($field_vdwdata[1][$j][0] eq $name1) and ($field_vdwdata[1][$j][1] eq $name2)
      or ($field_vdwdata[1][$j][0] eq $name2) and ($field_vdwdata[1][$j][1] eq $name1)) {
	push(@{$field_vdwdata[0]},\@{dclone \@{$field_vdwdata[1][$j]}});
	$exists = 1;
      }
    }
    if(not $exists) {
      print "calculating vdw parameters for $name1 $name2\n";
      @{$field_vdwdata[0][$field_numvdw[0]]}[0..2] = (uc($name1),uc($name2),"lj");
      @{$field_vdwdata[0][$field_numvdw[0]]}[3..4] = calc_vdw_params(uc($name1),uc($name2));
    }
    $field_numvdw[0]++;
  }
}

# sort indices lexically ascending
for($i=0;$i<$field_numvdw[0];$i++) {
  if($field_vdwdata[0][$i][0] gt $field_vdwdata[0][$i][1]) {
    $name1                   = $field_vdwdata[0][$i][1];
    $field_vdwdata[0][$i][1] = $field_vdwdata[0][$i][0];
    $field_vdwdata[0][$i][0] = $name1;
  }
}

@{$field_vdwdata[0]} = sort {$b->[2] cmp $a->[2] || $a->[0] cmp $b->[0] ||
                       $a->[1] cmp $b->[1]} @{$field_vdwdata[0]};

write_field_file($outputfilename,0);

exit 0;
