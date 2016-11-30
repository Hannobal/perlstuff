#!/usr/bin/perl
use Math::Trig;
use Switch;

use dlpoly_utility;
use hanno_utility;
#use Storable qw(dclone);


if($#ARGV<6) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file with surface\n";
  print "2. name of corresponding FIELD-file\n";
  print "3. name of mol2-file with molecule\n";
  print "4. name of output CONFIG-file\n";
  print "5. name of output FIELD-file\n";
  print "6. name of molecule to replace\n";
  print "7. molecule id (starting from 0)\n";
  exit 1;
}

$filenameincfg  = $ARGV[0];
$filenameinfld  = $ARGV[1];
$filenamemol2   = $ARGV[2];
$filenameoutcfg = $ARGV[3];
$filenameoutfld = $ARGV[4];
$repmolname     = $ARGV[5];
$repmolid       = $ARGV[6];

#load molecule
exit 1 if(read_mol2_file($filenamemol2,0)!=0);
$mol2name=$mol2_name[0];
$mol2name =~ s/\.\S*$//; # remove file extension
$mol2name = substr($mol2name,rindex($mol2name,'/')+1); # remove path (everything before last "/")

$mol2numox = 0;
for($a=0;$a<$mol2_numatoms[0];$a++) {
  if(uc($mol2_atomdata[0][$a][4]) eq "O" or uc($mol2_atomdata[0][$a][4]) eq "OH") {
    $mol2oxpos[$mol2numox][0]=$mol2_atomdata[0][$a][0];
    $mol2oxpos[$mol2numox][1]=$mol2_atomdata[0][$a][1];
    $mol2oxpos[$mol2numox][2]=$mol2_atomdata[0][$a][2];
    $mol2numox++;
  }
}

@mol2oxpos = sort { $a->[2] <=> $b->[2] } @mol2oxpos;
$c = sprintf("%.0f",$mol2_charge[0]);
if(not defined($bindingmode)) {
  if($c==0) {
    $bindingmode=0;
  }elsif($c==-1) {
    $bindingmode=1;
  }elsif($c==-2) {
    $bindingmode=2;
  }
}

# read FIELD file
exit 1 if(read_field_file($filenameinfld,0)!=0);
$foundmol    = -1;
# $foundsolv   = -1;
$foundrepmol = -1;
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /^$mol2name$/i) {
    $foundmol=$t;
#   } elsif($mol_name[0][$t] =~ /$solvname/i) {
#     $foundsolv=$t;
  } elsif($mol_name[0][$t] =~ /^$repmolname/i) {
    $foundrepmol=$t;
  }
}
die "error: could not find entry for oxygen in FIELD file $filenameinfld!\n" if($foundox==-1);
die "error: could not find entry for hydroxide ions in FIELD file $filenameinfld!\n" if($foundhydrox==-1);
die "error: could not find entry for molecule $mol2name in FIELD file $filenameinfld!\n" if($foundmol==-1);
die "error: could not find entry for molecule \"$repmolname\" in FIELD file $filenameinfld!\n" if($foundrepmol==-1);
# die "error: could not find entry for solvent molecule in FIELD file $filenameinfld!\n"
# if($foundsolv==-1 and $lsolv);

#read CONFIG file
exit 1 if(read_config_file($filenameincfg,0,0)!=0);

@pos=(0,0,99999);
for($a=0;$a<$mol_numatoms[0][$foundrepmol];$a++) {
  if($cdata[0][$foundrepmol][$repmolid][$a][2]<$pos[2]) {
    $pos[0]=$cdata[0][$foundrepmol][$repmolid][$a][0];
    $pos[1]=$cdata[0][$foundrepmol][$repmolid][$a][1];
    $pos[2]=$cdata[0][$foundrepmol][$repmolid][$a][2];
  }
}

splice(@{$cdata[0][$foundrepmol]},$repmolid,1);
$mol_numents[0][$foundrepmol]--;
@{$cdata[0][$foundmol][$mol_numents[0][$foundmol]]} = mol2_to_cdata(0);
move_mol(\@{$cdata[0][$foundmol][$mol_numents[0][$foundmol]]},\@pos);
$mol_numents[0][$foundmol]++;

write_config_file($filenameoutcfg,0,$config_title[0]);
write_field_file($filenameoutfld,0);
exit 0;
