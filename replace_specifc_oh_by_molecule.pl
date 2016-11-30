#!/usr/bin/perl
use Math::Trig;
use Switch;

use dlpoly_utility;
use hanno_utility;
use aloxsam_utility;
#use Storable qw(dclone);


if($#ARGV<5) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file with surface\n";
  print "2. name of corresponding FIELD-file\n";
  print "3. name of mol2-file with molecule\n";
  print "4. name of output CONFIG-file\n";
  print "5. name of output FIELD-file\n";
  print "6. list of molecule IDs for replacement\n";
  print "   you can add a rotation angle [deg.] separated with comma to each ID\n";
  print "   example: 10 12,120 15,60\n";
  exit 1;
}

$filenameincfg  = $ARGV[0];
$filenameinfld  = $ARGV[1];
$filenamemol2   = $ARGV[2];
$filenameoutcfg = $ARGV[3];
$filenameoutfld = $ARGV[4];
@ids = @ARGV[5..$#ARGV];
$surfdist=0;

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
$toa   = -1;
$toh   = -1;
$tal   = -1;
$tmol  = -1;
$tsolv = -1;
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /aluminum\s+free/i) {
    $tal=$t;
  } elsif($mol_name[0][$t] =~ /oxygen\s+free/i) {
    $toa=$t;
  } elsif($mol_name[0][$t] =~ /hydroxide/i) {
    $toh=$t;
  } elsif($mol_name[0][$t] =~ /^$mol2name$/i) {
    $tmol=$t;
  } elsif($mol_name[0][$t] =~ /$solvname/i) {
    $tsolv=$t;
  }
}
die "error: could not find entry for oxygen in FIELD file $filenameinfld!\n" if($toa==-1);
die "error: could not find entry for hydroxide ions in FIELD file $filenameinfld!\n" if($toh==-1);
die "error: could not find entry for molecule $mol2name in FIELD file $filenameinfld!\n" if($tmol==-1);
die "error: could not find entry for solvent molecule in FIELD file $filenameinfld!\n"
if($tsolv==-1 and $lsolv);

#read CONFIG file
exit 1 if(read_config_file($filenameincfg,0,0)!=0);
switch($periodic_key[0]) {
  case [1,2,6] {
  } case 0 {
    print "error: periodic conditions needed for calculation\n";
    exit 1;
  } else {
    print "sorry, this script is only implemented for orthogonal cells\n";
    exit 1;
  }
}

$zsurf = calc_zsurf(0,0,1);

@remoharr=();
$mmol = $mol_numents[0][$tmol];
foreach $id (@ids) {
  $id=~/(\d+),?(\S*)/;
  $moh   = $1;
  $angle = $2;
  if(not check_integer($moh)) {
    print "**** error: ID $moh is not an integer!\n"; exit 1;
  } elsif($m>=$mol_numents[0][$toh]) {
    print "**** error: ID $moh does not exist!\n"; exit 1;
  }
  if(length($angle)>0 and not check_real($angle)) {
    print "**** error: angle $angle for ID $moh is not a real number!\n"; exit 1;
  }
  @pos = @{$cdata[0][$toh][$moh][0]}[0..2];
  @{$cdata[0][$tmol][$mmol]} = mol2_to_cdata(0);
  if(length($angle)>0) {
    @rotmatrix = &gen_rot_matrix([0,0,1], $angle/180.0*pi);
    rotate_molecule(\@{$cdata[0][$tmol][$mmol]}, \@rotmatrix);
  }
  move_mol(\@{$cdata[0][$tmol][$mmol]},\@pos);
  $mol_numents[0][$tmol]++;
  $mmol++;
  push(@remoharr,$moh);
}

@remoharr = sort {$a <=> $b} @remoharr;
for($i=$#remoharr;$i>=0;$i--) {
  remove_mol_entity(0,0,$toh,$remoharr[$i]);
}
$i=$#remoharr+1;
print "removed $i hydroxide ions\n";


write_config_file($filenameoutcfg,0,$config_title[0]);
write_field_file($filenameoutfld,0);
exit 0;
