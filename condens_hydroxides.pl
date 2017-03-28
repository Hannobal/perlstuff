#!/usr/bin/perl

use Switch;
use hanno_utility;
use dlpoly_utility;
use Storable qw(dclone);

$lrmwater  = 0;
$maxohdist = 2.5;

if($#ARGV<3) {
  print "input format:\n";
  print " 1. name of CONFIG/REVCON-file with surface\n";
  print " 2. name of corresponding FIELD-file\n";
  print " 3. name of output CONFIG-file\n";
  print " 4. name name of output FIELD-file\n";
  print "-h <real>   maximum length of H-bond (default: $maxohdist)\n";
  print "-w          remove water (default: no)\n";
  exit 99;
}

$filenameincfg  = $ARGV[0];
$filenameinfld  = $ARGV[1];
$filenameoutcfg = $ARGV[2];
$filenameoutfld = $ARGV[3];

for($i=4;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-h/) {
      $i++;
      if(not check_real($ARGV[$i])) {
        print "**** error: maximum length of H-bond must be a real number!\n";
        exit 99;
      }
      $maxohdist=$ARGV[$i];
    } case (/^-w/) {
      $lrmwater=1;
    } else {
      print "**** error: unknown flag $ARGV[$i]\n";
      exit 99;
    }
  }
}

# read FIELD file
exit 99 if(read_field_file($filenameinfld,0)!=0);

#read CONFIG file
exit 99 if(read_config_file($filenameincfg,0,0)!=0);

our $toh  = -1;
our @tox  = ();
our $twat = -1;
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /water/i) {
    $twat = $t;
  } elsif($mol_name[0][$t] =~ /hydroxide/i) {
    $toh = $t;
  } elsif($mol_name[0][$t] =~ /oxygen/i) {
    push(@tox,$t);
  }
}

if($twat<0 and not $lrmwater) {
  print "**** error: did not find entry for water in FIELD file\n";
  $err=1;
}
if(not @tox) {
  print "**** error: did not find entry for oxygen in FIELD file\n";
  $err=1;
}
exit 99 if($err>0);

$changed=0;
for($m=0;$m<$mol_numents[0][$toh];$m++) {
  $dsqmin = 999;
  for($m2=0;$m2<$mol_numents[0][$toh];$m2++) {
    next if($m==$m2);
    $dsq=calc_dsq_orthocell(0,$toh,$m,1, $toh,$m2,0);
    if($dsq<$dsqmin) {
      $dsqmin=$dsq;
      $mmin=$m2;
    }
  }
  next if(sqrt($dsqmin)>$maxohdist);
  $changed++;
  # deprotonate hydroxide
  @proton=splice(@{$cdata[0][$toh][$m]},1,1);
  if($mmin>$m) {
    @tmp=splice(@{$cdata[0][$toh]},$mmin,1);
    push(@{$cdata[0][$tox[0]]},splice(@{$cdata[0][$toh]},$m,1));
  } else {
    push(@{$cdata[0][$tox[0]]},splice(@{$cdata[0][$toh]},$m,1));
    @tmp=splice(@{$cdata[0][$toh]},$mmin,1);
  }
  $mol_numents[0][$tox[0]]++;
  push(@{$tmp[0]},@proton);
  $mol_numents[0][$toh]-=2;
  if(not $lrmwater) {
    push(@{$cdata[0][$twat]},@tmp);
    $cdata[0][$twat][$mol_numents[0][$twat]][2][9]=$mol_atomdata[0][$twat][2][0];
    $mol_numents[0][$twat]++;
  }
#   last;
  $m--;
}

write_config_file($filenameoutcfg,0,$config_title[0],">",0);
write_field_file($filenameoutfld,0);
print "transferred $changed protons\n";
exit 0;