#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;
use aloxsam_utility;
use Switch;

$maxdist=3;
$maxdsqo=0;
$maxdsqoh=0;
$luseoh=0;
$luseox=0;

if($#ARGV<3) {
  print "input format:\n";
  print " 1. name of CONFIG/REVCON-file with surface\n";
  print " 2. name of corresponding FIELD-file\n";
  print " 3. name of output CONFIG-file\n";
  print " 4. name name of output FIELD-file\n";
  print "-d  <real>  cutoff for distance to cations (default: $maxdist)\n";
  print "-oh <real>  include hydroxides with specified cutoff for H-bonds\n";
  print "-o  <real>  include oxygens with specified cutoff for H-bonds\n";
  exit 1;
}

$filenameincfg  = $ARGV[0];
$filenameinfld  = $ARGV[1];
$filenameoutcfg = $ARGV[2];
$filenameoutfld = $ARGV[3];


for($i=4;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-d$/) {
      $i++;
      if(not check_real($ARGV[$i])) {
        print "**** error: maximum distance must be a real number!\n";
        exit 1;
      }
      $maxdist=$ARGV[$i];
    } case (/^-oh$/) {
      $luseoh=1;
      $i++;
      if(not check_real($ARGV[$i])) {
        print "**** error: maximum distance must be a real number!\n";
        exit 1;
      }
      $maxdistoh=$ARGV[$i];
    } case (/^-o$/) {
      $luseox=1;
      $i++;
      if(not check_real($ARGV[$i])) {
        print "**** error: maximum distance must be a real number!\n";
        exit 1;
      }
      $maxdisto=$ARGV[$i];
    } else {
      print "**** error: expected optional flag but got \"$ARGV[$i]\"!\n";
      exit 1;
    }
  }
}

# read FIELD file
exit 99 if(read_field_file($filenameinfld,0)!=0);

#read CONFIG file
exit 99 if(read_config_file($filenameincfg,0,0)!=0);

$zmin=999;
@tpas=();
@acidh=();
@acido=();
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /-PA/) {
    push(@tpas,$t);
    @{$acidh[$t]}=find_acidic_protons(0,$t);
    @{$acido[$t]}=find_acidic_oxygens(0,$t,1);
    for($m=0;$m<$mol_numents[0][$t];$m++) {
      for($a=0;$a<$mol_numatoms[0][$t];$a++) {
        $zmin=min($cdata[0][$t][$m][$a][2],$zmin);
      }
    }
  }
}
$zmin-=$maxdist+0.5;

@cations=();
@oxides=();
@hydroxides=();
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /(iron|aluminum|titanium)/i) {
    for($m=0;$m<$mol_numents[0][$t];$m++) {
      push(@cations,[$t,$m]) if($cdata[0][$t][$m][0][2]>$zmin)
    }
  } elsif($luseoh and $mol_name[0][$t] =~ /hydroxid/i) {
    for($m=0;$m<$mol_numents[0][$t];$m++) {
      push(@hydroxides,[$t,$m]) if($cdata[0][$t][$m][0][2]>$zmin)
    }
  } elsif($luseox and $mol_name[0][$t] =~ /(oxid|oxygen)/i) {
    for($m=0;$m<$mol_numents[0][$t];$m++) {
      push(@oxides,[$t,$m]) if($cdata[0][$t][$m][0][2]>$zmin)
    }
  }
}

$maxdsq=$maxdist*$maxdist;
$maxdsqo=$maxdisto*$maxdisto;
$maxdsqoh=$maxdistoh*$maxdistoh;
$removed=0;
foreach $t (@tpas) {
  loopmols:for($m=0;$m<$mol_numents[0][$t];$m++) {
    # check for bond between PA-oxygen and cation
    for($i=0;$i<@cations;$i++) {
      $to=$cations[$i][0];
      $mo=$cations[$i][1];
      foreach $a (@{$acido[$t]}) {
        $dsq=calc_dsq_orthocell(0,$t,$m,$a, $to,$mo,0);
        next loopmols if($dsq<$maxdsq);
      }
    }
    # check for bond between PA and hydroxide
    for($i=0;$i<@hydroxides;$i++) {
      $to=$hydroxides[$i][0];
      $mo=$hydroxides[$i][1];
      foreach $a (@{$acido[$t]}) {
        $dsq=calc_dsq_orthocell(0,$t,$m,$a, $to,$mo,1);
        next loopmols if($dsq<$maxdsqoh);
      }
      foreach $a (@{$acidh[$t]}) {
        $dsq=calc_dsq_orthocell(0,$t,$m,$a, $to,$mo,0);
        next loopmols if($dsq<$maxdsqoh);
      }
    }
    # check for bond between PA-hydrogen and oxide
    for($i=0;$i<@oxides;$i++) {
      $to=$oxides[$i][0];
      $mo=$oxides[$i][1];
      foreach $a (@{$acido[$t]}) {
        $dsq=calc_dsq_orthocell(0,$t,$m,$a, $to,$mo,0);
        next loopmols if($dsq<$maxdsqo);
      }
    }
    # no bond was found, so remove the molecule
    remove_mol_entity(0,0,$t,$m);
    $m--;
    $removed++;
  }
}

if($removed>0) {
  print "removed $removed molecules\n";
  write_config_file($filenameoutcfg,0,$config_title[0],">",0);
  write_field_file($filenameoutfld,0);
  exit 0;
} else {
  print "no changes made\n";
  exit 1;
}