#!/usr/bin/perl

use Switch;
use hanno_utility;
use dlpoly_utility;
use aloxsam_utility;
use Storable qw(dclone);

$minpadist = -1;
$maxohdist = 2.5;
$lrmwater  = 0;
$ldouble   = 0;
$postol    = 1;
@lattice   = ();
@latticecenter=(0,0);

if($#ARGV<3) {
  print "input format:\n";
  print " 1. name of CONFIG/REVCON-file with surface\n";
  print " 2. name of corresponding FIELD-file\n";
  print " 3. name of output CONFIG-file\n";
  print " 4. name name of output FIELD-file\n";
  print "-h <real>   maximum length of H-bond (default: $maxohdist)\n";
  print "-p <real>   minimum distance between PAs (recommended: 3.5)\n";
  print "-w          remove water (default: no)\n";
  print "-l <4*real> define a lattice\n";
  print "-c <2*real> define center for the lattice\n";
  print "-t <real>   tolerance for lattice positions (default: $postol)\n";
  exit;
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
    } case (/^-d/) {
      $ldouble=1;
    } case (/^-p/) {
      $i++;
      if(not check_real($ARGV[$i])) {
        print "**** error: minimum distance between PAs must be a real number!\n";
        exit 99;
      }
      $minpadist=$ARGV[$i];
    }  case (/^-t/) {
      $i++;
      if(not check_real($ARGV[$i])) {
        print "**** error: tolerance for lattice positions must be a real number!\n";
        exit 99;
      }
      $postol=$ARGV[$i];
    }  case (/^-t/) {
      if(not check_real(@ARGV[$i+1..$i+4])) {
        print "**** error: lattice vector components must be four real numbers!\n";
        exit 99;
      }
      $lattice[0][0]=$ARGV[$i+1];
      $lattice[0][1]=$ARGV[$i+2];
      $lattice[1][1]=$ARGV[$i+3];
      $lattice[1][0]=$ARGV[$i+4];
      $i+=4;
    }  case (/^-t/) {
      if(not check_real(@ARGV[$i+1..$i+2])) {
        print "**** error: lattice center must be two real numbers!\n";
        exit 99;
      }
      $latticecenter[1][1]=$ARGV[$i+3];
      $latticecenter[1][0]=$ARGV[$i+4];
      $i+=4;
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


our @tpa  = ();
our @tox  = ();
our $twat = -1;
our $toh  = -1;
$err  = 0;
our @acidh = ();

for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /water/i) {
    $twat = $t;
  } elsif($mol_name[0][$t] =~ /hydroxide/i) {
    $toh = $t;
  } elsif($mol_name[0][$t] =~ /oxygen/i) {
    push(@tox,$t);
  } elsif($mol_name[0][$t] =~ /-PA$/) {
    push(@tpa,[$t,-1,-1]);
  }
}

if($twat<0 and not $lrmwater) {
  print "**** error: did not find entry for water in FIELD file\n";
  $err=1;
}
if($toh<0) {
  print "**** error: did not find entry for hydroxide ions in FIELD file\n";
  $err=1;
}

if(not @tpa) {
  print "**** error: no PAs were found!\n";
  $err=1;
}

for($i=0;$i<@tpa;$i++) {
  for($t=0;$t<$field_nummols[0];$t++) {
    if($mol_name[0][$t] =~ /$mol_name[0][$tpa[$i][0]]-d/) {
      $tpa[$i][1] = $t;
    }
    if($mol_name[0][$t] =~ /$mol_name[0][$tpa[$i][0]]-2d/) {
      $tpa[$i][2] = $t;
    }
  }
  if($tpa[$i][1]<0) {
    print "**** error: could not find singly deprotonated form of $mol_name[0][$tpa[$i][0]]\n";
    $err=1;
    next;
  }
  if($tpa[$i][2]<0 and $ldouble) {
    print "**** error: could not find doubly deprotonated form of $mol_name[0][$tpa[$i][0]]\n";
    $err=1;
    next;
  }
  $t   = $tpa[$i][0];
  $td  = $tpa[$i][1];
  $t2d = $tpa[$i][2] if($ldouble);
}

exit 99 if($err>0);
$zsurf=9e20;
for($i=0;$i<@tpa;$i++) {
  foreach $t (@{$tpa[$i]}) {
    next if($t<0);
    @{$acidh[$t]}=find_acidic_protons(0,$t);
    @{$acido[$t]}=find_acidic_oxygens(0,$t);
    for($m=0;$m<$mol_numents[0][$t];$m++) {
      for($a=0;$a<$mol_numatoms[0][$t];$a++) {
        $zsurf=min($zsurf,$cdata[0][$t][$m][$a][2]);
      }
    }
  }
}
$zsurf-=$maxohdist+0.1;

@ocands=find_ocands($zsurf);

# do the replacements
$changed=0;
$err=0;
looppas : for($i=0;$i<@tpa;$i++) {
  if($ldouble) {
    $jmax=2;
  } else {
    $jmax=1;
  }
  for($j=0;$j<$jmax;$j++) {
    $t=$tpa[$i][$j];
    $td=$tpa[$i][$j+1];
    next looppas if($td<0);
    for($m=$mol_numents[0][$t]-1;$m>=0;$m--) {
      $dsqmin = 999;
      foreach $a (@{$acidh[$t]}) {
        foreach $to (@tox,$toh) {
          for($o=0;$o<@{$ocands[$to]};$o++) {
            $mo=$ocands[$to][$o];
            $dsq=calc_dsq_orthocell(0,$t,$m,$a, $to,$mo,0);
            if($dsq<$dsqmin) {
              $dsqmin=$dsq;
              $tomin=$to;
              $omin=$o;
              $amin=$a;
            }
          }
        }
      }
#       print "t $t m $m a $a to $tomin mo $momin dist ",sqrt($dsqmin),"\n";
      if(sqrt($dsqmin)>$maxohdist) {
#         print "no proton transfer possible for $mol_name[0][$t] number $m\n";
        next;
      }
      $changed++;
      # deprotonate acid
      @proton=splice(@{$cdata[0][$t][$m]},$amin,1);
      if($j==0) {
        # adjust position of O/OH
        if($amin==$acidh[$t][0]) {$aold=$acidh[$t][1];} else {$aold=$acidh[$t][0];}
        $at[0][0] = $mol_bondatoms[0][$t][$aold][0]; # oh to remain oh
        $at[1][0] = $mol_bondatoms[0][$t][$amin][0]; # oh to be made o
        $at[2][0] = $acido[$t][0];
        $at[3][0] = $aold;
        $at[0][1] = $mol_bondatoms[0][$td][$acidh[$td][0]][0];
        $at[1][1] = $acido[$td][0];
        $at[2][1] = $acido[$td][1];
        $at[3][1] = $acidh[$td][0];
        @at = sort { $b->[0] <=> $a->[0] } @at;
        for($k=0;$k<@at;$k++) {
          $at[$k][0]-- if($at[$k][0]>$amin); # proton has been spliced already
          $at[$k][2]=splice(@{$cdata[0][$t][$m]},$at[$k][0],1);
        }
        @at = sort { $a->[1] <=> $b->[1] } @at;
        for($k=0;$k<@at;$k++) {
          splice(@{$cdata[0][$t][$m]},$at[$k][1],0,[@{$at[$k][2]}]);
        }
      }
      push(@{$cdata[0][$td]},splice(@{$cdata[0][$t]},$m,1));
      $mol_numents[0][$t]--;
      $mol_numents[0][$td]++;
      # protonate oxygen/hydroxide
      $mo=$ocands[$tomin][$omin];
      @tmp=splice(@{$cdata[0][$tomin]},$mo,1);
      $mol_numents[0][$tomin]--;
      push(@{$tmp[0]},@proton);
      if($tomin==$toh) {
        if(not $lrmwater) {
          push(@{$cdata[0][$twat]},@tmp);
          $cdata[0][$twat][$mol_numents[0][$twat]][2][9]=$mol_atomdata[0][$twat][2][0];
          $mol_numents[0][$twat]++;
        }
      } else {
        $tmp[0][1][9]=$mol_atomdata[0][$toh][1][0];
        push(@{$cdata[0][$toh]},splice(@tmp,0,1));
        $mol_numents[0][$toh]++;
      }
#       remove_index_from_list($ocands[$tomin],$omin);
        @ocands=find_ocands($zsurf);
      push(@{$ocands[$toh]},$#{$cdata[$toh]}) if($tomin!=$toh);
    }
  }
}

write_config_file($filenameoutcfg,0,$config_title[0],">",0);
# write_config_file($filenameoutcfg,0);
write_field_file($filenameoutfld,0);
print "transferred $changed protons\n";
exit 0;

sub remove_index_from_list {
  my @arr=@{$_[0]};
  print "asdg $#arr\n";
  my $index=$_[1]; # index within the list, not the value stored in the array!
  for(my $i=0;$i<@arr;$i++) {
    $arr[$i]-- if($arr[$i]>$arr[$index]);
  }
  splice(@arr,$index,1);
}

sub find_ocands {
  my @lst=();
  foreach my $t (@tox,$toh) {
    for(my $m=0;$m<$mol_numents[0][$t];$m++) {
      push(@{$lst[$t]},$m) if($cdata[0][$t][$m][0][2]>$_[0]);
    }
  }
  return @lst;
}