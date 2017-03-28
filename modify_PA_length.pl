#!/usr/bin/perl

# still needs generalization. Currently only works for
# a bunch of C18-PA that are cut to C14/C10-PA to yield
# a 1:1:1 ratio

use dlpoly_utility;
use aloxsam_utility;
use Storable qw(dclone);

if($#ARGV<3) {
  print "input format:\n";
  print " 1. name of CONFIG/REVCON-file with surface\n";
  print " 2. name of corresponding FIELD-file\n";
  print " 3. name of output CONFIG-file\n";
  print " 4. name name of output FIELD-file\n";
  exit 1;
}

$filenameincfg  = $ARGV[0];
$filenameinfld  = $ARGV[1];
$filenameoutcfg = $ARGV[2];
$filenameoutfld = $ARGV[3];


# read FIELD file
exit 99 if(read_field_file($filenameinfld,0)!=0);

#read CONFIG file
exit 99 if(read_config_file($filenameincfg,0,0)!=0);

our @tpa  = (); # 1st index: protonation (PA, PA-d, PA-2d) 2nd: chain length
our @atid = ();

for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /C(\d+)-PA/) {
    $j=$1;
    if($mol_name[0][$t] =~ /-PA$/) {
      $i=0;
    } elsif($mol_name[0][$t] =~ /-PA-d$/) {
      $i=1;
    } elsif($mol_name[0][$t] =~ /-PA-2d$/) {
      $i=2;
    } else {
      print "PA $mol_name[0][$t] name has an unknown pattern!\n";
      exit 1;
    }
    $tpas[$i][$j]=$t;
    @{$atid[$t]}=get_atids(0,$t);
#     print "$mol_name[0][$t] ",join(" ",@{$atid[$t]}),"\n";
    if(not @{$atid[$t]}) {
      print "**** error analyzing the structure of $mol_name[0][$t]\n";
      exit 1;
    }
  }
}

# generate a 1:1:1 ratio by cutting a third of the C18 to C14 and C10
for($i=0;$i<3;$i++) {
  $ts=$tpas[$i][18];
  $jmax=int($mol_numents[0][$ts]/3.0);
  for($j=0;$j<$jmax;$j++) {
    if(cut_PA(0,0, $ts,int(rand($mol_numents[0][$ts])), $tpas[$i][14])!=0) {
      print "**** error cutting PA\n"; exit 1;
    }
    if(cut_PA(0,0, $ts,int(rand($mol_numents[0][$ts])), $tpas[$i][10])!=0) {
      print "**** error cutting PA\n"; exit 1;
    }
  }
}

write_config_file($filenameoutcfg,0,$config_title[0],">",0);
write_field_file($filenameoutfld,0);

sub get_atids {
  my $fi = $_[0];
  my $t  = $_[1];
  my @arr=();
  my($i,$j,$last,$curr,$next);
  # first is phosphorous
  $arr[0]=find_pa_phosphorous($fi,$t);
  # add o first
  foreach $i (@{$mol_bondatoms[$fi][$t][$arr[0]]}) {
    push(@arr,$i) if($mol_atomdata[$fi][$t][$i][0]=~/^O$/i);
  }
  # add OH groups
  foreach $i (@{$mol_bondatoms[$fi][$t][$arr[0]]}) {
    next if($mol_atomdata[$fi][$t][$i][0]!~/^OH$/i);
    push(@arr,$i);
    foreach $j (@{$mol_bondatoms[$fi][$t][$i]}) {
#       print "$mol_atomdata[$fi][$t][$i][0] $mol_atomdata[$fi][$t][$j][0]\n";
      push(@arr,$j) if($mol_atomdata[$fi][$t][$j][0]=~/^H/);
    }
  }
  $last=$arr[0];
  $curr=-1;
  foreach $i (@{$mol_bondatoms[$fi][$t][$arr[0]]}) {
    if($mol_atomdata[$fi][$t][$i][0]=~/^[CH]/i) {
      $curr=$i;
      last;
    }
  }
  # now add CH2 groups and CH3 group
  while($last!=$curr) {
    push(@arr,$curr);
    foreach $i (@{$mol_bondatoms[$fi][$t][$curr]}) {
      if($mol_atomdata[$fi][$t][$i][0]=~/^C/i) {
        $next=$i if($i!=$last);
      } elsif($mol_atomdata[$fi][$t][$i][0]=~/^H/i) {
        push(@arr,$i);
      } elsif($mol_atomdata[$fi][$t][$i][0]!~/^P5/i) {
        print "unexpected atom $mol_atomdata[$fi][$t][$i][0] in $mol_name[$fi][$t]\n";
        return undef;
      }
    }
    $last=$curr;
    $curr=$next;
  }
  return @arr;
}

sub cut_PA {
  my $fi = $_[0];
  my $ci = $_[1];
  my $ts = $_[2]; #source molecule-id
  my $ms = $_[3]; #source entity-id
  my $tt = $_[4]; #target molecule-id
  return 1 if(not (@{$atid[$ts]} and @{$atid[$tt]}));
  return 1 if($#{$atid[$ts]} < $#{$atid[$tt]});
  return 1 if($ms>=$mol_numents[$ci][$ts]);
  my $mt = $mol_numents[$fi][$tt];
  for(my $i=0;$i<@{$atid[$tt]};$i++) {
    @{$cdata[$ci][$tt][$mt][$atid[$tt][$i]]}
      = @{dclone \@{$cdata[$ci][$ts][$ms][$atid[$ts][$i]]}};
    $cdata[$ci][$tt][$mt][$atid[$tt][$i]][9]=$mol_atomdata[$ci][$tt][$atid[$tt][$i]][0];
  }
  remove_mol_entity($ci,$fi,$ts,$ms);
  $mol_numents[$ci][$tt]++;
  return 0;
}