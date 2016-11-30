#!/usr/bin/perl

# adds hydrogen atoms to O^2- ions in suitable positions on the 
# aluminum oxide surface of a spherical nanoparticle to neutralize the system.

use dlpoly_utility;
use hanno_utility;
use strict vars;

my($incfg,$infld,$outcfg,$outfld,$sphrad,$m,$t);
our($toh,$tal,$to,@alcands,@ocands,@alcheck,@ocheck,$rmradsq,$chkradsq,$dcoordsq,$totcharge);


if($#ARGV<4) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file\n";
  print "2. name of corresponding FIELD-file\n";
  print "3. name of output CONFIG-file\n";
  print "4. name of output FIELD-file\n";
  print "5. particle radius in A\n";
  exit;
}

$incfg=$ARGV[0];
$infld=$ARGV[1];
$outcfg=$ARGV[2];
$outfld=$ARGV[3];
$sphrad=$ARGV[4];

exit 1 if(read_field_file($infld,0)!=0);
exit 1 if(read_config_file($incfg,0,0)!=0);

$rmradsq  = ($sphrad-2)**2; # alox can be removed, oxygen protonated if outside of this radius
$chkradsq = ($sphrad-5)**2;  # only use atoms outside of this radius to check coordination number
$dcoordsq = 2.5**2;

$toh = -1;
$tal = -1;
$to  = -1;
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t]=~/hydrox/i) {
    $toh = $t;
  } elsif($mol_name[0][$t]=~/oxygen/i and not $mol_name[0][$t]=~/(frozen|fix)/i) {
    $to  = $t;
  } elsif($mol_name[0][$t]=~/aluminum/i and not $mol_name[0][$t]=~/(frozen|fix)/i) {
    $tal = $t;
  }
}

&error("cannot find aluminum, oxygen and/or hoydroxide in FIELD file!",2) if($toh<0 or $tal<0 or $to<0);

# calc the net charge
$totcharge=0;
for($t=0;$t<$field_nummols[0];$t++) {
  $totcharge+=$mol_numents[0][$t]*$mol_charge[0][$t];
}


&build_lists;
for($m=0;$m<@alcheck;$m++) {
  print "$alcheck[$m] ";
}
print "\n";
&calc_coordnr;

sub build_lists {
  my($m,$dsq);
  @alcands=(); # al atoms that can be removed
  @ocands=();  # o atoms that can be protonated
  @alcheck=(); # al atoms for coordination check of o atoms
  @ocheck=();  # o atoms for coordination check of al atoms
  for($m=0;$m<$mol_numents[0][$tal];$m++) {
    $dsq = $cdata[0][$tal][$m][0][0]*$cdata[0][$tal][$m][0][0]
	+ $cdata[0][$tal][$m][0][1]*$cdata[0][$tal][$m][0][1]
	+ $cdata[0][$tal][$m][0][2]*$cdata[0][$tal][$m][0][2];
    push(@alcands,[$m,0]) if($dsq>$rmradsq);
    push(@alcheck,$m) if($dsq>$chkradsq);
  }

  for($m=0;$m<$mol_numents[0][$to];$m++) {
    $dsq = $cdata[0][$to][$m][0][0]*$cdata[0][$to][$m][0][0]
	+ $cdata[0][$to][$m][0][1]*$cdata[0][$to][$m][0][1]
	+ $cdata[0][$to][$m][0][2]*$cdata[0][$to][$m][0][2];
    push(@ocands,[$m,0]) if($dsq>$rmradsq);
    push(@ocheck,$m) if($dsq>$chkradsq);
  }
}

sub calc_coordnr {
  my($i,$j,$m,$o,$dx,$dy,$dz);
  # for al atoms
  for($i=0;$i<@alcands;$i++) {
    $m=$alcands[$i][0];
    foreach $o (@ocheck) {
      $dx = $cdata[0][$to][$o][0][0] - $cdata[0][$tal][$m][0][0];
      $dy = $cdata[0][$to][$o][0][1] - $cdata[0][$tal][$m][0][1];
      $dz = $cdata[0][$to][$o][0][2] - $cdata[0][$tal][$m][0][2];
      if($dx*$dx+$dy*$dy+$dz*$dz<$dcoordsq) {
	$alcands[$i][1]++;
      }
    }
    for($o=0;$o<$mol_numents[0][$toh];$o++) {
      $dx = $cdata[0][$to][$o][0][0] - $cdata[0][$tal][$m][0][0];
      $dy = $cdata[0][$to][$o][0][1] - $cdata[0][$tal][$m][0][1];
      $dz = $cdata[0][$to][$o][0][2] - $cdata[0][$tal][$m][0][2];
      if($dx*$dx+$dy*$dy+$dz*$dz<$dcoordsq) {
	$alcands[$i][1]++;
      }
    }
#     if($alcands[$i][1]<6) {
#       print "$m ";
#     }
  }
  # for o atoms
  for($i=0;$i<@ocands;$i++) {
    $o=$ocands[$i][0];
    foreach $m (@alcheck) {
      $dx = $cdata[0][$tal][$m][0][0] - $cdata[0][$to][$o][0][0];
      $dy = $cdata[0][$tal][$m][0][1] - $cdata[0][$to][$o][0][1];
      $dz = $cdata[0][$tal][$m][0][2] - $cdata[0][$to][$o][0][2];
      if($dx*$dx+$dy*$dy+$dz*$dz<$dcoordsq) {
	$ocands[$i][1]++;
      }
    }
    if($ocands[$i][1]<3) {
      print "$o ";
    }
  }
}

sub error {
  print "**** error: $_[0]\n";
  if(defined($_[1])) {
    exit $_[1];
  } else {
    exit 99;
  }
}