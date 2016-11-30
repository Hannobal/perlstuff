#!/usr/bin/perl

# adds hydrogen atoms to O^2- ions in suitable positions on the 
# aluminum oxide surface to neutralize the system.

if($#ARGV<3) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file\n";
  print "2. name of corresponding FIELD-file\n";
  print "3. name of output CONFIG-file\n";
  print "4. name of output FIELD-file\n";
  print "5. optional: minimum distance from H+ to residue in A\n";
  exit;
}

if($ARGV[4]>0) {
  $distmin=$ARGV[4]
} else {
  $distmin=2; # minimum distance from H+ to residue
}

$zmin=7; # to obtain only the topmost atoms

open(CONFIG, "<", $ARGV[0]) or die "Can't open CONFIG-file:\n$!";
open(FIELD, "<", $ARGV[1]) or die "Can't open FIELD-file:\n$!";
open(OUTPUT_LOG, ">", "addh.xyz") or die "Can't create log-file:\n$!";

#read header of CONFIG-file
$title=<CONFIG>;
$_=<CONFIG>;
$config_header[0]=$_;
($config_key,$periodic_key) = /^\s*(\S+)\s+(\S+)/;
if($periodic_key!=0) { # read box size if necessary
  for($j=0;$j<=2;$j++) {
    $_=<CONFIG>;
    $config_header[$j+1]=$_;
    ($vec[$j][0], $vec[$j][1], $vec[$j][2]) = /\s*(\S*)\s*(\S*)\s*(\S*)/;
  }
  #define shift vectors for periodic conditions
  $shiftx[0]=0;
  $shiftx[1]=$vec[0][0]+$vec[1][0]+$vec[2][0];
  $shiftx[2]=-$shiftx[1];
  $shifty[0]=0;
  $shifty[1]=$vec[0][1]+$vec[1][1]+$vec[2][1];
  $shifty[2]=-$shifty[1];
}

# read FIELD file and rest of CONFIG file
$netcharge=0;
$foundox=-1;
$foundhydrox=-1;
while(<FIELD>) {
  if(/molecules/ or /MOLECULES/) { last; }
}
($field_nummols) = /^\s*\S+\s+(\S+)/;
#   print "found $field_nummols molecules\n";
if($field_nummols==0) {
  print "no molecules found in FIELD file\n";
  exit;
}

$i=0; #atom index
$tot_negchargesites=0;
for($t=0;$t<$field_nummols;$t++) {
  $j = 0; #number of atoms counted for the molecule
  $_=<FIELD>;
  $mol_name[$t]=$_;
  chomp($mol_name[$t]);
  $startatomindex[$t]=$i;
  if(/^\s*oxygen/i) {
    if(/free/i) {
      $foundox=$t;
    }
  } elsif(/^\s*hydroxide/i) {
    $foundhydrox=$t;
  }
  $_=<FIELD>;
  ($mol_numents[$t])=/^\s*\S+\s+(\S+)/;
  $_=<FIELD>;
  ($mol_numatoms[$t])=/^\s*\S+\s+(\S+)/;
#   print "found ".$mol_numatoms[$t]." in $t. molecule named ".$mol_name[$t]." with ".$mol_numents[$t]." entities\n";
  while(<FIELD>) {
    if(/FINISH/i or /CLOSE/i or /BONDS/i or /ANGLES/i or /DIHEDRALS/i) {
      if(not (/FINISH/i)) {
	while(<FIELD>) { 
	  if(/FINISH/i) { last; }
	}
      }
      last;
    }
    ($mol_atomname[$t][$j],$mol_atommass[$t][$j],$mol_atomcharge[$t][$j])=/^\s*(\S+)\s+(\S+)\s+(\S+)/;
#     printf "%8s %3.5f %3.5f\n", $mol_atomname[$t][$j],$mol_atommass[$t][$j],$mol_atomcharge[$t][$j];
    $netcharge += $mol_numents[$t]*$mol_atomcharge[$t][$j];
    $j++;
  }
  if($j!=$mol_numatoms[$t]) {
    print "error: ".$mol_numatoms[$t]." expected for molecule ".$mol_name[$t]." but ".$j." were found\n";
    exit;
  }
  # read atom coordinates for the molecule from the CONFIG file
  for($m=0;$m<$mol_numents[$t];$m++) {
    $mol_phosox=0;
    @phosoxposx=();
    @phosoxposy=();
    @phosoxposz=();
    for($a=0;$a<$mol_numatoms[$t];$a++) {
      $_=<CONFIG>;
      ($name[$i],$index[$i]) = /^(\S+)\s+(\d+)/;
      if($name[$i] ne $mol_atomname[$t][$a]) {
	print "error: atom name".$name[$i]." in CONFIG-file does not match\n";
	print "name ".$mol_atomname[$t][$a]." in FIELD file for molecule ".$mol_name[$t]."\n";
	exit;
      }
      $_=<CONFIG>;
      ($x[$i],$y[$i],$z[$i]) = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
      if($config_key>0) { #read velocity
	$_=<CONFIG>;
	($vx[$i], $vy[$i], $vz[$i]) = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
	if($config_key>1) { #read force
	  $_=<CONFIG>;
	  ($fx[$i], $fy[$i], $fz[$i]) = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
	}
      }
      if(uc($name[$i]) eq 'OA') {
	if($z[$i]>$zmin) {
	  push(@surfox,$i);
	}
      } elsif(uc($name[$i]) eq 'AL') {
	if($z[$i]>$zmin) {
	  push(@surfal,$i);
	}
      } elsif(uc($name[$i]) eq 'O') {
	if($z[$i]<14 and $mol_atomcharge[$t][$a]<-0.8) {
	  $phosoxposx[$mol_phosox] = $x[$i];
	  $phosoxposy[$mol_phosox] = $y[$i];
	  $phosoxposz[$mol_phosox] = $z[$i];
	  if($mol_phosox>0) {
	    $dmin=100;
	    $dz = $phosoxposz[0]-$phosoxposz[$mol_phosox];
	    for($k=0;$k<3;$k++) {
	      for($l=0;$l<3;$l++) {
		$dx = $phosoxposx[0]-$phosoxposx[$mol_phosox]+$shiftx[$k];
		$dy = $phosoxposy[0]-$phosoxposy[$mol_phosox]+$shifty[$l];
		$dist = sqrt($dx**2+$dy**2+$dz**2);
		if($dist<$dmin) {
		  $dmin=$dist;
		  @endvec=($shiftx[$k],$shifty[$l]);
		}
	      }
	    }
	    $phosoxposx[$mol_phosox] += $endvec[0];
	    $phosoxposy[$mol_phosox] += $endvec[1];
	  }
	  $mol_phosox++;
	}
      } elsif(not uc($name[$i]) eq 'AL' and $z[$i]<=11.0) {
	push(@compare,$i);
      }
#       printf "%3s %5d %3.2f %3.2f %3.2f\n",$name[$i],$i,$x[$i],$y[$i],$z[$i];
      $i++;
    }
    if($mol_phosox>0) {
      for($k=0;$k<@phosoxposx;$k++) {
	$negchargeposx[$tot_negchargesites] += $phosoxposx[$k];
	$negchargeposy[$tot_negchargesites] += $phosoxposy[$k];
	$negchargeposz[$tot_negchargesites] += $phosoxposz[$k];
      }
      $negchargeposx[$tot_negchargesites] /= $mol_phosox;
      $negchargeposy[$tot_negchargesites] /= $mol_phosox;
      $negchargeposz[$tot_negchargesites] /= $mol_phosox;
      $tot_negchargesites++;
    }
  }
  $endatomindex[$t]=$i;
}
close(CONFIG, FIELD);

if($foundox==-1) {
  print "error: could not find entry for oxygen in FIELD file!\n";
  exit;
}
if($foundhydrox==-1) {
  print "error: could not find entry for hydroxide ions in FIELD file!\n";
  exit;
}
$netcharge = sprintf "%.0f", $netcharge;
print "net charge of the system: $netcharge\n";
if($netcharge>=0) {
  print "error: cannot compensate charge by adding H+ ions!\n";
  exit;
}
$netcharge=abs($netcharge);

# for($j=0;$j<@compare;$j++) {
#   printf "%-3s %5d %10.2f %10.2f %10.2f\n",$name[$compare[$j]], $compare[$j], $x[$compare[$j]], $y[$compare[$j]], $z[$compare[$j]];
# }
# print OUTPUT_CONFIG "surfox \n";
#  for($j=0;$j<@surfox;$j++) {
#    printf OUTPUT_CONFIG "%-3s %5d %10.2f %10.2f %10.2f\n",$name[$surfox[$j]], $surfox[$j], $x[$surfox[$j]], $y[$surfox[$j]], $z[$surfox[$j]];
#  }
# print "$#surfox $#compare\n";

############# check distances for possible candidates ###################
print "checking for possible candidates\n";
for($i=0;$i<@surfox;$i++) {
  $possible=1;
  $a=$surfox[$i];
  $dmin=100;
  #calculate position of hydrogen
  for($j=0;$j<@surfal;$j++) {
    $dz = $z[$a]-$z[$surfal[$j]];
    for($k=0;$k<3;$k++) {
      for($l=0;$l<3;$l++) {
	$dx = $x[$a]-$x[$surfal[$j]]+$shiftx[$k];
	$dy = $y[$a]-$y[$surfal[$j]]+$shifty[$l];
	$dist = sqrt($dx**2+$dy**2+$dz**2);
	if($dist<$dmin) {
	  $dmin=$dist;
	  $nearestal=$surfal[$j];
	  $hx[$i]=$x[$a]+sgn($dx)*0.3;
	  $hy[$i]=$y[$a]+sgn($dy)*0.3;
	  $hz[$i]=$z[$a]+0.85;
	}
      }
    }
  }
  #check with environment
  for($j=0;$j<@compare;$j++) {
    $dz = $hz[$i]-$z[$compare[$j]];
    for($k=0;$k<3;$k++) {
      for($l=0;$l<3;$l++) {
	$dx = $hx[$i]-$x[$compare[$j]]+$shiftx[$k];
	$dy = $hy[$i]-$y[$compare[$j]]+$shiftx[$l];
	$dist = sqrt($dx**2+$dy**2+$dz**2);
	if($dist<$distmin) {
	  $possible=0;
	  last;
	}
      }
      if($possible==0) {
	last;
      }
    }
    if($possible==0) {
      last;
    }
  }
  if($possible) {
    push @candidates,$i;
  }
}

if(@candidates<$netcharge) {
  print "error:only ".($#candidates+1)." possible sites were
  found to compensate $netcharge charges!\n";
  exit;
} else {
  print ($#candidates+1);
  print " candidates found.\n";  
}

########## try to fit to phosphonate groups #############################
print "trying to fit to phosphonate groups.\n";
@negchargeposx2=@negchargeposx;
@negchargeposy2=@negchargeposy;
@negchargeposz2=@negchargeposz;
$netchargenew=$netcharge;
while($netchargenew>0 and @negchargeposx2>0) {
  $i=int(rand(@negchargeposx2));
  $dmin=100;
  for($j=0;$j<@candidates;$j++) {
    $dz = @negchargeposz2[$i]-$hz[$candidates[$j]];
    for($k=0;$k<3;$k++) {
      for($l=0;$l<3;$l++) {
	$dx = @negchargeposx2[$i]-$hx[$candidates[$j]]+$shiftx[$k];
	$dy = @negchargeposy2[$i]-$hy[$candidates[$j]]+$shifty[$l];
	$dist = sqrt($dx**2+$dy**2+$dz**2);
	if($dist<$dmin) {
	  $dmin=$dist;
	  $nearestcand=$j;
	}
      }
    }
  }
  push(@used,splice(@candidates,$nearestcand,1));
  splice(@negchargeposx2,$i,1);
  splice(@negchargeposy2,$i,1);
  splice(@negchargeposz2,$i,1);
  $netchargenew--;
}

########## randomly select candidates, remove the O and add the OH ######
while($netchargenew>0) {
  $j=int(rand(@candidates));
  push(@used,splice(@candidates,$j,1));
  $netchargenew--;
}
@used = sort {$a <=> $b} @used;


$asdf=$netcharge+@negchargeposx+@candidates;
print OUTPUT_LOG "$asdf\n";
for($i=0;$i<@negchargeposx;$i++) {
  printf OUTPUT_LOG "\nN   %12.5f %12.5f %12.5f", $negchargeposx[$i], $negchargeposy[$i], $negchargeposz[$i];
}
for($i=0;$i<@candidates;$i++) {
  printf OUTPUT_LOG "\nO   %12.5f %12.5f %12.5f", $x[$surfox[$candidates[$i]]],
  $y[$surfox[$candidates[$i]]],$z[$surfox[$candidates[$i]]];
}

############# write output CONFIG #######################################

open(OUTPUT_CONFIG, ">", $ARGV[2]) or die "Can't create output CONFIG-file:\n$!";
open(OUTPUT_FIELD, ">", $ARGV[3]) or die "Can't create output FIELD-file:\n$!";

print OUTPUT_CONFIG $title;
for $a (@config_header) {
  print OUTPUT_CONFIG $a;
}
$a=0;
$j=0;
for($t=0;$t<$field_nummols;$t++) {
  if($t==$foundox) {
    for($i=$startatomindex[$t];$i<$endatomindex[$t];$i++) {
      if($i==$surfox[$used[$j]]) {
	$j++;
	next;
      }
      $a++;
      printf OUTPUT_CONFIG "%-8s" ,$name[$i];
      printf OUTPUT_CONFIG "%10u\n",$a; #atom index
      printf OUTPUT_CONFIG " %19.11g %19.11g %19.11g\n",$x[$i],$y[$i],$z[$i];
      if($config_key>0) {
	printf OUTPUT_CONFIG " %19.11g %19.11g %19.11g\n",$vx[$i],$vy[$i],$vz[$i];
      }
      if($config_key>1) {
	printf OUTPUT_CONFIG " %19.11g %19.11g %19.11g\n",$fx[$i],$fy[$i],$fz[$i];
      }
    }
  } elsif($t==$foundhydrox) {
    for($i=$startatomindex[$t];$i<$endatomindex[$t];$i++) {
      $a++;
      printf OUTPUT_CONFIG "%-8s" ,$name[$i];
      printf OUTPUT_CONFIG "%10u\n",$a; #atom index
      printf OUTPUT_CONFIG " %19.11g %19.11g %19.11g\n",$x[$i],$y[$i],$z[$i];
      if($config_key>0) {
	printf OUTPUT_CONFIG " %19.11g %19.11g %19.11g\n",$vx[$i],$vy[$i],$vz[$i];
      }
      if($config_key>1) {
	printf OUTPUT_CONFIG " %19.11g %19.11g %19.11g\n",$fx[$i],$fy[$i],$fz[$i];
      }
    }
    for($i=0;$i<@used;$i++) {
      $a++;
      printf OUTPUT_CONFIG "%-8s" ,"OX";
      printf OUTPUT_CONFIG "%10u\n",$a; #atom index
      printf OUTPUT_CONFIG " %19.11g %19.11g %19.11g\n",$x[$surfox[$used[$i]]],$y[$surfox[$used[$i]]],$z[$surfox[$used[$i]]];
      if($config_key>0) {
	printf OUTPUT_CONFIG "%19.11g %19.11g %19.11g\n",$vx[$surfox[$used[$i]]],$vy[$surfox[$used[$i]]],$vz[$surfox[$used[$i]]];
      }
      if($config_key>1) {
	printf OUTPUT_CONFIG "%19.11g %19.11g %19.11g\n",$fx[$surfox[$used[$i]]],$fy[$surfox[$used[$i]]],$fz[$surfox[$used[$i]]];
      }
      $a++;
      printf OUTPUT_CONFIG "%-8s" ,"HG";
      printf OUTPUT_CONFIG "%10u\n",$a; #atom index
      printf OUTPUT_CONFIG "%19.11g %19.11g %19.11g\n",$hx[$used[$i]],$hy[$used[$i]],$hz[$used[$i]];
      printf OUTPUT_LOG "\nH %19.11f %19.11f %19.11f",$hx[$used[$i]],$hy[$used[$i]],$hz[$used[$i]];
      if($config_key>0) {
	printf OUTPUT_CONFIG "%19.11g %19.11g %19.11g\n",0,0,0;
      }
      if($config_key>1) {
	printf OUTPUT_CONFIG "%19.11g %19.11g %19.11g\n",0,0,0;
      }
    }
  } else {
    for($i=$startatomindex[$t];$i<$endatomindex[$t];$i++) {
      $a++;
      printf OUTPUT_CONFIG "%-8s" ,$name[$i];
      printf OUTPUT_CONFIG "%10u\n",$a; #atom index
      printf OUTPUT_CONFIG " %19.11g %19.11g %19.11g\n",$x[$i],$y[$i],$z[$i];
      if($config_key>0) {
	printf OUTPUT_CONFIG " %19.11g %19.11g %19.11g\n",$vx[$i],$vy[$i],$vz[$i];
      }
      if($config_key>1) {
	printf OUTPUT_CONFIG " %19.11g %19.11g %19.11g\n",$fx[$i],$fy[$i],$fz[$i];
      }
    }
  }
}

############# write output FIELD ########################################

open(FIELD, "<", $ARGV[1]) or die "Can't open FIELD-file:\n$!";
while(<FIELD>) {
  if(/^\s*oxygen/ or /^\s*OXYGEN/ or /^\s*Oxygen/ ) {
    print OUTPUT_FIELD $_;
    $_=<FIELD>;
    ($n) = /^\s*\S+\s+(\S+)/;
    $n -= $#used+1;
    print OUTPUT_FIELD "NUMMOLS $n\n";
  } elsif(/^\s*hydroxide/ or /^\s*HYDROXIDE/ or /^\s*Hydroxide/ ) {
    print OUTPUT_FIELD $_;
    $_=<FIELD>;
    ($n) = /^\s*\S+\s+(\S+)/;
    $n += $#used+1;
    print OUTPUT_FIELD "NUMMOLS $n\n";
  } else {
    print OUTPUT_FIELD $_;
    
  }
}
close(FIELD, OUTPUT_FIELD);

sub sgn {
  return ($_[0] <=> 0);
}