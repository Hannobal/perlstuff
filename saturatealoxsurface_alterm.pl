#!/usr/bin/perl
$pi=3.141592;
$doh=1;
$dohal=1.7;

$doo=2.74375;
$tolerance=0.2; # tolerance for surface oxygens
$ohperiodicvec[0][0]=$doo;
$ohperiodicvec[0][1]=0;
$ohperiodicvec[1][0]=$doo*cos($pi/3);
$ohperiodicvec[1][1]=$doo*sin($pi/3);

if($#ARGV<3) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file\n";
  print "2. name of corresponding FIELD-file\n";
  print "3. name of output CONFIG-file\n";
  print "4. name of output FIELD-file\n";
  print "[5. cutoff for surface atoms]\n";
  exit;
}

if($#ARGV>3) {
  $zmin=$ARGV[4]; # to obtain only the topmost atoms
} else {
  $zmin=7;
}

open(CONFIG, "<", $ARGV[0]) or die "Can't open CONFIG-file:\n$!";
open(FIELD, "<", $ARGV[1]) or die "Can't open FIELD-file:\n$!";
# open(OUTPUT_LOG, ">", "addh.xyz") or die "Can't create log-file:\n$!";

#read header of CONFIG-file
$title=<CONFIG>;
$_=<CONFIG>;
$config_header[0]=$_;
($config_key,$periodic_key) = /^\s*(\S+)\s+(\S+)/;
if($periodic_key!=0) { # read box size if necessary
  for($j=0;$j<=2;$j++) {
    $_=<CONFIG>;
    $config_header[$j+1]=$_;
    ($vec[$j][0], $vec[$j][1], $vec[$j][2]) = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
  }
}

# read FIELD file and rest of CONFIG file
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
    if(/FINISH/ or /finish/ or /CLOSE/ or /close/ or /BONDS/
       or /bonds/ or /ANGLES/ or /angles/ or /DIHEDRALS/ or /dihedrals/) {
      if(not (/FINISH/ or /finish/)) {
	while(<FIELD>) { 
	  if(/FINISH/ or /finish/) { last; }
	}
      }
      last;
    }
    ($mol_atomname[$t][$j],$mol_atommass[$t][$j],$mol_atomcharge[$t][$j])=/^\s*(\S+)\s+(\S+)\s+(\S+)/;
#     printf "%8s %3.5f %3.5f\n", $mol_atomname[$t][$j],$mol_atommass[$t][$j],$mol_atomcharge[$t][$j];
    $j++;
  }
  if($j!=$mol_numatoms[$t]) {
    print "error: ".$mol_numatoms[$t]." expected for molecule ".$mol_name[$t]." but ".$j." were found\n";
    exit;
  }
  # read atom coordinates for the molecule from the CONFIG file
  for($m=0;$m<$mol_numents[$t];$m++) {
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
      } elsif(not uc($name[$i]) eq 'AL' and $z[$i]<=11.0) {
	push(@compare,$i);
      }
#       printf "%3s %5d %3.2f %3.2f %3.2f\n",$name[$i],$i,$x[$i],$y[$i],$z[$i];
      $i++;
    }
  $endatomindex[$t]=$i;
  }
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

############# check distances for possible candidates ###################
print "calculating positions for hydroxides\n";
for($i=0;$i<@surfal;$i++) {
  $nearox=0;
  @nearestox=();
  @nearestoxvec=();
  $a=$surfal[$i];
  $dmin=100;
  #calculate position of hydrogen
  for($j=0;$j<@surfox;$j++) {
    $o=$surfox[$j];
    $dz = $z[$a]-$z[$o];
    for($k=-1;$k<2;$k++) {
      for($l=-1;$l<2;$l++) {
	$dx = $x[$a]-$x[$o]+$k*$vec[0][0]+$l*$vec[1][0];#-0.5*$ohperiodicvec[0][0]-0.5*$ohperiodicvec[1][0]; # UNCOMMENT THIS FOR REGULAR LATTICE INSTEAD OF RANDOM POS
	$dy = $y[$a]-$y[$o]+$k*$vec[0][1]+$l*$vec[1][1];#-0.5*$ohperiodicvec[0][1]-0.5*$ohperiodicvec[1][1]; # UNCOMMENT THIS FOR REGULAR LATTICE INSTEAD OF RANDOM POS
	$dist = sqrt($dx**2+$dy**2+$dz**2);
	if($dist<2) {
	  $nearestox[$nearox]=$j;
	  $nearestoxvecx[$nearox]=-$dx/$dist;
	  $nearestoxvecy[$nearox]=-$dy/$dist;
	  $nearestoxvecz[$nearox]=-$dz/$dist;
	  $nearox++;
	}
	if($dist<$dmin) {
	  $dist=$dmin;
	  $asdf=$nearox-1;
	}
      }
    }
  }
  # replace 1 oxygen by hydroxide
  $r=int(rand($nearox));
  # $r=$asdf; # UNCOMMENT THIS FOR REGULAR LATTICE INSTEAD OF RANDOM POS
  $o=$surfox[$nearestox[$r]];
  $newohox[@newohox]=$x[$o];
  $newohoy[@newohoy]=$y[$o];
  $newohoz[@newohoz]=$z[$o];
  $newohhx[@newohhx]=$x[$o]+$doh*cos($pi/3)*( sin($pi/6)*$nearestoxvecy[$r]+cos($pi/6)*$nearestoxvecx[$r]);
  $newohhy[@newohhy]=$y[$o]+$doh*cos($pi/3)*(-sin($pi/6)*$nearestoxvecx[$r]+cos($pi/6)*$nearestoxvecy[$r]);
  $newohhz[@newohhz]=$z[$o]+$doh*sin($pi/3);
  push(@used,splice(@surfox,$nearestox[$r],1));
  # add hydroxide above Al
  $angle=rand(2*$pi);
  $newohox[@newohox]=$x[$a];
  $newohoy[@newohoy]=$y[$a];
  $newohoz[@newohoz]=$z[$a]+$dohal;
  $newohhx[@newohhx]=$x[$a]+$doh*sin($pi/3)*cos($angle);
  $newohhy[@newohhy]=$y[$a]+$doh*sin($pi/3)*sin($angle);
  $newohhz[@newohhz]=$z[$a]+$dohal+$doh*cos($pi/3);
}

@used = sort {$a <=> $b} @used;
# $asdf=$numaddh+@candidates;
# print OUTPUT_LOG "$asdf\n";
# for($i=0;$i<@candidates;$i++) {
#   printf OUTPUT_LOG "\nO   %12.5f %12.5f %12.5f", $x[$surfox[$candidates[$i]]],
#   $y[$surfox[$candidates[$i]]],$z[$surfox[$candidates[$i]]];
# }

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
      if($i==$used[$j]) {
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
    for($i=0;$i<@newohhx;$i++) {
      $a++;
      printf OUTPUT_CONFIG "%-8s" ,"OX";
      printf OUTPUT_CONFIG "%10u\n",$a; #atom index
      printf OUTPUT_CONFIG " %19.11g %19.11g %19.11g\n",$newohox[$i],$newohoy[$i],$newohoz[$i];
      if($config_key>0) {
	printf OUTPUT_CONFIG "%19.11g %19.11g %19.11g\n",0,0,0;
      }
      if($config_key>1) {
	printf OUTPUT_CONFIG "%19.11g %19.11g %19.11g\n",0,0,0;
      }
      $a++;
      printf OUTPUT_CONFIG "%-8s" ,"HG";
      printf OUTPUT_CONFIG "%10u\n",$a; #atom index
      printf OUTPUT_CONFIG "%19.11g %19.11g %19.11g\n",$newohhx[$i],$newohhy[$i],$newohhz[$i];
#       printf OUTPUT_LOG "\nH %19.11f %19.11f %19.11f",$newohhx[$i],$newohhy[$i],$newohhz[$i];
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
  if(/^\s*oxygen/i and /free/i) {
    print OUTPUT_FIELD $_;
    $_=<FIELD>;
    ($n) = /^\s*\S+\s+(\S+)/;
    $n -= $#used+1;
    print OUTPUT_FIELD "NUMMOLS $n\n";
  } elsif(/^\s*hydroxide/i) {
    print OUTPUT_FIELD $_;
    $_=<FIELD>;
    ($n) = /^\s*\S+\s+(\S+)/;
    $n += $#newohhx+1;
    print OUTPUT_FIELD "NUMMOLS $n\n";
  } else {
    print OUTPUT_FIELD $_;
    
  }
}
close(FIELD, OUTPUT_FIELD);