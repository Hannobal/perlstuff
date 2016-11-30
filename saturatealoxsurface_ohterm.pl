#!/usr/bin/perl
$pi        = 3.141592;
$doh       = 1;
$dohal     = 1.7;
$tolerance = 0.2;
$bothsides = 0;
$thickness = 0;

use hanno_utility;

if($#ARGV<3) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file\n";
  print "2. name of corresponding FIELD-file\n";
  print "3. name of output CONFIG-file\n";
  print "4. name of output FIELD-file\n";
  print "[5. \"b\" for hydroxylation of both sides]\n";
  print "[6. thickness of the slab]\n";
  exit;
}


if($#ARGV==4) {
  if($ARGV[4] =~ /^b/i) {
    $bothsides = 1;
  } elsif(check_real($ARGV[4])) {
    $thickness=$ARGV[4]; # to obtain only the topmost atoms
  } else {
    print "**** error: \"$ARGV[4]\" is not a valid input at this point!\n";
    exit 1;
  }
} elsif($#ARGV>4) {
  if($ARGV[4] =~ /^b/i) {
    $bothsides = 1;
  } else {
    print "**** error: expected \"b\"!\n";
    exit 1;
  }
  if(check_real($ARGV[5])) {
    $thickness = $ARGV[5];
  } else {
    print "**** error: cutoff must be a real number!\n";
  }
}

open(CONFIG, "<", $ARGV[0]) or die "Can't open CONFIG-file:\n$!";
open(FIELD, "<", $ARGV[1]) or die "Can't open FIELD-file:\n$!";
# open(OUTPUT_LOG, ">", "addh.xyz") or die "Can't create log-file:\n$!";

#read header of CONFIG-file
$title=<CONFIG>;
$_=<CONFIG>;
$config_header[0]=$_;
($config_key,$periodic_key) = /^\s+(\S+)\s+(\S+)/;
if($periodic_key!=0) { # read box size if necessary
  for($j=0;$j<=2;$j++) {
    $_=<CONFIG>;
    $config_header[$j+1]=$_;
    ($vec[$j][0], $vec[$j][1], $vec[$j][2]) = /\s*(\S*)\s*(\S*)\s*(\S*)/;
  }
}

# read FIELD file and rest of CONFIG file
$foundox=-1;
$foundhydrox=-1;
$foundal=-1;
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
$alzmin = 9999;
$alzmax = -9999;
$oazmin = 9999;
$oazmax = -9999;
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
  } elsif(/^\s*aluminum/i ) {
    if(/free/i) {
      $foundal=$t;
    }
  }
  $_=<FIELD>;
  ($mol_numents[$t])=/^\s*\S+\s+(\S+)/;
  $_=<FIELD>;
  ($mol_numatoms[$t])=/^\s*\S+\s+(\S+)/;
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
  # read atom coordinates for the molecule from the CONFIG file# indices for cdata array:
  # 0 name
  # 1 x
  # 2 y
  # 3 z
  # 4 velocity x
  # 5 velocity y
  # 6 velocity z
  # 7 force x
  # 8 force y
  # 9 force z
  for($m=0;$m<$mol_numents[$t];$m++) {
    for($a=0;$a<$mol_numatoms[$t];$a++) {
      $_=<CONFIG>;
      ($cdata[$t][$m][$a][0]) = /^\s*(\S+)/;
      if($cdata[$t][$m][$a][0] ne $mol_atomname[$t][$a]) {
	print "error: atom name".$name[$t][$m][$a]." in CONFIG-file does not match\n";
	print "name ".$mol_atomname[$t][$a]." in FIELD file for molecule ".$mol_name[$t]."\n";
	exit;
      }
      $_=<CONFIG>;
      ($cdata[$t][$m][$a][1],$cdata[$t][$m][$a][2],$cdata[$t][$m][$a][3]) = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
      if($config_key>0) { #read velocity
	$_=<CONFIG>;
	($cdata[$t][$m][$a][4],$cdata[$t][$m][$a][5],$cdata[$t][$m][$a][6]) = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
	if($config_key>1) { #read force
	  $_=<CONFIG>;
	  ($cdata[$t][$m][$a][7],$cdata[$t][$m][$a][8],$cdata[$t][$m][$a][9]) = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
	}
      }
      if(uc($cdata[$t][$m][$a][0]) eq 'AL') {
	$alzmin = $cdata[$t][$m][$a][3] if($cdata[$t][$m][$a][3]<$alzmin);
	$alzmax = $cdata[$t][$m][$a][3] if($cdata[$t][$m][$a][3]>$alzmax);
      } elsif(uc($cdata[$t][$m][$a][0]) eq 'OA') {
	$oazmin = $cdata[$t][$m][$a][3] if($cdata[$t][$m][$a][3]<$oazmin);
	$oazmax = $cdata[$t][$m][$a][3] if($cdata[$t][$m][$a][3]>$oazmax);
      } else {
	print "**** error: script only works for plain AlO-surfaces! \n";
	exit 1;
      }
    }
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
if($foundal==-1) {
  print "error: could not find entry for aluminum ions in FIELD file!\n";
  exit;
}

############# actual replacement ########################################

if($oazmax-$oazmin < $thickness) {
  print "**** warning: slab is too thing for desired thickness ($thickness A).\n";
  $thickness = $oazmax-$oazmin;
  print "              new thickness set to $thickness A.\n";
}

if($thickness eq 0 and $bothsides) {
  $zmincutoff = $oazmin+$tolerance;
  $zmaxcutoff = $oazmax-$tolerance;
} elsif($thickness ne 0 and $bothsides) {
  $zmincutoff = $oazmin+$tolerance;
  $zmaxcutoff = $oazmin+$thickness-$tolerance;
} elsif($thickness ne 0 and not $bothsides) {
  $zmincutoff = $oazmin-999;
  $zmaxcutoff = $oazmin+$thickness-$tolerance;
} else {
  $zmincutoff = $oazmin-999;
  $zmaxcutoff = $oazmax-$tolerance;
}

$middle = ($zmincutoff+$zmaxcutoff)/2;
for($t=0;$t<$field_nummols;$t++) {
  for($m=0;$m<$mol_numents[$t];$m++) {
    for($a=0;$a<$mol_numatoms[$t];$a++) {
      if(uc($cdata[$t][$m][$a][0]) eq 'OA') {
	if($cdata[$t][$m][$a][3]>$zmaxcutoff or $cdata[$t][$m][$a][3]<$zmincutoff) {
	  $mol_numents[$t]--;
	  if(($cdata[$t][$m][$a][3]>$zmaxcutoff and $cdata[$t][$m][$a][3]<$zmaxcutoff+3*$tolerance)
	  or ($cdata[$t][$m][$a][3]<$zmincutoff and $cdata[$t][$m][$a][3]>$zmincutoff-3*$tolerance)) {
	    @{$cdata[$foundhydrox][$mol_numents[$foundhydrox]][0]} = @{$cdata[$t][$m][$a]};
	    @{$cdata[$foundhydrox][$mol_numents[$foundhydrox]][1]} = @{$cdata[$t][$m][$a]};
	    $cdata[$foundhydrox][$mol_numents[$foundhydrox]][0][0] = 'OX';
	    $cdata[$foundhydrox][$mol_numents[$foundhydrox]][1][0] = 'HG';
	    $cdata[$foundhydrox][$mol_numents[$foundhydrox]][1][3] += $doh * sgn($cdata[$t][$m][$a][3]-$middle);
	    $mol_numents[$foundhydrox]++;
	  }
	  splice(@{$cdata[$t]},$m,1);
	  $m--;
	}
      } elsif(uc($cdata[$t][$m][$a][0]) eq 'AL') {
	if($cdata[$t][$m][$a][3]>$zmaxcutoff or $cdata[$t][$m][$a][3]<$zmincutoff) {
	  $mol_numents[$t]--;
	  splice(@{$cdata[$t]},$m,1);
	  $m--;
	}
      }
    }
  }
}
############# write output CONFIG #######################################

open(OUTPUT_CONFIG, ">", $ARGV[2]) or die "Can't create output CONFIG-file:\n$!";
open(OUTPUT_FIELD, ">", $ARGV[3]) or die "Can't create output FIELD-file:\n$!";

print OUTPUT_CONFIG $title;
for $a (@config_header) {
  print OUTPUT_CONFIG $a;
}
$index=1;
for($t=0;$t<$field_nummols;$t++) {
  for($m=0;$m<$mol_numents[$t];$m++) {
    for($a=0;$a<$mol_numatoms[$t];$a++) {
      printf OUTPUT_CONFIG "%-8s" ,$cdata[$t][$m][$a][0];
      printf OUTPUT_CONFIG "%10u\n",$index; #atom index
      printf OUTPUT_CONFIG " %19.11g %19.11g %19.11g\n",$cdata[$t][$m][$a][1],$cdata[$t][$m][$a][2],$cdata[$t][$m][$a][3];
      if($config_key>0) {
	printf OUTPUT_CONFIG " %19.11g %19.11g %19.11g\n",$cdata[$t][$m][$a][4],$cdata[$t][$m][$a][5],$cdata[$t][$m][$a][6];
      }
      if($config_key>1) {
	printf OUTPUT_CONFIG " %19.11g %19.11g %19.11g\n",$cdata[$t][$m][$a][7],$cdata[$t][$m][$a][8],$cdata[$t][$m][$a][9];
      }
      $index++;
    }
  }
}

############# write output FIELD ########################################

open(FIELD, "<", $ARGV[1]) or die "Can't open FIELD-file:\n$!";
$t=0;
while(<FIELD>) {
  if(/NUMMOLS/) {
    print OUTPUT_FIELD "NUMMOLS $mol_numents[$t]\n";
    $t++;
  } else {
    print OUTPUT_FIELD $_;    
  }
}
close(FIELD, OUTPUT_FIELD);
