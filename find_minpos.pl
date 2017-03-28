#!/usr/bin/perl

use Switch;
use dlpoly_utility;
use hanno_utility;

if($#ARGV<1) {
  print "input format:\n";
  print " 1. name of CONFIG/HISTORY file\n";
  print " 2. name of corresponding FIELD file\n";
  print "[3. mode: \"mol\"/\"atoms\"\n]";
  print "[4. list of atom types/molecule names]\n";
  exit;
}

$cfgfilename = $ARGV[0];
$fldfilename = $ARGV[1];
exit 1 if(read_field_file($fldfilename,0) != 0);

if($#ARGV>1) {
  @types   = @ARGV[3..$#ARGV];
  if($ARGV[2] =~ /^atom/i) {
    $mode = 0;
    @types = @{$field_atomtypes[0]} if($#types==-1);
  } elsif($ARGV[2] =~ /^mol/i) {
    $mode = 1;
    @types = @{$mol_name[0]} if($#types==-1);
  } else {
    print "**** error: mode must be \"mol\" or \"atom\"\n";
    exit 1;
  }
} else {
  $mode = 1;
  @types = @{$mol_name[0]}
}

if($cfgfilename =~ /.cfg/ or $cfgfilename =~ /CONFIG/ or $cfgfilename =~ /REVCON/) {
  exit 1 if(read_config_file($cfgfilename,0,0) != 0);
}elsif($cfgfilename =~ /.dlpolyhist/ or $cfgfilename =~ /HISTORY/) {
  open($fhhist, '<', $cfgfilename) or die "**** error: can't open $cfgfilename:\n$!\n";
  exit 1 if(read_history_timestep($fhhist,0,0) != 0);
  close($fhhist);
} else {
  print "**** error: file format not recognized!\n";
}

if($mode==0) {
  for($a=0;$a<$field_numatomtypes[0];$a++) {
    @{$minpos{$field_atomtypes[0][$a]}} = ( 9e20, 9e20, 9e20);
    @{$maxpos{$field_atomtypes[0][$a]}} = (-9e20,-9e20,-9e20);
  }
  for($t=0;$t<@{$cdata[0]};$t++) {
    for($m=0;$m<@{$cdata[0][$t]};$m++) {
      for($a=0;$a<@{$cdata[0][$t][$m]};$a++) {
	for($c=0;$c<3;$c++) {
	  if($cdata[0][$t][$m][$a][$c] < $minpos{$cdata[0][$t][$m][$a][9]}[$c]) {
	    $minpos{$cdata[0][$t][$m][$a][9]}[$c] = $cdata[0][$t][$m][$a][$c];
	  }
	  if($cdata[0][$t][$m][$a][$c] > $maxpos{$cdata[0][$t][$m][$a][9]}[$c]) {
	    $maxpos{$cdata[0][$t][$m][$a][9]}[$c] = $cdata[0][$t][$m][$a][$c];
	  }
	} # end for c
      } # end for a
    } # end for m
  } # end for t

  @minpos_tot = ( 9e20, 9e20, 9e20);
  @maxpos_tot = (-9e20,-9e20,-9e20);
  printf "%-7s %10s %10s %10s %10s %10s %10s\n", 'name','min x','min y','min z','max x','max y','max z';
  for($a=0;$a<$field_numatomtypes[0];$a++) {
    for($b=0;$b<@types;$b++) {
      if($field_atomtypes[0][$a] eq $types[$b]) {
	if($minpos{$field_atomtypes[0][$a]}[0]==9e20) {
  # 	printf "%-7s", $field_atomtypes[0][$a];
  # 	printf " %14s %14s %14s",   '----------','----------','----------';
  # 	printf " %14s %14s %14s\n", '----------','----------','----------';
	} else {
	  printf "%-7s", $field_atomtypes[0][$a];
	  printf " %10.5f %10.5f %10.5f",   @{$minpos{$field_atomtypes[0][$a]}};
	  printf " %10.5f %10.5f %10.5f\n", @{$maxpos{$field_atomtypes[0][$a]}};
	  for($c=0;$c<3;$c++) {
	    $minpos_tot[$c] = $minpos{$field_atomtypes[0][$a]}[$c] if($minpos{$field_atomtypes[0][$a]}[$c]<$minpos_tot[$c]);
	    $maxpos_tot[$c] = $maxpos{$field_atomtypes[0][$a]}[$c] if($maxpos{$field_atomtypes[0][$a]}[$c]>$maxpos_tot[$c]);
	  }
	}
	last;
      }
    }
  }
  printf "%-7s", 'total';
  printf " %10.5f %10.5f %10.5f",   @minpos_tot;
  printf " %10.5f %10.5f %10.5f\n", @maxpos_tot;
} else {
  for($t=0;$t<$field_nummols[0];$t++) {
    $cycle=1;
    for($u=0;$u<@types;$u++) {
      if($mol_name[0][$t] eq $types[$u]) {
	$cycle=0; last;
      }
    }
    next if($cycle);
    if($mol_numents[0][$t]==0) {
      splice(@types,$u,1); next;
    };
    @{$minpos{$mol_name[0][$t]}} = ( 9e20, 9e20, 9e20);
    @{$maxpos{$mol_name[0][$t]}} = (-9e20,-9e20,-9e20);
    for($m=0;$m<@{$cdata[0][$t]};$m++) {
      for($a=0;$a<@{$cdata[0][$t][$m]};$a++) {
	for($c=0;$c<3;$c++) {
	  if($cdata[0][$t][$m][$a][$c] < $minpos{$mol_name[0][$t]}[$c]) {
	    $minpos{$mol_name[0][$t]}[$c] = $cdata[0][$t][$m][$a][$c];
	  }
	  if($cdata[0][$t][$m][$a][$c] > $maxpos{$mol_name[0][$t]}[$c]) {
	    $maxpos{$mol_name[0][$t]}[$c] = $cdata[0][$t][$m][$a][$c];
	  }
	} # end for c
      } # end for a
    } # end for m
  } # end for t
  @minpos_tot = ( 9e20, 9e20, 9e20);
  @maxpos_tot = (-9e20,-9e20,-9e20);
  printf "%-15s %10s %10s %10s %10s %10s %10s\n", 'name','min x','min y','min z','max x','max y','max z';
  for($u=0;$u<@types;$u++) {
    printf "%-15s", $types[$u];
    printf " %10.5f %10.5f %10.5f",   @{$minpos{$types[$u]}};
    printf " %10.5f %10.5f %10.5f\n", @{$maxpos{$types[$u]}};
    for($c=0;$c<3;$c++) {
      $minpos_tot[$c] = $minpos{$types[$u]}[$c] if($minpos{$types[$u]}[$c]<$minpos_tot[$c]);
      $maxpos_tot[$c] = $maxpos{$types[$u]}[$c] if($maxpos{$types[$u]}[$c]>$maxpos_tot[$c]);
    }
  }
  if($#types>0) {
    printf "%-15s", 'total';
    printf " %10.5f %10.5f %10.5f",   @minpos_tot;
    printf " %10.5f %10.5f %10.5f\n", @maxpos_tot;
  }
}
