#!/usr/bin/perl

$|=1;

#calculate the average position or center of mass of atoms over a trajectory

use hanno_utility;
use dlpoly_utility;
use Switch;


if($#ARGV<0) {
  print "input format:\n";
  print " output directory for the output files\n";
  print "-t n*<str>   list of molecules to analyze\n";
  print "-s <int>     start frame\n";
  print "-e <int>     end frame\n";
  print "-m           use center of mass instead of geometric average\n";
  exit;
}

$startframe    = 0;
$endframe      = 9e20;
@includetypes  = ();
$lincludetypes = 0;
$lcentermass   = 0;

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case(/^-s/i) {
      $i++;
      $startframe = $ARGV[$i];
      &error(1) if not check_integer($startframe);
    } case(/^-e/i) {
      $i++;
      $endframe = $ARGV[$i];
      &error(2) if not check_integer($endframe);
    } case(/^-t/i) {
      $i++;
      $lincludetypes=1;
      until($ARGV[$i]=~/^-/) {
	push(@includetypes,$ARGV[$i]);
	$i++;
      }
      $i--;
    } case(/^-m/i) {
      $lcentermass=1;
    } else {
      &error(10,$ARGV[$i]);
    }
  }
}

$histname = "HISTORY";
$fldname  = "FIELD";
$outdir   = $ARGV[0];

@directories=get_numeric_directories(".");
if(@directories) {
  @directories=get_numeric_directories(".",$startframe,$endframe);
}
push(@directories,".") if(-e "./$histname");
&error(8) if(not @directories);

exit 1 if(read_field_file("$directories[0]/$fldname",0)!=0);

$masstot=0;
$numatomstot=0;
if(not $lincludetypes) {
  for($t=0;$t<$field_nummols[0];$t++) {
    $masstot     += $mol_mass[0][$t]*$mol_numents[0][$t];
    $numatomstot += $mol_numents[0][$t]*$mol_numatoms[0][$t];
  }
  @includetypes = @{$mol_name[0]};
  @tlst         = (0..$#{$mol_name[0]});
} else {
  for($i=0;$i<@includetypes;$i++) {
    $tlst[$i] = -1;
    for($t=0;$t<$field_nummols[0];$t++) {
      if($mol_name[0][$t]=~/^$includetypes[$i]$/i) {
	$tlst[$i] = $t;
	$masstot+=$mol_mass[0][$t]*$mol_numents[0][$t];
	$numatomstot+=$mol_numents[0][$t]*$mol_numatoms[0][$t];
	last;
      }
    }
    &error(11,$includetypes[$i]) if($tlst[$i]==-1);
  }
}

mkdir($outdir) if(not -d $outdir);
if($lcentermass) {
  &error(4) if(not open($fhavpos,">","$outdir/AV_MASSCENTER"));
  print $fhavpos "# center of mass ";
} else {
  &error(4) if(not open($fhavpos,">","$outdir/AV_POSITION"));
  print $fhavpos "# average position ";
}
if($lincludetypes) {
  print $fhavpos "using molecules ",join(" ",@includetypes),"\n"
} else {
  print $fhavpos "using all molecules\n";
}

$firsthist=1;
$f=0;
foreach $dir (@directories) {
  exit 1 if(not open($fhhist,"<","$dir/$histname"));
  if($startframe>0 and $firsthist) {
    &error(5,$startframe) if(not find_history_timestep($fhhist,$startframe));
  }
  $firsthist=0;
  $err=0;
  print "\ranalyzing directory $dir";
  while($err==0) {
    $err=read_history_timestep($fhhist,0,0);
    last   if($err<0);
    exit 1 if($err>0);
    last if($frame_number[0]>$endframe);
    @center=(0,0,0);
    if($lcentermass) {
      foreach $t (@tlst) {
	for($m=0;$m<$mol_numents[0][$t];$m++) {
	  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	    $center[0] += $cdata[0][$t][$m][$a][0]*$mol_atomdata[0][$t][$a][1];
	    $center[1] += $cdata[0][$t][$m][$a][1]*$mol_atomdata[0][$t][$a][1];
	    $center[2] += $cdata[0][$t][$m][$a][2]*$mol_atomdata[0][$t][$a][1];
	  }
	}
      }
      $center[0] /= $masstot;
      $center[1] /= $masstot;
      $center[2] /= $masstot;
    } else {
      foreach $t (@tlst) {
	for($m=0;$m<$mol_numents[0][$t];$m++) {
	  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	    $center[0] += $cdata[0][$t][$m][$a][0];
	    $center[1] += $cdata[0][$t][$m][$a][1];
	    $center[2] += $cdata[0][$t][$m][$a][2];
	  }
	}
      }
      $center[0] /= $numatomstot;
      $center[1] /= $numatomstot;
      $center[2] /= $numatomstot;
    }
    if($frame_number[0]<$frame_number[1]) {
      $frame_offset+=$frame_number[1];
    }
    $frame_number[1]=$frame_number[0];
    $frame_number[0]+=$frame_offset;
    printf $fhavpos "%10u %17.10g %17.10g %17.10g\n",$frame_number[0],@center;
  }
}
print "\n";
    


sub error {
  if($_[0]==1) {
    print "**** error: start frame must be an integer number!\n";
  } elsif($_[0]==2) {
    print "**** error: end frame must be an integer number!\n";
  } elsif($_[0]==3) {
    print "**** error: carbon atom has $_[1] hydrogen atoms bound!\n";
  } elsif($_[0]==4) {
    print "**** error: could not open output file!\n";
  } elsif($_[0]==5) {
    print "**** error: did not find timestep $_[1]!\n";
  } elsif($_[0]==6) {
    print "**** error: histogram resolution must be greater than zero!\n";
  } elsif($_[0]==7) {
    print "**** error: gauche cutoff must be greater than zero!\n";
  } elsif($_[0]==8) {
    print "**** error: no suitable directories found!\n";
  } elsif($_[0]==9) {
    print "**** error: angular resolution must be divisor of 360!\n";
  } elsif($_[0]==10) {
    print "**** error: unknown flag $_[1]!\n";
  } elsif($_[0]==11) {
    print "**** error: molecule $_[1] was not found in FIELD data!\n";
  } else {
    print "**** unknown error!\n";
    exit 99;
  }
  exit $_[0]
}