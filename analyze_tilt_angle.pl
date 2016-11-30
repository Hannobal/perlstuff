#!/usr/bin/perl

# calculates the tilt angles, dihedral angles and chain lengths of carbon backbone of
# certain phosponic acids from DL_POLY HISTORY files and their develeopment in time

use Math::Trig;
use Cwd;
use Switch;
use hanno_utility;
use dlpoly_utility;
use aloxsam_utility;
use Storable qw(dclone);
@includetypes = ();

use Math::Trig;

$histrestiltdeg = 2;
$histresdihdeg  = 5;
$histreslen     = 2;
$startframe = 0;
$endframe   = 9e20;
$histname   = "HISTORY";
$fldname    = "FIELD";
$outdir     = $ARGV[0];

if($#ARGV<0) {
  print "input format:\n";
  print " 1. output directory for the output files\n";
  print "-rt <real>   histogram resolution for tilt angle (default: $histrestiltdeg deg)]\n";
  print "-rd <real>   histogram resolution for dihedrals (default: $histresdihdeg deg)]\n";
  print "-rl <real>   histogram resolution for chain lengths (default: $histreslen A)]\n";
  print "-t n*<str>   list of molecules to analyze\n";
  print "-s <int>     start frame\n";
  print "-e <int>     end frame\n";
  print "-lessoutput  do not print current frame\n";
  exit;
}

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    # general flags
    case (/^-rd/i) {
      $i++;
      $histresdihdeg = $ARGV[$i];
      &error(1,"dihedral histogram resolution must be >=0") if($histresdihdeg<=0);
    } case (/^-rl/i) {
      $i++;
      $histreslen = $ARGV[$i];
      &error(1,"chain length histogram resolution must be >=0")if($histreslen<=0);
    } case (/^-rt/i) {
      $i++;
      $histrestiltdeg = $ARGV[$i];
      &error(1,"tilt angle histogram resolution must be >=0") if($histrestiltdeg<=0);
    } case(/^-s/i) {
      $i++;
      $startframe = $ARGV[$i];
      &error(1,"start frame must be an integer") if not check_integer($startframe);
    } case(/^-e/i) {
      $i++;
      $endframe = $ARGV[$i];
      &error(1,"end frame must be an integer") if not check_integer($endframe);
    } case(/^-t/i) {
      $i++;
      until($ARGV[$i]=~/^-/ or $i>$#ARGV) {
	push(@includetypes,$ARGV[$i]);
	$i++;
      }
      $i--;
    } case(/^-lessoutput/i) {
      $lessoutput=1;
    } else {
      &error(10);
    }
  }
}

$histrestiltrad = $histrestiltdeg/180*pi;
$histresdihrad  = $histresdihdeg/180*pi;
@directories=get_numeric_directories(".");
if(@directories) {
  @directories=get_numeric_directories(".",$startframe,$endframe);
}
push(@directories,".") if(-e "./$histname");
&error(3,"no suitable directories found") if(not @directories);

exit 1 if(read_field_file("$directories[0]/$fldname",0)!=0);
@tpas=();
$numpas=0;
$numdihetot=0;
for($t=0;$t<$field_nummols[0];$t++) {
  next if(not $mol_name[0][$t]=~/-(PA|CA)-\d?d/);
  next if($mol_numents[0][$t]==0);
  if(@includetypes) {
    next if(not contains(@includetypes,$mol_name[0][$t]));
  }
  push(@tpas,$t);
  $numpas+=$mol_numents[0][$t];
  # remove dihedrals with atoms other than C3
  for($d=0;$d<$mol_numdihedrals[0][$t];$d++) {
      for($i=1;$i<5;$i++) {
	if(not $mol_atomdata[0][$t][$mol_dihedraldata[0][$t][$d][$i]][0]=~/^(C3|OS)/) {
	splice(@{$mol_dihedraldata[0][$t]},$d,1);
	$mol_numdihedrals[0][$t]--;
	$d--;
	last;
      }
    }
  }
  # remove doublets
  for($d1=0;$d1<$mol_numdihedrals[0][$t];$d1++) {
    remdih:for($d2=$d1+1;$d2<$mol_numdihedrals[0][$t];$d2++) {
      for($i=1;$i<5;$i++) {
	if($mol_dihedraldata[0][$t][$d1][$i]!=$mol_dihedraldata[0][$t][$d2][$i]){
	  next remdih;
	}
      }
      splice(@{$mol_dihedraldata[0][$t]},$d2,1);
      $mol_numdihedrals[0][$t]--;
      $d2--;
    }
  }
  $numdihetot+=$mol_numdihedrals[0][$t]*$mol_numents[0][$t];
  &error(4,"did not find index for P5 atom in molecule $mol_name[0][$t]") if(not defined $refcindex{$mol_name[0][$t]});
  &error(4,"did not find index for terminal C atom in molecule $mol_name[0][$t]") if(not defined $refp5index{$mol_name[0][$t]});
}
&error(4,"no molecules to analyze") if(not @tpas);
$tpamax=$#tpas+1;
if(@tpas>1) {
  push(@{$mol_name[0]},"all");
  push(@tpas,$#{$mol_name[0]});
  $iall=$#tpas;
}

mkdir($outdir) if(not -d $outdir);
for($i=0;$i<@tpas;$i++) {
  $t=$tpas[$i];
  &error(4) if(not open($fhavdihedral[$i],">","$outdir/AV_CCCC_DIHEDRAL_".$mol_name[0][$t]));
  printf {$fhavdihedral[$i]} "%10s %17s\n","#time [ps]","average dihedral angle [deg]";
  &error(4) if(not open($fhavheadtail[$i],">","$outdir/AV_head-tail_".$mol_name[0][$t]));
  printf {$fhavheadtail[$i]} "%10s %17s\n","#time [ps]","average head-tail-distance [A]";
  &error(4) if(not open($fhavheadtail_z[$i],">","$outdir/AV_head-tail_z_".$mol_name[0][$t]));
  printf {$fhavheadtail_z[$i]} "%10s %17s\n","#time [ps]","av. z-component of h-t-dist [A].";
  &error(4) if(not open($fhavheadtail_xy[$i],">","$outdir/AV_head-tail_xy_".$mol_name[0][$t]));
  printf {$fhavheadtail_xy[$i]} "%10s %17s\n","#time [ps]","av. xy-component of h-t-dist [A].";
  &error(4) if(not open($fhavtilt[$i],">","$outdir/AV_tiltangle_".$mol_name[0][$t]));
  printf {$fhavheadtail_xy[$i]} "%10s %17s\n","#time [ps]","av. tilt angle [deg].";
}

$firsthist=1;
$f=0;
if($lessoutput) {
  print "begin analysis...\n";
}
foreach $dir (@directories) {
  &error(2,"$dir/$histname") if(not open($fhhist,"<","$dir/$histname"));
  if($startframe>0 and $firsthist) {
    &error(5,$startframe) if(not find_history_timestep($fhhist,$startframe));
  }
  $firsthist=0;
  $err=0;
  while($err==0) {
    $err=read_history_timestep($fhhist,0,0);
    last if($err!=0);
    last if($frame_number[0]>$endframe);
    print "\ranalyzing timestep $frame_number[0]" unless($lessoutput);
    $avtiltall  = 0;
    $avdiheall  = 0;
    $avlenall   = 0;
    $avlenzall  = 0;
    $avlenxyall = 0;
    @tiltangles = ();
    @dihedrals  = ();
    for($i=0;$i<$tpamax;$i++) {
      $t=$tpas[$i];
      $j=0;
      $avdihedral[$i] = 0;
      $avlen[$i]      = 0;
      $avlenz[$i]     = 0;
      $avelnxy[$i]    = 0;
      $avtilt[$i]     = 0;
      $refc=$refCindex{$mol_name[0][$t]};
      $refp=$refp5index{$mol_name[0][$t]};
      $histdih   = $i*5;
      $histlen   = $i*5+1;
      $histlenz  = $i*5+2;
      $histlenxy = $i*5+3;
      $histtilt  = $i*5+4;
      for($m=0;$m<$mol_numents[0][$t];$m++) {
	remap_molecule(\@{$cdata[0][$t][$m]},[0,1],\@{$size[0]},10);
	# analyze didehdral angles
	for($d=0;$d<$mol_numdihedrals[0][$t];$d++) {
	  ($a1,$a2,$a3,$a4)=@{$mol_dihedraldata[0][$t][$d]}[1..4];
	  $dx  = ($cdata[0][$t][$m][$a2][0]+$cdata[0][$t][$m][$a3][0])/2;
	  $dy  = ($cdata[0][$t][$m][$a2][1]+$cdata[0][$t][$m][$a3][1])/2;
	  $dxy = sqrt($dx*$dx+$dy*$dy);
	  $dihedral[$i][$j]=abs(calc_dihedral_angle(0,0,$d,$t,$m));
	  $avdihedral[$i] += $dihedral[$i][$j];
	  histogram_add_one($histdih,$dihedral[$i][$j],$histresdihrad);
	  $j++;
	}
	# analyze head-tail-distance
	$dx = $cdata[0][$t][$m][$refc][0]-$cdata[0][$t][$m][$refp][0];
	$dy = $cdata[0][$t][$m][$refc][1]-$cdata[0][$t][$m][$refp][1];
	$dz[$i][$m]   = $cdata[0][$t][$m][$refc][2]-$cdata[0][$t][$m][$refp][2];
	$dxy[$i][$m]  = sqrt($dx*$dx+$dy*$dy);
	$dxyz[$i][$m] = sqrt($dx*$dx+$dy*$dy+$dz[$i][$m]*$dz[$i][$m]);
	histogram_add_one($histlen,$dxyz[$i][$m],$histreslen);
	histogram_add_one($histlenz,$dz[$i][$m],$histreslen);
	histogram_add_one($histlenxy,$dxy[$i][$m],$histreslen);
	$avlen[$i]   += $dxyz[$i][$m];
	$avlenz[$i]  += $dz[$i][$m];
	$avlenxy[$i] += $dxy[$i][$m];
	# analyze tilt angle
	$angle[$i][$m] = atan2($dxy[$i][$m],$dz[$i][$m]);
	histogram_add_one($histtilt,$angle[$i][$m],$histrestiltrad);
	$avtilt[$i] += $angle[$i][$m];
      }
      if(@tpas>1) {
	$avtiltall  += $avtilt[$i];
	$avdiheall  += $avdihedral[$i];
	$avlenall   += $avlen[$i];
	$avlenzall  += $avlenz[$i];
	$avlenxyall += $avlenxy[$i];
      }
      $avdihedral[$i] /= $mol_numents[0][$t]*$mol_numdihedrals[0][$t]*pi/180;
      $avlen[$i]      /= $mol_numents[0][$t];
      $avlenz[$i]     /= $mol_numents[0][$t];
      $avlenxy[$i]    /= $mol_numents[0][$t];
      $avtilt[$i]     /= $mol_numents[0][$t]*pi/180;
      printf {$fhavdihedral[$i]}    "%10u %17.10e\n",$frame_number[0],$avdihedral[$i];
      printf {$fhavheadtail[$i]}    "%10u %17.10e\n",$frame_number[0],$avlen[$i];
      printf {$fhavheadtail_z[$i]}  "%10u %17.10e\n",$frame_number[0],$avlenz[$i];
      printf {$fhavheadtail_xy[$i]} "%10u %17.10e\n",$frame_number[0],$avlenxy[$i];
      printf {$fhavtilt[$i]}        "%10u %17.10e\n",$frame_number[0],$avtilt[$i];
    }
    if(@tpas>1) {
      $avdiheall  /= $numdihetot*pi/180;
      $avlenall   /= $numpas;
      $avlenzall  /= $numpas;
      $avlenxyall /= $numpas;
      $avtiltall  /= $numpas;
      printf {$fhavdihedral[$iall]}    "%10u %17.10e\n",$frame_number[0],$avdiheall;
      printf {$fhavheadtail[$iall]}    "%10u %17.10e\n",$frame_number[0],$avlenall;
      printf {$fhavheadtail_z[$iall]}  "%10u %17.10e\n",$frame_number[0],$avlenzall;
      printf {$fhavheadtail_xy[$iall]} "%10u %17.10e\n",$frame_number[0],$avlenxyall;
      printf {$fhavtilt[$iall]} "%10u %17.10e\n",$frame_number[0],$avtiltall;
    }
    $f++;
  }
  close($fhhist);
}

if($lessoutput) {
  print "last frame analyzed: $frame_number[0]\n";
} else {
  print "\n";
}

for($i=0;$i<@tpas;$i++) {
  close($fhavdihedral[$i],$fhavheadtail[$i],$fhavheadtail_z[$i],$fhavheadtail_xy[$i]);
}
for($i=0;$i<$tpamax;$i++) {
  $t=$tpas[$i];
  $histdih   = $i*5;
  $histlen   = $i*5+1;
  $histlenz  = $i*5+2;
  $histlenxy = $i*5+3;
  $histtilt  = $i*5+4;
  &histogram_normalize_integral($histdih,$histresdihdeg);
  &histogram_normalize_integral($histlen,$histreslen);
  &histogram_normalize_integral($histlenz,$histreslen);
  &histogram_normalize_integral($histlenxy,$histreslen);
  &histogram_normalize_integral($histtilt,$histrestiltdeg);
  
  open($fhhist,">","$outdir/HIST_DIHEDRAL_".$mol_name[0][$t]);
  print $fhhist "# Histogram of CCCC dihedral angles [deg] normalized to integral=1\n";
  for($j=$histminindex[$histdih];$j<=$histmaxindex[$histdih];$j++) {
    printf $fhhist "%17.10g %17.10g %10u\n",$histresdihdeg*$j,$histogram[$histdih][$j],$histnum[$histdih][$j];
  }
  close($fhhist);
  
  open($fhhist,">","$outdir/HIST_H-T-DIST_".$mol_name[0][$t]);
  print $fhhist "# Histogram of head to tail length [A] normalized to integral=1\n";
  for($j=$histminindex[$histlen];$j<=$histmaxindex[$histlen];$j++) {
    printf $fhhist "%17.10g %17.10g %10u\n",$histreslen*$j,$histogram[$histlen][$j],$histnum[$histlen][$j];
  }
  close($fhhist);
  
  open($fhhist,">","$outdir/HIST_H-T-DIST_Z_".$mol_name[0][$t]);
  print $fhhist "# Histogram of z-component of head to tail length [A] normalized to integral=1\n";
  for($j=$histminindex[$histlenz];$j<=$histmaxindex[$histlenz];$j++) {
    printf $fhhist "%17.10g %17.10g %10u\n",$histreslen*$j,$histogram[$histlenz][$j],$histnum[$histlenz][$j];
  }
  close($fhhist);
  
  open($fhhist,">","$outdir/HIST_H-T-DIST_XY_".$mol_name[0][$t]);
  print $fhhist "# Histogram of xy-component of head to tail length [A] normalized to integral=1\n";
  for($j=$histminindex[$histlenxy];$j<=$histmaxindex[$histlenxy];$j++) {
    printf $fhhist "%17.10g %17.10g %10u\n",$histreslen*$j,$histogram[$histlenxy][$j],$histnum[$histlenxy][$j];
  }
  close($fhhist);
  
  open($fhhist,">","$outdir/HIST_TILT_".$mol_name[0][$t]);
  print $fhhist "# Histogram of tilt angle [deg] normalized to integral=1\n";
  for($j=$histminindex[$histtilt];$j<=$histmaxindex[$histtilt];$j++) {
    printf $fhhist "%17.10g %17.10g %10u\n",$histrestiltdeg*$j,$histogram[$histtilt][$j],$histnum[$histtilt][$j];
  }
  close($fhhist);
}
$numframes=$f+1;



sub error {
  if($_[0]==1) {
    print "**** error: $_[1] must be an integer number!\n";
  } elsif($_[0]==2) {
    print "**** error: could not open file $_[1]!\n";
  } elsif($_[0]==3) {
    print "**** error: carbon atom has $_[1] hydrogen atoms bound!\n";
  } elsif($_[0]==4) {
    print "**** error: could not open output file!\n";
  } elsif($_[0]==5) {
    print "**** error: did not find timestep $_[1]!\n";
  } elsif($_[0]==6) {
    print "**** error: histogram resolution for $_[1] must be greater than zero!\n";
  } elsif($_[0]==7) {
    print "**** error: gauche cutoff must be greater than zero!\n";
  } elsif($_[0]==8) {
    print "**** error: no suitable directories found!\n";
  } elsif($_[0]==9) {
    print "**** error: reference carbon/phosphorus atoms are not defined for molecule $_[1]!\n";
  } elsif($_[0]==10) {
    print "**** error: unknown flag $_[1]!\n";
  } else {
    print "**** unknown error!\n";
    exit 99;
  }
  print "**** error: $_[1]!\n" if(defined $_[1]);
  exit $_[0]
}