#!/usr/bin/perl

# calculates the tilt angles (to the surface normal => 0 deg for perpendicular
# molecules) of certain phosponic acids from DL_POLY HISTORY files and their
# develeopment in time for force_docking_z-type systems

use Math::Trig;
use Cwd;
use Switch;
use hanno_utility;
use dlpoly_utility;
use aloxsam_utility;

$histrestiltdeg  = 0.5;
$histresdihdeg   = 0.5;
$histreslen      = 0.5;
$gauchecutoffdeg = 30;
$startframe = 0;
$endframe   = 9e20;
$histname   = "HISTORY";
$fldname    = "FIELD";
$outdir     = $ARGV[0];
$subdir     = "run";
$lnumsubdir = 0;

if($#ARGV<0) {
  print "input format:\n";
  print " 1. <output-directory>\n";
  print "-gc <real>   cutoff for gauche defects (default: $gauchecutoffdeg deg)\n";
  print "-rt <real>   histogram resolution for tilt angle (default: $histrestiltdeg deg)]\n";
  print "-rd <real>   histogram resolution for dihedrals (default: $histresdihdeg deg)]\n";
  print "-rl <real>   histogram resolution for chain lengths (default: $histreslen A)]\n";
  print "-d  <str>    subdirectory to analyze (default: $subdir, \"num\" for last numeric directory)\n";
  print "-t n*<str>   list of molecules to analyze\n";
  print "-f <2*int>   start and end frame\n";
  print "-s <2*int>   start and end step\n";
  print "-lessoutput  do not print current frame\n";
  exit 1;
} 

$outdir = $ARGV[0];

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    # general flags
    case (/^-gc/i) {
      $i++;
      $gauchecutoffdeg = $ARGV[$i];
       &error(1,"gauche cutoff must be >=0") if($gauchecutoffdeg<=0);
    } case (/^-rd/i) {
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
    } case(/^-f/i) {
      $i++;
      $startframe = $ARGV[$i];
      &error(1,"start frame must be an integer") if not check_integer($startframe);
      $i++;
      $endframe = $ARGV[$i];
      &error(1,"end frame must be an integer") if not check_integer($endframe);
    } case(/^-s/i) {
      $i++;
      $startstep = $ARGV[$i];
      &error(1,"start step must be an integer") if not check_integer($startstep);
      $i++;
      $endstep = $ARGV[$i];
      &error(1,"end step must be an integer") if not check_integer($endstep);
    } case(/^-d/i) {
      $i++;
      $subdir=$ARGV[$i];
      if($subdir eq "num") {
	$lnumsubdir=1;
      } else {
	$lnumsubdir=0;
      }
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
      &error(2,"could not interpret flag $ARGV[$i]");
    }
  }
}

&error(1,"end frame must not be smaller than start frame") if($endframe<$startframe);
$histrestiltrad  = $histrestiltdeg/180*pi;
$histresdihrad   = $histresdihdeg/180*pi;
$gauchecutoffrad = $gauchecutoffdeg/180*pi;

############### find the relevant directories to read #########################
print "reading directories\n";
@directories=get_numeric_directories(".",$startstep,$endstep);
&error(3,"no suitable directories found") if(not @directories);
if($lnumsubdir) {
  @numdirs = get_numeric_directories($directories[0]);
  &error(3,"no numeric subdirectories found") if(@numdirs);
  if(-d "$directories[0]/min") {
    exit 1 if(read_field_file("$directories[0]/$numdirs[$#numdirs]/min/FIELD",0)!=0);
  } else {
    exit 1 if(read_field_file("$directories[0]/$numdirs[$#numdirs]/FIELD",0)!=0);
  }
} else {
  &error(3,"directory $directories[0]/$subdir does not exist") if(not -d "$directories[0]/$subdir");
  @numdirs = get_numeric_directories("$directories[0]/$subdir");
  if(@numdirs) {
    exit 1 if(read_field_file("$directories[0]/$subdir/$numdirs[0]/FIELD",0)!=0);
  } else {
    exit 1 if(read_field_file("$directories[0]/$subdir/FIELD",0)!=0);
  }
}

############### analyze the FIELD file ########################################

@tpas=();
for($t=0;$t<$field_nummols[0];$t++) {
  next if(not $mol_name[0][$t]=~/-(PA|CA)-\d?d/);
  next if($mol_numents[0][$t]==0);
  if(@includetypes) {
    next if(not contains(@includetypes,$mol_name[0][$t]));
  }
  push(@tpas,$t);
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
  &error(4,"did not find index for P5 atom in molecule $mol_name[0][$t]") if(not defined $refcindex{$mol_name[0][$t]});
  &error(4,"did not find index for terminal C atom in molecule $mol_name[0][$t]") if(not defined $refp5index{$mol_name[0][$t]});
}

&error(4,"no molecules to analyze") if(not @tpas);
$tpamax=$#tpas+1;
if(@tpas>1) {
  push(@{$mol_name[0]},"all");
  push(@tpas,$#{$mol_name[0]});
  $iall=$#tpas;
  $histdihall    = $iall*6;
  $histlenall    = $iall*6+1;
  $histlenzall   = $iall*6+2;
  $histlenxyall  = $iall*6+3;
  $histtiltall   = $iall*6+4;
  $histgaucheall = $iall*6+5;
}

############### read the HISTORY file from each directory #####################

mkdir($outdir) if(not -d $outdir);
for($i=0;$i<@tpas;$i++) {
  $t=$tpas[$i];
  &error(4) if(not open($fhavdihedral[$i],">","$outdir/AV_CCCC_DIHEDRAL_".$mol_name[0][$t]));
  printf {$fhavdihedral[$i]} "%10s %17s %17s\n","#time [ps]","av dih ang [deg]","stdev";
  &error(4) if(not open($fhavheadtail[$i],">","$outdir/AV_head-tail_".$mol_name[0][$t]));
  printf {$fhavheadtail[$i]} "%10s %17s %17s\n","#time [ps]","av head-tail [A]","stdev";
  &error(4) if(not open($fhavheadtail_z[$i],">","$outdir/AV_head-tail_z_".$mol_name[0][$t]));
  printf {$fhavheadtail_z[$i]} "%10s %17s %17s\n","#time [ps]","av z of h-t [A]","stdev";
  &error(4) if(not open($fhavheadtail_xy[$i],">","$outdir/AV_head-tail_xy_".$mol_name[0][$t]));
  printf {$fhavheadtail_xy[$i]} "%10s %17s %17s\n","#time [ps]","av xy of h-t [A]","stdev";
  &error(4) if(not open($fhavgauche[$i],">","$outdir/NUMGAUCHE_".$mol_name[0][$t]));
  printf {$fhavgauche[$i]} "%10s %17s %17s %17s\n","#time [ps]","density [nm-2]","percent of diheds", "# gauche";
  &error(4) if(not open($fhavtilt[$i],">","$outdir/AV_tiltangle_".$mol_name[0][$t]));
  printf {$fhavtilt[$i]} "%10s %17s %17s\n","#time [ps]","av tilt ang [deg]","stdev";
  &error(4) if(not open($fhavdihedral_tot[$i],">","$outdir/AV_CCCC_DIHEDRAL_".$mol_name[0][$t]."_TOT"));
  printf {$fhavdihedral_tot[$i]} "%10s %17s %17s\n","#time [ps]","av dih ang [deg]","stdev";
  &error(4) if(not open($fhavheadtail_tot[$i],">","$outdir/AV_head-tail_".$mol_name[0][$t]."_TOT"));
  printf {$fhavheadtail_tot[$i]} "%10s %17s %17s\n","#time [ps]","av head-tail [A]","stdev";
  &error(4) if(not open($fhavheadtail_z_tot[$i],">","$outdir/AV_head-tail_z_".$mol_name[0][$t]."_TOT"));
  printf {$fhavheadtail_z_tot[$i]} "%10s %17s %17s\n","#time [ps]","av z of h-t [A]","stdev";
  &error(4) if(not open($fhavheadtail_xy_tot[$i],">","$outdir/AV_head-tail_xy_".$mol_name[0][$t]."_TOT"));
  printf {$fhavheadtail_xy_tot[$i]} "%10s %17s %17s\n","#time [ps]","av xy of h-t [A]","stdev";
  &error(4) if(not open($fhavtilt_tot[$i],">","$outdir/AV_tiltangle_".$mol_name[0][$t]."_TOT"));
  printf {$fhavtilt_tot[$i]} "%10s %17s %17s\n","#time [ps]","av tilt ang [deg]","stdev";
  &error(4) if(not open($fhavgauche_tot[$i],">","$outdir/AV_NUMGAUCHE_".$mol_name[0][$t]."_TOT"));
  printf {$fhavgauche_tot[$i]} "%10s %17s %17s\n","#time [ps]","av \% of diheds","stdev";
  &error(4) if(not open($fhhistdih[$i],">","$outdir/HIST_DIHEDRAL_".$mol_name[0][$t]));
  print {$fhhistdih[$i]} "# Histogram of CCCC dihedral angles [deg] normalized to integral=1\n";
  &error(4) if(not open($fhhistlen[$i],">","$outdir/HIST_H-T-DIST_".$mol_name[0][$t]));
  print {$fhhistlen[$i]} "# Histogram of head to tail length [A] normalized to integral=1\n";
  &error(4) if(not open($fhhistlenz[$i],">","$outdir/HIST_H-T-DIST_Z_".$mol_name[0][$t]));
  print {$fhhistlenz[$i]} "# Histogram of z-component of head to tail length [A] normalized to integral=1\n";
  &error(4) if(not open($fhhistlenxy[$i],">","$outdir/HIST_H-T-DIST_XY_".$mol_name[0][$t]));
  print {$fhhistlenxy[$i]} "# Histogram of xy-component of head to tail length [A] normalized to integral=1\n";
  &error(4) if(not open($fhhisttilt[$i],">","$outdir/HIST_TILT_".$mol_name[0][$t]));
  print {$fhhisttilt[$i]} "# Histogram of tilt angle [deg] normalized to integral=1\n";
  &error(4) if(not open($fhhistgauche[$i],">","$outdir/HIST_NUMGAUCHE_".$mol_name[0][$t]));
  print {$fhhistgauche[$i]} "# Histogram of number of gauche defects normalized to integral=1\n";
}

if($lessoutput) {
  print "begin analysis...\n";
}

foreach $dir (@directories) {
  $nummols = $dir;
  if($lnumsubdir) {
    @numdirs = get_numeric_directories("$dir");
    $dir .= "/$numdirs[$#numdirs]" if(@numdirs>0);
    $dir .= "/min" if(-d "$dir/min");
  } else {
    $dir .= "/$subdir";
  }
  $firsthist=1;
  $f=0;
  $numdihetot = 0;
  $numpas     = 0;
  for($i=0;$i<$tpamax;$i++) {
    $numpas     += $mol_numents[0][$tpas[$i]];
    $numdihetot += $mol_numdihedrals[0][$tpas[$i]];
  }
  print "analyzing directory $dir\n";
  @directories2=get_numeric_directories($dir,$startframe,$endframe);
  push(@directories2,".") if(-e "$dir/HISTORY");
  if(not @directories2) {
    print "**** warning: did not find HISTORY files in $dir\n";
    next;
  }
  &error(6) if(0!=update_numents("$dir/$directories2[0]/FIELD",0));
  for($i=0;$i<@tpas;$i++) {
    print  {$fhavdihedral[$i]}    "# directory $dir\n";
    print  {$fhavheadtail[$i]}    "# directory $dir\n";
    print  {$fhavheadtail_z[$i]}  "# directory $dir\n";
    print  {$fhavheadtail_xy[$i]} "# directory $dir\n";
    print  {$fhavtilt[$i]}        "# directory $dir\n";
    $numgauche[$i][$f]=0;
  }
  foreach $dir2 (@directories2) {
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
      for($i=0;$i<$tpamax;$i++) {
	$t=$tpas[$i];
	$j=0;
	$refc=$refCindex{$mol_name[0][$t]};
	$refp=$refp5index{$mol_name[0][$t]};
	$histdih    = $i*6;
	$histlen    = $i*6+1;
	$histlenz   = $i*6+2;
	$histlenxy  = $i*6+3;
	$histtilt   = $i*6+4;
	$histgauche = $i*6+5;
	for($m=0;$m<$mol_numents[0][$t];$m++) {
	  remap_molecule(\@{$cdata[0][$t][$m]},[0,1],\@{$size[0]},10);
	  # analyze didehdral angles
	  for($d=0;$d<$mol_numdihedrals[0][$t];$d++) {
	    ($a1,$a2,$a3,$a4)=@{$mol_dihedraldata[0][$t][$d]}[1..4];
	    $dihedral[$i][$j]=abs(calc_dihedral_angle(0,0,$d,$t,$m));
	    histogram_add_one($histdih,   $dihedral[$i][$j],$histresdihrad);
	    histogram_add_one($histdihall,$dihedral[$i][$j],$histresdihrad) if(@tpas>1);
	    $j++;
	  }
	  # analyze head-tail-distance
	  $dx = $cdata[0][$t][$m][$refc][0]-$cdata[0][$t][$m][$refp][0];
	  $dy = $cdata[0][$t][$m][$refc][1]-$cdata[0][$t][$m][$refp][1];
	  $dz[$i][$m]   = $cdata[0][$t][$m][$refc][2]-$cdata[0][$t][$m][$refp][2];
	  $dxy[$i][$m]  = sqrt($dx*$dx+$dy*$dy);
	  $dxyz[$i][$m] = sqrt($dx*$dx+$dy*$dy+$dz[$i][$m]*$dz[$i][$m]);
	  histogram_add_one($histlen,   $dxyz[$i][$m],$histreslen);
	  histogram_add_one($histlenall,$dxyz[$i][$m],$histreslen) if(@tpas>1);
	  histogram_add_one($histlenz,   $dz[$i][$m],$histreslen);
	  histogram_add_one($histlenzall,$dz[$i][$m],$histreslen) if(@tpas>1);
	  histogram_add_one($histlenxy,   $dxy[$i][$m],$histreslen);
	  histogram_add_one($histlenxyall,$dxy[$i][$m],$histreslen) if(@tpas>1);
	  # analyze tilt angle
	  $tilt[$i][$m] = atan2($dxy[$i][$m],$dz[$i][$m]);
	  histogram_add_one($histtilt,   $tilt[$i][$m],$histrestiltrad);
	  histogram_add_one($histtiltall,$tilt[$i][$m],$histrestiltrad) if(@tpas>1);
	}
	
############### calculate averages over frame #################################

	($av,$stdev) = calc_stdev(@{$dihedral[$i]});
	printf {$fhavdihedral[$i]}    "%10u %17.10g %17.10g\n",$frame_number[0],$av/pi*180.0,$stdev/pi*180.0;
	push(@{$dihedraltot[$i]},@{$dihedral[$i]});
	
	($av,$stdev) = calc_stdev(@{$dxyz[$i]});
	printf {$fhavheadtail[$i]}    "%10u %17.10g %17.10g\n",$frame_number[0],$av,$stdev;
	push(@{$dxyztot[$i]},@{$dxyz[$i]});
	
	($av,$stdev) = calc_stdev(@{$dz[$i]});
	printf {$fhavheadtail_z[$i]}  "%10u %17.10g %17.10g\n",$frame_number[0],$av,$stdev;
	push(@{$dztot[$i]},@{$dz[$i]});
	
	($av,$stdev) = calc_stdev(@{$dxy[$i]});
	printf {$fhavheadtail_xy[$i]} "%10u %17.10g %17.10g\n",$frame_number[0],$av,$stdev;
	push(@{$dxytot[$i]},@{$dxy[$i]});
	
	($av,$stdev) = calc_stdev(@{$tilt[$i]});
	printf {$fhavtilt[$i]}        "%10u %17.10g %17.10g\n",$frame_number[0],$av/pi*180.0,$stdev/pi*180.0;
	push(@{$tilttot[$i]},@{$tilt[$i]});
      }
      
      if(@tpas>1) {
	($av,$stdev) = calc_stdev(@dihedral);
	printf {$fhavdihedral[$iall]}    "%10u %17.10g %17.10g\n",$frame_number[0],$av/pi*180.0,$stdev/pi*180.0;
	($av,$stdev) = calc_stdev(@dxyz);
	printf {$fhavheadtail[$iall]}    "%10u %17.10g %17.10g\n",$frame_number[0],$av,$stdev;
	($av,$stdev) = calc_stdev(@dz);
	printf {$fhavheadtail_z[$iall]}  "%10u %17.10g %17.10g\n",$frame_number[0],$av,$stdev;
	($av,$stdev) = calc_stdev(@dxy);
	printf {$fhavheadtail_xy[$iall]} "%10u %17.10g %17.10g\n",$frame_number[0],$av,$stdev;
	($av,$stdev) = calc_stdev(@tilt);
	printf {$fhavtilt[$iall]}        "%10u %17.10g %17.10g\n",$frame_number[0],$av/pi*180.0,$stdev/pi*180.0;
      }
      $f++;
    }
    close($fhhist);
    
    undef @dihedral;
    undef @dz;
    undef @dxy;
    undef @dxyz;
    undef @angle;
  }
  
############### print histograms ##############################################
  
  for($i=0;$i<@tpas;$i++) {
    print {$fhavdihedral[$i]} "\n\n";
    print {$fhavheadtail[$i]} "\n\n";
    print {$fhavheadtail_xy[$i]} "\n\n";
    print {$fhavheadtail_z[$i]} "\n\n";
    print {$fhavtilt[$i]} "\n\n";
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
    $histdih    = $i*6;
    $histlen    = $i*6+1;
    $histlenz   = $i*6+2;
    $histlenxy  = $i*6+3;
    $histtilt   = $i*6+4;
    $histgauche = $i*6+5;
    &histogram_normalize_integral($histdih,$histresdihdeg);
    &histogram_normalize_integral($histlen,$histreslen);
    &histogram_normalize_integral($histlenz,$histreslen);
    &histogram_normalize_integral($histlenxy,$histreslen);
    &histogram_normalize_integral($histtilt,$histrestiltdeg);
    &histogram_normalize_integral($histgauche,1);
    
    print {$fhhistdih[$i]} "# directory $dir\n";
    for($j=$histminindex[$histdih];$j<=$histmaxindex[$histdih];$j++) {
      printf {$fhhistdih[$i]} "%17.10g %17.10g %10u\n",$histresdihdeg*$j,$histogram[$histdih][$j],$histnum[$histdih][$j];
    }
    print {$fhhistdih[$i]} "\n\n";
    
    print {$fhhistlen[$i]} "# directory $dir\n";
    for($j=$histminindex[$histlen];$j<=$histmaxindex[$histlen];$j++) {
      printf {$fhhistlen[$i]} "%17.10g %17.10g %10u\n",$histreslen*$j,$histogram[$histlen][$j],$histnum[$histlen][$j];
    }
    print {$fhhistlen[$i]} "\n\n";
    
    print {$fhhistlenz[$i]} "# directory $dir\n";
    for($j=$histminindex[$histlenz];$j<=$histmaxindex[$histlenz];$j++) {
      printf {$fhhistlenz[$i]} "%17.10g %17.10g %10u\n",$histreslen*$j,$histogram[$histlenz][$j],$histnum[$histlenz][$j];
    }
    print {$fhhistlenz[$i]} "\n\n";
    
    print {$fhhistlenxy[$i]} "# directory $dir\n";
    for($j=$histminindex[$histlenxy];$j<=$histmaxindex[$histlenxy];$j++) {
      printf {$fhhistlenxy[$i]} "%17.10g %17.10g %10u\n",$histreslen*$j,$histogram[$histlenxy][$j],$histnum[$histlenxy][$j];
    }
    print {$fhhistlenxy[$i]} "\n\n";
    
    print {$fhhisttilt[$i]} "# directory $dir\n";
    for($j=$histminindex[$histtilt];$j<=$histmaxindex[$histtilt];$j++) {
      printf {$fhhisttilt[$i]} "%17.10g %17.10g %10u\n",$histrestiltdeg*$j,$histogram[$histtilt][$j],$histnum[$histtilt][$j];
    }
    print {$fhhisttilt[$i]} "\n\n";
    
    print {$fhhistgauche[$i]} "# directory $dir\n";
    for($j=$histminindex[$histgauche];$j<=$histmaxindex[$histgauche];$j++) {
      printf {$fhhistgauche[$i]} "%17.10g %17.10g %10u\n",$histresgauchedeg*$j,$histogram[$histgauche][$j],$histnum[$histgauche][$j];
    }
    print {$fhhistgauche[$i]} "\n\n";
    
############### calculate averages over trajectory ############################
    
    ($av,$stdev) = calc_stdev(@{$dihedraltot[$i]});
    printf {$fhavdihedral_tot[$i]}    "%10u %17.10g %17.10g\n",$nummols,$av/pi*180.0,$stdev/pi*180.0;
    ($av,$stdev) = calc_stdev(@{$dxyztot[$i]});
    printf {$fhavheadtail_tot[$i]}    "%10u %17.10g %17.10g\n",$nummols,$av,$stdev;
    ($av,$stdev) = calc_stdev(@{$dztot[$i]});
    printf {$fhavheadtail_z_tot[$i]}  "%10u %17.10g %17.10g\n",$nummols,$av,$stdev;
    ($av,$stdev) = calc_stdev(@{$dxytot[$i]});
    printf {$fhavheadtail_xy_tot[$i]} "%10u %17.10g %17.10g\n",$nummols,$av,$stdev;
    ($av,$stdev) = calc_stdev(@{$tilttot[$i]});
    printf {$fhavtilt_tot[$i]}        "%10u %17.10g %17.10g\n",$nummols,$av/pi*180.0,$stdev/pi*180.0;
  }
  
  if(@tpas>1) {
    ($av,$stdev) = calc_stdev(@dihedraltot);
    printf {$fhavdihedral_tot[$iall]}    "%10u %17.10g %17.10g\n",$frame_number[0],$av/pi*180.0,$stdev/pi*180.0;
    ($av,$stdev) = calc_stdev(@dxyztot);
    printf {$fhavheadtail[$iall]}    "%10u %17.10g %17.10g\n",$frame_number[0],$av,$stdev;
    ($av,$stdev) = calc_stdev(@dztot);
    printf {$fhavheadtail_z_tot[$iall]}  "%10u %17.10g %17.10g\n",$frame_number[0],$av,$stdev;
    ($av,$stdev) = calc_stdev(@dxytot);
    printf {$fhavheadtail_xy_tot[$iall]} "%10u %17.10g %17.10g\n",$frame_number[0],$av,$stdev;
    ($av,$stdev) = calc_stdev(@tilttot);
    printf {$fhavtilt_tot[$iall]}        "%10u %17.10g %17.10g\n",$frame_number[0],$av/pi*180.0,$stdev/pi*180.0;
    
    &histogram_normalize_integral($histdihall,$histresdihdeg);
    &histogram_normalize_integral($histlenall,$histreslen);
    &histogram_normalize_integral($histlenzall,$histreslen);
    &histogram_normalize_integral($histlenxyall,$histreslen);
    &histogram_normalize_integral($histtiltall,$histrestiltdeg);
    &histogram_normalize_integral($histgaucheall,1);
    
    print {$fhhistdih[$iall]} "# directory $dir\n";
    for($j=$histminindex[$histdihall];$j<=$histmaxindex[$histdihall];$j++) {
      printf {$fhhistdih[$iall]} "%17.10g %17.10g %10u\n",$histresdihdeg*$j,$histogram[$histdihall][$j],$histnum[$histdihall][$j];
    }
    print {$fhhistdih[$iall]} "\n\n";
    
    print {$fhhistlen[$iall]} "# directory $dir\n";
    for($j=$histminindex[$histlenall];$j<=$histmaxindex[$histlenall];$j++) {
      printf {$fhhistlen[$iall]} "%17.10g %17.10g %10u\n",$histreslen*$j,$histogram[$histlenall][$j],$histnum[$histlenall][$j];
    }
    print {$fhhistlen[$iall]} "\n\n";
    
    print {$fhhistlenz[$iall]} "# directory $dir\n";
    for($j=$histminindex[$histlenzall];$j<=$histmaxindex[$histlenzall];$j++) {
      printf {$fhhistlenz[$iall]} "%17.10g %17.10g %10u\n",$histreslen*$j,$histogram[$histlenzall][$j],$histnum[$histlenzall][$j];
    }
    print {$fhhistlenz[$iall]} "\n\n";
    
    print {$fhhistlenxy[$iall]} "# directory $dir\n";
    for($j=$histminindex[$histlenxyall];$j<=$histmaxindex[$histlenxyall];$j++) {
      printf {$fhhistlenxy[$iall]} "%17.10g %17.10g %10u\n",$histreslen*$j,$histogram[$histlenxyall][$j],$histnum[$histlenxyall][$j];
    }
    print {$fhhistlenxy[$iall]} "\n\n";
    
    print {$fhhisttilt[$iall]} "# directory $dir\n";
    for($j=$histminindex[$histtiltall];$j<=$histmaxindex[$histtiltall];$j++) {
      printf {$fhhisttilt[$iall]} "%17.10g %17.10g %10u\n",$histrestiltdeg*$j,$histogram[$histtiltall][$j],$histnum[$histtiltall][$j];
    }
    print {$fhhisttilt[$iall]} "\n\n";
    
    histogram_clear($histdihall,$histlenall,$histlenzall,$histlenxyall,$histtiltall);
  }
  
  histogram_clear($histdih,$histlen,$histlenz,$histlenxy,$histtilt);
  
  undef @dihedraltot;
  undef @dztot;
  undef @dxytot;
  undef @dxyztot;
  undef @tilttot;
}


sub update_numents {
  # arguments: filename, fi
  my($fh,$t,@new);
  return 1 if(not open($fh,"<",$_[0]));
  $t=0;
  while($_=<$fh>) {
    if(/NUMMOLS\s+(\d+)/i) {
      $new[$t]=$1;
      $t++;
    }
  }
  return 1 if($field_nummols[$_[1]] != $t);
  @{$mol_numents[$_[1]]}=@new;
  return 0;
}

sub error {
  if(defined $_[1]) {
    print "**** error $_[1]\n";
  } else {
    print "**** unknown error!\n";
  }
  exit $_[0]
}