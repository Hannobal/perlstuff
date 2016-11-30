#!/usr/bin/perl

use hanno_utility;
use dlpoly_utility;
use Switch;
use Storable qw(dclone);
use Math::Trig;

$|=1;

$lessoutput = 0;
$langular   = 0;
$lspatial   = 0;
@includetypes = ();
$histname = "HISTORY";
$fldname  = "FIELD";
$outdir   = $ARGV[0];
$gauchecutoffdeg = 30;
$histresdihdeg   = 5;
$histresangdeg   = 10;
$histresxy       = 2;
$startframe      = 0;
$endframe        = 9e20;
$radiuscutoff    = 9999;

if($#ARGV<0) {
  print "input format:\n";
  print " 1. output directory for the output files\n";
  print "-gc <real>   cutoff for gauche defects (default: $gauchecutoffdeg deg)\n";
  print "-rd <real>   histogram resolution for dihedrals (default: $histresdihdeg deg)\n";
  print "-rs <real>   spatial histogram resolution (default: no spatial analysis)\n";
  print "-ra <real>   angular histogram resolution (default: no angular analysis)\n";
  print "-c <real     cutoff around center (default: $radiuscutoff A)\n";
  print "-t n*<str>   list of molecules to analyze\n";
  print "-s <int>     start frame\n";
  print "-e <int>     end frame\n";
  print "-lessoutput  do not print current frame\n";
  exit;
}

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    # general flags
    case (/^-gc/i) {
      $i++;
      $gauchecutoffdeg = $ARGV[$i];
      &error(7) if($gauchecutoffdeg<=0);
    } case (/^-rd/i) {
      $i++;
      $histresdihdeg = $ARGV[$i];
      &error(6) if($histresdihdeg<=0);
    } case (/^-rs/i) {
      $i++;
      $histresxy = $ARGV[$i];
      $lspatial=1;
      &error(6) if($histresxy<=0);
    } case (/^-ra/i) {
      $i++;
      $histresangdeg = $ARGV[$i];
      &error(6) if($histresangdeg<=0);
      &error(9) if(360%$histresangdeg!=0);
      $langular=1;
    } case (/^-c/i) {
      $i++;
      $radiuscutoff = $ARGV[$i];
      &error(6) if($radiuscutoff<=0);
    } case(/^-s/i) {
      $i++;
      $startframe = $ARGV[$i];
      &error(1) if not check_integer($startframe);
    } case(/^-e/i) {
      $i++;
      $endframe = $ARGV[$i];
      &error(2) if not check_integer($endframe);
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
      &error(10,$ARGV[$i]);
    }
  }
}

$gauchecutoffrad = $gauchecutoffdeg/180*pi;
$histresdihrad   = $histresdihdeg/180*pi;
$histresangrad   = $histresangdeg/180*pi;

@directories=get_numeric_directories(".");
if(@directories) {
  @directories=get_numeric_directories(".",$startframe,$endframe);
}
push(@directories,".") if(-e "./$histname");
&error(8) if(not @directories);

exit 1 if(read_field_file("$directories[0]/$fldname",0)!=0);

# find the molecule types to analyze
if(not @includetypes) {
  @includetypes = @{$mol_name[0]};
  @tlst         = (0..$#{$mol_name[0]});
} else {
  for($i=0;$i<@includetypes;$i++) {
    $tlst[$i] = -1;
    for($t=0;$t<$field_nummols[0];$t++) {
      if($mol_name[0][$t]=~/^$includetypes[$i]$/i) {
	$tlst[$i] = $t;
	last;
      }
    }
    &error(11,$includetypes[$i]) if($tlst[$i]==-1);
  }
}

# remove dihedrals with atoms other than C3
foreach $t (@tlst) {
  for($d=0;$d<$mol_numdihedrals[0][$t];$d++) {
      for($i=1;$i<5;$i++) {
	if(not $mol_atomdata[0][$t][$mol_dihedraldata[0][$t][$d][$i]][0]=~/C3/) {
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
}

for($i=0;$i<@includetypes;$i++) {
  if($mol_numdihedrals[0][$tlst[$i]]==0) {
    print "**** warning: molecule $mol_name[0][$tlst[$i]] has no carbon backbone to analyze!\n";
    splice(@includetypes,$i,1);
    splice(@tlst,$i,1);
    $i--;
  }
}

mkdir($outdir) if(not -d $outdir);
&error(4) if(not open($fhavdihedral,">","$outdir/AV_CCCC_DIHEDRAL"));
print $fhavdihedral "# average dihedral angle in degrees within $radiuscutoff A from center\n";
print $fhavdihedral "# molecules used for analysis: ",join(" ",@includetypes),"\n";
printf $fhavdihedral "#%9s %17s %17s\n","frame","dihedral","stdev";
&error(4) if(not open($fhnumgauche,">","$outdir/NUMGAUCHE"));
print $fhnumgauche "# percentage of dihedral angles <= $gauchecutoffdeg (cutoff=$radiuscutoff A)\n";
print $fhnumgauche "# molecules used for analysis: ",join(" ",@includetypes),"\n";
printf $fhnumgauche "%10s %17s %17s %10s\n","frame","density [nm-2]","percent of dihedrals","# gauche";

$firsthist=1;
$f=0;
if($lessoutput) {
  print "begin analysis...\n";
} else {
  print "\n" 
}
foreach $dir (@directories) {
  exit 1 if(not open($fhhist,"<","$dir/$histname"));
  if($startframe>0 and $firsthist) {
    &error(5,$startframe) if(not find_history_timestep($fhhist,$startframe));
  }
  $firsthist=0;
  $err=0;
  while($err==0) {
    $err=read_history_timestep($fhhist,0,0);
    last   if($err<0);
    exit 1 if($err>0);
    last if($frame_number[0]>$endframe);
    print "\ranalyzing timestep $frame_number[0]" unless($lessoutput);
    $fn[$f]=$frame_number[0];
    $numgauche[$f] = 0;
    $numdihtotcut = 0;
    undef @valdihedral;
    $minsize = min($size[0][0],$size[0][1],$radiuscutoff);
    $i=0;
    foreach $t (@tlst) {
      for($m=0;$m<$mol_numents[0][$t];$m++) {
  #       remap_molecule(\@{$cdata[0][$t][$m]},[0,1],\@{$size[0]},10);
	for($d=0;$d<$mol_numdihedrals[0][$t];$d++) {
	  ($a1,$a2,$a3,$a4)=@{$mol_dihedraldata[0][$t][$d]}[1..4];
	  @{$pos[0]}=@{$cdata[0][$t][$m][$a1]}[0..2];
	  @{$pos[1]}=@{$cdata[0][$t][$m][$a2]}[0..2];
	  @{$pos[2]}=@{$cdata[0][$t][$m][$a3]}[0..2];
	  @{$pos[3]}=@{$cdata[0][$t][$m][$a4]}[0..2];
	  #remap if necessary
	  @{$pos[0]} = unwrap(\@{$pos[0]},\@{$pos[1]});
	  @{$pos[2]} = unwrap(\@{$pos[2]},\@{$pos[1]});
	  @{$pos[3]} = unwrap(\@{$pos[3]},\@{$pos[1]});
	  # analyze
	  $dx=($pos[1][0]+$pos[2][0])/2;
	  $dy=($pos[1][1]+$pos[2][1])/2;
	  $dxy = sqrt($dx*$dx+$dy*$dy);
	  $value=abs(calc_dihedral(\@pos));
	  if($lspatial) {
	    histogram_add_value(0,$dxy,$histresxy,$value);
	    histogram_add_one(5,$dxy,$histresxy);
	  }
	  if($value<$gauchecutoffrad) {
	    histogram_add_one(2,$dxy,$histresxy) if($lspatial);
	    histogram_add_one(7,$d,1); #histogram of gauche defects vs. position in chain
	  }
	  if($dxy<$radiuscutoff) {
	    $numdihtotcut++;
	    if($value<$gauchecutoffrad) {
	      $numgauche[$f]++;
	    }
	    $valdihedral[$i]=$value;
	    histogram_add_one(1,$value,$histresdihrad);
	    $i++;
	  }
	  if($langular) {
	    if($dxy<=$minsize) {
	      $angle=atan2($dy,$dx);
	      $angle+=2*pi if($angle<0);
	      histogram_add_one(6,$angle,$histresangrad);
	      histogram_add_value(3,$angle,$histresangrad,$value);
	      if($value<$gauchecutoffrad) {
		histogram_add_one(4,$angle,$histresangrad);
	      }
	    }
	  }
	}
      }
    }
    $areatot = calc_area($radiuscutoff,0,@{$size[0]})/100;
    ($value,$stdev) = calc_stdev(@valdihedral);
    printf $fhavdihedral "%10u %17.10g %17.10g\n",$fn[$f],$value*180/pi,$stdev*180/pi;
    printf $fhnumgauche  "%10u %17.10g %17.10g %10u\n",$fn[$f],$numgauche[$f]/$areatot,$numgauche[$f]/$numdihtotcut*100,$numgauche[$f];
    $f++;
  }
  close($fhhist);
}

if($lessoutput) {
  print "last frame analyzed: $frame_number[0]\n";
} else {
  print "\n";
}

close($fhnumgauche,$fhavdihedral);
print "\n";
$numframes=$f+1;
# &histogram_normalize_individual(0,3);     # av. dihedral vs. radius/angle
&histogram_normalize_integral(1,$histresdihrad); # overall dihedral histogram
# &histogram_normalize_value(2,$numframes); # radial distribution of gauche defects
# &histogram_normalize_value(4,$numframes); # angular distribution of gauche defects
&histogram_normalize_integral(7,1); # histogram of gauche defect position in chain

open($fhhist,">","$outdir/HISTOGRAM_DIHEDRAL");
print $fhhist "# histogram of all dihedral angles [deg] (radial cutoff=$radiuscutoff A)\n";
print $fhhist "# averaged over $numframes snapshots\n";
print $fhhist "# molecules used for analysis: ",join(" ",@includetypes),"\n";
for($i=$histminindex[1];$i<=$histmaxindex[1];$i++) {
  printf $fhhist "%17.10g %17.10g\n",$histresdihdeg*$i,$histogram[1][$i];
}
close($fhhist);

open($fhhist,">","$outdir/HISTOGRAM_GAUCHE_POS");
print $fhhist "# position of gauche defects (cutoff=$gauchecutoffdeg deg) in chain\n";
print $fhhist "# averaged over $numframes snapshots\n";
print $fhhist "# molecules used for analysis: ",join(" ",@includetypes),"\n";
printf $fhhist "%10s %17s\n","position","percent of dihedrals";
for($i=$histminindex[7];$i<=$histmaxindex[7];$i++) {
  printf $fhhist "%10u %17.10g\n",$i,$histogram[7][$i]*100;
}
close($fhhist);

if($lspatial) {
  open($fhhist,">","$outdir/AV_DIHEDRAL_RADIAL_DISTRIB");
  print $fhhist "# average CCCC dihedral angle [deg] as function of distance from center [A]\n";
  print $fhhist "# averaged over $numframes snapshots\n";
  print $fhhist "# molecules used for analysis: ",join(" ",@includetypes),"\n";
  printf $fhhist "#%16s %17s %17s %10s\n","radius [A]","value [deg]","stdev [deg]","# values";
  for($i=$histminindex[0];$i<=$histmaxindex[0];$i++) {
    ($av,$stdev) = calc_stdev(@{$histvals[0][$i]});
    printf $fhhist "%17.10g %17.10g %17.10g %10u\n",$histresxy*$i,$av/pi*180,$stdev/pi*180,$histnum[0][$i];
  }
  close($fhhist);

  open($fhhist,">","$outdir/NUMGAUCHE_RADIAL_DISTRIB");
  print $fhhist "# number of gauche defects (cutoff=$gauchecutoffdeg deg) per frame as function of distance from center [A]\n";
  print $fhhist "# averaged over $numframes snapshots\n";
  print $fhhist "# molecules used for analysis: ",join(" ",@includetypes),"\n";
  printf $fhhist "%17s %17s %17s\n","radius [A]","density [nm-2]","percent of dihedrals";
  for($i=$histminindex[2];$i<=$histmaxindex[2];$i++) {
    $area = calc_area(($i+1)*$histresxy,$i*$histresxy,@{$size[0]})/100; # in nm^2
    next if($area==0);
    printf $fhhist "%17.10g %17.10g %17.10g\n",$histresxy*$i,$histogram[2][$i]/($numframes*$area),$histogram[2][$i]/$histogram[5][$i]*100;
  }
  close($fhhist);
}

if($langular) {
  open($fhhist,">","$outdir/AV_DIHEDRAL_ANGULAR_DISTRIB");
  print $fhhist "# average CCCC dihedral angle [deg] within  $minsize A of the center as polar plot [deg]\n";
  print $fhhist "# averaged over $numframes snapshots\n";
  print $fhhist "# molecules used for analysis: ",join(" ",@includetypes),"\n";
  printf $fhhist "#%16s %17s %17s %10s\n","angle [deg]","value [deg]","stdev [deg]","# values";
  for($i=$histminindex[3];$i<=$histmaxindex[3];$i++) {
    ($av,$stdev) = calc_stdev(@{$histvals[3][$i]});
    next if($area==0);
    printf $fhhist "%17.10g %17.10g %10u\n",$histresangdeg*$i,$av/pi*180,$stdev/pi*180,$histnum[3][$i];
  }
  ($av,$stdev) = calc_stdev(@{$histvals[3][0]});
  printf $fhhist "%17.10g %17.10g %10u\n",360,$av/pi*180,$stdev/pi*180,$histnum[3][0];
  close($fhhist);

  open($fhhist,">","$outdir/NUMGAUCHE_ANGULAR_DISTRIB");
  print $fhhist "# number of gauche defects (cutoff=$gauchecutoffdeg deg) per frame within $minsize A of the center as polar plot [deg]\n";
  print $fhhist "# averaged over $numframes snapshots\n";
  print $fhhist "# molecules used for analysis: ",join(" ",@includetypes),"\n";
  printf $fhhist "%17s %17s %17s\n","angle [deg]","density [nm-2]","percent of dihedrals";
  $area = pi*$minsize*$minsize*$histresangdeg/360.0/100;
  for($i=$histminindex[4];$i<=$histmaxindex[4];$i++) {
    next if($area==0);
    printf $fhhist "%17.10g %17.10g %17.10g\n",$histresangdeg*$i,$histogram[4][$i]/($numframes*$area),$histogram[4][$i]/$histogram[6][$i]*100;
  }
  printf $fhhist "%17.10g %17.10g %17.10g\n",360,$histogram[4][0]/($numframes*$area),$histogram[4][0]/$histogram[6][0]*100;
  close($fhhist);
}

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

sub unwrap {
  my @pos1 = @{$_[0]};
  my @pos2 = @{$_[1]};
  my($d,$c);
  for($c=0;$c<2;$c++) {
    $d = $pos1[$c] - $pos2[$c];
    if($d > $size[0][$c]) {
      $pos1[$c] -= 2 * $size[0][$c];
    } elsif($d <-$size[0][$c]) {
      $pos1[$c] += 2 * $size[0][$c];
    }
  }
  return @pos1;
}

sub calc_dihedral {
  my @pos  = @{$_[0]};
  my (@vec,@vecc,@vecn1,@vecn2);
  $vecc[0] = $pos[2][0]-$pos[1][0];
  $vecc[1] = $pos[2][1]-$pos[1][1];
  $vecc[2] = $pos[2][2]-$pos[1][2];
  $vec[0]  = $pos[0][0]-$pos[1][0];
  $vec[1]  = $pos[0][1]-$pos[1][1];
  $vec[2]  = $pos[0][2]-$pos[1][2];
  @vecn1=vector_product(\@vec,\@vecc);
  @vecn1=normalize_vector(\@vecn1);
  $vec[0]  = $pos[3][0]-$pos[2][0];
  $vec[1]  = $pos[3][1]-$pos[2][1];
  $vec[2]  = $pos[3][2]-$pos[2][2];
  @vecn2=vector_product(\@vec,\@vecc);
  @vecn2=normalize_vector(\@vecn2);
  return acos(dot_product(\@vecn1,\@vecn2));
}

sub calc_area {
  my(@syssize,$router,$rinner,$area,$halfdiag,$h);
  $router = $_[0];
  $rinner = $_[1];
  @syssize = @_[2..3];
  $halfdiag = sqrt($syssize[0]*$syssize[0]+$syssize[1]*$syssize[1]);
  return 0 if($rinner>$halfdiag);
  $router=$halfdiag if($router>$halfdiag);
  $area = pi*($router*$router-$rinner*$rinner);
  # check x-axis with outer radius
  $h = $router-$syssize[0];
  $area -= 2*segment_area($router,$h) if($h>0);
  # check y-axis with outer radius
  $h = $router-$syssize[1];
  $area -= 2*segment_area($router,$h) if($h>0);
  # check x-axis with inner radius
  $h = $rinner-$syssize[0];
  $area += 2*segment_area($rinner,$h) if($h>0);
  # check y-axis with inner radius
  $h = $rinner-$syssize[1];
  $area += 2*segment_area($rinner,$h) if($h>0);
  return $area;
}

sub segment_area {
  my $r = $_[0];
  my $h = $_[1];
  return $r*$r*acos(1-$h/$r)-sqrt(2*$r*$h-$h*$h)*($r-$h);
}