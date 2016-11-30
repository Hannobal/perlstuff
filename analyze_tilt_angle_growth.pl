#!/usr/bin/perl

# calculates the tilt angles (to the surface normal => 0 deg for perpendicular
# molecules) of certain phosponic acids from DL_POLY HISTORY files and their
# develeopment in time for force_docking_z-type systems

use Math::Trig;
use Cwd;
use Switch;
use hanno_utility;
use dlpoly_utility;

if($#ARGV<0) {
  print "input format:\n";
  print " 1. <output-directory>\n";
  print "[2. histogram resolution (default: 0.5 deg)]\n";
  print "[3. start and end step]\n";
  exit 1;
} 

$outdir = $ARGV[0];

if($#ARGV>0) {
  $histres=$ARGV[1];
  if($histres<=0) {
    print "**** error: histogram resolution must be greater than 0!\n";
    exit 1;
  }
} else {
  $histres=0.5;
}

if ($#ARGV>1) {
  $startstep = $ARGV[2];
} else {
  $startstep = 1;
}

if ($#ARGV>2) {
  $endstep = $ARGV[3];
} else {
  $endstep = 9e20;
}

check_integer($startstep) or die "**** error: $startstep is not an integer number!\n";
check_integer($endstep) or die "**** error: $endstep is not an integer number!\n";

if($endstep<$startstep) {
  print "**** error: end frame must not be smaller than start frame\n";
  exit;
}

%refp5index = (
  # neutral PAs
  "C10-PA"     => 10,
  "C14-PA"     => 10,
  "C17-PA"     => 17,
  "C1-PA"      =>  1,
  "C60-C18-PA" => 10,
  #single deprotonated PAs
  "BTBT-C12-PA-d"    => 10,
  "C10-PA-d"         => 10,
  "C12-C60-C18-PA-d" => 10,
  "C14-PA-d"         => 10,
  "C17-PA-d"         => 17,
  "C1-PA-d"          =>  1,
  "C2-4T-C12-PA-d"   => 10,
  "C60-C18-PA-d"     => 10,
  "F15-C18-PA-d"     => 10,
  "PHDA-d"           => 10,
  #double deprotonated PAs
  "C10-PA-2d" => 10
);

%refCindex = (
  # neutral PAs
  "C10-PA"     => 0,
  "C14-PA"     => 38,
  "C17-PA"     => 15,
  "C1-PA"      => 0,
  "C60-C18-PA" => 52,
  #single deprotonated PAs
  "BTBT-C12-PA-d"    => 36,
  "C10-PA-d"         => 0,
  "C12-C60-C18-PA-d" => 52,
  "C14-PA-d"         => 38,
  "C17-PA-d"         => 15,
  "C1-PA-d"          => 0,
  "C2-4T-C12-PA-d"   => 36,
  "C60-C18-PA-d"     => 51, # !!! was 52 before moving HG to bottom of atom list
  "F15-C18-PA-d"     => 52,
  "PHDA-d"           => 48,
  #double deprotonated PAs
  "C10-PA-2d" => 0
);

$maxhistindex=int(90/$histres);

############### find the relevant directories to read #########################
print "reading directories\n";
@directories=get_numeric_directories(".",$startstep,$endstep);

############### read the HISTORY file from each directory #####################

foreach $dir (@directories) {
  @subdirs = get_numeric_directories("$dir");
  $dir .= "/$subdirs[$#subdirs]" if(@subdirs>0);
  $dir .= "/min" if(-d "$dir/min");
  if(not -e "$dir/HISTORY") {
    print "**** warning: did not find file $dir/HISTORY\n";
    next;
  }
  print "analyzing $dir/HISTORY\n";
  # read FIELD file
  exit 1 if(read_field_file("$dir/FIELD",0)!=0);
  $d   = int($dir);
  
  # read HISTORY file
  if(not open($fhhist,"<","$dir/HISTORY")) {
    print "**** error: could not open file $dir/HISTORY: $!\n";
    exit 1;
  }
  $err = 0;
  $f   = 0;
  $minangle[$d] =  9e20;
  $maxangle[$d] = -9e20;
  
  $numpads{"all"}[$d] = 0;
  for($t=0;$t<$field_nummols[0];$t++) {
    next unless($mol_name[0][$t] =~ /PA-2?d/);
    $numpads{$t}[$d]     = $mol_numents[0][$t];
    $numpads{"all"}[$d] += $mol_numents[0][$t];
    $tav_angle{$t}[$d]  = 0;
    $tav_length{$t}[$d] = 0;
    $tav_len_xy{$t}[$d] = 0;
    $tav_len_z{$t}[$d]  = 0;
  }
  $tav_angle{"all"}[$d]  = 0;
  $tav_length{"all"}[$d] = 0;
  $tav_len_xy{"all"}[$d] = 0;
  $tav_len_z{"all"}[$d]  = 0;
  
  while($err==0) {
    $err=read_history_timestep($fhhist,0,0);
    last if($err<0);
    exit 1 if($err>0);
    # analyze timestep
    if($periodic_key[0]==6) { # for mirroring
      @remapaxes = (0,1);
    } elsif($periodic_key[0]>0) {
      @remapaxes = (0,1,2);
    }
    for($t=0;$t<$field_nummols[0];$t++) {
      next unless($mol_name[0][$t] =~ /PA-2?d/);
      
      # initialize single molceule histogram
      for($i=0;$i<=$maxhistindex;$i++) {
	$histogram{$t}[$d][$f][$i]=0;
      }
      
      $av_angle{$t}[$d][$f]  = 0;
      $av_length{$t}[$d][$f] = 0;
      $av_len_xy{$t}[$d][$f] = 0;
      $av_len_z{$t}[$d][$f]  = 0;
      
      # calculate angle and fill histogram data
      $refC  = $refCindex{$mol_name[0][$t]};
      $refp5 = $refp5index{$mol_name[0][$t]};
      for($m=0;$m<$mol_numents[0][$t];$m++) {
	remap_molecule(\@{$cdata[0][$t][$m]},\@remapaxes,\@{$size[0]}) if($periodic_key[0]>0);
	$dx = $cdata[0][$t][$m][$refC][0]-$cdata[0][$t][$m][$refp5][0];
	$dy = $cdata[0][$t][$m][$refC][1]-$cdata[0][$t][$m][$refp5][1];
	$dz = $cdata[0][$t][$m][$refC][2]-$cdata[0][$t][$m][$refp5][2];
	$len_xy{$t}[$m][$d][$f]  = sqrt($dx*$dx+$dy*$dy);
	$len_z{$t}[$m][$d][$f]   = $dz;
	$length{$t}[$m][$d][$f]  = sqrt($dx*$dx+$dy*$dy+$dz*$dz);
	$angle{$t}[$m][$d][$f] = atan2($len_xy{$t}[$m][$d][$f],$dz)/pi*180;
	$minangle[$d] = $angle{$t}[$m][$d][$f] if($angle{$t}[$m][$d][$f]<$minangle[$d]);
	$maxangle[$d] = $angle{$t}[$m][$d][$f] if($angle{$t}[$m][$d][$f]>$maxngle[$d]);
	$histogramindex = int($angle{$t}[$m][$d][$f]/$histres);
	$av_angle{$t}[$d][$f]    += $angle{$t}[$m][$d][$f];
	$av_angle{"all"}[$d][$f] += $angle{$t}[$m][$d][$f];
	$av_length{$t}[$d][$f]   += $length{$t}[$m][$d][$f];
	$av_len_xy{$t}[$d][$f]   += $len_xy{$t}[$m][$d][$f];
	$av_len_z{$t}[$d][$f]    += $len_z{$t}[$m][$d][$f];
	$tav_angle{$t}[$d]     += $angle{$t}[$m][$d][$f];
	$tav_angle{"all"}[$d]  += $angle{$t}[$m][$d][$f];
	$tav_length{$t}[$d]    += $length{$t}[$m][$d][$f];
	$tav_len_xy{$t}[$d]    += $len_xy{$t}[$m][$d][$f];
	$tav_len_z{$t}[$d]     += $len_z{$t}[$m][$d][$f];
	$tav_length{"all"}[$d] += $length{$t}[$m][$d][$f];
	$tav_len_xy{"all"}[$d] += $len_xy{$t}[$m][$d][$f];
	$tav_len_z{"all"}[$d]  += $len_z{$t}[$m][$d][$f];
	$histogram{$t}[$d][$f][$histogramindex]++;
	$tavhistogram{$t}[$d][$histogramindex]++;
      }# done for $m
      if($numpads{$t}[$d]>0) {
	$av_angle{$t}[$d][$f]    /= $numpads{$t}[$d];
	$av_length{$t}[$d][$f]   /= $numpads{$t}[$d];
	$av_len_xy{$t}[$d][$f]   /= $numpads{$t}[$d];
	$av_len_z{$t}[$d][$f]    /= $numpads{$t}[$d];
      }
    }
    $av_angle{"all"}[$d][$f]  /= $numpads{"all"}[$d];
    $f++;
  }
  close($fhhist);
  $totframes[$d]=$f;
  $minhistindex[$d] = int($minangle[$d]/$histres)-1;
  $minhistindex[$d] = 0 if($minhistindex[$d]<0);
  $maxhistindex[$d] = int($maxangle[$d]/$histres)+1;
  foreach $key ( keys %tav_angle ) {
    next if($numpads{$key}[$d]==0);
    $tav_angle{$key}[$d]  /= $totframes[$d]*$numpads{$key}[$d];
    $tav_length{$key}[$d] /= $totframes[$d]*$numpads{$key}[$d];
    $tav_len_xy{$key}[$d] /= $totframes[$d]*$numpads{$key}[$d];
    $tav_len_z{$key}[$d]  /= $totframes[$d]*$numpads{$key}[$d];
    for($i=$minhistindex[$d];$i<=$maxhistindex[$d];$i++) {
      $tavhistogram{$key}[$d][$i] /= $totframes[$d];
    }
    print "$key ",$tav_angle{$key}[$d],"\n";
  }
}

############### generate the output files #####################################

# !!! histograms and average angles per frame are calculated but not outputted

if(not -d $outdir) {
  mkdir $outdir or die "cannot create output-directory: $!";
}

foreach $key ( keys %tav_angle ) {
  if(check_integer($key)) {$name=$mol_name[0][$key]; $numents=$mol_numents[0][$key]} else {$name='all';$numents=$numpads}
  open(TAV_ANGLE, ">","$outdir/tav_angle_$name");
  print  TAV_ANGLE  "# time and entity averaged tilt angle for $name\n";
  printf TAV_ANGLE  "#%12s %16s %12s %12s\n", 'growthstep', 'tilt angle [deg]','# entities','# frames';
  foreach $dir (@directories) {
    $d=int($dir);
    if(defined( $av_angle{$key}[$d])) {
      printf  TAV_ANGLE " %12u %16.10f %12u %12u\n",$d,$tav_angle{$key}[$d],$numpads{$key}[$d],$totframes[$d];
    }
  }
  close(TAV_ANGLE);
}


foreach $key ( keys %tav_length ) {
  if(check_integer($key)) {$name=$mol_name[0][$key]; $numents=$mol_numents[0][$key]} else {$name='all';$numents=$numpads}
  open(TAV_LENGTH, ">","$outdir/tav_length_$name");
  open(TAV_LENGTH_Z, ">","$outdir/tav_z_length_$name");
  open(TAV_LENGTH_XY, ">","$outdir/tav_xy_length_$name");
  print  TAV_LENGTH  "# time and entity averaged chain length between atoms ".$refPindex{$name}.
                     " and ".$refCindex{$name}."for $name\n";
  printf TAV_LENGTH  "#%12s %16s %12s\n", 'growthstep', 'chain length [A]','# entities';
  print  TAV_LENGTH_Z  "# time and entity averaged chain length in z-direction between atoms ".$refPindex{$name}.
                     " and ".$refCindex{$name}."for $name\n";
  printf TAV_LENGTH_Z  "#%12s %16s %12s\n", 'growthstep', 'chain length [A]','# entities';
  print  TAV_LENGTH_XY  "# time and entity averaged chain length in xy-direction between atoms ".$refPindex{$name}.
                     " and ".$refCindex{$name}."for $name\n";
  printf TAV_LENGTH_XY  "#%12s %16s %12s\n", 'growthstep', 'chain length [A]','# entities';
  foreach $dir (@directories) {
    $d=int($dir);
    if(defined( $av_angle{$key}[$d])) {
      printf  TAV_LENGTH " %12u %16.10f %16.10f %16.10f %12u\n",$d,$tav_length{$key}[$d],$numpads{$key}[$d],
        $tav_len_xy{$key}[$d],$numpads{$key}[$d], $tav_len_z{$key}[$d],$numpads{$key}[$d];
      printf  TAV_LENGTH_Z " %12u %16.10f %12u\n",$d,$tav_len_z{$key}[$d],$numpads{$key}[$d];
      printf  TAV_LENGTH_XY " %12u %16.10f %12u\n",$d,$tav_len_xy{$key}[$d],$numpads{$key}[$d];
    }
  }
  close(TAV_LENGTH);
  close(TAV_LENGTH_Z);
  close(TAV_LENGTH_XY);
}