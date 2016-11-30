#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;

if($#ARGV<1) {
  print "input format:\n";
  print "1. <output-directory>\n";
  print "2. split job? <y/n/h>\n";
  print "3. [histogram resolution (default: 0.5 A)]\n";
  print "4. [begin frame/temperature]\n";
  print "5. [end frame/temperature]\n";
  exit 1;
} 

$outdir = $ARGV[0];

if($#ARGV>1) {
  $histstep=$ARGV[2];
  if($histstep<=0) {
    print "**** error: histogram resolution must be greater than 0!\n";
    exit 1;
  }
} else {
  $histstep=0.5;
}

if ($#ARGV>2) {
  $startframe = $ARGV[3];
} else {
  $startframe = 0;
}

if ($#ARGV>3) {
  $endframe = $ARGV[4];
} else {
  $endframe = 9e20;
}

$nextpercent = 5;

if(not $ARGV[1] =~ /^[ynh]/i) {
  print "**** error: did not recognize input for \"split job?\"!";
  exit 1;
} elsif($ARGV[1] =~ /^y/i) {
  $splitjob = 1;
  $heatup   = 0;
} elsif($ARGV[1] =~ /^h/i) {
  $heatup   = 1;
  $splitjob = 1;
} else {
  $heatup   = 0;
  $splitjob = 0;
}

check_integer($startframe) or die "error: startframe $startframe is not an integer number!\n";
check_integer($endframe)   or die "error: endframe $endframe is not an integer number!\n";

if($endframe<$startframe) {
  print "**** error: end frame must not be smaller than start frame\n";
  exit;
}

@directories = get_numeric_directories(".",$startframe,$endframe);

############### read the HISTORY file from each directory #####################
$f         = 0;
$starttime = 0;
$mincenterhistindex   = 9e20;
$maxcenterhistindex   = -9e20;
$minbottomhistindex   = 9e20;
$maxbottomhistindex   = -9e20;
loopdir:foreach $dir (@directories) {
  print "reading directory $dir\n" if $splitjob;
  $err = &read_field_file("$dir/FIELD",0);
  if($err>0) {
    exit 1;
  }
  # create array with numbers for molecule types that contain C60
  @c60mols  = ();
  @c60atoms = ();
  $numfullatomstot = 0;
  $numfulls        = 0;
  for($t=0;$t<$field_nummols[0];$t++) {
    if($mol_name[0][$t] =~ /C60/i) {
      $numfulls += $mol_numents[0][$t];
      push(@c60mols,$t);
      $numfullatoms[$t] = 0;
      for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	if(uc($mol_atomdata[0][$t][$a][0]) eq "CA" or 
	   uc($mol_atomdata[0][$t][$a][0]) eq "CX") {
	  $numfullatomstot += $mol_numents[0][$t];
	  $numfullatoms[$t]++;
	  push(@{$c60atoms[$t]},$a);
	}
      }
    }
  }
  open(CONTROL, "<", "$dir/CONTROL") or print "file CONTROL was not found in $dir\n";
  # read the time step interval out of the CONTROL file
  while(<CONTROL>) {
    if(/\s*timestep/) {
      ($timestep) = /\s*timestep\s+(\S+)/;
    } elsif(/traj/) {
      ($trajinterval) = /\s*traj\S*\s+\S+\s+(\S+)/;
    }
  }
  close(CONTROL);
  ######################## read HISTORY file ##################################
  if(not open($fhhist, "<", "$dir/HISTORY")) {
    print "**** warning: file HISTORY was not found in directory $dir\n";
    next loopdir;
  }
  while(read_history_timestep($fhhist,0,0)==0) {
#     print "evaluating frame $f\n";
    $frame[$f]=$frame_number[0];
    if($frame_numatoms[0] != $field_numatoms[0]) {
      print "**** error in HISOTRY file: $field_numatoms[0] expected, $frame_numatoms[0] given in HISTORY!\n";
      exit 1;
    }
#     print TEST "$numatoms\n\n";
    last if($frame[$f]>$endframe and not $heatup);
    $time[$f]=$frame[$f]*$timestep;
    $time[$f]+=$starttime if $time[$f]<=$starttime;
    # determine surface z-position
    $surfzmax=-9e20;
    for($t=0;$t<$field_nummols[0];$t++) {
      next if(not $mol_name[0][$t] =~ /aluminum/i);
      for($m=0;$m<$mol_numents[0][$t];$m++) {
	$surfzmax = $cdata[0][$t][$m][0][2] if($surfzmax<$cdata[0][$t][$m][0][2]);
      }
    }
    # remap fullerene molecules 
#     for($t=0;$t<$field_nummols[0];$t++) {
#       next if(not $mol_name[0][$t] =~ /C60/i); # fullerene molecules only!
#       for($m=0;$m<$mol_numents[0][$t];$m++) {
# 	remap_molecule(\@{$cdata[0][$t][$m]},[0,1],\@{$size[0]});
#       }
#     }
    
############### evaluate the fullerene positions ##############################
    if($frame[$f]>=$startframe or $heatup) {
      $tot_full_mindist_bottom[$f] = 9e20;
      $b = 0; #buckyball index
      foreach $t (values @c60mols) {
	for($m=0;$m<$mol_numents[0][$t];$m++) {
	  $fullbottom[$f][$b]=9e20;
	  foreach $a (values @{$c60atoms[$t]}) {
# 	    $fullcenter[$f][$b][0] += $cdata[0][$t][$m][$a][0];
# 	    $fullcenter[$f][$b][1] += $cdata[0][$t][$m][$a][1];
	    $fullcenter[$f][$b][2] += $cdata[0][$t][$m][$a][2];
	    if($cdata[0][$t][$m][$a][2]<$fullbottom[$f][$b]) {
	      $fullbottom[$f][$b] = $cdata[0][$t][$m][$a][2];
	    }
	  }
# 	  $fullcenter[$f][$b][0] /= $numfullatoms[$t];
# 	  $fullcenter[$f][$b][1] /= $numfullatoms[$t];
	  $fullcenter[$f][$b][2] /= $numfullatoms[$t];
	  $fullcenter[$f][$b][2] -= $surfzmax;
	  $fullbottom[$f][$b]    -= $surfzmax;
	  $av_full_mindist_center[$f] += $fullcenter[$f][$b][2];
	  $av_full_mindist_bottom[$f] += $fullbottom[$f][$b];
	  if($tot_full_mindist_bottom[$f]>$fullbottom[$f][$b]) {
	    $tot_full_mindist_bottom[$f] = $fullbottom[$f][$b];
	  }
	  # analyze minimum and maximum values for histograms
	  $centerhistindex = int($fullcenter[$f][$b][2]/$histstep);
	  $histogramcenter[$f]{$centerhistindex}++;
	  $bottomhistindex = int($fullbottom[$f][$b]/$histstep);
	  $histogrambottom[$f]{$bottomhistindex}++;
	  $maxcenterhistindex = $centerhistindex if($centerhistindex > $maxcenterhistindex);
	  $mincenterhistindex = $centerhistindex if($centerhistindex < $mincenterhistindex);
	  $maxbottomhistindex = $bottomhistindex if($bottomhistindex > $maxbottomhistindex);
	  $minbottomhistindex = $bottomhistindex if($bottomhistindex < $minbottomhistindex);
	  $b++;
	}
      }
      $av_full_mindist_center[$f] /= $numfulls;
      $av_full_mindist_bottom[$f] /= $numfulls;
      $f++;
    }
  }
  close($fhhist);
  if($err>1) {
    exit 1;
  }
  $starttime=$time[$f-1];
}

$totframes = $f;

############### generate the output files #####################################
if(not -d $outdir) {
  mkdir $outdir or die "cannot create output-directory: $!";
}
#open(FULLERENES, ">", "FUL_DISTANCE") or die "cannot open file FUL_DISTANCE: $!";
open(FULL_AV_DIST_BOTTOM, ">", "$outdir/full_av_dist_bottom") or die "cannot open file full_av_dist_bottom: $!";
open(FULL_AV_DIST_CENTER, ">", "$outdir/full_av_dist_center") or die "cannot open file full_av_dist_center: $!";
open(FULL_DIST_MIN, ">", "$outdir/full_dist_min") or die "cannot open file full_dist_min: $!";
open(FULL_HIST_DIST_CENTER, ">", "$outdir/full_hist_dist_center") or die "cannot open file full_hist_dist_center: $!";
open(FULL_HIST_DIST_BOTTOM, ">", "$outdir/full_hist_dist_bottom") or die "cannot open file full_hist_dist_bottom: $!";
open(FULL_SINGLE_DIST_CENTER, ">", "$outdir/full_single_dist_center") or die "cannot open file full_single_dist_center: $!";
open(FULL_SINGLE_DIST_BOTTOM, ">", "$outdir/full_single_dist_bottom") or die "cannot open file full_single_dist_bottom: $!";

printf FULL_AV_DIST_CENTER "#%12s %53s\n", 'time [ps]', 'average distance of fullerene center to surface [A]';
printf FULL_AV_DIST_BOTTOM    "#%12s %53s\n", 'time [ps]', 'average distance of fullerene bottom to surface [A]';
printf FULL_DIST_MIN       "#%12s %53s\n", 'time [ps]', 'minimum distance of all fullerene bottoms to surface [A]';

for($f=0;$f<$totframes;$f++) {
  printf FULL_AV_DIST_CENTER "%13.4f %13.10f\n", $time[$f], $av_full_mindist_center[$f];
  printf FULL_AV_DIST_BOTTOM    "%13.4f %13.10f\n", $time[$f], $av_full_mindist_bottom[$f];
  printf FULL_DIST_MIN       "%13.4f %13.10f\n", $time[$f], $tot_full_mindist_bottom[$f];
  print FULL_HIST_DIST_CENTER "timestep $f ($time[$f] ps)\n";
  print FULL_HIST_DIST_BOTTOM "timestep $f ($time[$f] ps)\n";
  for($i=$mincenterhistindex;$i<=$maxcenterhistindex;$i++) {
    printf FULL_HIST_DIST_CENTER "%13.4f %13u\n", $i*$histstep, $histogramcenter[$f]{$i};
  }
  for($i=$minbottomhistindex;$i<=$maxbottomhistindex;$i++) {
    printf FULL_HIST_DIST_BOTTOM "%13.4f %13u\n", $i*$histstep, $histogrambottom[$f]{$i};
  }
  print FULL_HIST_DIST_CENTER "\n\n";
  print FULL_HIST_DIST_BOTTOM "\n\n";
}

print FULL_SINGLE_DIST_CENTER "#single fullerenes $b\n";
print FULL_SINGLE_DIST_BOTTOM "#single fullerenes $b\n";
for($b=0;$b<$numfulls;$b++) {
  print FULL_SINGLE_DIST_CENTER "#fullerene number $b\n";
  print FULL_SINGLE_DIST_BOTTOM "#fullerene number $b\n";
  printf FULL_SINGLE_DIST_CENTER "#%12s %53s\n", 'time [ps]', 'average distance of fullerene center to surface [A]';
  printf FULL_SINGLE_DIST_BOTTOM "#%12s %53s\n", 'time [ps]', 'average distance of fullerene bottom to surface [A]';
  for($f=0;$f<$totframes;$f++) {
  printf FULL_SINGLE_DIST_CENTER "%13.4f %13.10f\n", $time[$f], $fullcenter[$f][$b];
  printf FULL_SINGLE_DIST_BOTTOM "%13.4f %13.10f\n", $time[$f], $fullbottom[$f][$b];
  }
  print FULL_SINGLE_DIST_CENTER "\n\n";
  print FULL_SINGLE_DIST_BOTTOM "\n\n";
}
