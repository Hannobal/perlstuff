#!/usr/bin/perl

# calculates the tilt angles of certain phosponic acids from DL_POLY HISTORY
# files and their develeopment in time

use Math::Trig;
use Cwd;
use Switch;
use hanno_utility;
use dlpoly_utility;
# use aloxsam_utility;

if($#ARGV<1) {
  print "input format:\n";
  print "1. <output-directory>\n";
  print "2. split job? <y/n/h>\n";
  print "3. [begin frame]\n";
  print "4. [end frame]\n";
  print "5. [step frame]\n";
  exit 1;
} 

$outdir = $ARGV[0];
if(not -d $outdir) {
  mkdir $outdir or die "cannot create output-directory: $!";
}

if ($#ARGV>1) {
  $startframe = $ARGV[2];
} else {
  $startframe = 0;
}

if ($#ARGV>2) {
  $endframe = $ARGV[3];
} else {
  $endframe = 9e20;
}

if ($#ARGV>3) {
  $stepframe = $ARGV[4];
} else {
  $stepframe = 1;
}

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

check_integer($startframe) or die "**** error: $startframe is not an integer number!\n";
check_integer($endframe) or die "**** error: $endframe is not an integer number!\n";

if($endframe<$startframe) {
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
  "C60-C18-PA-d"     => 52,
  "F15-C18-PA-d"     => 52,
  "PHDA-d"           => 48,
  #double deprotonated PAs
  "C10-PA-2d" => 0
);

############### find the relevant directories to read #########################
opendir(DIR, ".");
if(not $splitjob) {
  @directories=(".")
} else {
  @files = readdir(DIR);
  closedir(DIR);
  # find the relevant folders
  $i=0;
  foreach $file (@files) {
    if (check_integer($file)) {
      $files2[$i]=$file;
      $i++;
    }
  }
  @files = sort {$a <=> $b} @files2;
  # find index of start directory
  if($heatup) {
    $startdirindex = 0;
    $enddirindex   = $#files;
  } else {
    $startdirindex=-1;
    for($i=0;$i<=$#files;$i++) {
      if($files[$i]>=$startframe) {
	$startdirindex = $i;
	last;
      }
    }
    if($startdirindex==-1) {
      print "start frame is larger than simulation frame";
      exit 1;
    }
    # find index of end directory
    for($i=0;$i<=$#files;$i++) {
      if(($files[$i]>=$endframe) or ($i==$#files)) {
	$enddirindex = $i;
	last;
      }
    }
  }
  @directories=@files[$startdirindex..$enddirindex];
}

$endframe=$directories[$enddirindex] if($endframe==9e20);
$numframes=$endframe-$startframe;

############### read the HISTORY file from each directory #####################
$f         = 0;
$starttime = 0;
$minangle = 9e20;
$maxangle = -9e20;
foreach $dir (@directories) {
  print "reading directory $dir\n" if $splitjob;
  exit 1 if(read_field_file("$dir/FIELD",0)!=0);
  
  ######################## read HISTORY file ##################################
  if(not open($fhhist, "<", "$dir/HISTORY")) {
    print "**** warning: file HISTORY was not found in directory $dir\n";
    next;
  }
  $err=read_history_timestep($fhhist,0,0);
  while($err==0) {
    exit 1 if ($err>0);
    last if($frame_number[0]>$endframe);
    if(($frame_number[0] % $stepframe) == 0) {
#       print "analyzing frame $frame_number[0]\n";
      $zmin = $zsurf-0.5;
      $zmax = $zsurf+0.7;
      for($t=0;$t<$field_nummols[0];$t++) {
	if($mol_name[0][$t] =~ /alumini?um/i
	or $mol_name[0][$t] =~ /oxygen/i
	or $mol_name[0][$t] =~ /hydroxide/i) {
	  @{$cdata[0][$t]} = ();
	  $frame_numatoms[0] -= $mol_numents[0][$t]*$mol_numatoms[0][$t];
	} elsif($mol_name[0][$t] =~ /PA-d/) {
	  $p5=$refp5index{$mol_name[0][$t]};
	  mol:for($m=0;$m<$mol_numents[0][$t];$m++) {
	    remap_molecule(\@{$cdata[0][$t][$m]},[0,1],\@{$size[0]});
	    # determine lowest oxygen atom in molecule
	    $minpos=9e20;
	    for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	      if($mol_atomdata[0][$t][$a][0] eq 'O' and $cdata[0][$t][$m][$a][2]<$minpos) {
		$minpos=$cdata[0][$t][$m][$a][2];
		$amin=$a;
	      }
	      # add a proton
	    } # end for a
	      $a=$mol_numatoms[0][$t];
	      $cdata[0][$t][$m][$a][9]       = 'HG';
	      @rot_matrix=gen_rot_matrix([0,0,1],rand(2*pi));
	      @vec = rotate_vector([0,-0.939692620785908,-0.342020143325669],\@rot_matrix);
	      $cdata[0][$t][$m][$a][0] = $vec[0]+$cdata[0][$t][$m][$amin][0];
	      $cdata[0][$t][$m][$a][1] = $vec[1]+$cdata[0][$t][$m][$amin][1];
	      $cdata[0][$t][$m][$a][2] = $vec[2]+$cdata[0][$t][$m][$amin][2];
	  } # end for m
	  $frame_numatoms[0] += $mol_numents[0][$t];
	} # end if PA-d
      } #end for t
      $number = sprintf("%09u",$frame_number[0]);
      open($fhxyz,">","$outdir/$number.xyz");
      write_xyz_timestep($fhxyz,0,"");
      close($fhxyz);
      print "$number\n";
    } # end if frame_number
    $err=read_history_timestep($fhhist,0,0);
  } # end while read_history_timestep
  close($fhhist);
}

$totframes=$f;