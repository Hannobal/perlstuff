#!/usr/bin/perl

use hanno_utility;
use dlpoly_utility;
use Switch;
use Storable qw(dclone);
use Math::Trig;

$lessoutput = 0;
@includetypes = ();
$outdir       = $ARGV[0];
$anglecutdeg  = 20;
$distcut      = 3.0;
$startframe   = 0;
$endframe     = 9e20;


if($#ARGV<0) {
  print "input format:\n";
  print " 1. output directory for the output files\n";
  print "-ca <real>   angular cutoff (default: $anglecutdeg deg)\n";
  print "-cr <real>   distance cutoff (default: $distcut A)\n";
  print "-t n*<str>   list of molecules to analyze\n";
  print "-s <int>     start frame\n";
  print "-e <int>     end frame\n";
  print "-lessoutput  do not print current frame\n";
  exit;
}

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    # general flags
    case (/^-ca/i) {
      $i++;
      $anglecutdeg = $ARGV[$i];
      &error(4,"angular cutoff must be a real number >0 <180") if($anglecutdeg<=0 or $anglecutdeg>180);
    } case (/^-cr/i) {
      $i++;
      $distcut = $ARGV[$i];
      &error(4,"cutoff radius must be a real number >=0") if($distcut<=0);
    } case(/^-s/i) {
      $i++;
      $startframe=$ARGV[$i];
      &error(1,"start frame must be an integer") if not check_integer($startframe);
      &error(1) if not check_integer($startframe);
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
      &error(2,"could not interpret flag $ARGV[$i]");
    }
  }
}
print "$startframe $endframe\n";
&error(1,"end frame must not be smaller than start frame") if($endframe<$startframe);
$anglecutrad = (180-$anglecutdeg)/180*pi;
$distcutsq   = $distcut*$distcut;

############### find the relevant directories to read #########################
@directories=get_numeric_directories(".",$startframe,$endframe);
push(@directories,".") if(-e "./HISTORY");
&error(3,"no suitable directories found") if(not @directories);

############### analyze the FIELD file ########################################

exit 1 if(read_field_file("$directories[0]/FIELD",0)!=0);

@donors=();
@acceptors=();
@donormols = ();
@accmols   = ();
for($t=0;$t<$field_nummols[0];$t++) {
  if(@includetypes) {
    next if(not contains(@includetypes,$mol_name[0][$t]));
  }
  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
    if($mol_atomdata[0][$t][$a][0]=~/^(O|N|S)/i) {
      push(@acceptors,[$t,$a]);
      push(@accmols,$t) if(not contains(@accmols,$t));
      for $h (@{$mol_bondatoms[0][$t][$a]}) {
	if($mol_atomdata[0][$t][$h][0]=~/^H/i) {
	  push(@donors,[$t,$a,$h]);
	  push(@donormols,$t) if(not contains(@donormols,$t));
	}
      }
    }
  }
}

if(not @acceptors or not @donors) {
  &error(3,"no acceptors and/or donors found!");
}

print "using $distcut A and $anglecutdeg deg as cutoff\n";
print "found the following donors:\n";
for($d=0;$d<@donors;$d++) {
  $td=$donors[$d][0];
  $ad=$donors[$d][1];
  $ah=$donors[$d][2];
  printf "%-20s %4u (%2s) %4u (%2s)\n",$mol_name[0][$td],$ad,$mol_atomdata[0][$td][$ad][0],$ah,$mol_atomdata[0][$td][$ah][0];
}

print "found the following acceptors:\n";
for($c=0;$c<@acceptors;$c++) {
  $tc=$acceptors[$c][0];
  $ac=$acceptors[$c][1];
  printf "%-20s %4u (%2s)\n",$mol_name[0][$tc],$ac,$mol_atomdata[0][$tc][$ac][0];
}

############### read the HISTORY file from each directory #####################

mkdir($outdir) if(not -d $outdir);
open($fhout,">","$outdir/HBONDS");
print $fhout "# Hydrogen bonds donors-acceptors\n";
printf $fhout "# columns (donor-acceptor):\n#%3u %s\n#%3u %s\n", 1,"frame",2,"all (total number of H-Bonds)";
$i=3;
foreach $td (@donormols) {
  foreach $tc (@accmols) {
    printf $fhout "#%3u %20s -- %-20s\n", $i,$mol_name[0][$td],$mol_name[0][$tc];
    $i++;
  }  
}
print $fhout "\n";

$firsthist = 1;
$f         = 0;
foreach $dir (@directories) {
  print "analyzing directory $dir\n" if($lessoutput);
  &error(2,"cannot open $dir/$histname") if(not open($fhhist,"<","$dir/HISTORY"));
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
    $numall = 0;
    foreach $td (@donormols) {
      foreach $tc (@accmols) {
	$numhbonds[$td][$tc]=0;
      }
    }
    for($d=0;$d<@donors;$d++) {
      $td=$donors[$d][0];
      $ad=$donors[$d][1];
      $ah=$donors[$d][2];
      for($c=0;$c<@acceptors;$c++) {
	$tc=$acceptors[$c][0];
	$ac=$acceptors[$c][1];
	for($md=0;$md<$mol_numents[0][$td];$md++) {
	  for($mc=0;$mc<$mol_numents[0][$tc];$mc++) {
	    next if($td==$tc and $mc==$md and $ac==$ad);
# 	    print "$md $ad $ah  $mc $ac  ";
	    $dsq = calc_dsq_orthocell(0,$td,$md,$ad, $tc,$mc,$ac);
# 	    print " ",sqrt($dsq);
	    next if($dsq>$distcutsq);
	    $angle = calc_angle_orthocell(0, $td,$md,$ad, $td,$md,$ah, $tc,$mc,$ac);
# 	    print " ",$angle/pi*180;
	    next if($angle<$anglecutrad);
# 	    if($dsq>$distcutsq or $angle<$anglecutrad) {
# 	      print "  no hbond\n"
# 	    } else {
# 	      print "  hbond\n"
# 	    }
	    $numhbonds[$td][$tc]++;
	    $numall++;
	  }
	}
      }
    }
    # print output
    printf $fhout "%10u %7u",$frame_number[0],$numall;
    foreach $td (@donormols) {
      foreach $tc (@accmols) {
	printf $fhout " %7u",$numhbonds[$td][$tc];
      }  
    }
    print $fhout "\n";
    $f++;
  }
  close($fhhist);
}
print "\n" unless($lessoutput);



close($fhout);

sub error {
  if(defined $_[1]) {
    print "**** error $_[1]\n";
  } else {
    print "**** unknown error!\n";
  }
  exit $_[0]
}
(END)
