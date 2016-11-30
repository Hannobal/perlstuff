#!/usr/bin/perl

#calculates the centers of fullerene cages and outputs them as xyz trajectory file

use dlpoly_utility;
use aloxsam_utility;
use hanno_utility;
use Switch;

$outdir="analysis";
$startframe=0;
$endframe=9.0e20;
$histres=0.1;
$histname="HISTORY";
$fldname="FIELD";

for($i=0;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case("-h") {
      print "input format:\n",
      "-o <str>    output directory\n",
      "-s <int>    start frame\n",
      "-e <int>    end frame\n";
      &error(1);
    } case("-o") {
      $i++;
      $outdir = $ARGV[$i];
    } case("-s") {
      $i++;
      $startframe = $ARGV[$i];
      &error(1,"start frame must be an integer number!") if not check_integer($startframe);
    } case("-e") {
      $i++;
      $endframe = $ARGV[$i];
      &error(1,"end frame must be an integer number!") if not check_integer($endframe);
    } else {
      &error(1,"Unknown flag $ARGV[$i]!");
    }
  }
}


# find directories
@directories=get_numeric_directories(".");
@directories=get_numeric_directories(".",$startframe,$endframe);
push(@directories,".") if(-e "./$histname");
&error(2,"No suitable directories found") if(not @directories);

# find ions and surface
$dir = $directories[0];
&error(3) if(read_field_file("$directories[0]/$fldname",0)!=0);
@tfull=();
@numfullatoms=();
$numfulltot=0;
for($t=0;$t<$field_nummols[0];$t++) {
  next if(not $mol_name[0][$t]=~/(C60)/i);
  push(@tfull,$t);
  $numfulltot+=$mol_numents[0][$t];
  $numfullatoms[$t]=0;
  findfullatoms:for($a=0;$a<$mol_numatoms[0][$t];$a++) {
    next  if(not $mol_atomdata[0][$t][$a][0] =~ /(^CA$|^CX$)/i);
    foreach $b (@{$mol_bondatoms[0][$t][$a]}) {
      #skip atoms that are bound to e.g. C3 or OS (bingel atom)
      next findfullatoms if(not $mol_atomdata[0][$t][$b][0] =~ /(^CA$|^CX$)/i);
    }
    $numfullatoms[$t]++;
    push(@{$afull[$t]},$a);
  }
  print "molecule $mol_name[0][$t] contains $numfullatoms[$t] C60 atoms\n";
}
&error(4, "Did not find fullerenes in FIELD file!") if(not @tfull);
&error(5) if(read_config_file("$directories[0]/CONFIG",0,0)!=0);
undef @cdata;


mkdir($outdir) if(not -d $outdir);
open($fhxyz,">","$outdir/FULL_POS.xyz");
#analyze HISTORY files
$firsthist=1;
if($lessoutput) {
  print "begin analysis...\n";
}
foreach $dir (@directories) {
  &error(7,"Could not open HISTOY file in directory $dir!") if(not open($fhhist,"<","$dir/$histname"));
  if($startframe>0 and $firsthist) {
    &error(6,"Did not find frame $startframe in file $dir/$histname!")
    if(not find_history_timestep($fhhist,$startframe));
  }
  $firsthist=0;
  $err=0;
  while($err==0) {
    $err=read_history_timestep($fhhist,0,0);
    last   if($err<0);
    exit 1 if($err>0);
    last if($frame_number[0]>$endframe);
    print "\ranalyzing timestep $frame_number[0]" unless($lessoutput);
    print $fhxyz "$numfulltot\nframe $frame_number[0]\n";
    foreach $t (@tfull) {
      for($m=0;$m<$mol_numents[0][$t];$m++) {
	@pos=(0,0,0);
	remap_molecule(\@{$cdata[0][$t][$m]},[0,1],\@{$size[0]},$afull[$t][0]);
	foreach $a (@{$afull[$t]}) {
	  $pos[0] += $cdata[0][$t][$m][$a][0];
	  $pos[1] += $cdata[0][$t][$m][$a][1];
	  $pos[2] += $cdata[0][$t][$m][$a][2];
	}
	$pos[0] /= $numfullatoms[$t];
	$pos[1] /= $numfullatoms[$t];
	$pos[2] /= $numfullatoms[$t];
	if($pos[0]>$size[0][0]) { $pos[0] -= 2*$size[0][0] }
	elsif($pos[0]<-$size[0][0]) { $pos[0] += 2*$size[0][0] }
	if($pos[1]>$size[0][1]) { $pos[1] -= 2*$size[0][1] }
	elsif($pos[1]<-$size[0][1]) { $pos[1] += 2*$size[0][1] }
	printf $fhxyz "X   %10.5f %10.5f %10.5f\n",@pos;
      }
    }
  }
  close($fhhist);
}
close($fhxyz);

if($lessoutput) {
  print "last frame analyzed: $frame_number[0]\n";
} else {
  print "\n";
}

sub error {
  print "**** error: $_[1]\n" if(defined($_[1]));
  exit $_[0];
}