#!/usr/bin/perl
use Switch;
use Storable qw(dclone);

use dlpoly_utility;
use hanno_utility;

# repeats the content of a DL_POLY HISTORY file in x- and y- direction
# and cuts a surface with a certain slab thickness

if($#ARGV<6) {
  print "input format:\n";
  print " 1. dl_poly HISTORY-file\n";
  print " 2. corresponding dl_poly FIELD-file\n";
  print " 3. output-HISTORY with ending\n";
  print " 4. output-FIELD with ending\n";
  print " 5. repetitions in direction of vector a (integer)\n";
  print " 6. repetitions in direction of vector b (integer)\n";
  print " 7. repetitions in direction of vector c (integer)\n";
  print "[start and end frame]\n";
  exit;
}

$startframe = 0;
$endframe   = 9e20;
$filenameinhist  = $ARGV[0];
$filenameinfld   = $ARGV[1];
$filenameouthist = $ARGV[2];
$filenameoutfld  = $ARGV[3];
$rep[0]=$ARGV[4];
$rep[1]=$ARGV[5];
$rep[2]=$ARGV[6];
$factor=$rep[0]*$rep[1]*$rep[2];
if($#ARGV>6) {
  $startframe = $ARGV[7];
  if(not check_integer($startframe)) {
    print "**** error: start frame must be integer!\n"; exit 1;
  }
}
if($#ARGV>7) {
  $endframe = $ARGV[8];
  if(not check_integer($endframe)) {
    print "**** error: end frame must be integer!\n"; exit 1;
  }
}

if(not check_integer(@rep) or $rep[0]<1 or $rep[1]<1 or $rep[2]<1) {
  print "**** error: repetitions must be integer numbers >= 1\n";
  exit 1;
}

exit 1 if(read_field_file($filenameinfld,0) != 0);
if(not open($fhinhist,"<",$filenameinhist)) {
  print "**** error: input file $filenameinhist could not be opened!\n";
  exit 1;
}

if(not open($fhouthist,">",$filenameouthist)) {
  print "**** error: output file $filenameouthist could not be opened!\n";
  exit 1;
}

while(read_history_timestep($fhinhist,0,0)==0) {
  next if($frame_number[0]<$startframe);
  last if($frame_number[0]>$endframe);
  if($periodic_key[0]==0) {
    print "**** error: periodic conditions required for calculation\n";
    exit 1;
  } elsif(not ($periodic_key[0]==1 or $periodic_key[0]==2 or $periodic_key[0]==6)) {
    print "**** error:, this script only works for orthogonal cells!\n";
    exit 1;
  } elsif($periodic_key[0]==6) { # for mirroring
    @remapaxes = (0,1);
  } else {
    @remapaxes = (0,1,2);
  }
  print "processing frame $frame_number[0]\n";
  for($t=0;$t<$field_nummols[0];$t++) {
    for($m=0;$m<$mol_numents[0][$t];$m++) {
      remap_molecule(\@{$cdata[0][$t][$m]},\@remapaxes,\@{$size[0]});
    }
  }

  $dx = -(($rep[0]-1)*$cell[0][0][0]+($rep[1]-1)*$cell[0][1][0]+($rep[2]-1)*$cell[0][2][0])/2;
  $dy = -(($rep[0]-1)*$cell[0][0][1]+($rep[1]-1)*$cell[0][1][1]+($rep[2]-1)*$cell[0][2][1])/2;
  $dz = -(($rep[0]-1)*$cell[0][0][2]+($rep[1]-1)*$cell[0][1][2]+($rep[2]-1)*$cell[0][2][2])/2;

  @{$cdata[1]} = ();
  for($t=0;$t<$field_nummols[0];$t++) {
    @{$cdata[1][$t]}=();
    $newmolid=0;
    for($m=0;$m<$mol_numents[0][$t];$m++) {
      for($x=0;$x<$rep[0];$x++) {
	for($y=0;$y<$rep[1];$y++) {
	  for($z=0;$z<$rep[2];$z++) {
	    @{$cdata[1][$t][$newmolid]} = @{dclone \@{$cdata[0][$t][$m]}};
	    for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	      $cdata[1][$t][$newmolid][$a][0] += $x*$cell[0][0][0]+$y*$cell[0][1][0]+$z*$cell[0][2][0]+$dx;
	      $cdata[1][$t][$newmolid][$a][1] += $x*$cell[0][0][1]+$y*$cell[0][1][1]+$z*$cell[0][2][1]+$dy;
	      $cdata[1][$t][$newmolid][$a][2] += $x*$cell[0][0][2]+$y*$cell[0][1][2]+$z*$cell[0][2][2]+$dz;
	    } # end for $a
	    $newmolid++;
	  } # end for $z
	} # end for $y
      } # end for $x
    } # end for $m
  } # end for $t
  @{$cell[1]} = @{$cell[0]};
  for($c=0;$c<3;$c++) {
    $cell[1][$c][0] *= $rep[$c];
    $cell[1][$c][1] *= $rep[$c];
    $cell[1][$c][2] *= $rep[$c];
  }

  $periodic_key[1] = $periodic_key[0];
  $config_key[1]   = $config_key[0];
  $frame_numatoms[1] = $frame_numatoms[0]*$factor;
  write_history_timestep($fhouthist,0,1);
}
close($fhinhist,$fhouthist);

for($t=0;$t<$field_nummols[0];$t++) {
  $mol_numents[0][$t] *= $factor;
}

write_field_file($filenameoutfld,0);