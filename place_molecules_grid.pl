#!/usr/bin/perl
use Math::Trig;
use Switch;

use dlpoly_utility;
use hanno_utility;
#use Storable qw(dclone);

$doo=2.755333333;
$minmoldist= $doo+1;
$tolerance=0.2;
$periodicvec[0][0]=$doo;
$periodicvec[0][1]=0;
$periodicvec[1][0]=$doo*cos(pi/3);
$periodicvec[1][1]=$doo*sin(pi/3);
@islandcenter=(0,0);
@holecenter=(0,0);

if($#ARGV<8) {
  print "input format:\n";
  print " 1. name of CONFIG/REVCON-file with surface\n";
  print " 2. name of corresponding FIELD-file\n";
  print " 3. name of mol2-file with molecule\n";
  print " 4. name of output CONFIG-file\n";
  print " 5. name of output FIELD-file\n";
  print " 6. x-component of vector 1\n";
  print " 7. y-component of vector 1\n";
  print " 8. x-component of vector 2\n";
  print " 9. y-component of vector 2\n";
  print "10. z-position\n";
  print "-gc <2*real> define grid origin\n";
  print "-i  <real>  create island with defined radius\n";
  print "-h  <real>  create hole with defined radius\n";
  print "-ic <2*real> center of the island\n";
  print "-s  <str>   name of solvent (overlapping molecules will be removed)\n";
  print "-sv <2*real> shift vector \n";
  exit 1;
}

$filenameincfg  = $ARGV[0];
$filenameinfld  = $ARGV[1];
$filenamemol2   = $ARGV[2];
$filenameoutcfg = $ARGV[3];
$filenameoutfld = $ARGV[4];
$molperiodicvec[0][0]=$ARGV[5];
$molperiodicvec[0][1]=$ARGV[6];
$molperiodicvec[1][0]=$ARGV[7];
$molperiodicvec[1][1]=$ARGV[8];
$xpos=0;
$ypos=0;
$zpos=$ARGV[9];
$molshiftvec[0]=0;
$molshiftvec[1]=0;

for($i=10;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-sv/) {
      $i++;
      $molshiftvec[0]=$ARGV[$i];
      $i++;
      $molshiftvec[1]=$ARGV[$i];
      if(not check_real(@molshiftvec)) {
	print "**** error: argument for flag -sv must be two real numbers!\n";
	exit 1;
      }
    } case (/^-s/) {
      $i++;
      $solvname = $ARGV[$i];
      $lsolv    = 1;
    } case (/^-hc/) {
      $i++;
      if(not check_real(@ARGV[$i..$i+1])) {
	print "**** error: center for island must be real numbers!\n";
	exit 1;
      }
      @islandcenter=@ARGV[$i..$i+1];
      $i++;
    } case (/^-h/) {
      $i++;
      if(not check_real($ARGV[$i])) {
	print "**** error: hole radius must be a real number!\n";
	exit 1;
      }
      $holersq=$ARGV[$i]*$ARGV[$i];
    } case (/^-gc/) {
      $i++;
      if(not check_real(@ARGV[$i..$i+1])) {
	print "**** error: origin of grid must be real numbers!\n";
	exit 1;
      }
      $xpos = $ARGV[$i];
      $i++;
      $ypos = $ARGV[$i];
    } case (/^-ic/) {
      $i++;
      if(not check_real(@ARGV[$i..$i+1])) {
	print "**** error: center for island must be real numbers!\n";
	exit 1;
      }
      @holecenter=@ARGV[$i..$i+1];
      $i++;
    } case (/^-i/) {
      $i++;
      if(not check_real($ARGV[$i])) {
	print "**** error: island radius must be a real number!\n";
	exit 1;
      }
      $islandrsq=$ARGV[$i]*$ARGV[$i];
    } else {
      print "**** error: expected optional flag but got \"$ARGV[$i]\"!\n";
      exit 1;
    }
  }
}

if($molperiodicvec[0][0]==$molperiodicvec[1][0] and $molperiodicvec[0][1]==$molperiodicvec[1][1]) {
  print "error: vectors must not be equal\n";
  exit 1;
}

# indices for cdata array:
# 0 name
# 1 x
# 2 y
# 3 z
# 4 velocity x
# 5 velocity y
# 6 velocity z
# 7 force x
# 8 force y
# 9 force z

#load molecule
exit 1 if(read_mol2_file($filenamemol2,0)!=0);
$mol2name=$mol2_name[0];
$mol2name =~ s/\.\S*$//; # remove file extension
$mol2name = substr($mol2name,rindex($mol2name,'/')+1); # remove path (everything before last "/")

# read FIELD file
exit 1 if(read_field_file($filenameinfld,0)!=0);
$foundmol    = -1;
$foundsolv   = -1;
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /^$mol2name$/i) {
    $foundmol=$t;
  } elsif($mol_name[0][$t] =~ /$solvname/i) {
    $foundsolv=$t;
  }
}

die "error: could not find entry for molecule $mol2name in FIELD file $filenameinfld!\n" if($foundmol==-1);
die "error: could not find entry for solvent molecule in FIELD file $filenameinfld!\n"
if($foundsolv==-1 and $lsolv);

#read CONFIG file
exit 1 if(read_config_file($filenameincfg,0,0)!=0);
switch($periodic_key[0]) {
  case [1,2,6] {
  } case 0 {
    print "error: periodic conditions needed for calculation\n";
    exit 1;
  } else {
    print "sorry, this script is only implemented for orthogonal cells\n";
    exit 1;
  }
}

$imaxx = 100;
$imaxy = 100;
$i=0;

$pos[2]=$zpos;
gridposy : for($iy=-$imaxy;$iy<=$imaxy;$iy++) {
  gridposx : for($ix=-$imaxx;$ix<=$imaxx;$ix++) {
    # determine position according to grid
    $pos[0] = $ix*$molperiodicvec[0][0]+$iy*$molperiodicvec[1][0]+$molshiftvec[0]+$xpos;
    $pos[1] = $ix*$molperiodicvec[0][1]+$iy*$molperiodicvec[1][1]+$molshiftvec[1]+$ypos;
    next if($pos[0] < -$size[0][0] or $pos[0] > $size[0][0] or $pos[1] < -$size[0][1] or $pos[1] > $size[0][1]);
    if(defined($islandrsq)) {
      $dx=$pos[0]-$islandcenter[0];
      $dy=$pos[1]-$islandcenter[1];
      next if(($dx*$dx+$dy*$dy)>$islandrsq);
    }
    if(defined($holersq)) {
      $dx=$pos[0]-$holecenter[0];
      $dy=$pos[1]-$holecenter[1];
      next if(($dx*$dx+$dy*$dy)<$holersq);
    }
    # quickly check with molecules in neighboring cells to see if the molecule fits
    for($j=0;$j<@finalpos;$j++) {
      for($x=-1;$x<2;$x++) {
	for($y=-1;$y<2;$y++) {
	  $dx=$pos[0]-$finalpos[$j][0]-$x*$cell[0][0][0]-$y*$cell[0][1][0];
	  $dy=$pos[1]-$finalpos[$j][1]-$x*$cell[0][0][1]-$y*$cell[0][1][1];
	  $dist=sqrt($dx*$dx+$dy*$dy);
	  if($dist<3) {
	    next gridposx
	  }
	}
      }
    }
    # add a molecule at the new position
    push(@{$finalpos},[@pos]);
    $m=@{$cdata[0][$foundmol]};
    @{$cdata[0][$foundmol][$m]} = mol2_to_cdata(0);
    move_mol(\@{$cdata[0][$foundmol][$m]},\@pos);
    $field_numatoms[0] += $mol_numatoms[0][$foundmol];
    $mol_numents[0][$foundmol]++;
  }
}

write_config_file($filenameoutcfg,0,$config_title[0]);
write_field_file($filenameoutfld,0);
exit 0;
