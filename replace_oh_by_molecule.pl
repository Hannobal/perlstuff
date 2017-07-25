#!/usr/bin/perl
use Math::Trig;
use Switch;

use dlpoly_utility;
use hanno_utility;
#use Storable qw(dclone);

$tolerance=0.5;
@islandcenter=(0,0);
@holecenter=(0,0);

#          x         
#         A          
#        /           
#     b/             
#     /              
#   /                
#  /        a        
# x---------------->x

if($#ARGV<8) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file with surface\n";
  print "2. name of corresponding FIELD-file\n";
  print "3. name of mol2-file with molecule\n";
  print "4. name of output CONFIG-file\n";
  print "5. name of output FIELD-file\n";
  print "6. reconst. vector 1 factor for a\n";
  print "7. reconst. vector 1 factor for b\n";
  print "8. reconst. vector 2 factor for a\n";
  print "9. reconst. vector 2 factor for b\n";
  print "-b  [0123]   zero/mono/bi/tri-dentate binding\n";
  print "-i  <real>   create circular island with defined radius\n";
  print "-r  <2*real> create rectangular island with defined width and height\n";
  print "-h  <real>   create hole with defined radius\n";
  print "-ic <2*real> center of the island\n";
  print "-d  <real>   additional distance to surface\n";
  print "-t  <real>   tolerance for OH positions (default: $tolerance A)\n";
  print "-s  <str>    name of solvent (overlapping molecules will be removed)\n";
  print "-sa <real>   shift vector in multiples of a\n";
  print "-sb <real>   shift vector in multiples of b\n";
  print "-c  <2*real> center of molecule\n";
  print "-f  <str>    surface (alumina0001,mag111)\n";
  exit 1;
}

$filenameincfg  = $ARGV[0];
$filenameinfld  = $ARGV[1];
$filenamemol2   = $ARGV[2];
$filenameoutcfg = $ARGV[3];
$filenameoutfld = $ARGV[4];
$factor[0][0]=$ARGV[5]; # first index: new vector, second index: old vector
$factor[0][1]=$ARGV[6];
$factor[1][0]=$ARGV[7];
$factor[1][1]=$ARGV[8];
@molshiftvec=(0,0);
$surfdist=0;
$system="alumina0001";

for($i=9;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-sa/) {
      $i++;
      $molshiftvec[0]=$ARGV[$i];
      if(not check_real($molshiftvec[0])) {
	print "**** error: argument for flag -sa must be a real number!\n";
	exit 1;
      }
    } case (/^-sb/) {
      $i++;
      $molshiftvec[1]=$ARGV[$i];
      if(not check_real($molshiftvec[1])) {
	print "**** error: argument for flag -sb must be a real number!\n";
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
      @holecenter=@ARGV[$i..$i+1];
      $i++;
    } case (/^-h/) {
      $i++;
      if(not check_real($ARGV[$i])) {
	print "**** error: hole radius must be a real number!\n";
	exit 1;
      }
      $holersq=$ARGV[$i]*$ARGV[$i];
    } case (/^-ic/) {
      $i++;
      if(not check_real(@ARGV[$i..$i+1])) {
	print "**** error: center for island must be real numbers!\n";
	exit 1;
      }
      @islandcenter=@ARGV[$i..$i+1];
      $i++;
    } case (/^-i/) {
      $i++;
      if(not check_real($ARGV[$i])) {
	print "**** error: island radius must be a real number!\n";
	exit 1;
      }
      $islandrsq=$ARGV[$i]*$ARGV[$i];
    } case (/^-r/) {
      $i++;
      if(not check_real($ARGV[$i])) {
	print "**** error: argument for flag -t must be a real number!\n";
	exit 1;
      }
      $islandwidth=$ARGV[$i]/2.0;
      $i++;
      if(not check_real($ARGV[$i])) {
	print "**** error: argument for flag -t must be a real number!\n";
	exit 1;
      }
      $islandheight=$ARGV[$i]/2.0;
    } case (/^-b/) {
      $i++;
      if(not $ARGV[$i] =~ /^[0123]+$/) {
	print "**** error: argument for flag -b must be 0,1,2 or 3!\n";
	exit 1;
      }
      $bindingmode=$ARGV[$i];
    } case (/^-d/) {
      $i++;
      $surfdist=$ARGV[$i];
      if(not check_real($surfdist)) {
	print "**** error: argument for flag -d must be a real number!\n";
	exit 1;
      }
    } case (/^-t/) {
      $i++;
      $tolerance=$ARGV[$i];
      if(not check_real($tolerance)) {
	print "**** error: argument for flag -t must be a real number!\n";
	exit 1;
      }
    } case (/^-c/) {
      $i++;
      if(not check_real(@ARGV[$i..$i+1])) {
	print "**** error: center must be real numbers!\n";
	exit 1;
      }
      @center=@ARGV[$i..$i+1];
      $i++;
    } case (/^-f/) {
      $i++;
      $system=$ARGV[$i];
    } else {
      print "**** error: expected optional flag but got \"$ARGV[$i]\"!\n";
      exit 1;
    }
  }
}

if($factor[0][0]==$factor[1][0] and $factor[0][1]==$factor[1][1]) {
  print "error: vectors must not be equal\n";
  exit 1;
}

if($system eq "alumina0001") {
  $doo=2.75552;
  $ohperiodicvec[0][0]=$doo;
  $ohperiodicvec[0][1]=0;
  $ohperiodicvec[1][0]=$doo*cos(pi/3);
  $ohperiodicvec[1][1]=$doo*sin(pi/3);
} elsif($system eq "mag111") {
  $doo=2.96861104411750;
  $ohperiodicvec[1][0]=0;
  $ohperiodicvec[1][1]=$doo;
  $ohperiodicvec[0][0]=$doo*cos(pi/6);
  $ohperiodicvec[0][1]=$doo*sin(pi/6);
} else {
  die "error: unkown surface type!\n";
}
$minmoldist= $doo+1;

$molperiodicvec[0][0]=$factor[0][0]*$ohperiodicvec[0][0]+$factor[0][1]*$ohperiodicvec[1][0];
$molperiodicvec[0][1]=$factor[0][0]*$ohperiodicvec[0][1]+$factor[0][1]*$ohperiodicvec[1][1];
$molperiodicvec[1][0]=$factor[1][0]*$ohperiodicvec[0][0]+$factor[1][1]*$ohperiodicvec[1][0];
$molperiodicvec[1][1]=$factor[1][0]*$ohperiodicvec[0][1]+$factor[1][1]*$ohperiodicvec[1][1];

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

$mol2numox = 0;
for($a=0;$a<$mol2_numatoms[0];$a++) {
  if(uc($mol2_atomdata[0][$a][4]) eq "O" or uc($mol2_atomdata[0][$a][4]) eq "OH") {
    $mol2oxpos[$mol2numox][0]=$mol2_atomdata[0][$a][0];
    $mol2oxpos[$mol2numox][1]=$mol2_atomdata[0][$a][1];
    $mol2oxpos[$mol2numox][2]=$mol2_atomdata[0][$a][2];
    $mol2numox++;
  }
}

@mol2oxpos = sort { $a->[2] <=> $b->[2] } @mol2oxpos;
$c = sprintf("%.0f",$mol2_charge[0]);
if(not defined($bindingmode)) {
  if($c==0) {
    $bindingmode=0;
  }elsif($c==-1) {
    $bindingmode=1;
  }elsif($c==-2) {
    $bindingmode=2;
  }
}

# read FIELD file
exit 1 if(read_field_file($filenameinfld,0)!=0);
$foundhydrox = -1;
$foundmol    = -1;
$foundsolv   = -1;
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /hydroxide/i) {
    $foundhydrox=$t;
  } elsif($mol_name[0][$t] =~ /^$mol2name$/i) {
    $foundmol=$t;
  } elsif($mol_name[0][$t] =~ /$solvname/i) {
    $foundsolv=$t;
  }
}
die "error: could not find entry for hydroxide ions in FIELD file $filenameinfld!\n" if($foundhydrox==-1);
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

# find maximum z value of surface
$zsurf=-9e20;
$xmin=9e20; # to find hydroxide in lower left corner
$ymin=9e20; # to find hydroxide in lower left corner
$t=$foundhydrox;
for($m=0;$m<$mol_numents[0][$t];$m++) {
  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
    if($cdata[0][$t][$m][$a][9]=~/^O/i) {
      if($cdata[0][$t][$m][$a][2]>$zsurf) {
	$zsurf=$cdata[0][$t][$m][$a][2];
      }
      if($cdata[0][$t][$m][$a][0]<$xmin) {
	$xmin=$cdata[0][$t][$m][$a][0];
      }
      if($cdata[0][$t][$m][$a][1]<$ymin) {
	$ymin=$cdata[0][$t][$m][$a][1];
      }
    }
  }
}

# create array of candidates
$zmintol=$zsurf-$tolerance;
$zmaxtol=$zsurf+$tolerance;
$candindex=0;
for($m=0;$m<$mol_numents[0][$t];$m++) {
  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
    if(($cdata[0][$t][$m][$a][9]=~/^O/i)
	and $cdata[0][$t][$m][$a][2]>$zmintol
	and $cdata[0][$t][$m][$a][2]<$zmaxtol) {
      if($cdata[0][$t][$m][$a][0]<$xmin) {
	$xmin=$cdata[0][$t][$m][$a][0];
      }
      if($cdata[0][$t][$m][$a][2]<$ymin) {
	$ymin=$cdata[0][$t][$m][$a][2];
      }
      $candidate[$candindex][0]=$cdata[0][$t][$m][$a][0];
      $candidate[$candindex][1]=$cdata[0][$t][$m][$a][1];
      $candidate[$candindex][2]=$cdata[0][$t][$m][$a][2];
      $candidate[$candindex][3]=$t;
      $candidate[$candindex][4]=$m;
      $candindex++;
    }
  }
}

print "found ",($#candidate+1)," candiates\n";

if(not @center) {
  @center=($xmin,$ymin);
}

$borderxmin = $xmin-0.2;
$borderxmax = $xmin+2*$size[0]-0.2;
$borderymin = $ymin-0.2;
$borderymax = $ymin+2*$size[1]-0.2;

$imaxx = int(2*$size[0][0]/$doo+2)+2;
$imaxy = int(2*$size[0][1]/$ohperiodicvec[1][1])+2;
$i=0;

@remoharr=();
@finalpos=();
gridposy : for($iy=-$imaxy;$iy<=$imaxy;$iy++) {
  gridposx : for($ix=-$imaxx;$ix<=$imaxx;$ix++) {
    # determine position according to grid
    $pos[0] = $center[0]+$ix*$molperiodicvec[0][0]+$iy*$molperiodicvec[1][0]
	    +$molshiftvec[0]*$ohperiodicvec[0][0]+$molshiftvec[1]*$ohperiodicvec[1][0];
    $pos[1] = $center[1]+$ix*$molperiodicvec[0][1]+$iy*$molperiodicvec[1][1]
	    +$molshiftvec[0]*$ohperiodicvec[0][1]+$molshiftvec[1]*$ohperiodicvec[1][1];
    next if($pos[0] < -$size[0][0] or $pos[0] > $size[0][0] or $pos[1] < -$size[0][1] or $pos[1] > $size[0][1]);
    if(defined($islandrsq)) {
      $dx=$pos[0]-$islandcenter[0];
      $dy=$pos[1]-$islandcenter[1];
      next if(($dx*$dx+$dy*$dy)>$islandrsq);
    } elsif(defined($islandwidth)) {
      next if($pos[0]-$islandcenter[0] > $islandwidth);
      next if($pos[0]-$islandcenter[0] < -$islandwidth);
      next if($pos[1]-$islandcenter[1] > $islandheight);
      next if($pos[1]-$islandcenter[1] < -$islandheight);
    }
    if(defined($holersq)) {
      $dx=$pos[0]-$holecenter[0];
      $dy=$pos[1]-$holecenter[1];
      next if(($dx*$dx+$dy*$dy)<$holersq);
    }
    # find nearest hydroxide and remove it
    $dmin=9e20;
    for($j=0;$j<=$#candidate;$j++) {
      for($x=-1;$x<2;$x++) {
	for($y=-1;$y<2;$y++) {
	  $dx=$pos[0]+$mol2oxpos[0][0]-$candidate[$j][0]-$x*$cell[0][0][0]-$y*$cell[0][1][0];
	  $dy=$pos[1]+$mol2oxpos[0][1]-$candidate[$j][1]-$x*$cell[0][0][1]-$y*$cell[0][1][1];
	  $dist=sqrt($dx**2+$dy**2);
	  if($dist<$dmin) {
	    $dmin=$dist;
	    $remoh=$j;
	  }
	}
      }
    }
    next if($dmin>1.5);
    $pos[0] = $candidate[$remoh][0];
    $pos[1] = $candidate[$remoh][1];
    $pos[2] = $candidate[$remoh][2];
    push(@remoharr,splice(@candidate,$remoh,1)) if($bindingmode>0);
    # replace additional hydroxides if bindingmode>1
    for($b=1;$b<$bindingmode;$b++) {
      $dmin=9e20;
      for($j=0;$j<=$#candidate;$j++) {
	for($x=-1;$x<2;$x++) {
	  for($y=-1;$y<2;$y++) {
	    $dx=$pos[0]+$mol2oxpos[$b][0]-$candidate[$j][0]-$x*$cell[0][0][0]-$y*$cell[0][1][0];
	    $dy=$pos[1]+$mol2oxpos[$b][1]-$candidate[$j][1]-$x*$cell[0][0][1]-$y*$cell[0][1][1];
	    $dist=sqrt($dx**2+$dy**2);
	    if($dist<$dmin) {
	      $dmin=$dist;
	      $remoh=$j;
	    }
	  }
	}
      }
      $pos[2] += $candidate[$remoh][2];
      push(@remoharr,splice(@candidate,$remoh,1));
    }
    if($bindingmode>0) {
      $pos[2]/=$bindingmode;
    }
    $pos[2]+=$surfdist;
    # quickly check with molecules in neighboring cells to see if the molecule fits
    for($j=0;$j<@finalpos;$j++) {
      for($x=-1;$x<2;$x++) {
	for($y=-1;$y<2;$y++) {
	  $dx=$pos[0]-$finalpos[$j][0]-$x*$cell[0][0][0]-$y*$cell[0][1][0];
	  $dy=$pos[1]-$finalpos[$j][1]-$x*$cell[0][0][1]-$y*$cell[0][1][1];
	  $dist=sqrt($dx*$dx+$dy*$dy);
	  if($dist<3) {
	    # undo changes
	    for($a=0;$a<$bindingmode;$a++) {
	      push(@candidate,pop(@remoharr));
	    }
	    next gridposx
	  }
	}
      }
    }
    # add a molecule at the new position
    @{$finalpos[@finalpos]} = @pos;
    $m=@{$cdata[0][$foundmol]};
    @{$cdata[0][$foundmol][$m]} = mol2_to_cdata(0);
    move_mol(\@{$cdata[0][$foundmol][$m]},\@pos);
    $field_numatoms[0] += $mol_numatoms[0][$foundmol];
    $frame_numatoms[0] += $mol_numatoms[0][$foundmol];
    $mol_numents[0][$foundmol]++;
  }
}
@remoharr = sort {$a->[4] <=> $b->[4]} @remoharr;
for($i=$#remoharr;$i>=0;$i--) {
  remove_mol_entity(0,0,$remoharr[$i][3],$remoharr[$i][4]);
}
$i=$#remoharr+1;
print "removed $i hydroxide ions\n";
write_config_file($filenameoutcfg,0,$config_title[0]);
write_field_file($filenameoutfld,0);
exit 0;
