#!/usr/bin/perl

use Switch;
use dlpoly_utility;
use hanno_utility;
use Storable qw(dclone);
use Math::Trig;

if($#ARGV<5) {
  print " 1. name of CONFIG-file with unit cell\n";
  print " 2. name of corresponding FIELD file\n";
  print " 3. name of output CONFIG-file\n";
  print " 4. name of output FIELD-file\n";
  print " 5. height of nanoparticle\n";
  print " 6. inner radius of nanoparticle\n";
  print "[7. rotation angle of hexagonal prism]\n";
  print "[8. shift center of hexagonal prism]\n";
  exit 1;
}

$inangle = 0;
@center  = (0,0,0);

$incfg     = $ARGV[0];
$infld     = $ARGV[1];
$outcfg    = $ARGV[2];
$outfld    = $ARGV[3];
$inheight  = $ARGV[4];
$inradius  = $ARGV[5];
$inangle   = $ARGV[6] if($#ARGV>5);
$center[0] = $ARGV[7] if($#ARGV>6);
$center[1] = $ARGV[8] if($#ARGV>7);
$center[2] = $ARGV[9] if($#ARGV>8);

error(1) if(not check_real($inheight, $inradius, $inangle));
error(2,$infld) if(read_field_file($infld,0)    != 0);
error(2,$incfg) if(read_config_file($incfg,0,0) != 0);

$inangle  *= pi/180.0;
$inheight /= 2.0;
$outradius = $inradius*2.0/sqrt(3.0);
# generate boundary planes for nanoparticle
for($i=0;$i<6;$i++) {
  $angle = pi*$i/3.0;
  $normvec[$i][0] = cos($angle);
  $normvec[$i][1] = sin($angle);
  $normvec[$i][2] = 0;
  $planepoint[$i][0] = $normvec[$i][0]*$inradius+$center[0];
  $planepoint[$i][1] = $normvec[$i][1]*$inradius+$center[1];
  $planepoint[$i][2] = $center[2];
}
@{$normvec[6]}    = (0,0, 1);
@{$normvec[7]}    = (0,0,-1);
@{$planepoint[6]} = ($center[0],$center[1], $inheight+$center[2]);
@{$planepoint[7]} = ($center[0],$center[1],-$inheight+$center[2]);

open(ASDF,">","box.xyz");
print ASDF "8\n\n";
for($i=0;$i<@normvec;$i++) {
  print ASDF "H ",join(" ",@{$planepoint[$i]}),"   ",join(" ",@{$normvec[$i]}),"\n";
}
close(ASDF);

# center the system

@minpos = ( 1e9, 1e9, 1e9);
@maxpos = (-1e9,-1e9,-1e9);
for($t=0;$t<$field_nummols[0];$t++) {
  for($m=0;$m<$mol_numents[0][$t];$m++) {
    for($a=0;$a<$mol_numatoms[0][$t];$a++) {
      for($c=0;$c<3;$c++) {
	$minpos[$c] = $cdata[0][$t][$m][$a][$c] if($cdata[0][$t][$m][$a][$c] < $minpos[$c]);
	$maxpos[$c] = $cdata[0][$t][$m][$a][$c] if($cdata[0][$t][$m][$a][$c] > $maxpos[$c]);
      }
    }
    remap_molecule(\@{$cdata[0][$t][$m]},\@remapaxes,\@{$size[0]});
  }
}

for($c=0;$c<3;$c++) {
  $dist = ($minpos[$c]+$maxpos[$c])/2;
  for($t=0;$t<$field_nummols[0];$t++) {
    for($m=0;$m<$mol_numents[0][$t];$m++) {
      for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	$cdata[0][$t][$m][$a][$c]-=$dist;
      }
    }
  }
}


# enlarge system

$max[0] = int($outradius/$size[0][1])+1;
$max[1] = int($outradius/$size[0][2])+1;
$max[2] = int($inheight/$size[0][2])+1;
$periodic_key[1] = 0;
copy_field(0,1);
@failed=(0,0,0,0,0,0,0,0);
for($t=0;$t<$field_nummols[0];$t++) {
  @{$cdata[1][$t]}=();
  $mol_numents[1][$t]=0;
  $newmolid=0;
  for($m=0;$m<$mol_numents[0][$t];$m++) {
    for($x=-$max[0];$x<=$max[0];$x++) {
      for($y=-$max[1];$y<=$max[1];$y++) {
	loopmol:for($z=-$max[2];$z<=$max[2];$z++) {
	  @{$cdata[1][$t][$newmolid]} = @{dclone \@{$cdata[0][$t][$m]}};
	  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	    $cdata[1][$t][$newmolid][$a][0] += $x*$cell[0][0][0]+$y*$cell[0][1][0]+$z*$cell[0][2][0];
	    $cdata[1][$t][$newmolid][$a][1] += $x*$cell[0][0][1]+$y*$cell[0][1][1]+$z*$cell[0][2][1];
	    $cdata[1][$t][$newmolid][$a][2] += $x*$cell[0][0][2]+$y*$cell[0][1][2]+$z*$cell[0][2][2];
	  }
	  $mol_numents[1][$t]++;
	  # check boundaries
	  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	    for($i=0;$i<@normvec;$i++) {
	      $diffvec[0] = $cdata[1][$t][$newmolid][$a][0] - $planepoint[$i][0];
	      $diffvec[1] = $cdata[1][$t][$newmolid][$a][1] - $planepoint[$i][1];
	      $diffvec[2] = $cdata[1][$t][$newmolid][$a][2] - $planepoint[$i][2];
	      $dx   = $diffvec[0]*$normvec[$i][0];
	      $dy   = $diffvec[1]*$normvec[$i][1];
	      $dz   = $diffvec[2]*$normvec[$i][2];
	      $dist = $dx+$dy+$dz;
# 	      print "$newmolid $dx $dy $dz\n";
	      if($dist>0) {
		$failed[$i]++;
		remove_mol_entity(1,1,$t,$newmolid);
		next loopmol;
	      }
	    }
	  }
	  $newmolid++;
	}
      }
    }
  }
}
print join(" ",@failed),"\n";
write_field_file($outfld,1);
write_config_file($outcfg,1,"alox nanoparticle",">",1);

sub error {
  if($_[0]==1) {
    print "**** error: height and radius must be real numbers!\n";
  } elsif($_[0]==2) {
    print "**** error: could not open file $_[1]!\n";
  } elsif($_[0]==4) {
    print "**** error: could not open output file!\n";
  } else {
    print "**** unknown error!\n";
    exit 99;
  }
  exit $_[0]
}
