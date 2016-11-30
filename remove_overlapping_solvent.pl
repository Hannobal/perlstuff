#!/usr/bin/perl
use Math::Trig;
use Switch;

use dlpoly_utility;
use hanno_utility;

#use Storable qw(dclone);

if($#ARGV<5) {
  print "input format:\n";
  print "1. name of CONFIG/REVCON-file with surface\n";
  print "2. name of corresponding FIELD-file\n";
  print "3. name of output CONFIG-file\n";
  print "4. name of output FIELD-file\n";
  print "5. name of solvent molecule\n";
  print "6. minimum distance between atoms\n";
  print "[list of molecule names to check]\n";
  exit 1;
}

$filenameincfg  = $ARGV[0];
$filenameinfld  = $ARGV[1];
$filenameoutcfg = $ARGV[2];
$filenameoutfld = $ARGV[3];
$solvname       = $ARGV[4];
$dmin           = $ARGV[5];
@includetype    = @ARGV[6..$#ARGV];
die "**** error: minimum distance must be a real number!\n" if(not check_real($dmin));
$dminsq         = $dmin*$dmin;

# read FIELD file
exit 1 if(read_field_file($filenameinfld,0)!=0);
$foundsolv   = -1;
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /$solvname/i) {
    $s=$t;
  }
}
die "**** error: could not find entry for solvent molecule in FIELD file $filenameinfld!\n"
if($s==-1 and $lsolv);

for($i=0;$i<@includetype;$i++) {
  $found=0;
  if(not contains(@{$mol_name[0]},$includetype[$i])) {
    die "**** error: could not find entry for molecule named $includetype[$i] ",
    "in FIELD file $filenameinfld!\n";
  }
}

#read CONFIG file
exit 1 if(read_config_file($filenameincfg,0,0)!=0);

# remove solvent
switch($periodic_key[0]) {
  case [0] {
    for($t=0;$t<$field_nummols[0];$t++) {
      next if ($t==$s);
      if(@includetype) {
	next if(not contains(@includetype,$mol_name[0][$t]));
      }
      refmol : for($m=0;$m<$mol_numents[0][$t];$m++) {
	refat : for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	  solvmol : for($n=0;$n<$mol_numents[0][$s];$n++) {
	    solvat : for($b=0;$b<$mol_numatoms[0][$s];$b++) {
	      $dx=$cdata[0][$t][$m][$a][0]-$cdata[0][$s][$n][$b][0];
	      $dy=$cdata[0][$t][$m][$a][1]-$cdata[0][$s][$n][$b][1];
	      $dz=$cdata[0][$t][$m][$a][2]-$cdata[0][$s][$n][$b][2];
	      $dsq = $dx*$dx+$dy*$dy+$dz*$dz;
	      if($dsq<$dminsq) {
		remove_mol_entity(0,0,$s,$n);
		$n--;
		next solvmol;
	      }
	    }
	  }
	}
      }
    }
  } case [1,2] {
    for($t=0;$t<$field_nummols[0];$t++) {
      next if ($t==$s);
      if(@includetype) {
	next if(not contains(@includetype,$mol_name[0][$t]));
      }
      refmol : for($m=0;$m<$mol_numents[0][$t];$m++) {
	refat : for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	  solvmol : for($n=0;$n<$mol_numents[0][$s];$n++) {
	    solvat : for($b=0;$b<$mol_numatoms[0][$s];$b++) {
	      for($x=-1;$x<2;$x++) {
		$sx = $x*$cell[0][0][0];
		$dx=$cdata[0][$t][$m][$a][0]-$cdata[0][$s][$n][$b][0]-$sx;
		for($y=-1;$y<2;$y++) {
		  $sy = $y*$cell[0][1][1];
		  $dy=$cdata[0][$t][$m][$a][1]-$cdata[0][$s][$n][$b][1]-$sy;
		  for($z=-1;$z<2;$z++) {
		    $sz = $z*$cell[0][2][2];
		    $dz=$cdata[0][$t][$m][$a][2]-$cdata[0][$s][$n][$b][2]-$sz;
		    $dsq = $dx*$dx+$dy*$dy+$dz*$dz;
		    if($dsq<$dminsq) {
		      remove_mol_entity(0,0,$s,$n);
		      $n--;
		      next solvmol;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  } case [6] {
    for($t=0;$t<$field_nummols[0];$t++) {
      next if ($t==$s);
      if(@includetype) {
	next if(not contains(@includetype,$mol_name[0][$t]));
      }
      print "analyzing $mol_name[0][$t]\n";
      refmol : for($m=0;$m<$mol_numents[0][$t];$m++) {
	refat : for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	  solvmol : for($n=0;$n<$mol_numents[0][$s];$n++) {
	    solvat : for($b=0;$b<$mol_numatoms[0][$s];$b++) {
	      $dz=$cdata[0][$t][$m][$a][2]-$cdata[0][$s][$n][$b][2];
	      next if($dz*$dz>$dminsq);
	      for($x=-1;$x<2;$x++) {
		$sx = $x*$cell[0][0][0];
		$dx=$cdata[0][$t][$m][$a][0]-$cdata[0][$s][$n][$b][0]-$sx;
		for($y=-1;$y<2;$y++) {
		  $sy = $y*$cell[0][1][1];
		  $dy=$cdata[0][$t][$m][$a][1]-$cdata[0][$s][$n][$b][1]-$sy;
		  $dsq = $dx*$dx+$dy*$dy+$dz*$dz;
		  if($dsq<$dminsq) {
		    remove_mol_entity(0,0,$s,$n);
		    $n--;
		    next solvmol;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  } else {
    print "sorry, this script is only implemented for orthogonal cells\n";
    exit 1;
  }
}

write_config_file($filenameoutcfg,0,$config_title[0]);
write_field_file($filenameoutfld,0);
exit 0;
