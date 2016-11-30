#!/usr/bin/perl

use Math::Trig;
use dlpoly_utility;
use hanno_utility;

if($#ARGV<7) {
  print "input format:\n";
  print "1. name of input CONFIG file\n";
  print "2. name of input FIELD file\n";
  print "3. name of output CONFIG file\n";
  print "4. name of output FIELD file\n";
  print "5. name of molecule to replace\n";
  print "6. name of mol2 file with new molecule\n";
  print "7. percentage of molecules to replace (e.g. 10)\n";
  print "8. minimum distance between two atoms\n";
  exit 1;
}

$incfg      = $ARGV[0];
$infld      = $ARGV[1];
$outcfg     = $ARGV[2];
$outfld     = $ARGV[3];
$replacemol = $ARGV[4];
$inmol2     = $ARGV[5];
$percentage = $ARGV[6];
$mindistsq  = $ARGV[7]*$ARGV[7];
@unitvector = ([1,0,0],[0,1,0],[0,0,1]);

exit 1 if(read_field_file($infld,0)!=0);
exit 1 if(read_config_file($incfg,0,0)!=0);
exit 1 if(read_mol2_file($inmol2,0)!=0);

$ta = -1;
$tb = -1;
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] eq $replacemol) {
    $ta=$t
  }elsif($mol_name[0][$t] eq $mol2_name[0]) {
    $tb=$t
  }
}

if($tb<0 or $ta<0) {
  print "**** error: did not find entry for molecules in input FIELD file\n";
  exit 1;
}

$numreplace = int($mol_numents[0][$ta]*$percentage/100);

print "replacing $numreplace molecules\n";
for($i=0;$i<$numreplace;$i++) {
  print "replacing moleucle $i of $numreplace\n";
  # randomly select molecule and remap if necessary
  $m=int(rand($mol_numents[0][$ta]));
  if($periodic_key[0]==6) {
    remap_molecule(\@{$cdata[0][$ta][$m]},[0,1],\@{$size[0]});
  } elsif($periodic_key[0]>0) {
    remap_molecule(\@{$cdata[0][$ta][$m]},[0,1,2],\@{$size[0]});
  }
  # determine molecule position
  for($a=0;$a<$mol_numatoms[0][$ta];$a++) {
    $pos[0] += $cdata[0][$ta][$m][$a][0];
    $pos[1] += $cdata[0][$ta][$m][$a][1];
    $pos[2] += $cdata[0][$ta][$m][$a][2];
  }
  $pos[0] /= $mol_numatoms[0][$ta];
  $pos[1] /= $mol_numatoms[0][$ta];
  $pos[2] /= $mol_numatoms[0][$ta];
  # remove the molecule
  @removedmol=splice(@{$cdata[0][$ta]},$m,1);
  $mol_numents[0][$ta]--;
  # add new molecule
  $mb=$mol_numents[0][$tb];
  @{$cdata[0][$tb][$mb]} = mol2_to_cdata(0);
  $mol_numents[0][$tb]++;
  # move molecule to position
  move_mol(\@{$cdata[0][$tb][$mb]}, \@pos);
  # randomly rotate new molecule until it fits
  while($fit==0) {
    for($c=0;$c<3;$c++) {
      $angle = rand(2)*pi;
      @rotmatrix = gen_rot_matrix($unitvector[$c], $angle);
      rotate_molecule(\@{$cdata[0][$tb][$mb]}, \@rotmatrix);
    }
    $fit=1;
    # check fit
    if($periodic_key[0]==6) {
      checkloop1:for($t=0;$t<$field_nummols[0];$t++) {
	for($m=0;$m<$mol_numents[0][$t];$m++) {
	  next if($t==$tb and $m==$mb);
	  for($ab=0;$ab<$mol_numatoms[0][$tb];$ab++) {
	    for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	      for($x=0;$x<2;$x++) {
		for($y=0;$y<2;$y++) {
		  $pervec[0] = $x*$cell[0][0][0]+$y*$cell[0][1][0]+$z*$cell[0][2][0];
		  $pervec[1] = $x*$cell[0][0][1]+$y*$cell[0][1][1]+$z*$cell[0][2][1];
		  $dx = $cdata[0][$tb][$mb][$ab][0]-$cdata[0][$t][$m][$a][0]+$pervec[0];
		  $dy = $cdata[0][$tb][$mb][$ab][1]-$cdata[0][$t][$m][$a][1]+$pervec[1];
		  $dz = $cdata[0][$tb][$mb][$ab][2]-$cdata[0][$t][$m][$a][2]+$pervec[2];
		  $distsq = $dx*$dx+$dy*$dy+$dz*$dz;
		  if($distsq<$mindistsq) {
		    $fit = 0;
		    print "$t $m $a $distsq\n";
		    last checkloop1;
		  } # end if $distsq
		} # end for $y
	      } # end for $x
	    } # end for $a
	  } # end for $ab
	} # end for $m
      } # end for $t
    } elsif($periodic_key[0]>0) {
      checkloop2:for($t=0;$t<$field_nummols[0];$t++) {
	for($m=0;$m<$mol_numents[0][$t];$m++) {
	  next if($t==$tb and $m==$mb);
	  for($x=0;$x<2;$x++) {
	    for($y=0;$y<2;$y++) {
	      for($z=0;$z<2;$z++) {
		$pervec[0] = $x*$cell[0][0][0]+$y*$cell[0][1][0]+$z*$cell[0][2][0];
		$pervec[1] = $x*$cell[0][0][1]+$y*$cell[0][1][1]+$z*$cell[0][2][1];
		$pervec[2] = $x*$cell[0][0][2]+$y*$cell[0][1][2]+$z*$cell[0][2][2];
		for($ab=0;$ab<$mol_numatoms[0][$tb];$ab++) {
		  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
		    $dx = $cdata[0][$tb][$mb][$ab][0]-$cdata[0][$t][$m][$a][0]+$pervec[0];
		    $dy = $cdata[0][$tb][$mb][$ab][1]-$cdata[0][$t][$m][$a][1]+$pervec[1];
		    $dz = $cdata[0][$tb][$mb][$ab][2]-$cdata[0][$t][$m][$a][2]+$pervec[2];
		    $distsq = $dx*$dx+$dy*$dy+$dz*$dz;
		    if($distsq<$mindistsq) {
		      $fit = 0;
		      print "$distsq\n";
		      last checkloop2;
		    } # end if $distsq
		  } # end for $a
		} # end for $ab
	      } # end for $z
	    } # end for $y
	  } # end for $x
	} # end for $m
      } # end for $t
    } else {
      checkloop3:for($t=0;$t<$field_nummols[0];$t++) {
	for($m=0;$m<$mol_numents[0][$t];$m++) {
	  next if($t==$tb and $m==$mb);
	  for($ab=0;$ab<$mol_numatoms[0][$tb];$ab++) {
	    for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	      $dx = $cdata[0][$tb][$mb][$ab][0]-$cdata[0][$t][$m][$a][0];
	      $dy = $cdata[0][$tb][$mb][$ab][1]-$cdata[0][$t][$m][$a][1];
	      $dz = $cdata[0][$tb][$mb][$ab][2]-$cdata[0][$t][$m][$a][2];
	      $distsq = $dx*$dx+$dy*$dy+$dz*$dz;
	      if($distsq<$mindistsq) {
		$fit = 0;
		print "$distsq\n";
		last checkloop3;
	      } # end if $distsq
	    } # end for $a
	  } # end for $ab
	} # end for $m
      } # end for $t
    }
  } # end while
}

write_field_file($outfld,0);
write_config_file($outcfg,0,$config_title[0],">",0);