#!/usr/bin/perl

# use strict;
use dlpoly_utility;
use hanno_utility;
use Math::Trig;
use List::Util 'shuffle';
use Storable qw(dclone);

$tolerance = 0.2;
$surfdist  = 4;
$maxdistsq = 1;
$maxiter   = 10;
$force     = 5000;

if($#ARGV<4) {
  print "input format:\n";
  print "1. CONFIG file\n";
  print "2. FIELD file\n";
  print "3. minimum distance between two binding sites\n";
  print "4. list of molecule names and ratio\n";
  print "   example: C10-PA:C14-PA:C60-C18-PA 1:1:3\n";
  exit 1;
}

$incfgname=$ARGV[0];
$infldname=$ARGV[1];
$mindistsq=$ARGV[2]*$ARGV[2];
@molnames=split(':',$ARGV[3]);
@molratio=split(':',$ARGV[4]);
if($mindistsq==0) {
  print STDERR "**** error: minimum distance must be a real number != 0\n";
  exit 1;
}
if(@molnames != @molratio) {
  print STDERR "**** error: number of molecule ratios (",$#molratio+1,") does not equal ";
  print "number of molecule names (",$#molnames,")\n";
  exit 1;
}

$err=0;
$err=1 if(read_field_file($infldname,0)!=0);
$err=1 if(read_config_file($incfgname,0,0)!=0);
exit 1 if $err;

for($i=0;$i<@molnames;$i++) {
  $err=1 if(read_mol2_file("$molnames[$i].mol2",$i)!=0);
  # determine molecule boundaries
  @{$mol2_minpos[$i]} = ( 9e20,  9e20,  9e20);
  @{$mol2_maxpos[$i]} = (-9e20, -9e20, -9e20);
  for($a=0;$a<$mol2_numatoms[$i];$a++) {
    for($c=0;$c<3;$c++) {
      $mol2_maxpos[$i][$c]=$mol2_atomdata[$i][$a][$c] if($mol2_atomdata[$i][$a][$c]>$mol2_maxpos[$i][$c]);
      $mol2_minpos[$i][$c]=$mol2_atomdata[$i][$a][$c] if($mol2_atomdata[$i][$a][$c]<$mol2_minpos[$i][$c]);
    }
  }
  $tm[$i]=-1;
  for($t=0;$t<$field_nummols[0];$t++) {
    if($mol_name[0][$t] =~ /$molnames[$i]/ and not ($mol_name[0][$t] =~ /frozen/i or $mol_name[0][$t] =~ /fix/i)) {
      $tm[$i]=$t;
      last;
    }
  }
  if($tm[$i]<0) {
    print STDERR "**** error: molecule named $molnames[$i] was not found in FIELD file!\n";
    $err=1;
  }
}
$toh=-1;
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /hydroxide/i) {
    $toh=$t;
  }
}
if($toh<0) {
  print STDERR "**** error: entry for hydroxide ions was not found in FIELD file!\n";
  $err=1;
}
exit 1 if $err;

# build array for adding new molecules
for($i=0;$i<@molnames;$i++) {
  for($j=0;$j<$molratio[$i];$j++) {
    push(@addmolorder,$i)
  }
}
@addmolorder = shuffle(@addmolorder);
print STDOUT "order for adding molecules:\n@addmolorder\n";


# determine surface z-position
$zsurf=-9e20;
for($m=0;$m<$mol_numents[0][$toh];$m++) {
  $zsurf=$cdata[0][$toh][$m][0][2] if($zsurf<$cdata[0][$toh][$m][0][2]);
}
$zmintol=$zsurf-$tolerance;
$zmaxtol=$zsurf+$tolerance;

print STDOUT "considering Hydroxides between z=$zmintol and $zmaxtol\n";

for($m=0;$m<$mol_numents[0][$toh];$m++) {
  if($cdata[0][$toh][$m][0][2]>$zmintol and $cdata[0][$toh][$m][0][2]<$zmaxtol) {
    push(@candidates,$m);
  }
}

print STDOUT "found ",$#candidates+1," candidates\n";

$nextmol=0;  # the index for the addmolorder-array
$success=1; # whether the iteration over the previous candidate worked
$canditer=0;
loopcand : while(@candidates) {
  $canditer++;
  $cdirname = sprintf("%04d",$canditer);
  mkdir $cdirname;
  $nextmol=($nextmol+1)%@addmolorder if($success);
  $n=$addmolorder[$nextmol];
  $c=rand(@candidates);
  $m=$candidates[$c];
  splice(@candidates,$c,1);
  open(CANDS,">","candidates.xyz");
  print CANDS $#candidates+2,"\n";
  for($i=0;$i<@candidates;$i++) {
    printf CANDS "\nH %10.5f %10.5f %10.5f %5d", @{$cdata[0][$toh][$candidates[$i]][0]}[0..2],$candidates[$i];
  }
    printf CANDS "\nH %10.5f %10.5f %10.5f", @currpos;
  close(CANDS);
  @currpos = @{$cdata[0][$toh][$m][0]}[0..2];
  @curroh  = splice(@{$cdata[0][$toh]},$m,1); # ! now values in @candidates are wrong !
  $mol_numents[0][$toh]--;
  print STDOUT "trying to replace OH at position @currpos with $molnames[$n]\n";
  # determine maximum z-position
  $zmax=-9e20;
  for($t=0;$t<$field_nummols[0];$t++) {
    for($m2=0;$m2<$mol_numents[0][$t];$m2++) {
      for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	$zmax=$cdata[0][$t][$m2][$a][2] if($zmax<$cdata[0][$t][$m2][$a][2]);
      }
    }
  }
  # add new molecule
  if($success) {
    unshift(@{$cdata[0][$tm[$n]]},undef);
    $mol_numents[0][$tm[$n]]++;
  }
  $iter=0;
  $success=0;
  loopiter : while($iter<$maxiter) {
    $iter++;
    $idirname = sprintf("%03d",$iter);
    print STDOUT "begin iteration $iter\n";
    @{$cdata[0][$tm[$n]][0]} = mol2_to_cdata($n);
    # rotate around z-axis
    @rotmatrix = &gen_rot_matrix([0,0,1], rand(2*pi));
    rotate_molecule(\@{$cdata[0][$tm[$n]][0]}, \@rotmatrix);
    # shift to random position
    $shiftz = $zmax + $surfdist - $mol2_minpos[2];
    @shiftvec = ($currpos[0],$currpos[1],$shiftz);
    move_mol(\@{$cdata[0][$tm[$n]][0]},\@shiftvec);
    # add distance restraint (+1.5 for oxygen-phosphorus distance)
    $index = 11;
    for($t=0;$t<$tm[$n];$t++) {
      $index += $mol_numents[0][$t]*$mol_numatoms[0][$t];
    }
    $field_numextern[0]=1;
    @{$field_externdata[0][0]} = ('fcon',$index,@currpos[0..1],$currpos[2]+1.5,$force);
    # write dlpoly input files
    write_config_file('CONFIG',0,$config_title[0]);
    write_field_file('FIELD',0);
    ###### run dlpoly ######
    # $test = `mpirun.openmpi /raid/home/dietrich/bin/dl_class_1.9_2dnpt/execute/DLPOLY.X`;
    $test = `mpirun.openmpi -n 3 /home/dietrich/work/dl_class_1.9_2dnpt/execute/DLPOLY.X`;
    # analyze the output
    exit 1 if(read_config_file('REVCON',0,1)!=0);
    mkdir "$cdirname/$idirname";
    `mv candidates.xyz FIELD CONFIG HISTORY REVCON REVIVE STATIS RDFDAT ZDNDAT CFGMIN $cdirname/$idirname`;
    loopchecksuccess : for($a=0;$a<$mol_numatoms[0][$tm[$n]];$a++) {
      if($mol_atomdata[0][$tm[$n]][$a][0] eq 'O') {
	if($mol_atomdata[0][$tm[$n]][$mol_bondatoms[0][$tm[$n]][$a][0]][0] eq 'P5') {
	  for($x=-1;$x<2;$x++) {
	    for($y=-1;$y<2;$y++) {
	      $dx = $currpos[0] - $cdata[1][$tm[$n]][0][$a][0] + $x*$cell[0][0][0]+$y*$cell[0][0][1];
	      $dy = $currpos[1] - $cdata[1][$tm[$n]][0][$a][1] + $x*$cell[0][1][0]+$y*$cell[0][1][1];
	      $dz = $currpos[2] - $cdata[1][$tm[$n]][0][$a][2];
	      $distsq = $dx*$dx + $dy*$dy + $dz*$dz;
	      if($distsq < $maxdistsq) {
		$success=1;
		print STDOUT "success!\n\n";
		last loopchecksuccess;
	      } # end distance check
	    } # end for $y
	  } # end for $x
	} # end if bonded to P5
      } # end if 'O'
    } # end for $a
    if($success) {
      # copy REVCON data to input
      exit 1 if(read_config_file("$cdirname/$idirname/REVCON",0,0)!=0);
      # correct indices of @candidates array and remove obsolete ones
      for($i=0;$i<@candidates;$i++) {
	$candidates[$i]-- if($candidates[$i]>$m);
	print STDERR "**** error: $cdata[0][$toh][$candidates[$i]][0][9] (cand $i, $candidate[$i]) is no oxygen!\n"
	    if($cdata[0][$toh][$candidates[$i]][0][9] ne 'OX');
	for($x=-1;$x<2;$x++) {
	  for($y=-1;$y<2;$y++) {
	    $dx=$cdata[0][$toh][$candidates[$i]][0][0] - $currpos[0] + $x*$cell[0][0][0]+$y*$cell[0][0][1];
	    $dy=$cdata[0][$toh][$candidates[$i]][0][1] - $currpos[1] + $x*$cell[0][1][0]+$y*$cell[0][1][1];
	    $dz=$cdata[0][$toh][$candidates[$i]][0][2] - $currpos[2];
	    $distsq = $dx*$dx + $dy*$dy + $dz*$dz;
	    if($distsq<12) {
	      print "removing candidate at @{$cdata[0][$toh][$candidates[$i]][0]}[0..2]\n";
	      splice(@candidates,$i,1);
	      $i--;
	    }
	  } # end for $y
	} # end for $x
      } # end for @candidates
      next loopcand;
    }
  } # end loopiter
  print STDOUT "no success with this candidate! removing it from list\n";
  # reinsert hydroxide
  splice(@{$cdata[0][$toh]},$m,0,@curroh);
}

exit 0;
