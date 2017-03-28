#!/usr/bin/perl
use Math::Trig;
use Switch;

use hanno_utility;
use dlpoly_utility;
use aloxsam_utility;
use Storable qw(dclone);

# use Storable qw(dclone);

if($#ARGV<3) {
  print "input format:\n";
  print " 1. name of CONFIG/REVCON-file with surface\n";
  print " 2. name of corresponding FIELD-file\n";
  print " 3. name of output CONFIG-file\n";
  print " 4. name of output FIELD-file\n";
  print " 5. minimum distance between PAs (recommended: 3.5)\n";
  print "[6. remove water? (y/n)] (default: n)\n";
  print "[7. double deprotonation(y/n)] (default: n)\n";
  print "[8. drag force strength (default=0)]\n";
  print "[9. tolerance for end position (default=-1)]\n";
  print "[10. lattice vectors (4*int) only for single deprot.)]\n";
  print "[11. lattice center (2*real) only for single deprot.)]\n";
  exit 99;
}

$filenameincfg  = $ARGV[0];
$filenameinfld  = $ARGV[1];
$filenameoutcfg = $ARGV[2];
$filenameoutfld = $ARGV[3];
$minpadist      = $ARGV[4];
$maxohdist      = 4;
if(not check_real($minpadist)) {
  print "**** error: minimum distance between PAs must be a real number!\n";
  exit 99;
}

if($#ARGV>4) {
  $lohrem       = $ARGV[5];
  if($lohrem =~ /^y/i) {
    $lohrem=1;
  } elsif($lohrem =~ /^n/i) {
    $lohrem=0;
  } else {
    print "**** error: could not interpret input for 'remove water?'\n";
    exit 99;
  }
} else {
  $lohrem = 0;
}

if($#ARGV>5) {
  $ldouble       = $ARGV[6];
  if($ldouble =~ /^y/i) {
    $ldouble=1;
  } elsif($ldouble =~ /^n/i) {
    $ldouble=0;
  } else {
    print "**** error: could not interpret input for 'double deprotonation?'\n";
    exit 99;
  }
} else {
  $ldouble = 0;
}

if($#ARGV>6) {
  $force = $ARGV[7];
  if(check_real($force)) {
    if(not $lohrem or $ldouble) {
      print "**** error: double deprotonation must be switched off and water removal on for drag force!\n";
      exit 99;
    }
    $ldrag=1;
  } else {
    print "**** error: could not interpret input for drag force\n";
    exit 99;
  }
} else {
  $ldrag = 0;
}

if($#ARGV>7) {
  $tol = $ARGV[8];
  if(not check_real($tol)) {
    print "**** error: could not interpret input for tolerance\n";
    exit 99;
  }
} else {
  $tol=-1;
}

@gridvec=();
if($#ARGV>8) {
  if($#ARGV<12) {
    print "**** error: lattice vectors must be four numbers!\n";
    exit 99;
  }
  @gridvec = @ARGV[9..12];
  if(not check_real(@gridvec)) {
    print "**** error: could not interpret input for lattice vectors\n";
    exit 99;
  }
}
@gridcenter=(0,0);
if($#ARGV>12) {
  if($#ARGV<14) {
    print "**** error: lattice center must be two numbers!\n";
    exit 99;
  }
  @gridcenter = @ARGV[13..14];
  if(not check_real(@gridcenter)) {
    print "**** error: could not interpret input for lattice center\n";
    exit 99;
  }
}

if($ldouble and $mindistsq>0) {
  print "**** warning: double deprotonation might not occur due to minimum PA-PA distance!\n";
}
if($ldouble and @gridvec) {
  print "**** error: lattice vectors may only be specified in single deprotonation mode!\n";
  exit 99;
}

# read FIELD file
exit 99 if(read_field_file($filenameinfld,0)!=0);

#read CONFIG file
exit 99 if(read_config_file($filenameincfg,0,0)!=0);

$toh      = -1;
$twat     = -1;
@foundmol = ();
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /water/i) {
    $twat = $t;
  } elsif($mol_name[0][$t] =~ /hydroxide/i) {
    $toh = $t;
  } elsif($mol_name[0][$t] =~ /-PA$/) {
    push(@foundmol,[$t,-1,-1]);
  }
}

if(@gridvec){
  @candidates = find_hydroxide_candidates(0,0,1,$minpadist,4,@gridvec,@gridcenter);
} else {
  @candidates = find_hydroxide_candidates(0,0,1,$minpadist,4);
}

if(not @candidates) {
  print "**** error: no candidates were found!\n";
  exit 99;
}

if($twat<0 and not $lohrem) {
  print "**** error: did not find entry for water in FIELD file\n";
  $err=1;
}
if($toh<0) {
  print "**** error: did not find entry for hydroxide ions in FIELD file\n";
  $err=1;
}

if(not @foundmol) {
  print "**** error: no PAs were found!\n";
  $err=1;
}

exit 99 if($err>0);

# find corresponding deprotonated forms of PAs and atom indices for oxygen
for($i=0;$i<@foundmol;$i++) {
  for($t=0;$t<$field_nummols[0];$t++) {
    if($mol_name[0][$t] =~ /$mol_name[0][$foundmol[$i][0]]-d/) {
      $foundmol[$i][1] = $t;
    }
    if($mol_name[0][$t] =~ /$mol_name[0][$foundmol[$i][0]]-2d/) {
      $foundmol[$i][2] = $t;
    }
  }
  if($foundmol[$i][1]<0) {
    print "**** error: could not find singly deprotonated form of $mol_name[0][$foundmol[$i][0]]\n";
    $err=1;
    next;
  }
  if($foundmol[$i][2]<0 and $ldouble) {
    print "**** error: could not find doubly deprotonated form of $mol_name[0][$foundmol[$i][0]]\n";
    $err=1;
    next;
  }
  $t   = $foundmol[$i][0];
  $td  = $foundmol[$i][1];
  $t2d = $foundmol[$i][2] if($ldouble);
  # find acidic hydrogens for PA
  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
    if($mol_atomdata[0][$t][$a][0] eq 'P5') {
      for $b (@{$mol_bondatoms[0][$t][$a]}) {
	for $c (@{$mol_bondatoms[0][$t][$b]}) {
	  if($mol_atomdata[0][$t][$c][0]=~/(HG|HO)/) {
	    push(@{$pah[$i]},$c);
	  }
	}
      }
    }
  }
}

exit 99 if($err>0);

$zsurf=-9e20;
$hydroxo=-1;
$hydroxh=-1;
for($a=0;$a<$mol_numatoms[0][$toh];$a++) {
#   print "$mol_atomdata[0][$toh][$a][0]\n";
  $hydroxo=$a if($mol_atomdata[0][$toh][$a][0] =~ /^O/i);
  $hydroxh=$a if($mol_atomdata[0][$toh][$a][0] =~ /^H/i);
}
for($m=0;$m<$mol_numents[0][$toh];$m++) {
  $zsurf=$cdata[0][$toh][$m][$hydroxo][2] if($zsurf<$cdata[0][$toh][$m][$hydroxo][2]);
}

# do the replacements

$changed=0;
$err=0;
for($i=0;$i<@foundmol;$i++) {
  $t   = $foundmol[$i][0];
  $td  = $foundmol[$i][1];
  $t2d = $foundmol[$i][2] if($ldouble);
  ############ first deprotonation ##########################
  $mold=-1;
  loopmol : for($m=0;$m<$mol_numents[0][$t];$m++) {
    $mindistsq=9e20;
    $mold++;
    loopatoms : foreach $a (@{$pah[$i]}) {
      # quickly check the distance
      if($cdata[0][$t][$m][$a][2]>$zsurf+$maxohdist+3) {
        next loopatoms;
      }
      #find the nearest hydroxide oxygen
      for($x=-1;$x<2;$x++) {
        for($y=-1;$y<2;$y++) {
          $shiftx = $x*$cell[0][0][0] + $y*$cell[0][1][0];
          $shifty = $x*$cell[0][0][1] + $y*$cell[0][1][1];
          for($j=0;$j<@candidates;$j++) {
            $n  = $candidates[$j];
            $dx = $cdata[0][$t][$m][$a][0] - $cdata[0][$toh][$n][$hydroxo][0] + $shiftx;
            $dy = $cdata[0][$t][$m][$a][1] - $cdata[0][$toh][$n][$hydroxo][1] + $shifty;
            $dz = $cdata[0][$t][$m][$a][2] - $cdata[0][$toh][$n][$hydroxo][2];
            $distsq = $dx*$dx + $dy*$dy + $dz*$dz;
            if($distsq<$mindistsq) {
              $mindistsq = $distsq;
              $nearest   = $j;
              $ah        = $a;
            }
          } # end for $y
        } # end for $x
      } # end for $n
    } #end for $a
    if(sqrt($mindistsq)>$maxohdist) {
      print "**** error: no suitable hydroxide found in vicinity of $mol_name[0][$t] number $mold!\n";
      next;
    }
    $changed++;
    @{$water[0]}=@{dclone \@{$cdata[0][$toh][$candidates[$nearest]]}};
    push(@remohs,splice(@candidates,$nearest,1));
    push(@{$water[0]},splice(@{$cdata[0][$t][$m]},$ah,1)); # shift proton
    if($lohrem==0) {
      $water[0][0][9] = 'OW'; # change atom type of water oxygen
      push(@{$cdata[0][$twat]},splice(@water,0,1)); # add water to cdata array
      $last=$#{$cdata[0][$twat]};
      $mol_numents[0][$twat]++;
    }
    # change the name of the oxygen
    $o = $mol_bondatoms[0][$t][$ah][0];
    $cdata[0][$t][$m][$o][9] = 'O';
    # sort the oxygens
    if($mol_atomdata[0][$td][$o][0] ne 'O') {
      if($mol_name[0][$t] eq 'C60-C18-PA') {
	# exchange positions 11 & 13
	splice(@{$cdata[0][$t][$m]},11,0,splice(@{$cdata[0][$t][$m]},13,1));
      } else {
	# exchange positions 12 & 13
	splice(@{$cdata[0][$t][$m]},12,0,splice(@{$cdata[0][$t][$m]},13,1));
      }
      $o=13;
    }
    push(@{$cdata[0][$td]},splice(@{$cdata[0][$t]},$m,1)); # move molecule to deprotonated ones
    $mol_numents[0][$t]--;
    $mol_numents[0][$td]++;
    $m--;
    ##### add drag force
    if($ldrag) {
      $field_numextern[0]=1;
      # determine atom index
      $index=$o+1;
      for($t2=0;$t2<$td;$t2++) {
	$index+= $mol_numents[0][$t2]*$mol_numatoms[0][$t2]
      }
      $index+=($mol_numents[0][$td]-1)*$mol_numatoms[0][$td];
      @{$field_externdata[0][0]} = ('fcon 6',$index,@{$water[0][0]}[0..2],$tol,$force);
    }
  } #end for $m
  ############ second deprotonation ##########################
  next if(not $ldouble);
  for($m=0;$m<$mol_numents[0][$td];$m++) {
    $mindistsq=9e20;
    for($a=0;$a<$mol_numatoms[0][$td];$a++) {
      if($cdata[0][$td][$m][$a][9] =~/HG/) {
	# quickly check the distance
	if($cdata[0][$td][$m][$a][2]>$zsurf+$maxohdist) {
# 	  print "**** warning: proton is too far from surface!\n";
	  next;
	}
	#find the nearest hydroxide oxygen
	for($x=-1;$x<2;$x++) {
	  for($y=-1;$y<2;$y++) {
	    $shiftx = $x*$cell[0][0][0] + $y*$cell[0][1][0];
	    $shifty = $x*$cell[0][0][1] + $y*$cell[0][1][1];
	    for($i=0;$i<@candidates;$i++) {
	      $n  = $candidates[$i];
	      $dx = $cdata[0][$td][$m][$a][0] - $cdata[0][$toh][$n][$hydroxo][0] + $shiftx;
	      $dy = $cdata[0][$td][$m][$a][1] - $cdata[0][$toh][$n][$hydroxo][1] + $shifty;
	      $dz = $cdata[0][$td][$m][$a][2] - $cdata[0][$toh][$n][$hydroxo][2];
	      $distsq = $dx*$dx + $dy*$dy + $dz*$dz;
	      if($distsq<$mindistsq) {
		$mindistsq = $distsq;
		$nearest   = $i;
		$ah        = $a;
	      }
	    } # end for $y
	  } # end for $x
	} # end for $n
      } # end if HG
    } #end for $a
    next if(sqrt($mindistsq)>$maxohdist);
    $changed++;
    @{$water[0]}=@{dclone \@{$cdata[0][$toh][$candidates[$nearest]]}};
    push(@remohs,splice(@candidates,$nearest,1));
    push(@{$water[0]},splice(@{$cdata[0][$td][$m]},$ah,1)); # shift proton
    if($lohrem==0) {
      $water[0][0][9] = 'OW'; # change atom type of water oxygen
      push(@{$cdata[0][$twat]},splice(@water,0,1)); # add water to cdata array
      $last=$#{$cdata[0][$twat]};
      $mol_numents[0][$twat]++;
    }
    # change the name of the oxygen
    $o = $mol_bondatoms[0][$td][$ah][0];
    $cdata[0][$td][$m][$o][9] = 'O';
    push(@{$cdata[0][$t2d]},splice(@{$cdata[0][$td]},$m,1)); # move molecule to deprotonated ones
    $mol_numents[0][$td]--;
    $mol_numents[0][$t2d]++;
    $m--;
  } #end for $m
} #end for $i

exit 1 if($err);
if($changed>0) {
  # remove hydroxides
  @remohs = sort {$a <=> $b} @remohs;
  for($i=$#remohs;$i>=0;$i--) {
    $n  = $remohs[$i];
    splice(@{$cdata[0][$toh]},$n,1);
    $mol_numents[0][$toh]--;
    splice(@remohs,$i,1);
  }
  write_config_file($filenameoutcfg,0,$config_title[0]);
  write_field_file($filenameoutfld,0);
  print "transfered $changed protons\n";
  exit 0;
} else {
  print "**** error: no changes made!\n";
  exit 1;
}
