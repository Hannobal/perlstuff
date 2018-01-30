package lammps_utility;

use 5.014002;
use strict;
# use warnings;
use POSIX qw(ceil);
use Math::Trig;
use Storable qw(dclone);

use hanno_utility;
use dlpoly_utility;

require Exporter;
our @ISA = qw(Exporter);
our %EXPORT_TAGS = ( 'all' => [ qw() ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw(@lmp_numstyles @lmp_styles @lmp_coeffs $lmp_units @lmp_fixes
@lmp_variables @lmp_units @lmp_atomtypes @ldata @ldata_names

read_lammps_data read_lammpstrj_timestep read_dcd_timestep find_last_dcd_timestep
close_dcd_file read_logfile open_dcd_write write_dcd_timestep skip_dcd_timestep
goto_last_dcd_timestep goto_dcd_timestep);

sub read_lammps_data {
  my $filename = $_[0];
  my $fi       = $_[1];
  my $ci       = $_[2];
  my($filehandle,$keyword,@bonds,@linedata,@atoms,@vel,@angles,@dihedrals,@impropers,
  @masses,@bdtypes,@angtypes,@dihtypes,@imptypes,
  $natoms,$nbonds,$nangles,$ndihedrals,$nimpropers,$nattypes,$nbdtypes,$nangtypes,
  $ndihtypes,$nimptypes,$nmax,$arrptr);
  if((not defined($fi)) or $fi<0) {
    print "**** error: index \$fi must not be smaller than zero in sub read_lammps_data\n";
    return -1;
  }
  
  if((not defined($ci)) or $ci<0) {
    print "**** error: index \$ci must not be smaller than zero in sub read_lammps_data\n";
    return -1;
  }
  
  if(not open($filehandle, "<", $filename)) {
    print "**** error: Can't open LAMMPS data file \"$filename\": $!\n";
    return 1;
  }
  
  clear_field_data($fi);
  undef @{$cdata[$ci]};
  undef @{$cell[$ci]};
  undef @{$size[$ci]};
  $field_filename[$fi]=$filename;
  
  ############## read header from data file ##########################
  
  $field_title[$fi] = <$filehandle>;
  $field_title[$fi] =~ s/(^\s+|\s+$)//g;
  $config_title[$ci] = $field_title[$fi];
  @{$cell[$ci]} = ([0,0,0],[0,0,0],[0,0,0]);
  while(<$filehandle>) {
    $_ =~ s/(^\s+|#.*$|\s+$)//g;
    next if(length($_)==0);
    @linedata = split(/\s+/,$_);
    last if(not check_real($linedata[0]));
    if(/atoms/) {
      $natoms=$linedata[0];
    } elsif(/bonds/) {
      $nbonds=$linedata[0];
    } elsif(/angles/) {
      $nangles=$linedata[0];
    } elsif(/dihedrals/) {
      $ndihedrals=$linedata[0];
    } elsif(/impropers/) {
      $nimpropers=$linedata[0];
    } elsif(/atom\s+types/) {
      $nattypes=$linedata[0];
    } elsif(/bond\s+types/) {
      $nbdtypes=$linedata[0];
    } elsif(/angle\s+types/) {
      $nangtypes=$linedata[0];
    } elsif(/dihedral\s+types/) {
      $ndihtypes=$linedata[0];
    } elsif(/xlo\s+xhi/) {
      $cell[$ci][0][0]=$linedata[1]-$linedata[0];
    } elsif(/ylo\s+yhi/) {
      $cell[$ci][1][1]=$linedata[1]-$linedata[0];
    } elsif(/zlo\s+zhi/) {
      $cell[$ci][2][2]=$linedata[1]-$linedata[0];
    } elsif(/xy\s+xz\s+yz/) {
      $cell[$ci][1][0]=$linedata[0];
      $cell[$ci][2][0]=$linedata[1];
      $cell[$ci][2][1]=$linedata[2];
    } else {
      print "**** error: unknown line in file $filename!\n";
      close($filehandle); return 1;
    }
  }
  
  ############## read rest of data file ##########################
  do {
    $_ =~ s/(^\s+|#.*$|\s+$)//g;
    my $keyword = $_;
    $_ = <$filehandle>; $_ =~ s/(^\s+|#.*$|\s+$)//g;
    if(length($_)!=0) {
      print "**** error: empty line expected at beginning of section $keyword!\n";
      close($filehandle); return 1;
    }
    if($keyword eq "Masses") {
      $arrptr=\@masses;
      $nmax=$nattypes;
    } elsif($keyword eq "Bond Coeffs") {
      $arrptr=\@bdtypes;
      $nmax=$nbdtypes;
    } elsif($keyword eq "Angle Coeffs") {
      $arrptr=\@angtypes;
      $nmax=$nangtypes;
    } elsif($keyword eq "Dihedral Coeffs") {
      $arrptr=\@dihtypes;
      $nmax=$ndihtypes;
    } elsif($keyword eq "Improper Coeffs") {
      $arrptr=\@imptypes;
      $nmax=$nimptypes;
    } elsif($keyword eq "Atoms") {
      $arrptr=\@atoms;
      $nmax=$natoms;
    } elsif($keyword eq "Velocities") {
      $arrptr=\@vel;
      $nmax=$natoms;
    } elsif($keyword eq "Bonds") {
      $arrptr=\@bonds;
      $nmax=$nbonds;
    } elsif($keyword eq "Angles") {
      $arrptr=\@angles;
      $nmax=$nangles;
    } elsif($keyword eq "Dihedrals") {
      $arrptr=\@dihedrals;
      $nmax=$ndihedrals;
    } elsif($keyword eq "Impropers") {
      $arrptr=\@impropers;
      $nmax=$nimpropers;
    } else {
      print "**** error: unknown keyword \"$keyword\"!\n";
      close($filehandle); return 1;
    }
    for(my $i=0;$i<$nmax;$i++) {
      $_ = <$filehandle>; $_ =~ s/(^\s+|#.*$|\s+$)//g;
      @linedata = split(/\s+/,$_);
      my $index = splice(@linedata,0,1)-1;
      @{${$arrptr}[$index]} = @linedata;
#       push(@{$arrptr},[@linedata]);
    }
    while(<$filehandle>) {
      $_ =~ s/(^\s+|#.*$|\s+$)//g;
      last if(length($_)!=0);
    }
  } while(not eof($filehandle));
  close($filehandle);
  # maybe we should check for sanity here...
  
  ############## adjust atom and type indices ####################
  for(my $i=0;$i<@bonds;$i++) {
    $bonds[$i][0]--;
    $bonds[$i][1]--;
    $bonds[$i][2]--;
  }
  @bonds = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @bonds;
  
  for(my $i=0;$i<@angles;$i++) {
    $angles[$i][0]--;
    $angles[$i][1]--;
    $angles[$i][2]--;
    $angles[$i][3]--;
  }
  @angles = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3] } @angles;
  
  for(my $i=0;$i<@dihedrals;$i++) {
    $dihedrals[$i][0]--;
    $dihedrals[$i][1]--;
    $dihedrals[$i][2]--;
    $dihedrals[$i][3]--;
    $dihedrals[$i][4]--;
  }
  @dihedrals = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3]  || $a->[4] <=> $b->[4] } @dihedrals;
  
  for(my $i=0;$i<@impropers;$i++) {
    $impropers[$i][0]--;
    $impropers[$i][1]--;
    $impropers[$i][2]--;
    $impropers[$i][3]--;
    $impropers[$i][4]--;
  }
  @impropers = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3]  || $a->[4] <=> $b->[4] } @impropers;
  
  ############## build the bond map ##############################
  my @atinmol=();
  my @molsat=();
  my @atbnd=();
  for(my $i=0;$i<@bonds;$i++) {
    my $at1 = $bonds[$i][1];
    my $at2 = $bonds[$i][2];
    push(@{$atbnd[$at1]},$i);
    push(@{$atbnd[$at2]},$i);
    if(defined($atinmol[$at1])) {
      if(defined($atinmol[$at2])) {
        if($atinmol[$at1]!=$atinmol[$at2]) {
          # we have to join the fragments
          my $tmp = $atinmol[$at2];
          foreach $a (@{$molsat[$tmp]}) {
            $atinmol[$a]=$atinmol[$at1];
          }
          push(@{$molsat[$atinmol[$at1]]},splice(@{$molsat[$tmp]},0));
        } # otherwise it's a ring
      } else {
        # put the second atom in the same mol as the first
        push(@{$molsat[$atinmol[$at1]]},$at2);
        $atinmol[$at2] = $atinmol[$at1];
      }
    } elsif(defined($atinmol[$at2])) {
      # put the first atom in the same mol as the second
      push(@{$molsat[$atinmol[$at2]]},$at1);
      $atinmol[$at1] = $atinmol[$at2];
    } else {
      # create a new molecule and pack both atoms in it
      push(@molsat,[$at1,$at2]);
      $atinmol[$at1] = $#molsat;
      $atinmol[$at2] = $#molsat;
    }
  }
  undef @atinmol;
  # delete empty mols
  for(my $i=$#molsat;$i>=0;$i--) {
    splice(@molsat,$i,1) if($molsat[$i]=="");
  }
  @molsat = sort {$a->[0] <=> $b->[0]} @molsat;
  ############## sort into mol types #############################
  my @molintype=();
  my @typemols=();
  itermol:for(my $i=0;$i<@molsat;$i++) {
    # compare with molecules already checked
    cmpmol:for(my $j=$i-1;$j>=0;$j--) {
      next if($#{$molsat[$i]} != $#{$molsat[$j]}); # number of atoms
      # compare atom types
      for(my $k=0;$k<@{$molsat[$i]};$k++) {
        my $a = $molsat[$i][$k];
        my $b = $molsat[$j][$k];
        if($atoms[$a][1]!=$atoms[$b][1]) {
          next cmpmol;
        }
      }
      # found a match:
      $molintype[$i]=$molintype[$j];
      push(@{$typemols[$molintype[$j]]},$i);
      next itermol;
    }
    # no match found:
    push(@{$typemols[$#typemols+1]},$i);
    $molintype[$i]=$#typemols;
  }

  ############## fill field data #################################
  $field_nummols[$fi] = $#typemols+1;
  for(my $t=0;$t<@typemols;$t++) {
    $mol_numents[$fi][$t]  = $#{$typemols[$t]};
    $mol_numatoms[$fi][$t] = $#{$molsat[$typemols[$t][0]]};
    if($cell[$ci][1][0]!=0 or $cell[$ci][2][0]!=0 or $cell[$ci][2][1]!=0) {
      $periodic_key[$ci]=3;
    } else {
      $periodic_key[$ci]=2;
    }
    for(my $m=0;$m<@{$typemols[$t]};$m++) {
      my $mold=$typemols[$t][$m];
      for(my $a=0;$a<@{$molsat[$mold]};$a++) {
        my $aold=$molsat[$mold][$a];
        my $type=$atoms[$aold][1];
        $cdata[$ci][$t][$m][$a][0] = $atoms[$aold][3];
        $cdata[$ci][$t][$m][$a][1] = $atoms[$aold][4];
        $cdata[$ci][$t][$m][$a][2] = $atoms[$aold][5];
        if(@vel) {
          $cdata[$ci][$t][$m][$a][3] = $vel[$aold][0];
          $cdata[$ci][$t][$m][$a][4] = $vel[$aold][1];
          $cdata[$ci][$t][$m][$a][5] = $vel[$aold][2];
        }
        $cdata[$ci][$t][$m][$a][9]  = $type;
        $cdata[$ci][$t][$m][$a][11] = $aold;
        $cdata[$ci][$t][$m][$a][12] = $masses[$type];
        $cdata[$ci][$t][$m][$a][13] = $atoms[$aold][2];
      }
    }
    $mol_charge[$fi][$t]=0;
    $mol_mass[$fi][$t]=0;
    for(my $a=0;$a<@{$molsat[0]};$a++) {
      my $aold=$molsat[0][$a];
      my $type=$atoms[$molsat[0][$a]][1];
      $mol_atomdata[$fi][$t][$a][0] = $type;
      $mol_atomdata[$fi][$t][$a][1] = $masses[$type];
      $mol_atomdata[$fi][$t][$a][2] = $atoms[$aold][2];
      $mol_charge[$fi][$t] += $mol_atomdata[$fi][$t][$a][2];
      $mol_mass[$fi][$t]   += $mol_atomdata[$fi][$t][$a][1];
    }
  }
}

sub print2darr {
  my @arr = @{$_[0]};  
  for(my $i=0;$i<@arr;$i++) {
    print $#{$arr[$i]}+1," atoms: ";
    print join(" ",@{$arr[$i]}),"\n";
  }
}

our(@lmp_numstyles,@lmp_styles,@lmp_coeffs,$lmp_units,@lmp_fixes,%dcd_natoms,
%dcd_title,%dcd_timestep,%dcd_extra_block,%dcd_is_charmm,@lmp_atomtypes,
%dcd_frame,%dcd_startframe,%dcd_endframe,%dcd_step,%dcd_has_4dims,%dcd_nframes,,
%dcd_timestep_length,%dcd_startbyte,@ldata,@ldata_names,$fmass,$ftime,$flength,
$fenergy,$fpress,$fcharge,$fdipole);
# indices for arrays:
#   styles:  atomtypes:
# 0 bond     name/element
# 1 angle    mass
# 2 dihedral
# 3 improper
# 4 pair
# 5 coulomb

sub read_lammpstrj_timestep {

  my $filehandle = $_[0];
  my $fi         = $_[1]; #index of corresponding FIELD file
  my $ci         = $_[2];
  my($i,$j,$t,$m,$a,$molindex,@linedata);
  
  return 2 if($ci<0);
  
  if(not $_=<$filehandle>) {return -1} #end of file
  if(not /TIMESTEP/) {
    print "error reading lammpstrj file on line $.: \"ITEM: TIMESTEP\" expected but found:\n$_";
    return 2;
  }
  
  #reset values
  @{$cdata[$ci]} = ();
  @{$cell[$ci]}  = ();
  @{$size[$ci]}  = ();
  @{$minpos[$ci]}=( 9e20, 9e20, 9e20);
  @{$maxpos[$ci]}=(-9e20,-9e20,-9e20);
  undef $periodic_key[$ci];
  $frame_numatoms[$ci]=-1;
  
  $frame_number[$ci]=<$filehandle>; $frame_number[$ci]=~s/\s+$//;
  # print "frame $frame_number[$ci]\n";
  $molindex = 0;
  $periodic_key[$ci]=0;
  while($_=<$filehandle>) {
    if(not /ITEM/) {
      print "error reading lammpstrj file on line $.: expected ITEM but found:\n$_";
      return 2;
    }
    if(/NUMBER OF ATOMS/) {
      $frame_numatoms[$ci]=<$filehandle>; $frame_numatoms[$ci]=~s/\s+$//;
      if($fi>=0) {
	if($frame_numatoms[$ci]!=$field_numatoms[$fi]) {
	  print "**** error: number of atoms in field information ($field_numatoms[$fi])".
	  "does not macht number of atoms in lammpstrj file ($frame_numatoms[$ci])!\n";
	  return 2;
	}
      }
      # print "numatoms $frame_numatoms[$ci]\n";
    } elsif(/BOX BOUNDS/) {
      if(/pp pp pp/) {
	$periodic_key[$ci]=2;
      } elsif(/pp pp ff/) {
	$periodic_key[$ci]=6;
      } else {
	print "**** error: unknown box: $_";
	return 2;
      }
      # print "box $periodic_key[$ci]\n";
      if(/pp pp pp/ or /pp pp ff/) {
	for($i=0;$i<3;$i++) {
	  $_=<$filehandle>; @linedata = split(/\s+/, $_);
	  $cell[$ci][$i][$i] = $linedata[1]-$linedata[0];
	  # print "$cell[$ci][$i][$i]\n";
	}
      }
    } elsif(/ATOMS/) {
      # print "reading atoms\n";
      if($frame_numatoms[$ci]<0) {
	print "**** error: number of atoms not defined yet!\n";
	return 2;
      }
      @linedata = split(/\s+/, $_);
      my @indices=();
      for($i=2;$i<@linedata;$i++) {
	if($linedata[$i] =~ /^x/) {
	  push(@indices,0);
	} elsif($linedata[$i] =~ /^y/) {
	  push(@indices,1);
	} elsif($linedata[$i] =~ /^z/) {
	  push(@indices,2);
	} elsif($linedata[$i] =~ /^vx/) {
	  push(@indices,3);
	} elsif($linedata[$i] =~ /^vy/) {
	  push(@indices,4);
	} elsif($linedata[$i] =~ /^vz/) {
	  push(@indices,5);
	} elsif($linedata[$i] =~ /^fx/) {
	  push(@indices,6);
	} elsif($linedata[$i] =~ /^fy/) {
	  push(@indices,7);
	} elsif($linedata[$i] =~ /^fz/) {
	  push(@indices,8);
	} elsif($linedata[$i] =~ /^(element|type)/) {
	  push(@indices,9);
	} elsif($linedata[$i] =~ /^mol/) {
	  push(@indices,10);
	} elsif($linedata[$i] =~ /^id/) {
	  push(@indices,11);
	} elsif($linedata[$i] =~ /^mass/) {
	  push(@indices,12);
	} elsif($linedata[$i] =~ /^q/) {
	  push(@indices,13);
	} else {
	  push(@indices,-1);
	}
      }
      # print "indices: ",join(" ",@indices),"\n";
      if($fi>=0) {  # read using FIELD information
	for($t=0;$t<$field_nummols[$fi];$t++) {
	  @{$cdata[$ci][$t]}=();
	  for($m=0;$m<$mol_numents[$fi][$t];$m++) {
	    @{$cdata[$ci][$t][$m]}=();
	    for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
	      @{$cdata[$ci][$t][$m][$a]}=();
	      $_=<$filehandle>; @linedata = split(/\s+/, $_);
	      for($j=0;$j<@indices;$j++) {
		next if($indices[$j]<0);
		$cdata[$ci][$t][$m][$a][$indices[$j]]=$linedata[$j];
	      }
	      $minpos[$ci][0] = $cdata[$ci][$t][$m][$a][0] if($cdata[$ci][$t][$m][$a][0]<$minpos[$ci][0]);
	      $minpos[$ci][1] = $cdata[$ci][$t][$m][$a][1] if($cdata[$ci][$t][$m][$a][1]<$minpos[$ci][1]);
	      $minpos[$ci][2] = $cdata[$ci][$t][$m][$a][2] if($cdata[$ci][$t][$m][$a][2]<$minpos[$ci][2]);
	      $maxpos[$ci][0] = $cdata[$ci][$t][$m][$a][0] if($cdata[$ci][$t][$m][$a][0]>$maxpos[$ci][0]);
	      $maxpos[$ci][1] = $cdata[$ci][$t][$m][$a][1] if($cdata[$ci][$t][$m][$a][1]>$maxpos[$ci][1]);
	      $maxpos[$ci][2] = $cdata[$ci][$t][$m][$a][2] if($cdata[$ci][$t][$m][$a][2]>$maxpos[$ci][2]);
	    }
	  }
	}
      } else { # read w/o FIELD information
	for($i=0;$i<$frame_numatoms[$ci];$i++) {
	  $_=<$filehandle>; @linedata = split(/\s+/, $_);
	  $t=0;
	  $m=0;
	  $a=$i;
	  for($j=0;$j<@indices;$j++) {
	    next if($indices[$j]<0);
	    $cdata[$ci][$t][$m][$a][$indices[$j]]=$linedata[$j];
	  }
	  $minpos[$ci][0] = $cdata[$ci][$t][$m][$a][0] if($cdata[$ci][$t][$m][$a][0]<$minpos[$ci][0]);
	  $minpos[$ci][1] = $cdata[$ci][$t][$m][$a][1] if($cdata[$ci][$t][$m][$a][1]<$minpos[$ci][1]);
	  $minpos[$ci][2] = $cdata[$ci][$t][$m][$a][2] if($cdata[$ci][$t][$m][$a][2]<$minpos[$ci][2]);
	  $maxpos[$ci][0] = $cdata[$ci][$t][$m][$a][0] if($cdata[$ci][$t][$m][$a][0]>$maxpos[$ci][0]);
	  $maxpos[$ci][1] = $cdata[$ci][$t][$m][$a][1] if($cdata[$ci][$t][$m][$a][1]>$maxpos[$ci][1]);
	  $maxpos[$ci][2] = $cdata[$ci][$t][$m][$a][2] if($cdata[$ci][$t][$m][$a][2]>$maxpos[$ci][2]);
	}
      }
    } elsif(/TIMESTEP/) {
      seek($filehandle, -length($_), 1);
      last;
    } else {
      print "**** error: unknown item: $_";
      return 2;
    }
  }
  return 0;
}

sub clear_dcd {
  my $fh = $_[0];
  delete $dcd_natoms{$fh};
  delete $dcd_nframes{$fh};
  delete $dcd_timestep{$fh};
  delete $dcd_extra_block{$fh};
  delete $dcd_has_4dims{$fh};
  delete $dcd_is_charmm{$fh};
  delete $dcd_startframe{$fh};
  delete $dcd_step{$fh};
  delete $dcd_frame{$fh};
  delete $dcd_title{$fh};
  delete $dcd_timestep_length{$fh};
}

sub read_fortran_record {
  my $fh = $_[0];
  my($tmp,$content,$len);
  return undef unless(read($fh,$tmp, 4));
  $len = unpack('l',$tmp);
  return undef unless(read($fh, $content, $len));
  return undef unless(read($fh,$tmp, 4));
  return undef if($len != unpack('l',$tmp));
  return $content;
}

sub read_dcd_header {
  my $fh = $_[0];
  my($content,@data);
  unless(defined($fh)) {
    print "**** error: no DCD file in handle in call of read_dcd_header!\n";
    return 3;
  }
  $content = read_fortran_record($fh);
  unless(defined($content)) {
    print "**** error reading DCD header!\n";
    return 3;
  }
  @data=unpack('c4i9f1i10',$content);
#   for(my $i=0;$i<@data;$i++) {
#     print "$i $data[$i]\n";
#   }
#   print join(" ",@data),"\n";
  $content = read_fortran_record($fh);
  my $ntitle=unpack('l',$content);
  for(my $i=0;$i<$ntitle;$i++) {
    $dcd_title{$fh}[$i] = unpack('@'.(4+$i*80).'A80',$content);
  }
  $content = read_fortran_record($fh);
  unless(defined($content)) {
    clear_dcd($fh);
    print "**** error reading DCD header!\n";
    return 1;
  }
  $dcd_natoms{$fh}=unpack('l',$content);
  $dcd_startframe{$fh} = $data[5];
  $dcd_frame{$fh}      = $data[5]; # same as startframe
  $dcd_step{$fh}       = $data[6];
  $dcd_endframe{$fh}   = $data[7];
  $dcd_timestep{$fh}   = $data[13];
  if($data[23]!=0) {
    $dcd_extra_block{$fh} = $data[14];
    $dcd_has_4dims{$fh}   = $data[15];
    $dcd_is_charmm{$fh}   = 1;
  } else {
    $dcd_extra_block{$fh} = 0;
    $dcd_has_4dims{$fh}   = 0;
    $dcd_is_charmm{$fh}   = 0;
  }
  if($dcd_has_4dims{$fh}) {
    $dcd_timestep_length{$fh}=(4*$dcd_natoms{$fh}+8)*4;
  } else {
    $dcd_timestep_length{$fh}=(4*$dcd_natoms{$fh}+8)*3;
  }
  $dcd_timestep_length{$fh}+=56 if($dcd_extra_block{$fh});
  $dcd_step{$fh}=1 if($dcd_step{$fh}<=0);
  if($dcd_endframe{$fh}<$dcd_nframes{$fh}*$dcd_step{$fh}) {
    $dcd_endframe{$fh}=$dcd_nframes{$fh}*$dcd_step{$fh};
  }
  $dcd_startbyte{$fh} = tell($fh);
  return 0;
}

sub goto_dcd_timestep {
  # first argument:  filehandle
  # second argument: frame number (<0 => jump to end)
  unless(defined($_[0])) {
    print "**** error: no DCD file in handle in call of goto_dcd_timestep!\n";
    return 1;
  }
  unless(defined($_[1])) {
    print "**** error: frame number not defined in call of goto_dcd_timestep\n";
    return 1;
  }
  unless(defined($dcd_natoms{$_[0]})) {
    return 1 unless(read_dcd_header($_[0])==0);
  }
  return goto_last_dcd_timestep($_[0]) if($_[1]<0);
  if($_[1]<$dcd_startframe{$_[0]}) {
    print "**** error: desired frame $_[1] smaller than start frame of DCD ($dcd_startframe{$_[0]})!\n";
    return 1;
  }
  if($_[1]>$dcd_endframe{$_[0]}) {
    print "**** error: desired frame $_[1] larger than end frame of DCD ($dcd_endframe{$_[0]})!\n";
    return 1;
  }
  if(($_[1]-$dcd_startframe{$_[0]})%$dcd_step{$_[0]} != 0) {
    print "**** error: desired frame $_[1] does not fit frame pattern\n",
          "            of DCD (start $dcd_startframe{$_[0]} step $dcd_step{$_[0]})!\n";
    return 1;
  }
  seek($_[0],$dcd_startbyte{$_[0]}+$dcd_timestep_length{$_[0]}*($_[1]-$dcd_startframe{$_[0]})/$dcd_step{$_[0]},0);
  $dcd_frame{$_[0]} = $_[1];
  return 0;
}

sub skip_dcd_timestep {
  # first argument:  filehandle
  # second argument: number of frames to skip (optional)
  my $nskip=$_[1];
  $nskip=1 unless(defined($nskip));
  unless(defined($dcd_natoms{$_[0]})) {
    return 1 unless(read_dcd_header($_[0])==0);
  }
  seek($_[0],$nskip*$dcd_timestep_length{$_[0]},1);
  $dcd_frame{$_[0]}+=$dcd_step{$_[0]};
  return 0;
}

sub goto_last_dcd_timestep {
  unless(defined($dcd_natoms{$_[0]})) {
    return 1 unless(read_dcd_header($_[0])==0);
  }
  seek($_[0],-$dcd_timestep_length{$_[0]},2);
  $dcd_frame{$_[0]}+=$dcd_endframe{$_[0]}-$dcd_step{$_[0]};
  return 0;
}

sub read_dcd_timestep {
  my $fh=$_[0];
  my $fi=$_[1];
  my $ci=$_[2];
  my($content,@data,$cmax,$str);
  return 2 if($ci<0);
  return 3 unless(defined($fh));
  unless(defined($dcd_natoms{$fh})) {
    return 4 unless(read_dcd_header($fh)==0);
  }
  if($fi>=0 and $field_numatoms[$fi]!=$dcd_natoms{$fh}) {
    print "**** error: number of atoms in FIELD file does not match number of atoms in dcd file!\n";
    return 1;
  }
  $content = read_fortran_record($fh);
  return -1 unless(defined $content); # probably the end of the file
  #reset values
  @{$cdata[$ci]} = ();
  @{$cell[$ci]}  = ();
  @{$size[$ci]}  = ();
  @{$minpos[$ci]}=( 9e20, 9e20, 9e20);
  @{$maxpos[$ci]}=(-9e20,-9e20,-9e20);
  undef $frame_number[$ci];
  undef $frame_numatoms[$ci];
  undef $periodic_key[$ci];
  if($dcd_extra_block{$fh}) {
    @data=unpack("d6",$content); # cell data
    @{$cell[$ci]}=calc_cell_vecs($data[0],$data[2],$data[5],acos($data[4]),acos($data[3]),acos($data[1]));
    $size[$ci][0] = ($cell[$ci][0][0] + $cell[$ci][1][0] + $cell[$ci][2][0])/2;
    $size[$ci][1] = ($cell[$ci][0][1] + $cell[$ci][1][1] + $cell[$ci][2][1])/2;
    $size[$ci][2] = ($cell[$ci][0][2] + $cell[$ci][1][2] + $cell[$ci][2][2])/2;
#     print "cell: ",join(" ",@data),"\n";
  }
  $cmax=3;
  $cmax=4 if($dcd_has_4dims{$fh});
  $str="f$dcd_natoms{$fh}";
  if($fi<0) {
    for(my $c=0;$c<$cmax;$c++) {
      $content = read_fortran_record($fh);
      @data=unpack($str,$content);
      for(my $i=0;$i<$dcd_natoms{$fh};$i++) {
        $cdata[$ci][0][0][$i][11] = $i;
        $cdata[$ci][0][0][$i][$c] = $data[$i];
        $minpos[$ci][$c] = $cdata[$ci][0][0][$i][$c] if($cdata[$ci][0][0][$i][$c]<$minpos[$ci][$c]);
        $maxpos[$ci][$c] = $cdata[$ci][0][0][$i][$c] if($cdata[$ci][0][0][$i][$c]>$maxpos[$ci][$c]);
      }
    }
  } else {
    @data=();
    for(my $c=0;$c<$cmax;$c++) {
      $content = read_fortran_record($fh);
      push(@data,[unpack($str,$content)]);
    }
    my $i=0;
    my $molindex = 0;
    for(my $t=0;$t<$field_nummols[$fi];$t++) {
      @{$cdata[$ci][$t]}=();
      for(my $m=0;$m<$mol_numents[$fi][$t];$m++) {
        @{$cdata[$ci][$t][$m]}=();
        for(my $a=0;$a<$mol_numatoms[$fi][$t];$a++) {
          @{$cdata[$ci][$t][$m][$a]}=();
          for(my $c=0;$c<$cmax;$c++) {
            $cdata[$ci][$t][$m][$a][$c] = $data[$c][$i];
          }
          $cdata[$ci][$t][$m][$a][9]  = $mol_atomdata[$fi][$t][$a][0];
          $cdata[$ci][$t][$m][$a][10] = $molindex;
          $cdata[$ci][$t][$m][$a][11] = $i;
          $cdata[$ci][$t][$m][$a][12] = $mol_atomdata[$fi][$t][$a][1];
          $cdata[$ci][$t][$m][$a][13] = $mol_atomdata[$fi][$t][$a][2];
          $i++;
        }
        $molindex++;
      }
    }
  }
  $frame_numatoms[$ci]=$dcd_natoms{$fh};
  $frame_number[$ci]=$dcd_frame{$fh};
  $dcd_frame{$fh}+=$dcd_step{$fh};
  return 0;
}

sub close_dcd_file {
  clear_dcd($_[0]);
  close($_[0]);
}

sub find_lammps_coeff {
  # call example: $found=find_lammps_style(\@data,\@{$lmp_coeffs[0][0]});
  my @data=@{$_[0]};
  my @arr=@{$_[1]};
  return -1 if(not @data);
  return -1 if(not @arr);
  looparr:for(my $i=0;$i<@arr;$i++) {
    next if(@{$arr[$i]}!=@data);
    next if($data[0] ne $arr[$i][0]); # the style for bond/angle/...
    for(my $j=1;$j<@data;$j++) {
      next looparr if($arr[$i][$j] != $data[$j]);
    }
    return $i;
  }
  return -1;
}

sub find_lammps_atomtype {
  # call example: find_lammps_atomtype("HW",0);
  for(my $i=0;$i<@{$lmp_atomtypes[$_[1]]};$i++) {
    return $i if($_[0] eq $lmp_atomtypes[$_[1]][$i][0]);
  }
  return -1;
}


sub convert_bond_style {
  # convert from dlpoly to lammps
  my @res=();
  if($_[0]=~/^harm/) {
    $res[0] = "harmonic";
    $res[1] = $_[3]*0.5*$fenergy/($flength*$flength);
    $res[2] = $_[4]*$flength;
  } elsif($_[0]=~/^quar/) {
    $res[0] = "class2";
    $res[1] = $_[4]*$flength;
    $res[2] = $_[3]*0.5*$fenergy;
    $res[3] = $_[5]*$fenergy/3.0;
    $res[3] = $_[6]*0.25*$fenergy;
  } else {
    die "**** error: bond potential \"$_[0]\" is not implemented (yet)!\n";
  }
  return @res
}

sub convert_angle_style {
  # convert from dlpoly to lammps
  my @res=();
  if($_[0]=~/^harm/) {
    $res[0] = "harmonic";
    $res[1] = 0.5*$fenergy*$_[4];
    $res[2] = $_[5];
  } elsif($_[0]=~/^hcos/) {
    $res[0] = "cosine/squared";
    $res[1] = 0.5*$fenergy*$_[4];
    $res[2] = $_[5];
  } else {
    die "**** error: angle potential \"$_[0]\" is not implemented (yet)!\n";
  }
  return @res;
}

sub convert_dihedral_style {
  # convert from dlpoly to lammps
  my @res=();
  if($_[0]=~/^harm/) {
    $res[0] = "harmonic";
    $res[1] = 0.5*$fenergy*$_[5];
    $res[2] = $_[6];
  } elsif($_[0]=~/^cos/) {
    $res[0] = "charmm";
    $res[1] = $_[5]*$fenergy;
    $res[2] = $_[7];
    $res[3] = $_[6];
    $res[4] = 0;
  } else {
    die "**** error: dihedral potential \"$_[0]\" is not implemented (yet)!\n";
  }
  return @res;
}

sub convert_vdw {
  # call_example: @data=convert_vdw(0,@{$vdw_data[0][3]});
  my @res=();
  my($k,$l);
  $k = find_lammps_atomtype($_[0]);
  $k = find_lammps_atomtype($_[1]);
  if($_[0]=~/^lj/) {
    
  } else {
    die "**** error: vdw potential \"$_[0]\" is not implemented (yet)!\n";
  }
  return @res;
}

sub read_logfile {
  my $filename = $_[0];
  my $li       = $_[1]; # index for array of field-data
  my($fh,$line,$i,@data,$first,@split,$continue);
  if($li<0) {
    print "**** error: index \$li must not be smaller than zero in sub read_field_file\n";
    return 1;
  }
  if(not open($fh, "<", $filename)) {
    print "**** error: Can't open lammps log file \"$filename\": $!\n";
    return 1;
  }
  undef @{$ldata[$li]};
  @{$ldata[$li]} = ();
  $first=1;
  $continue=0;
  while(<$fh>) {
    $_=~s/(^\s+|\s+$)//g;
    $continue=1 if(/minimiz/i);
    if(/^[0-9\s\-\+eE.]+$/ and /[0-9]/) {
      if($continue) {
        $continue=2;
        next;
      }
      if($first) {
        $first=0;
        @{$ldata_names[$li]}=split(/\s+/,$line);
      }
      @split = split(/\s+/,$_);
      next if($#split!=$#{$ldata_names[$li]});
      push(@{$ldata[$li]},[@split]);
    } else {
      $continue=0 if($continue==2);
      $line=$_;
    }
  }
  close($fh);
  return 0;
}

1;

sub write_dcd_timestep {
  my($fh,$fi,$ci,$packed,$nbytes,$n);
  ($fh,$fi,$ci) = @_;
  return if(not defined $dcd_natoms{$fh});
  $nbytes = $dcd_natoms{$fh}*4;
  for(my $c=0;$c<3;$c++) {
    print $fh pack('l<',$nbytes);
    $n=0;
    for(my $t=0;$t<@{$cdata[$ci]};$t++) {
      for(my $m=0;$m<@{$cdata[$ci][$t]};$m++) {
        for(my $a=0;$a<@{$cdata[$ci][$t][$m]};$a++) {
          $n++;
          print $fh pack('f<',$cdata[0][$t][$m][$a][$c]);
        }
      }
    }
    print $fh pack('l<',$nbytes);
    if($n!=$dcd_natoms{$fh}) {
      print "**** error: number of atoms ($n) does not match dcd file ($dcd_natoms{$fh})\n";
    }
  }
}

sub open_dcd_write {
  my($fh,$filename,$natoms,$timestep,$startframe,$endframe,$stepframe,$has_4dims,
  $extra_block,$is_charmm,$title);
  ($filename,$natoms,$title,$extra_block,$has_4dims,$timestep) = @_;
  $extra_block = 0 if(not defined($extra_block));
  $has_4dims   = 0 if(not defined($has_4dims));
  $timestep    = 1 if(not defined($timestep));
  $title       = "" if(not defined($title));
  $startframe  = 0;
  $stepframe   = 1;
  
  if(not open($fh, ">:raw", $filename)) {
    print "**** error: Can't open DCD file \"$filename\": $!\n";
    return undef;
  }
  
  if($extra_block or $has_4dims) {
    $is_charmm=24;
    $dcd_is_charmm{$fh}=1;
  }
  $dcd_frame{$fh}       = 0;
  $dcd_is_charmm{$fh}   = 1;
  $dcd_natoms{$fh}      = $natoms;
  $dcd_timestep{$fh}    = $timestep;
  $dcd_extra_block{$fh} = $extra_block;
  $dcd_has_4dims{$fh}   = $has_4dims;
  $dcd_step{$fh}        = $stepframe;
  $dcd_startframe{$fh}  = $startframe;
  $dcd_title{$fh}[0]    = $title;
  $dcd_title{$fh}[1]    = "written by perl";
  
  print $fh pack('(lA[4]l[9]fl[11])<',
    84,
    "CORD",
    $dcd_frame{$fh},  # number of frames
    $dcd_startframe{$fh},
    $dcd_step{$fh},
    0,0,0,0,0,0, # six zeros
    $dcd_timestep{$fh},
    $dcd_extra_block{$fh},
    $dcd_has_4dims{$fh},
    0,0,0,0,0,0,0, # seven zeros,
    $is_charmm,
    84);
  print $fh pack('(llA[80]A[80]l)<',
    164,2,$dcd_title{$fh}[0],$dcd_title{$fh}[1],164);
  print $fh pack('(lll)<',
  4,$dcd_natoms{$fh},4);
  
  return $fh;
}

__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

lammps_utility - Perl extension for blah blah blah

=head1 SYNOPSIS

  use lammps_utility;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for lammps_utility, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Hanno Dietrich, E<lt>dietrich@E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Hanno Dietrich

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.2 or,
at your option, any later version of Perl 5 you may have available.


=cut