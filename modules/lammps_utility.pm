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

read_lammpstrj_timestep read_dcd_timestep close_dcd_file
read_logfile);

our(@lmp_numstyles,@lmp_styles,@lmp_coeffs,$lmp_units,@lmp_fixes,%dcd_natoms,
%dcd_title,%dcd_timestep,%dcd_extra_block,%dcd_is_charmm,@lmp_atomtypes,
%dcd_frame,%dcd_startframe,%dcd_step,@ldata,@ldata_names,
$fmass,$ftime,$flength,$fenergy,$fpress,$fcharge,$fdipole);
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

sub read_dcd_header {
  print "reading header\n";
  my $fh = $_[0];
  my($bytes,$bytes_read,$is_charmm,$extra_block,$has_4dims);
  $bytes_read = read($fh, $bytes, 88);
  return 1 unless($bytes_read==88);
  return 1 unless(unpack('l',$bytes)==84); # first entry must be "84"
  my $totframes  = unpack('@8l',$bytes);
  my $startframe = unpack('@12l',$bytes);
  my $stepframe  = unpack('@16l',$bytes);
  my $endframe = unpack('@20l',$bytes);
  my $timestep = unpack('@44f',$bytes);
  # check if it's a CHARMM file
  if(unpack('@84l',$bytes)!=0) {
    $is_charmm   = 1;
    $extra_block = 1 if(unpack('@48l',$bytes)!=0);
    $has_4dims   = 1 if(unpack('@52l',$bytes)!=0);
  } else {
    my $dcd_is_charmm=0;
  }
  # second header block
  $bytes_read = read($fh, $bytes, 4);
  return 1 unless(unpack('l',$bytes)==84); # first entry must be "84"
  $bytes_read = read($fh, $bytes, 4); # number of bytes 
  $bytes_read = read($fh, $bytes, 4);
  my $ntitle = unpack('l',$bytes),"\n";
  my @title=();
  for(my $i=0;$i<$ntitle;$i++) {
    $bytes_read = read($fh, $bytes, 80);
    $title[$i]=unpack('A80',$bytes),"\n";
  }
  $bytes_read = read($fh, $bytes, 4);
  $bytes_read = read($fh, $bytes, 4);
  return 1 unless(unpack('l',$bytes)==4); # here must be "4"
  $bytes_read = read($fh, $bytes, 4);
  my $natoms=unpack('l',$bytes);
  $bytes_read = read($fh, $bytes, 4);
  return 1 unless(unpack('l',$bytes)==4); # here must be "4"
  # copy data
  $dcd_natoms{$fh}=$natoms;
  $dcd_timestep{$fh}=$timestep;
  $dcd_extra_block{$fh}=$extra_block;
  $dcd_is_charmm{$fh}=$is_charmm;
  $dcd_startframe{$fh}=$startframe;
  $dcd_step{$fh}=$stepframe;
  for(my $i=0;$i<$ntitle;$i++) {
    $dcd_title{$fh}[$i]=$title[$i];
  }
  $dcd_frame{$fh}=0;
  print "... done\n";
  return 0;
}

sub read_dcd_timestep {
  my $fh=$_[0];
  my $fi=$_[1];
  my $ci=$_[2];
  my($bytes,$bytes_read);
  return 2 if($ci<0);
  return 3 unless(defined($fh));
  
  if(not defined $dcd_natoms{$fh}) {
    if(read_dcd_header($fh)!=0) {
      print "**** error reading header of dcd file!\n";
      return 4;
    }
  }
  if($fi>=0 and $field_numatoms[$fi]!=$dcd_natoms{$fh}) {
    print "**** error: number of atoms in FIELD file does not match number of atoms in dcd file!\n";
    return 1;
  }
  
  # check whether we're at the end of the file
  my $nbytes=-1;
  my $currpos=tell($fh);
  $bytes_read = read($fh, $bytes, 4);
  $nbytes = unpack('l',$bytes) if($nbytes<0);
  return -1 unless($bytes_read==4);
  seek($fh,$currpos,0);
  
  #reset values
  @{$cdata[$ci]} = ();
  @{$cell[$ci]}  = ();
  @{$size[$ci]}  = ();
  @{$minpos[$ci]}=( 9e20, 9e20, 9e20);
  @{$maxpos[$ci]}=(-9e20,-9e20,-9e20);
  undef $frame_number[0];
  undef $frame_numatoms[0];
  
  if($dcd_extra_block{$fh}) {
    $bytes_read = read($fh, $bytes, 4);
    unless(unpack('l',$bytes)==48) {# here must be "48"
      print "**** error: expected \"48\" in extra block of dcd file!";
      return 1;
    }
    $bytes_read = read($fh, $bytes, 48);
    my @d=unpack("d d d d d d",$bytes);
    $bytes_read = read($fh, $bytes, 4);
    unless(unpack('l',$bytes)==48) {# here must be "48"
      print "**** error: expected \"48\" in extra block of dcd file!";
      return 1;
    }
    @{$cell[$ci]}=calc_cell_vecs($d[0],$d[2],$d[5],acos($d[4]),acos($d[3]),acos($d[1]));
    $size[$ci][0] = ($cell[$ci][0][0] + $cell[$ci][1][0] + $cell[$ci][2][0])/2;
    $size[$ci][1] = ($cell[$ci][0][1] + $cell[$ci][1][1] + $cell[$ci][2][1])/2;
    $size[$ci][2] = ($cell[$ci][0][2] + $cell[$ci][1][2] + $cell[$ci][2][2])/2;
  }
  
  
  if($fi<0) {
    for(my $c=0;$c<3;$c++) {
      $bytes_read = read($fh, $bytes, 4);
      $nbytes = unpack('l',$bytes);
      unless($nbytes/4==$dcd_natoms{$fh}) { # must be 4*N bytes for x/y/z arrays
        print "**** error: unexpected end of dcd file!\n";
        return 1;
      }
      $bytes_read = read($fh, $bytes, $nbytes);
      unless($bytes_read==$nbytes) {
        print "**** error: unexpected end of dcd file!\n";
        return 1;
      }
      for(my $i=0;$i<$dcd_natoms{$fh};$i++) {
        $cdata[$ci][0][0][$i][$c] = unpack("\@".($i*4)."f",$bytes);
        $minpos[$ci][$c] = $cdata[$ci][0][0][0][$c] if($cdata[$ci][0][0][0][$c]<$minpos[$ci][$c]);
        $maxpos[$ci][$c] = $cdata[$ci][0][0][0][$c] if($cdata[$ci][0][0][0][$c]>$maxpos[$ci][$c]);
      }
      $bytes_read = read($fh, $bytes, 4);
      unless($nbytes/4==$dcd_natoms{$fh}) {
        print "**** error: unexpected end of dcd file!\n";
        return 1;
      }
    }
  } else {
    my $molindex = 0;
    for(my $t=0;$t<$field_nummols[$fi];$t++) {
      @{$cdata[$ci][$t]}=();
      for(my $m=0;$m<$mol_numents[$fi][$t];$m++) {
        @{$cdata[$ci][$t][$m]}=();
        for(my $a=0;$a<$mol_numatoms[$fi][$t];$a++) {
          $cdata[$ci][$t][$m][$a][9]  = $mol_atomdata[$fi][$t][$a][0];
          $cdata[$ci][$t][$m][$a][10] = $molindex;
          $cdata[$ci][$t][$m][$a][12] = $mol_atomdata[$fi][$t][$a][1];
          $cdata[$ci][$t][$m][$a][13] = $mol_atomdata[$fi][$t][$a][2];
        }
        $molindex++;
      }
    }
    for(my $c=0;$c<3;$c++) {
      $bytes_read = read($fh, $bytes, 4);
      $nbytes = unpack('l',$bytes);
      unless($nbytes/4==$dcd_natoms{$fh}) { # must be 4*N bytes for x/y/z arrays
        print "**** error: unexpected end of dcd file!\n";
        return 1;
      }
      $bytes_read = read($fh, $bytes, $nbytes);
      unless($bytes_read==$nbytes) {
        print "**** error: unexpected end of dcd file!\n";
        return 1;
      }
      my $i=0;
      for(my $t=0;$t<$field_nummols[$fi];$t++) {
        for(my $m=0;$m<$mol_numents[$fi][$t];$m++) {
          for(my $a=0;$a<$mol_numatoms[$fi][$t];$a++) {
            $cdata[$ci][$t][$m][$a][$c] = unpack("\@".($i*4)."f",$bytes);
            $minpos[$ci][$c] = $cdata[$ci][$t][$m][$a][$c] if($cdata[$ci][$t][$m][$a][$c]<$minpos[$ci][$c]);
            $maxpos[$ci][$c] = $cdata[$ci][$t][$m][$a][$c] if($cdata[$ci][$t][$m][$a][$c]>$maxpos[$ci][$c]);
            $i++;
          }
        }
      }
      $bytes_read = read($fh, $bytes, 4);
      unless($nbytes/4==$dcd_natoms{$fh}) {
        print "**** error: unexpected end of dcd file!\n";
        return 1;
      }
    }
  }
  $frame_number[$ci]=$dcd_startframe{$fh}+$dcd_frame{$fh}*$dcd_step{$fh};
  if($dcd_extra_block{$fh}) {
    if(check_orthogonal($ci)) {
      $periodic_key[$ci]=2;
    } else {
      $periodic_key[$ci]=3;
    }
  } else {
    $periodic_key[$ci]=0;
  }
  $config_key[$ci]=0;
  $dcd_frame{$fh}++;
  return 0;
}

sub close_dcd_file {
  undef $dcd_extra_block{$_[0]};
  undef $dcd_frame{$_[0]};
  undef $dcd_is_charmm{$_[0]};
  undef $dcd_natoms{$_[0]};
  undef $dcd_startframe{$_[0]};
  undef $dcd_step{$_[0]};
  undef $dcd_timestep{$_[0]};
  undef $dcd_title{$_[0]};
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
  my($fh,$line,$i,@data,$first,@split);
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
  while(<$fh>) {
    $_=~s/(^\s+|\s+$)//g;
    if(/^[0-9\s\-\+eE.]+$/ and /[0-9]/) {
      if($first) {
        $first=0;
        @{$ldata_names[$li]}=split(/\s+/,$line);
      }
      @split = split(/\s+/,$_);
      next if($#split!=$#{$ldata_names[$li]});
      push(@{$ldata[$li]},[@split]);
    } else {
      $line=$_;
    }
  }
  close($fh);
  return 0;
}

1;

__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

dlpoly_utility - Perl extension for blah blah blah

=head1 SYNOPSIS

  use dlpoly_utility;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for dlpoly_utility, created by h2xs. It looks like the
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