package dlpoly_utility;

use 5.014002;
use strict;
# use warnings;
use POSIX qw(ceil);
use Math::Trig;
use Storable qw(dclone);

use hanno_utility;

require Exporter;
our @ISA = qw(Exporter);
our %EXPORT_TAGS = ( 'all' => [ qw() ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw(@cell @size @config_key @config_title @periodic_key @cdata @frame_number @frame_numatoms @frame_timestep
    @field_numheader @field_nummols @field_numatoms @field_numvdw @field_numtbp @field_numfbp @field_nummetal
    @field_numextern @field_units @field_neutral @field_title @field_header @field_vdwdata @field_tbpdata
    @field_fbpdata @field_metaldata @field_repsdata @field_externdata @field_numatomtypes @field_atomtypes @field_filename
    @cell @size @config_key @config_title @periodic_key @cdata @frame_number @frame_numatoms @frame_timestep
    @field_numheader @field_nummols @field_numatoms @field_numvdw @field_numtbp @field_numfbp @field_numextern
    @field_units @field_neutral @field_title @field_header @field_vdwdata @field_tbpdata @field_fbpdata
    @field_externdata @field_numatomtypes @field_atomtypes @field_filename
    @mol_name @mol_charge @mol_mass @mol_numents @mol_numatoms @mol_numbonds @mol_numconstraints
    @mol_numangles @mol_numdihedrals @mol_numinversions @mol_numrigid @mol_numtether @mol_atomdata @mol_bonddata
    @mol_bondatoms @mol_constraintdata @mol_angledata @mol_dihedraldata @mol_inversiondata @mol_rigiddata
    @mol_tetherdata @minpos @maxpos @sdata @sdata_names @rdfdata @rdfnames @zdndata @zdnnames
    %element_names @thermodata @thermokeys
    
    @control_allpairs @control_cap @control_closetime @control_collect @control_coulmethod @control_coulprm
    @control_cutglob @control_cutprim @control_delr @control_densvar @control_distan @control_ensemble
    @control_ensembleprm @control_eps @control_equil @control_integrator @control_jobtime @control_minimmethod
    @control_minimprm @control_mult @control_nofic @control_nolink @control_novdw @control_optimmethod
    @control_optimprm @control_pressure @control_print @control_quaternion @control_rdf @control_rdfinterval
    @control_rdfwidth @control_reactmethod @control_reactprm @control_regauss @control_restartmethod
    @control_rlxtol @control_rvdw @control_scale @control_shake @control_shells @control_stack @control_stats
    @control_steps @control_temp @control_thermostat @control_timestep @control_title
    @control_traj @control_zdenprm
    
    read_field_file write_field_file read_config_file copy_config copy_field find_history_timestep
    skip_history_timestep read_history_timestep move_mol move_mol2 add_mol_velocity remap_molecule
    rotate_molecule write_config_file write_history_timestep remove_mol_type remove_mol_entity
    read_statis_timestep read_xyz_timestep write_xyz_timestep connection_depth calc_dihedral_angle
    read_control_file calc_dihedral_angle2 calc_temperature scale_temperature calc_rmsd calc_minmaxpos
    enlarge_system get_atom_list rotate_dihedral remap_molecule2 calc_dsq calc_dsq_orthocell
    calc_dvec_orthocell calc_dvec_orthocell_vec calc_center_of_mass_orthocell calc_dipole_moment
    cut_truncated_octahedron calc_cell_abc calc_cell_vecs print_statis_data_header print_statis_data
    calc_angle calc_angle_orthocell rotate_cell_vmd read_rdfdat_file print_rdfdat_data read_zdndat_file
    print_zdndat_data check_orthogonal clear_field_data find_molname cut_octahedron
    
    @mol2_name @mol2_numatoms @mol2_numbonds @mol2_numsubst @mol2_numfeat @mol2_numsets @mol2_moltype 
    @mol2_chargemethod @mol2_atomdata @mol2_bonddata @mol2_substdata @mol2_atomtypes @mol2_statusbits 
    @mol2_comment @mol2_charge
    
    read_mol2_file mol2_to_cdata write_mol2_file
    
    %gaff_atoms %gaff_bonds %gaff_angles %gaff_dihedrals %gaff_vdwparam
    
    read_gaff_file read_frcmod_file calc_vdw_params
    
    read_lammps_thermo print_thermo read_lammpstrj_timestep
    
    lj_ab_to_sigeps);

our $VERSION = '0.01';

our(@cell, @size, @config_key, @config_title, @periodic_key, @cdata, @frame_number, @frame_numatoms, @frame_timestep,
    @field_numheader, @field_nummols, @field_numatoms, @field_numvdw, @field_numtbp, @field_numfbp, @field_nummetal,
    @field_numextern, @field_units, @field_neutral, @field_title, @field_header, @field_vdwdata, @field_tbpdata,
    @field_fbpdata, @field_metaldata, @field_repsdata, @field_externdata, @field_numatomtypes, @field_atomtypes,
    @field_filename);

our(@mol_name, @mol_charge, @mol_mass, @mol_numents, @mol_numatoms, @mol_numbonds, @mol_numconstraints,
    @mol_numangles, @mol_numdihedrals, @mol_numinversions, @mol_numrigid, @mol_numtether, @mol_atomdata,
    @mol_bonddata, @mol_bondatoms, @mol_constraintdata, @mol_angledata, @mol_dihedraldata, @mol_inversiondata, 
    @mol_rigiddata, @mol_tetherdata, @minpos, @maxpos);

our(@sdata, @rdfdata, @rdfnames, @zdndata, @zdnnames);

our(@control_allpairs, @control_cap, @control_closetime, @control_collect, @control_coulmethod, @control_coulprm,
    @control_cutglob, @control_cutprim, @control_delr, @control_densvar, @control_distan, @control_ensemble,
    @control_ensembleprm, @control_eps, @control_equil, @control_integrator, @control_jobtime, @control_minimmethod,
    @control_minimprm, @control_mult, @control_nofic, @control_nolink, @control_novdw, @control_optimmethod,
    @control_optimprm, @control_pressure, @control_print, @control_quaternion, @control_rdf, @control_rdfinterval,
    @control_rdfwidth, @control_reactmethod, @control_reactprm, @control_regauss, @control_restartmethod,
    @control_rlxtol, @control_rvdw, @control_scale, @control_shake, @control_shells, @control_stack, @control_stats,
    @control_steps, @control_temp, @control_thermostat, @control_timestep, @control_title,
    @control_traj, @control_zdenprm);

our(@mol2_name,@mol2_numatoms,@mol2_numbonds,@mol2_numsubst,@mol2_numfeat,@mol2_numsets,@mol2_moltype,
    @mol2_chargemethod,@mol2_atomdata,@mol2_bonddata,@mol2_substdata,@mol2_atomtypes,@mol2_statusbits,
    @mol2_comment,@mol2_charge);
    
our(%gaff_atoms, %gaff_bonds, %gaff_angles, %gaff_dihedrals, %gaff_vdwparam);

our(@thermodata,@thermokeys);

our %element_names=(
  "C"  => "C",
  "C1" => "C",
  "C2" => "C",
  "C3" => "C",
  "CA" => "C",
  "CP" => "C",
  "CQ" => "C",
  "CC" => "C",
  "CD" => "C",
  "CE" => "C",
  "CF" => "C",
  "CG" => "C",
  "CH" => "C",
  "CX" => "C",
  "CY" => "C",
  "CU" => "C",
  "CV" => "C",
  "CZ" => "C",
  "H"  => "H",
  "H1" => "H",
  "H2" => "H",
  "H3" => "H",
  "H4" => "H",
  "H5" => "H",
  "HA" => "H",
  "HC" => "H",
  "HN" => "H",
  "HO" => "H",
  "HG" => "H",
  "HP" => "H",
  "HS" => "H",
  "HW" => "H",
  "HX" => "H",
  "F"  => "F",
  "CL" => "Cl",
  "BR" => "Br",
  "I"  => "I",
  "N"  => "N",
  "N1" => "N",
  "N2" => "N",
  "N3" => "N",
  "N4" => "N",
  "NA" => "N",
  "NB" => "N",
  "NC" => "N",
  "ND" => "N",
  "NE" => "N",
  "NF" => "N",
  "NH" => "N",
  "NO" => "N",
  "O"  => "O",
  "OH" => "O",
  "OS" => "O",
  "OW" => "O",
  "P"  => "P",
  "P2" => "P",
  "P3" => "P",
  "P4" => "P",
  "P5" => "P",
  "PB" => "P",
  "PC" => "P",
  "PD" => "P",
  "PE" => "P",
  "PF" => "P",
  "PX" => "P",
  "PY" => "P",
  "S"  => "S",
  "S2" => "S",
  "S4" => "S",
  "S6" => "S",
  "SH" => "S",
  "SS" => "S",
  "SX" => "S",
  "SY" => "S",
  "AL" => "Al",
  "OA" => "O",
  "OX" => "O",
  "XX" => "X",
  "NA+" => "Na"
);

#                      0      1       2     3      4      5       6      7       8     9      10     11     12     13    14     15
our @sdata_names = ('FRAME','TIME','ETOT','TTOT','ECFG','EVDW','ECOUL','EBND','EANG','EDIH','ETET','HTOT','TTOT','TROT','VOL','PRESS');

sub clear_field_data {
  my $fi       = $_[0]; # index for array of field-data
  if((not defined($fi)) or $fi<0) {
    print "**** error: index \$fi must not be smaller than zero in sub clear_field_data\n";
    return -1;
  }
  
  $field_numheader[$fi]    = 0;
  $field_nummols[$fi]      = 0;
  $field_numatoms[$fi]     = 0;
  $field_numatomtypes[$fi] = 0;
  $field_numvdw[$fi]       = 0;
  $field_numtbp[$fi]       = 0;
  $field_numfbp[$fi]       = 0;
  $field_nummetal[$fi]     = 0;
  $field_numextern[$fi]    = 0;
  $field_neutral[$fi]      = 0;
  $field_units[$fi]        = undef;
  undef @{$field_header[$fi]};
  undef @{$mol_name[$fi]};
  undef @{$mol_mass[$fi]};
  undef @{$mol_charge[$fi]};
  undef @{$mol_numents[$fi]};
  undef @{$mol_numatoms[$fi]};
  undef @{$mol_numbonds[$fi]};
  undef @{$mol_numconstraints[$fi]};
  undef @{$mol_numangles[$fi]};
  undef @{$mol_numdihedrals[$fi]};
  undef @{$mol_numinversions[$fi]};
  undef @{$mol_numrigid[$fi]};
  undef @{$mol_numtether[$fi]};
  undef @{$mol_atomdata[$fi]};
  undef @{$mol_bonddata[$fi]};
  undef @{$mol_bondatoms[$fi]};
  undef @{$mol_constraintdata[$fi]};
  undef @{$mol_angledata[$fi]};
  undef @{$mol_dihedraldata[$fi]};
  undef @{$mol_inversiondata[$fi]};
  undef @{$mol_rigiddata[$fi]};
  undef @{$mol_tetherdata[$fi]};
  undef @{$field_vdwdata[$fi]};
  undef @{$field_tbpdata[$fi]};
  undef @{$field_fbpdata[$fi]};
  undef @{$field_metaldata[$fi]};
  undef @{$field_repsdata[$fi]};
  undef @{$field_externdata[$fi]};
  undef @{$field_atomtypes[$fi]};
  
}

sub read_field_file {
  my $filename = $_[0];
  my $fi       = $_[1]; # index for array of field-data
  my($filehandle,$keyword,$t,$i,$j,$numentries,$finished,$type, $addtype);
  if(not defined($fi) or $fi<0) {
    print "**** error: index \$fi must not be smaller than zero in sub read_field_file\n";
    return -1;
  }
  $field_filename[$fi]=$filename;
  if(not open($filehandle, "<", $filename)) {
    print "**** error: Can't open FIELD file \"$filename\": $!\n";
    return 1;
  }
  
  clear_field_data($fi);
  my @linedata;
  
  ############## read header from FIELD file ##########################
  $field_title[$fi] = <$filehandle>;
  $field_title[$fi] =~ s/^\s+//;
  $field_title[$fi] =~ s/\s+$//;
  while(<$filehandle>) {
    $field_header[$fi][$field_numheader[$fi]] = $_;
    $field_numheader[$fi]++;
    if(/UNITS/i) {
      ($field_units[$fi]) = /UNITS\s+(\S+)/i;
    } elsif(/NEUT/i) {
      $field_neutral[$fi]   = 1;
    } elsif(/MOLECULES/i) {
      ($field_nummols[$fi]) = /(\d+)/;
      last;
    }
  }
  if($field_nummols[$fi]==0) {
    print "**** error: no molecules found in FIELD file $filename.\n";
    close($filehandle);
    return 1;
  }

  ############## read molecule data from FIELD file ###################
  #keywords SHELL and PMF are not implemented yet
  for(my $t=0;$t<$field_nummols[$fi];$t++) {
    $numentries                  = 0;
    $mol_numents[$fi][$t]        = 0;
    $mol_numatoms[$fi][$t]       = 0;
    $mol_numbonds[$fi][$t]       = 0;
    $mol_numconstraints[$fi][$t] = 0;
    $mol_numangles[$fi][$t]      = 0;
    $mol_numdihedrals[$fi][$t]   = 0;
    $mol_numinversions[$fi][$t]  = 0;
    $mol_numrigid[$fi][$t]       = 0;
    $mol_numtether[$fi][$t]      = 0;
    $mol_mass[$fi][$t]           = 0;
    $mol_charge[$fi][$t]         = 0;
    @{$mol_atomdata[$fi][$t]}       = ();
    @{$mol_bonddata[$fi][$t]}       = ();
    @{$mol_bondatoms[$fi][$t]}      = ();
    @{$mol_constraintdata[$fi][$t]} = ();
    @{$mol_angledata[$fi][$t]}      = ();
    @{$mol_dihedraldata[$fi][$t]}   = ();
    @{$mol_inversiondata[$fi][$t]}  = ();
    @{$mol_rigiddata[$fi][$t]}      = ();
    @{$mol_tetherdata[$fi][$t]}     = ();
    $finished = 0;
    $_ = <$filehandle>;
    $mol_name[$fi][$t] = $_;
    $mol_name[$fi][$t] =~ s/^\s+//;
    $mol_name[$fi][$t] =~ s/\s+$//;
    $_ = <$filehandle>;
    if(not /NUMMOLS/i) {
      print "**** error: expected entry \"NUMMOLS\" after molecule name $mol_name[$fi][$t] in FIELD file $filename\n";
      close($filehandle);
      return 1;
    }
    ($mol_numents[$fi][$t]) = /(\d+)/;
    while(<$filehandle>) {
      if(($keyword)=/^\s*(ATOMS|SHELL|BONDS|CONSTRAINTS|ANGLES|DIHEDRALS|INVERSIONS|RIGID|TETH|FINISH)/i) {
	($numentries)=/$keyword\s+(\d+)/;
	$keyword=uc($keyword);
      } else {
	print "**** error: unknown keyword in FIELD file $filename on line $.: $keyword\n $_ missing \"FINISH\"?\n";
	close($filehandle);
	return 1;
      }
      if(/FINISH/i) {
	$finished = 1;
	last;
      }
      for($i=0;$i<$numentries;$i++) {
	if(not $_=<$filehandle>) {
	  print "**** error: unexpected end of FIELD file $filename on line $. in section $keyword of molecule $mol_name[$fi][$t]\n";
	  return 1;
	}
	$_ =~ s/^\s+//; $_ =~ s/\s+$//;
	@linedata = split(/\s+/, $_);
	if($keyword=~/ATOMS/) {
	  # index    0    1    2      3    4    5
	  # format:  name mass charge nrep ifrz igrp
	  $mol_mass[$fi][$t]   += $linedata[1];
	  $mol_charge[$fi][$t] += $linedata[2];
	  for($j=0;$j<$linedata[3];$j++) {
	    @{$mol_atomdata[$fi][$t][$mol_numatoms[$fi][$t]]} = @linedata;
	    $addtype=1;
	    foreach $type (@{$field_atomtypes[$fi]}) {
	      if($type eq $linedata[0]) {$addtype=0; last;}
	    }
	    if($addtype) {
	      push(@{$field_atomtypes[$fi]},$linedata[0]);
	      $field_numatomtypes[$fi]++;
	    }
	    $mol_numatoms[$fi][$t]++;
	  }
	  $i+=$linedata[3]-1;
	} elsif($keyword=~/BONDS/) {
	  # index    0  1  2  3  4  5  6
	  # format:  k a1 a2 v1 v2 v3 v4
	  $linedata[1]--;
	  $linedata[2]--;
	  @{$mol_bonddata[$fi][$t][$mol_numbonds[$fi][$t]]} = @linedata;
	  push (@{$mol_bondatoms[$fi][$t][$linedata[1]]}, $linedata[2]);
	  push (@{$mol_bondatoms[$fi][$t][$linedata[2]]}, $linedata[1]);
	  $mol_numbonds[$fi][$t]++;
	} elsif($keyword=~/CONSTRAINTS/) {
	  $linedata[0]--;
	  $linedata[1]--;
	  # index    0  1  2
	  # format:  a1 a2 v1
	  @{$mol_constraintdata[$fi][$t][$mol_numconstraints[$fi][$t]]} = @linedata;
	  push (@{$mol_bondatoms[$fi][$t][$linedata[0]]}, $linedata[1]);
	  push (@{$mol_bondatoms[$fi][$t][$linedata[1]]}, $linedata[0]);
	  $mol_numconstraints[$fi][$t]++;
	} elsif($keyword=~/ANGLES/) {
	  $linedata[1]--;
	  $linedata[2]--;
	  $linedata[3]--;
	  # index    0  1  2  3  4  5  6  7  8  9
	  # format:  k a1 a2 a3 v1 v2 v3 v4 v5 v6
	  @{$mol_angledata[$fi][$t][$mol_numangles[$fi][$t]]} = @linedata;
	  $mol_numangles[$fi][$t]++;
	} elsif($keyword=~/DIHEDRALS/) {
	  $linedata[1]--;
	  $linedata[2]--;
	  $linedata[3]--;
	  $linedata[4]--;
	  # index    0  1  2  3  4  5  6  7  8  9
	  # format:  k a1 a2 a3 a4 v1 v2 v3 v4 v5
	  @{$mol_dihedraldata[$fi][$t][$mol_numdihedrals[$fi][$t]]} = @linedata;
	  $mol_numdihedrals[$fi][$t]++;
	} elsif($keyword=~/INVERSIONS/) {
	  $linedata[1]--;
	  $linedata[2]--;
	  $linedata[3]--;
	  $linedata[4]--;
	  # index    0  1  2  3  4  5  6
	  # format:  k a1 a2 a3 a4 v1 v2
	  @{$mol_inversiondata[$fi][$t][$mol_numinversions[$fi][$t]]} = @linedata;
	  $mol_numinversions[$fi][$t]++;
	} elsif($keyword=~/RIGID/) {
	  # index    0       1-16
	  # format:  #at  at1-atn
	  @{$mol_rigiddata[$fi][$t][$mol_numrigid[$fi][$t]]} = ($linedata[0]);
	  for($j=1;$j<@linedata;$j++) {
	    push(@{$mol_rigiddata[$fi][$t][$mol_numrigid[$fi][$t]]},$linedata[$j]-1);
	  }
	  $mol_numrigid[$fi][$t]++;
	} elsif($keyword=~/TETH/) {
	  $linedata[1]--;
	  # index    0  1  2  3  4  5
	  # format:  k a1 v1 v2 v3 v4
	  @{$mol_tetherdata[$fi][$t][$mol_numtether[$fi][$t]]} = @linedata;
	  $mol_numtether[$fi][$t]++;
	} # end if
      } # end for $i=$numentries
    } # end while(<$filehandle>)
    if(not $finished) {
      print "**** error reading FIELD file $filename on line $.!\n";
      close($filehandle);
      return 1;
    }
    $field_numatoms[$fi] += $mol_numents[$fi][$t]*$mol_numatoms[$fi][$t];
  } # end for $t
  
  ############## read nonbonding VDW data from FIELD file #############
  #Metal potentials are not implemented yet
  while(<$filehandle>) {
    if(($keyword)=/^\s*(VDW|TBP|FBP|METAL|EXTERN|REPSPHERE|CLOSE)/i) {
      ($numentries)=/$keyword\s+(\d+)/;
      $numentries=1 if(not defined($numentries));
      $keyword=uc($keyword);
    } else {
      print "**** error: unknown keyword in FIELD file $filename on line $.: $keyword\n $_ missing \"CLOSE\"?\n";
      close($filehandle);
      return 1;
    }
    if(/CLOSE/i) {
      $finished = 1;
      last;
    }
    for($i=0;$i<$numentries;$i++) {
      if(not $_=<$filehandle>) {
	print "**** error: unexpected end of FIELD file $filename on line $. in section $keyword\n";
	return 1;
      }
      $_ =~ s/^\s+//; $_ =~ s/\s+$//;
      @linedata = split(/\s+/, $_);
      if($keyword=~/VDW/) {
	# index    0  1  2  3  4  5  6  7
	# format: a1 a2  k v1 v2 v3 v4 v5
	@{$field_vdwdata[$fi][$i]} = @linedata;
	$field_numvdw[$fi]++;
      } elsif($keyword=~/TBP/) {
	# index    0  1  2  3  4  5  6  7  8
	# format: a1 a2 a3  k v1 v2 v3 v4 v5
	@{$field_tbpdata[$fi][$i]} = @linedata;
	$field_numtbp[$fi]++;
      } elsif($keyword=~/FBP/) {
	# index    0  1  2  3  4  5  6  7
	# format: a1 a2 a3 a4  k v1 v2 v3
	@{$field_fbpdata[$fi][$i]} = @linedata;
	$field_numfbp[$fi]++;
      } elsif($keyword=~/METAL/) {
	# index    0  1  2  3  4  5  6  7  8  9
	# format: a1 a2  k v1 v2 v3 v4 v5 v6 v7
	@{$field_metaldata[$fi][$i]} = @linedata;
	$field_nummetal[$fi]++;
      } elsif($keyword=~/EXTERN/) {
	# index    0  1  2  3  4  5
	# format:  k v1 v2 v3 v4 v5
	$field_externdata[$fi][$i][0] = $linedata[0]; # key
	if(not $_=<$filehandle>) {
	  print "**** error: unexpected end of FIELD file $filename on line $. in section $keyword\n";
	  return 1;
	}
	$_ =~ s/^\s+//; $_ =~ s/\s+$//;
	push(@{$field_externdata[$fi][$i]},split(/\s+/,$_));
	$field_numextern[$fi]++;
      } elsif($keyword=~/REPSPHERE/) {
	# index    0  1  2  3  4  5
	# format:  k v1 v2 v3 v4 v5
	@{$field_repsdata[$fi][$i]} = @linedata;
      } else {
	print "**** error: unknown keyword in FIELD file $filename on line $.: $keyword\n";
	close($filehandle);
	return 1;
      } # end if
    } # end for $i=$numentries
  } # end while(<$filehandle>)
  close($filehandle);
  return 0;
} # end subroutine read_field_file


sub write_field_file {
  my $filename = $_[0];
  my $fi       = $_[1]; # index for array of field-data
  return 1 if($fi<0);
  my($filehandle,$t,$i,$j);
  if(not open($filehandle, ">", $filename)) {
    print "**** error: Can't open FIELD file \"$filename\": $!\n";
    return 1;
  }
  print $filehandle "$field_title[$fi]\n";
  print $filehandle "UNITS $field_units[$fi]\n" if(defined($field_units[$fi]));
  print $filehandle "NEUT\n" if($field_neutral[$fi]);
  print $filehandle "MOLECULES $field_nummols[$fi]\n";
  for($t=0;$t<$field_nummols[$fi];$t++) {
    print $filehandle "$mol_name[$fi][$t]\n";
    print $filehandle "NUMMOLS $mol_numents[$fi][$t]\n";
    print $filehandle "ATOMS $mol_numatoms[$fi][$t]\n";
    for($i=0;$i<$mol_numatoms[$fi][$t];$i++) {
      printf $filehandle "%-8s %9.4f %12.8f %5u %5u %5u\n", @{$mol_atomdata[$fi][$t][$i]}[0..5];
      $i+=$mol_atomdata[$fi][$t][$i][3]-1
    }
    if($mol_numbonds[$fi][$t]>0) {
      print $filehandle "BONDS $mol_numbonds[$fi][$t]\n";
      for($i=0;$i<$mol_numbonds[$fi][$t];$i++) {
	printf $filehandle "%-5s %6u %6u %12.6f %12.6f\n",$mol_bonddata[$fi][$t][$i][0],
	$mol_bonddata[$fi][$t][$i][1]+1,$mol_bonddata[$fi][$t][$i][2]+1, @{$mol_bonddata[$fi][$t][$i]}[3..4];
      }
    }
    if($mol_numconstraints[$fi][$t]>0) {
      print $filehandle "CONSTRAINTS $mol_numconstraints[$fi][$t]\n";
      for($i=0;$i<$mol_numconstraints[$fi][$t];$i++) {
	printf $filehandle "%6u %6u %12.6f\n",$mol_constraintdata[$fi][$t][$i][0]+1,
	$mol_constraintdata[$fi][$t][$i][1]+1,$mol_constraintdata[$fi][$t][$i][2];
      }
    }
    if($mol_numangles[$fi][$t]>0) {
      print $filehandle "ANGLES $mol_numangles[$fi][$t]\n";
      for($i=0;$i<$mol_numangles[$fi][$t];$i++) {
	printf $filehandle "%-5s %6u %6u %6u %12.5f %12.5f\n", $mol_angledata[$fi][$t][$i][0],
	$mol_angledata[$fi][$t][$i][1]+1,$mol_angledata[$fi][$t][$i][2]+1,
	$mol_angledata[$fi][$t][$i][3]+1,@{$mol_angledata[$fi][$t][$i]}[4..5];
      }
    }
    if($mol_numdihedrals[$fi][$t]>0) {
      print $filehandle "DIHEDRALS $mol_numdihedrals[$fi][$t]\n";
      for($i=0;$i<$mol_numdihedrals[$fi][$t];$i++) {
	printf $filehandle "%-5s %6u %6u %6u %6u %12.5f %12.5f %8.3f %8.3f %8.3f\n",
	$mol_dihedraldata[$fi][$t][$i][0], $mol_dihedraldata[$fi][$t][$i][1]+1,$mol_dihedraldata[$fi][$t][$i][2]+1,
	$mol_dihedraldata[$fi][$t][$i][3]+1,$mol_dihedraldata[$fi][$t][$i][4]+1,
	@{$mol_dihedraldata[$fi][$t][$i]}[5..9];
      }
    }
    if($mol_numinversions[$fi][$t]>0) {
      print $filehandle "INVERSIONS $mol_numinversions[$fi][$t]\n";
      for($i=0;$i<$mol_numinversions[$fi][$t];$i++) {
	printf $filehandle "%-5s %6u %6u %6u %6u %12.5f %12.5f %5.3f %5.3f %5.3f\n",
	$mol_inversiondata[$fi][$t][$i][0], $mol_inversiondata[$fi][$t][$i][1]+1,$mol_inversiondata[$fi][$t][$i][2]+1,
	$mol_inversiondata[$fi][$t][$i][3]+1,$mol_inversiondata[$fi][$t][$i][4]+1,
	@{$mol_inversiondata[$fi][$t][$i]}[5..6];
      }
    }
    if($mol_numrigid[$fi][$t]>0) {
      print $filehandle "RIGID $mol_numrigid[$fi][$t]\n";
      for($i=0;$i<$mol_numrigid[$fi][$t];$i++) {
	printf $filehandle "%-6u", $mol_rigiddata[$fi][$t][$i][0];
	for($j=1;$j<@{$mol_rigiddata[$fi][$t][$i]};$j++) {
	  printf $filehandle " %6u", $mol_rigiddata[$fi][$t][$i][$j]+1;
	}
	print $filehandle "\n";
      }
    }
    if($mol_numtether[$fi][$t]>0) {
      print $filehandle "TETH $mol_numtether[$fi][$t]\n";
      for($i=0;$i<$mol_numtether[$fi][$t];$i++) {
	printf $filehandle "%-5s %6u", $mol_tetherdata[$fi][$t][$i][0], $mol_tetherdata[$fi][$t][$i][1]+1;
	for($j=2;$j<@{$mol_tetherdata[$fi][$t][$i]};$j++) {
	  printf $filehandle " %12.5f", $mol_tetherdata[$fi][$t][$i][$j];
	}
	print $filehandle "\n";
      }
    }
    print $filehandle "FINISH\n";
  }
  if($field_numvdw[$fi]>0) {
    print $filehandle "VDW $field_numvdw[$fi]\n";
    for($i=0;$i<$field_numvdw[$fi];$i++) {
      printf $filehandle "%-7s %-7s %5s",@{$field_vdwdata[$fi][$i]}[0..2];
      for($j=3;$j<@{$field_vdwdata[$fi][$i]};$j++) {
	printf $filehandle " %15.9G",$field_vdwdata[$fi][$i][$j];
      }
      print $filehandle "\n";
    }
  }
  if($field_numtbp[$fi]>0) {
    print $filehandle "TBP $field_numtbp[$fi]\n";
    for($i=0;$i<$field_numtbp[$fi];$i++) {
      printf $filehandle "%-7s %-7s %-7s %5s",@{$field_tbpdata[$fi][$i]}[0..3];
      for($j=4;$j<@{$field_tbpdata[$fi][$i]};$j++) {
	printf $filehandle " %12.5f",$field_tbpdata[$fi][$i][$j];
      }
      print $filehandle "\n";
    }
  }
  if($field_numfbp[$fi]>0) {
    print $filehandle "FBP $field_numfbp[$fi]\n";
    for($i=0;$i<$field_numfbp[$fi];$i++) {
      printf $filehandle "%-7s %-7s %-7s %-7s %5s",@{$field_fbpdata[$fi][$i]}[0..4];
      for($j=5;$j<@{$field_fbpdata[$fi][$i]};$j++) {
	printf $filehandle " %12.5f",$field_fbpdata[$fi][$i][$j];
      }
      print $filehandle "\n";
    }
  }
  if($field_nummetal[$fi]>0) {
    print $filehandle "METAL $field_nummetal[$fi]\n";
    for($i=0;$i<$field_nummetal[$fi];$i++) {
      printf $filehandle "%-7s %-7s %5s",@{$field_metaldata[$fi][$i]}[0..2];
      for($j=3;$j<@{$field_metaldata[$fi][$i]};$j++) {
	printf $filehandle " %12.5f",$field_metaldata[$fi][$i][$j];
      }
      print $filehandle "\n";
    }
  }
  if($field_numextern[$fi]>0) {
    print $filehandle "EXTERN $field_numextern[$fi]\n";
    for($i=0;$i<$field_numextern[$fi];$i++) {
      print $filehandle "$field_externdata[$fi][$i][0] $#{$field_externdata[$fi][$i]}\n";
      for($j=1;$j<@{$field_externdata[$fi][$i]};$j++) {
	printf $filehandle "%14.7g",$field_externdata[$fi][$i][$j];
      }
      print $filehandle "\n";
    }
  }
  print $filehandle "CLOSE";
  close($filehandle);
  return 0;
}

sub find_molname {
  # find a molecule name in a FIELD file
  my $fi       = $_[0]; # index for array of field-data
  my $name     = $_[1]; # index for array of field-data
  return -1 if($fi<0);
  return -1 if(not defined($field_nummols[$fi]));
  for(my $t=0;$t<$field_nummols[$fi];$t++) {
    return $t if($mol_name[$fi][$t]=~/^$name$/i)
  }
  return -1;
}

sub read_control_file {
  my $filename = $_[0];
  my $ci       = $_[1]; #index of CONTROL file
  my($filehandle,@linedata,$i);
  
  return 1 if($ci<0);
  
  #reset values
  $control_allpairs[$ci] = undef;
  $control_cap[$ci] = undef;
  $control_closetime[$ci] = undef;
  $control_collect[$ci] = undef;
  $control_coulmethod[$ci] = undef;
  @{$control_coulprm[$ci]} = ();
  $control_cutglob[$ci] = undef;
  $control_cutprim[$ci] = undef;
  $control_delr[$ci] = undef;
  $control_densvar[$ci] = undef;
  $control_distan[$ci] = undef;
  $control_ensemble[$ci] = undef;
  @{$control_ensembleprm[$ci]} = ();
  $control_eps[$ci] = undef;
  $control_equil[$ci] = undef;
  $control_integrator[$ci] = undef;
  $control_jobtime[$ci] = undef;
  $control_minimmethod[$ci] = undef;
  @{$control_minimprm[$ci]} = ();
  $control_mult[$ci] = undef;
  $control_nofic[$ci] = undef;
  $control_nolink[$ci] = undef;
  $control_novdw[$ci] = undef;
  $control_optimmethod[$ci] = undef;
  $control_optimprm[$ci] = undef;
  $control_pressure[$ci] = undef;
  $control_print[$ci] = undef;
  $control_quaternion[$ci] = undef;
  $control_rdf[$ci] = undef;
  $control_rdfinterval[$ci] = undef;
  $control_rdfwidth[$ci] = undef;
  $control_reactmethod[$ci] = undef;
  $control_reactprm[$ci] = undef;
  $control_regauss[$ci] = undef;
  $control_restartmethod[$ci] = undef;
  $control_rlxtol[$ci] = undef;
  $control_rvdw[$ci] = undef;
  $control_scale[$ci] = undef;
  $control_shake[$ci] = undef;
  $control_shells[$ci] = undef;
  $control_stack[$ci] = undef;
  $control_stats[$ci] = undef;
  $control_steps[$ci] = undef;
  $control_temp[$ci] = undef;
  $control_thermostat[$ci] = undef;
  $control_timestep[$ci] = undef;
  $control_title[$ci] = undef;
  @{$control_traj[$ci]} = ();
  @{$control_zdenprm[$ci]} = ();
  
  if(not open($filehandle, "<", $filename)) {
    print "**** error: Can't open CONFIG file \"$filename\": $!\n";
    return 1;
  }
  #read header of CONFIG file
  $control_title[$ci]=<$filehandle>;
  while(<$filehandle>) {
    $_ =~ s/^\s+//; $_ =~ s/\s+$//;
    next if(/^#/);
    next if(length($_)==0);
    @linedata = split(/\s+/, $_);
    for($i=0;$i<@linedata;$i++) {
      if($linedata[$i]=~/^#/) {
	splice(@linedata,$i,$#linedata-$i+1);
	last;
      }
    }
    if($linedata[0]=~/^all/i and $linedata[1]=~/^pairs/i) {
      $control_allpairs[$ci] = 1;
    } elsif($linedata[0]=~/^cap/i) {
      $control_cap[$ci] = $linedata[1];
      $control_cap[$ci] = 1000 if($#linedata<1);
    } elsif($linedata[0]=~/^close/i and $linedata[1]=~/^time/i) {
      $control_closetime[$ci] = $linedata[2];
    } elsif($linedata[0]=~/^collect/i) {
      $control_collect[$ci] = 1;
    } elsif($linedata[0]=~/^coul/i) {
      $control_coulmethod[$ci] = 'coul';
      @{$control_coulprm[$i]}  = ();
    } elsif($linedata[0]=~/^cut/i) {
      $control_cutglob[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^densvar/i) {
      $control_densvar[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^distan/i) {
      $control_distan[$ci] = 1;
    } elsif($linedata[0]=~/^delr/i) {
      $control_delr[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^ensemble/i) {
      if($linedata[1]=~/^nvt/i) {
	$control_ensemble[$ci]       = 'nvt';
	$control_thermostat[$ci]     = $linedata[2];
	@{$control_ensembleprm[$ci]} = subarray(\@linedata,3,$#linedata);
      } elsif($linedata[1]=~/^npt/i) {
	$control_ensemble[$ci]       = 'npt';
	$control_thermostat[$ci]     = $linedata[2];
	@{$control_ensembleprm[$ci]} = subarray(\@linedata,3,$#linedata);
      } elsif($linedata[1]=~/^nst/i) {
	$control_ensemble[$ci]       = 'nst';
	$control_thermostat[$ci]     = $linedata[2];
	@{$control_ensembleprm[$ci]} = subarray(\@linedata,3,$#linedata);
      } else {
	$control_ensemble[$ci]     = $linedata[1];
      }
    } elsif($linedata[0]=~/^eps/i) {
      $control_eps[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^equil/i) {
      $control_equil[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^ewald/i) { # precision or sum
      $control_coulmethod[$ci]  = 'ewald '.$linedata[1];
      @{$control_coulprm[$ci]} = subarray(\@linedata,2,$#linedata);
    } elsif($linedata[0]=~/^finish/i) {
      last;
    } elsif($linedata[0]=~/^hke/i) { # precision or sum
      $control_coulmethod[$ci]  = 'hke '.$linedata[1];
      @{$control_coulprm[$ci]} = subarray(\@linedata,2,$#linedata);
    } elsif($linedata[0]=~/^integrator/i) {
      $control_integrator[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^job/i and $linedata[1]=~/^time/i) {
      $control_jobtime[$ci] = $linedata[2];
    } elsif($linedata[0]=~/^minim/i) { # energy force or position
      $control_minimmethod[$ci] = $linedata[1];
      @{$control_minimprm[$ci]} = subarray(\@linedata,2,$#linedata);
    } elsif($linedata[0]=~/^mult/i) {
      $control_mult[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^no/i) {
      if($linedata[1]=~/^elec/i) {
	$control_coulmethod[$ci]    = undef;
	@{$control_coulmethod[$ci]} = ();
      } elsif($linedata[1]=~/^fic/i) {
	$control_nofic[$ci]  = $linedata[2];
      } elsif($linedata[1]=~/^link/i) {
	$control_nolink[$ci] = 1;
      } elsif($linedata[1]=~/^vdw/i) {
	$control_novdw[$ci]  = 1;
      }
    } elsif($linedata[0]=~/^optim/i) { # energy force or position
      $control_optimmethod[$ci] = $linedata[1];
      $control_optimprm[$ci]    = $linedata[2];
    } elsif($linedata[0]=~/^pres/i) {
      $control_pressure[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^prim/i) {
      $control_cutprim[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^print/i) {
      if($linedata[1]=~/^rdf/i) {
	$control_rdf[$ci] = 1;
      } else {
        $_=~/(\d+)/;
	$control_print[$ci] = $linedata[1];
      }
    } elsif($linedata[0]=~/^put/i) { # electrostatic shells on cores
      $control_shells[$ci] = 1;
    } elsif($linedata[0]=~/^quaternion/i) {
      $control_quaternion[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^rdf/i) {
      $control_rdf[$ci]    = 1;
      $control_rdfinterval[$ci] = $linedata[1];
      $control_rdfwidth[$ci]    = $linedata[2];
    } elsif($linedata[0]=~/^reaction/i) { # precision damped or undef
      if($#linedata>0) {
	$control_reactmethod[$ci] = 'reaction '.$linedata[1];
	$control_reactprm[$ci]    = $linedata[2];
      }else{
	$control_reactmethod[$ci] = 'reaction';
	$control_reactprm[$ci]    = undef;
      }
    } elsif($linedata[0]=~/^regauss/i) {
      $control_regauss[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^restart/i) { # precision damped or undef
      if($#linedata>0) {
	$control_restartmethod[$ci] = 'restart '.$linedata[1];
      }else{
	$control_restartmethod[$ci] = 'restart';
      }
    } elsif($linedata[0]=~/^rlxtol/i) {
      $control_rlxtol[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^rvdw/i) {
      $control_rvdw[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^scale/i) {
      $control_scale[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^shake/i) {
      $control_shake[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^shift/i) { # precision or sum
      $control_coulmethod[$ci]  = 'shift '.$linedata[1];
      @{$control_coulprm[$ci]} = subarray(\@linedata,2,$#linedata);
    } elsif($linedata[0]=~/^spme/i) { # precision or sum
      $control_coulmethod[$ci]  = 'spme '.$linedata[1];
      @{$control_coulprm[$ci]} = subarray(\@linedata,2,$#linedata);
    } elsif($linedata[0]=~/^stack/i) {
      $_=~/(\d+)/;
      $control_stack[$ci] = $1;
    } elsif($linedata[0]=~/^stats/i) {
      $_=~/(\d+)/;
      $control_stats[$ci] = $1;
    } elsif($linedata[0]=~/^steps/i) {
      $_=~/(\d+)/;
      $control_steps[$ci] = $1;
    } elsif($linedata[0]=~/^temp/i) {
      $_=~/(\d+)/;
      $control_temp[$ci] = $1;
    } elsif($linedata[0]=~/^traj/i) {
      @{$control_traj[$ci]} = subarray(\@linedata,1,3);
    } elsif($linedata[0]=~/^timestep/i) {
      $control_timestep[$ci] = $linedata[1];
    } elsif($linedata[0]=~/^zden/i) {
      @{$control_zdenprm[$ci]} = subarray(\@linedata,1,$#linedata);
    } elsif($linedata[0]=~/^zero/i) {
      $control_temp[$ci] = 0;
    } else {
      print "unknown directive in line $_\n";
    }
  }
} # end sub read_control_file

sub read_config_file {
  # opens and reads data from CONFIG file either using ($fi>=0) or disregarding ($fi<0)
  # information from previously read FIELD file and writes content to cdata array
  # example call:
  # $err = read_config_file("testfile.cfg",-1,0);
  # indices of cdata array
  # 0 1 2 3  4  5  6  7  8  9    10       11    12   13
  # x y z vx vy vz fx fy fz name molindex index mass charge
  my $filename = $_[0];
  my $fi       = $_[1]; #index of corresponding FIELD file
  my $ci       = $_[2];
  my($i,$j,$t,$m,$a,$molindex, $filehandle);
  
  return 1 if($ci<0);
  
  #reset values
  undef @{$cdata[$ci]};
  undef @{$cell[$ci]};
  undef @{$size[$ci]};
  $molindex      = 0;
  
  if(not open($filehandle, "<", $filename)) {
    print "**** error: Can't open CONFIG file \"$filename\": $!\n";
    return 1;
  }
  #read header of CONFIG file
  $config_title[$ci]=<$filehandle>;
  $_=<$filehandle>;
  ($config_key[$ci],$periodic_key[$ci]) = /^\s+(\S+)\s+(\S+)/;
  if($periodic_key[$ci]>0) {
    for($j=0;$j<3;$j++) {
      $_=<$filehandle>;
      ($cell[$ci][$j][0], $cell[$ci][$j][1], $cell[$ci][$j][2]) = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
    }
    $size[$ci][0] = ($cell[$ci][0][0] + $cell[$ci][1][0] + $cell[$ci][2][0])/2;
    $size[$ci][1] = ($cell[$ci][0][1] + $cell[$ci][1][1] + $cell[$ci][2][1])/2;
    $size[$ci][2] = ($cell[$ci][0][2] + $cell[$ci][1][2] + $cell[$ci][2][2])/2;
  }
  
  @{$minpos[$ci]}=( 9e20, 9e20, 9e20);
  @{$maxpos[$ci]}=(-9e20,-9e20,-9e20);
  $frame_numatoms[$ci] = 0;
  if($fi>=0) { # read using FIELD information
    for($t=0;$t<$field_nummols[$fi];$t++) {
      @{$cdata[$ci][$t]}=();
      for($m=0;$m<$mol_numents[$fi][$t];$m++) {
	@{$cdata[$ci][$t][$m]}=();
	for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
	@{$cdata[$ci][$t][$m][$a]}=();
	  if(not $_=<$filehandle>) {print "**** unexpected end of CONFIG file (t=$t,m=$m,a=$a,i=1)!\n"; return 1}
	  ($cdata[$ci][$t][$m][$a][9],$cdata[$ci][$t][$m][$a][11]) = /^\s*(\S+)\s*(\S*)/;
	  if($cdata[$ci][$t][$m][$a][9] ne $mol_atomdata[$fi][$t][$a][0]) {
	    print "**** error: atom name $cdata[$ci][$t][$m][$a][9] on line $. in CONFIG file '$filename' does not match name ",
	    "$mol_atomdata[$fi][$t][$a][0] in FIELD data for molecule $mol_name[$fi][$t] from file ",
	    "'$field_filename[$fi]'\n";
	    exit 1;
	  }
	  if(not $_=<$filehandle>) {print "**** unexpected end of CONFIG file (t=$t,m=$m,a=$a,i=2)!\n"; return 1}
	  ($cdata[$ci][$t][$m][$a][0],$cdata[$ci][$t][$m][$a][1],$cdata[$ci][$t][$m][$a][2])
	   = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
	  if($config_key[$ci]>0) { #read velocity
	    if(not $_=<$filehandle>) {print "**** unexpected end of CONFIG file (t=$t,m=$m,a=$a,i=3)!\n"; return 1}
	    ($cdata[$ci][$t][$m][$a][3],$cdata[$ci][$t][$m][$a][4],$cdata[$ci][$t][$m][$a][5])
	    = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
	    if($config_key[$ci]>1) { #read force
	      if(not $_=<$filehandle>) {print "**** unexpected end of CONFIG file (t=$t,m=$m,a=$a,i=4)!\n"; return 1}
	      ($cdata[$ci][$t][$m][$a][6],$cdata[$ci][$t][$m][$a][7],$cdata[$ci][$t][$m][$a][8])
	      = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
	    }
	  }
	  $cdata[$ci][$t][$m][$a][10]=$molindex;
	  $cdata[$ci][$t][$m][$a][12]=$mol_atomdata[$fi][$t][$a][1];
	  $cdata[$ci][$t][$m][$a][13]=$mol_atomdata[$fi][$t][$a][2];
	  $minpos[$ci][0] = $cdata[$ci][$t][$m][$a][0] if($cdata[$ci][$t][$m][$a][0]<$minpos[$ci][0]);
	  $minpos[$ci][1] = $cdata[$ci][$t][$m][$a][1] if($cdata[$ci][$t][$m][$a][1]<$minpos[$ci][1]);
	  $minpos[$ci][2] = $cdata[$ci][$t][$m][$a][2] if($cdata[$ci][$t][$m][$a][2]<$minpos[$ci][2]);
	  $maxpos[$ci][0] = $cdata[$ci][$t][$m][$a][0] if($cdata[$ci][$t][$m][$a][0]>$maxpos[$ci][0]);
	  $maxpos[$ci][1] = $cdata[$ci][$t][$m][$a][1] if($cdata[$ci][$t][$m][$a][1]>$maxpos[$ci][1]);
	  $maxpos[$ci][2] = $cdata[$ci][$t][$m][$a][2] if($cdata[$ci][$t][$m][$a][2]>$maxpos[$ci][2]);
	}
	$molindex++;
      }
    }
    $frame_numatoms[$ci] = $field_numatoms[$fi];
    # check whether the end of the file was reached
    if($_=<$filehandle>) {
      chomp($_);
      if(length($_)>0) {
	print "**** warning: CONFIG file $filename seems to be longer than content\n";
	print "              of FIELD file #$fi lets expect!\n";
      }
    }
  
  } else { # read w/o FIELD information
    $t=0; $m=0; $a=0;
    while(<$filehandle>) {
      $_ =~ s/^\s+//; $_ =~ s/\s+$//;
      last if(length($_)==0);
      ($cdata[$ci][$t][$m][$a][9],$cdata[$ci][$t][$m][$a][11]) = /^\s*(\S+)\s*(\S*)/;
      if(not $_=<$filehandle>) {print "**** unexpected end of CONFIG file (t=$t,m=$m,a=$a,i=2) on line $.!\n"; return 1}
      ($cdata[$ci][$t][$m][$a][0],$cdata[$ci][$t][$m][$a][1],$cdata[$ci][$t][$m][$a][2])
	= /^\s*(\S*)\s*(\S*)\s*(\S*)/;
      if($config_key[$ci]>0) { #read velocity
	if(not $_=<$filehandle>) {print "**** unexpected end of CONFIG file (t=$t,m=$m,a=$a,i=3) on line $.!\n"; return 1}
	($cdata[$ci][$t][$m][$a][3],$cdata[$ci][$t][$m][$a][4],$cdata[$ci][$t][$m][$a][5])
	= /^\s*(\S*)\s*(\S*)\s*(\S*)/;
	if($config_key[$ci]>1) { #read force
	  if(not $_=<$filehandle>) {print "**** unexpected end of CONFIG file (t=$t,m=$m,a=$a,i=4) on line $.!\n"; return 1}
	  ($cdata[$ci][$t][$m][$a][6],$cdata[$ci][$t][$m][$a][7],$cdata[$ci][$t][$m][$a][8])
	  = /^\s*(\S*)\s*(\S*)\s*(\S*)/;
	}
      }
      $minpos[$ci][0] = $cdata[$ci][$t][$m][$a][0] if($cdata[$ci][$t][$m][$a][0]<$minpos[$ci][0]);
      $minpos[$ci][1] = $cdata[$ci][$t][$m][$a][1] if($cdata[$ci][$t][$m][$a][1]<$minpos[$ci][1]);
      $minpos[$ci][2] = $cdata[$ci][$t][$m][$a][2] if($cdata[$ci][$t][$m][$a][2]<$minpos[$ci][2]);
      $maxpos[$ci][0] = $cdata[$ci][$t][$m][$a][0] if($cdata[$ci][$t][$m][$a][0]>$maxpos[$ci][0]);
      $maxpos[$ci][1] = $cdata[$ci][$t][$m][$a][1] if($cdata[$ci][$t][$m][$a][1]>$maxpos[$ci][1]);
      $maxpos[$ci][2] = $cdata[$ci][$t][$m][$a][2] if($cdata[$ci][$t][$m][$a][2]>$maxpos[$ci][2]);
    $a++;
    } # end while(<$filehandle>)
    $frame_numatoms[$ci]=$a;
  } # end if
  close($filehandle);
  return 0;
} # end subroutine read_config_file


sub copy_config {
  my $src       = $_[0];
  my $trg       = $_[1];
  return if($src==$trg);
  @{$cdata[$trg]}       = ();
  @{$cdata[$trg]}       = @{dclone \@{$cdata[$src]}};
  $config_key[$trg]     = $config_key[$src];
  $config_title[$trg]   = $config_title[$src];
  $periodic_key[$trg]   = $periodic_key[$src];
  $frame_numatoms[$trg] = $frame_numatoms[$src];
  if($periodic_key[$src] > 0) {
    @{$cell[$trg]}      = @{dclone \@{$cell[$src]}};
  }
}

sub copy_field {
  my $src       = $_[0];
  my $trg       = $_[1];
  return if($src==$trg);
  $field_numheader[$trg]    = $field_numheader[$src];
  $field_nummols[$trg]      = $field_nummols[$src];
  $field_numatoms[$trg]     = $field_numatoms[$src];
  $field_numatomtypes[$trg] = $field_numatomtypes[$src];
  $field_numvdw[$trg]       = $field_numvdw[$src];
  $field_numtbp[$trg]       = $field_numtbp[$src];
  $field_numfbp[$trg]       = $field_numfbp[$src];
  $field_nummetal[$trg]     = $field_nummetal[$src];
  $field_numextern[$trg]    = $field_numextern[$src];
  $field_neutral[$trg]      = $field_neutral[$src];
  $field_title[$trg]        = $field_title[$src];
  $field_units[$trg]           = $field_units[$src];
  @{$field_header[$trg]}       = @{$field_header[$src]};
  @{$mol_name[$trg]}           = @{$mol_name[$src]};
  @{$mol_mass[$trg]}           = @{$mol_mass[$src]};
  @{$mol_charge[$trg]}         = @{$mol_charge[$src]};
  @{$mol_numents[$trg]}        = @{$mol_numents[$src]};
  @{$mol_numatoms[$trg]}       = @{$mol_numatoms[$src]};
  @{$mol_numbonds[$trg]}       = @{$mol_numbonds[$src]};
  @{$mol_numconstraints[$trg]} = @{$mol_numconstraints[$src]};
  @{$mol_numangles[$trg]}      = @{$mol_numangles[$src]};
  @{$mol_numdihedrals[$trg]}   = @{$mol_numdihedrals[$src]};
  @{$mol_numinversions[$trg]}  = @{$mol_numinversions[$src]};
  @{$mol_numtether[$trg]}      = @{$mol_numtether[$src]};
  @{$mol_atomdata[$trg]}       = @{dclone \@{$mol_atomdata[$src]}};
  @{$mol_bonddata[$trg]}       = @{dclone \@{$mol_bonddata[$src]}};
  @{$mol_bondatoms[$trg]}      = @{dclone \@{$mol_bondatoms[$src]}};
  @{$mol_constraintdata[$trg]} = @{dclone \@{$mol_constraintdata[$src]}};
  @{$mol_angledata[$trg]}      = @{dclone \@{$mol_angledata[$src]}};
  @{$mol_dihedraldata[$trg]}   = @{dclone \@{$mol_dihedraldata[$src]}};
  @{$mol_inversiondata[$trg]}  = @{dclone \@{$mol_inversiondata[$src]}};
  @{$mol_tetherdata[$trg]}     = @{dclone \@{$mol_tetherdata[$src]}};
  @{$field_vdwdata[$trg]}      = @{dclone \@{$field_vdwdata[$src]}};
  @{$field_tbpdata[$trg]}      = @{dclone \@{$field_tbpdata[$src]}};
  @{$field_fbpdata[$trg]}      = @{dclone \@{$field_fbpdata[$src]}};
  @{$field_metaldata[$trg]}    = @{dclone \@{$field_metaldata[$src]}};
  @{$field_externdata[$trg]}   = @{dclone \@{$field_externdata[$src]}};
  @{$field_atomtypes[$trg]}    = @{$field_atomtypes[$src]};
}

sub find_history_timestep {
  # searches already opened HISTORY file for a certain frame, returns 1 if the frame
  # was found and leaves the cursor in the file at the position right before the frame
  # so that it can be read with subroutine read_history_timestep
  my $filehandle  = $_[0];
  my $targetframe = $_[1]; # <0 for last frame
  my $found       = 0;
  my $line_number = 0;
  my($frame_number,$i,$line_number_last_frame);
  while(<$filehandle>) {
    if(/timestep/) {
      $line_number_last_frame = $line_number;
      ($frame_number) = /timestep\s+(\S+)/;
      if($frame_number == $targetframe) {
	$found = 1;
	seek($filehandle,-length($_), 1);
	last;
      } elsif($frame_number > $targetframe) {
	print "**** error: target frame was not found in HISTORY file\n";
	return $found;
      }
    }
    $line_number++;
  }
  if($targetframe<0 and defined($line_number_last_frame)) {
    $found = 1;
    seek($filehandle,0,0); # rewind file
    for($i=0;$i<$line_number_last_frame;$i++) {$_=<$filehandle>;}
  }
  return $found;
} # end subroutine find_timestep


sub skip_history_timestep {
  my $filehandle  = $_[0];
  my($imax,$i,$frame_number,$frame_numatoms,$config_key,$periodic_key,$frame_timestep);
  if(not $_=<$filehandle>) {return -1} #end of file
  if(not /timestep/) {
    $_=<$filehandle>; $_=<$filehandle>;
    if( not /timestep/) {
      print "error reading HISTORY file on line $.: \"timestep\" expected but found:\n$_";
      return 1;
    }
  }
  ($frame_number, $frame_numatoms, $config_key, $periodic_key, $frame_timestep) =
  /timestep\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
  $imax = ($periodic_key>0)*3+(2+$config_key)*$frame_numatoms;
  for($i=0;$i<$imax;$i++) {
    if(not $_=<$filehandle>) {return -1} #end of file
  }
  return 0;
} # end subroutine skip_history_timestep


sub read_history_timestep {
  # reads the next timestep from an already opened HISTORY file and writes output to cdata array
  
  my $filehandle = $_[0];
  my $fi         = $_[1]; #index of corresponding FIELD file
  my $ci         = $_[2];
  my($i,$j,$t,$m,$a,$molindex);
  
  return 2 if($ci<0);
  
  if(not $_=<$filehandle>) {return -1} #end of file
  if(not /timestep/) {
    $_=<$filehandle>; $_=<$filehandle>;
    if( not /timestep/) {
      print "error reading HISTORY file on line $.: \"timestep\" expected but found:\n$_";
      return 2;
    }
  }
  
  #reset values
  @{$cdata[$ci]} = ();
  @{$cell[$ci]}  = ();
  @{$size[$ci]}  = ();
  @{$minpos[$ci]}=( 9e20, 9e20, 9e20);
  @{$maxpos[$ci]}=(-9e20,-9e20,-9e20);
  $molindex      = 0;
  
  ($frame_number[$ci], $frame_numatoms[$ci], $config_key[$ci], $periodic_key[$ci], $frame_timestep[$ci]) =
  /timestep\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
  if($periodic_key[$ci]>0) {
    for($j=0;$j<3;$j++) {
      $_=<$filehandle>; $_ =~ s/^\s+//; $_ =~ s/\s+$//;
      ($cell[$ci][$j][0], $cell[$ci][$j][1], $cell[$ci][$j][2]) = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
    }
    $size[$ci][0] = ($cell[$ci][0][0] + $cell[$ci][1][0] + $cell[$ci][2][0])/2;
    $size[$ci][1] = ($cell[$ci][0][1] + $cell[$ci][1][1] + $cell[$ci][2][1])/2;
    $size[$ci][2] = ($cell[$ci][0][2] + $cell[$ci][1][2] + $cell[$ci][2][2])/2;
  }
  if($fi>=0) {  # read using FIELD information
    for($t=0;$t<$field_nummols[$fi];$t++) {
      @{$cdata[$ci][$t]}=();
      for($m=0;$m<$mol_numents[$fi][$t];$m++) {
	@{$cdata[$ci][$t][$m]}=();
	for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
	  @{$cdata[$ci][$t][$m][$a]}=();
	  if(not $_=<$filehandle>) {
	  print "**** unexpected end of HISTORY file (t=$t,m=$m,a=$a,i=1)!\n"; return 1}
	  $_ =~ s/^\s+//; $_ =~ s/\s+$//;
	  ($cdata[$ci][$t][$m][$a][9],$cdata[$ci][$t][$m][$a][11], $cdata[$ci][$t][$m][$a][12], $cdata[$ci][$t][$m][$a][13])
	  = split(/\s+/,$_);
	  if($cdata[$ci][$t][$m][$a][9] ne $mol_atomdata[$fi][$t][$a][0]) {
	    print "**** error: atom name $cdata[$ci][$t][$m][$a][9] in HISTORY file on line $. does not match name ",
	    "$mol_atomdata[$fi][$t][$a][0] in FIELD data for molecule $mol_name[$fi][$t] from file ",
	    "'$field_filename[$fi]'\n";
	    exit 1;
	  }
	  if(not $_=<$filehandle>) {
	  print "**** unexpected end of HISTORY file on line $. (t=$t,m=$m,a=$a,i=2)!\n"; return 1}
	  $_ =~ s/^\s+//; $_ =~ s/\s+$//;
	  @{$cdata[$ci][$t][$m][$a]}[0..2] = split(/\s+/,$_);
	  if($config_key[$ci]>0) { #read velocity
	    if(not $_=<$filehandle>) {
	    print "**** unexpected end of HISTORY file on line $. (t=$t,m=$m,a=$a,i=3)!\n"; return 1}
	    $_ =~ s/^\s+//; $_ =~ s/\s+$//;
	    @{$cdata[$ci][$t][$m][$a]}[3..5] = split(/\s+/,$_);
	    if($config_key[$ci]>1) { #read force
	      if(not $_=<$filehandle>) {
	      print "**** unexpected end of HISTORY file on line $. (t=$t,m=$m,a=$a,i=4)!\n"; return 1}
	      $_ =~ s/^\s+//; $_ =~ s/\s+$//;
	      @{$cdata[$ci][$t][$m][$a]}[6..8] = split(/\s+/,$_);
	    }
	  }
	  $cdata[$ci][$t][$m][$a][10]=$molindex;
	  $minpos[$ci][0] = $cdata[$ci][$t][$m][$a][0] if($cdata[$ci][$t][$m][$a][0]<$minpos[$ci][0]);
	  $minpos[$ci][1] = $cdata[$ci][$t][$m][$a][1] if($cdata[$ci][$t][$m][$a][1]<$minpos[$ci][1]);
	  $minpos[$ci][2] = $cdata[$ci][$t][$m][$a][2] if($cdata[$ci][$t][$m][$a][2]<$minpos[$ci][2]);
	  $maxpos[$ci][0] = $cdata[$ci][$t][$m][$a][0] if($cdata[$ci][$t][$m][$a][0]>$maxpos[$ci][0]);
	  $maxpos[$ci][1] = $cdata[$ci][$t][$m][$a][1] if($cdata[$ci][$t][$m][$a][1]>$maxpos[$ci][1]);
	  $maxpos[$ci][2] = $cdata[$ci][$t][$m][$a][2] if($cdata[$ci][$t][$m][$a][2]>$maxpos[$ci][2]);
	}
	$molindex++;
      }
    }

  } else { # read w/o FIELD information
    $t=0; $m=0;
    for($a=0;$a<$frame_numatoms[$ci];$a++) {
      if(not $_=<$filehandle>) {
      print "**** unexpected end of HISTORY file (t=$t,m=$m,a=$a,i=1)!\n"; return 1}
      $_ =~ s/^\s+//; $_ =~ s/\s+$//;
      ($cdata[$ci][$t][$m][$a][9],$cdata[$ci][$t][$m][$a][11], $cdata[$ci][$t][$m][$a][12], $cdata[$ci][$t][$m][$a][13])
      = split(/\s+/,$_);
      if(not $_=<$filehandle>) {print "**** unexpected end of HISTORY file on line $. (t=$t,m=$m,a=$a,i=2)!\n"; return 1}
      $_ =~ s/^\s+//; $_ =~ s/\s+$//;
      @{$cdata[$ci][$t][$m][$a]}[0..2] = split(/\s+/,$_);
      if($config_key[$ci]>0) { #read velocity
	if(not $_=<$filehandle>) {print "**** unexpected end of HISTORY file on line $. (t=$t,m=$m,a=$a,i=3)!\n"; return 1}
	$_ =~ s/^\s+//; $_ =~ s/\s+$//;
	@{$cdata[$ci][$t][$m][$a]}[3..5] = split(/\s+/,$_);
	if($config_key[$ci]>1) { #read force
	  if(not $_=<$filehandle>) {print "**** unexpected end of HISTORY file on line $. (t=$t,m=$m,a=$a,i=4)!\n"; return 1}
	  $_ =~ s/^\s+//; $_ =~ s/\s+$//;
	  @{$cdata[$ci][$t][$m][$a]}[6..8] = split(/\s+/,$_);
	}
      }
      $minpos[$ci][0] = $cdata[$ci][$t][$m][$a][0] if($cdata[$ci][$t][$m][$a][0]<$minpos[$ci][0]);
      $minpos[$ci][1] = $cdata[$ci][$t][$m][$a][1] if($cdata[$ci][$t][$m][$a][1]<$minpos[$ci][1]);
      $minpos[$ci][2] = $cdata[$ci][$t][$m][$a][2] if($cdata[$ci][$t][$m][$a][2]<$minpos[$ci][2]);
      $maxpos[$ci][0] = $cdata[$ci][$t][$m][$a][0] if($cdata[$ci][$t][$m][$a][0]>$maxpos[$ci][0]);
      $maxpos[$ci][1] = $cdata[$ci][$t][$m][$a][1] if($cdata[$ci][$t][$m][$a][1]>$maxpos[$ci][1]);
      $maxpos[$ci][2] = $cdata[$ci][$t][$m][$a][2] if($cdata[$ci][$t][$m][$a][2]>$maxpos[$ci][2]);
    }
  }
  return 0;
} # end subroutine read_history_timestep


sub move_mol {
  # call example:
  # move_mol(\@{$cdata[0][$t][$m]},\@shiftvec);
  my @atoms = @{$_[0]};
  my @shiftvec = @{$_[1]};
  for($a=0;$a<@atoms;$a++) {
    $atoms[$a][0] += $shiftvec[0];
    $atoms[$a][1] += $shiftvec[1];
    $atoms[$a][2] += $shiftvec[2];
  }
}

sub move_mol2 {
  # call example:
  # move_mol2($ci,$t,$m,\@shiftvec);
  my $ci = $_[0];
  my $t  = $_[1];
  my $m  = $_[2];
  my @shiftvec = @{$_[3]};
  for($a=0;$a<@{$cdata[$ci][$t][$m]};$a++) {
    $cdata[$ci][$t][$m][$a][0] += $shiftvec[0];
    $cdata[$ci][$t][$m][$a][1] += $shiftvec[1];
    $cdata[$ci][$t][$m][$a][2] += $shiftvec[2];
  }
}


sub add_mol_velocity {
  # call example:
  # add_mol_velocity(\@{$cdata[0][$t][$m]},\@shiftvec);
  my @atoms = @{$_[0]};
  my @velvec = @{$_[1]};
  for($a=0;$a<@atoms;$a++) {
    $atoms[$a][3] += $velvec[0];
    $atoms[$a][4] += $velvec[1];
    $atoms[$a][5] += $velvec[2];
  }
}

sub check_orthogonal {
  my $ci = $_[0];
  return 0 if(abs($cell[$ci][0][1])>1e-14);
  return 0 if(abs($cell[$ci][0][2])>1e-14);
  return 0 if(abs($cell[$ci][1][0])>1e-14);
  return 0 if(abs($cell[$ci][1][2])>1e-14);
  return 0 if(abs($cell[$ci][2][0])>1e-14);
  return 0 if(abs($cell[$ci][2][1])>1e-14);
  return 1;
}

sub remap_molecule {
  #!!! only works for orthogonal cells. For non-orthogonal cells use remap_molecule2
  my($a,$c,$a_old,@d);
  # call example:
  # remap_molecule(\@{$cdata[0][$t][$m]},\@remapaxes,\@{$size[0]});
  my @atoms = @{$_[0]};
  # format of atoms array must be:
  # index 1: atom index
  # index 2: 0  1  2   3   4   5   6   7   8     9
  #          x  y  z  vx  vy  vz  fx  fy  fz  name
  # (cf. structure of @cdata in subroutine read_config_file)
  
  my @axis  = @{$_[1]};
  # array @axis may contain values
  # 0 = x-axis
  # 1 = y-axis
  # 2 = z-axis
  
  my @subsize = @{$_[2]};
  
  my $startindex = $_[3]; # optional
  
  if(defined($startindex)) {
    for($a=$startindex+1;$a<@atoms;$a++) {
      $a_old = $a-1;
      foreach $c (@axis) {
	$d[$c] = $atoms[$a][$c]-$atoms[$a_old][$c];
	if(abs($d[$c])>$subsize[$c]) {
	  $atoms[$a][$c] -= 2 * ($d[$c]<=>0) * $subsize[$c];
	}
      } # end for $c
    } # end for $a
    for($a=$startindex-1;$a>-1;$a--) {
      $a_old = $a+1;
      foreach $c (@axis) {
	$d[$c] = $atoms[$a][$c]-$atoms[$a_old][$c];
	if(abs($d[$c])>$subsize[$c]) {
	  $atoms[$a][$c] -= 2 * ($d[$c]<=>0) * $subsize[$c];
	}
      } # end for $c
    } # end for $a
  } else {
    for($a=1;$a<@atoms;$a++) {
      $a_old = $a-1;
      foreach $c (@axis) {
	$d[$c] = $atoms[$a][$c]-$atoms[$a_old][$c];
	if(abs($d[$c])>$subsize[$c]) {
	  $atoms[$a][$c] -= 2 * ($d[$c]<=>0) * $subsize[$c];
	}
      } # end for $c
    } # end for $a
  }
} # end subroutine remap_molecule


sub remap_molecule2 {
  my($a,$c,$a_old,@d,@d2,$dsq,$dminsq,$x,$y,$z,@min);
  # call example:
  # remap_molecule2($fi,$ci,$t,$m,[$aref]);
  my $fi = $_[0];
  my $ci = $_[1];
  my $t  = $_[2];
  my $m  = $_[3];
  my $startindex = $_[4]; # optional
  if($periodic_key[$ci]==0) {
    # do nothing if system is not periodic
  } elsif($periodic_key[$ci]==1 or $periodic_key[$ci]==2 or $periodic_key[$ci]==3) {
    for($a=1;$a<$mol_numatoms[$fi][$t];$a++) {
      $a_old = $a-1;
      $d[0] = $cdata[$ci][$t][$m][$a][0] - $cdata[$ci][$t][$m][$a_old][0];
      $d[1] = $cdata[$ci][$t][$m][$a][1] - $cdata[$ci][$t][$m][$a_old][1];
      $d[2] = $cdata[$ci][$t][$m][$a][2] - $cdata[$ci][$t][$m][$a_old][2];
      $dminsq = $d[0]*$d[0]+$d[1]*$d[1]+$d[2]*$d[2];
      @min=(0,0,0);
      for($x=-1;$x<2;$x++) {
	for($y=-1;$y<2;$y++) {
	  for($z=-1;$z<2;$z++) {
	    next if($x==0 and $y==0 and $z==0);
	    $d2[0] = $d[0] + $x*$cell[$ci][0][0] + $y*$cell[$ci][1][0] + $z*$cell[$ci][2][0];
	    $d2[1] = $d[1] + $x*$cell[$ci][0][1] + $y*$cell[$ci][1][1] + $z*$cell[$ci][2][1];
	    $d2[2] = $d[2] + $x*$cell[$ci][0][2] + $y*$cell[$ci][1][2] + $z*$cell[$ci][2][2];
	    $dsq = $d2[0]*$d2[0]+$d2[1]*$d2[1]+$d2[2]*$d2[2];
	    if($dsq<$dminsq) {
	      $dminsq = $dsq;
	      @min=($x,$y,$z);
	    }
	  }
	}
      }
      $cdata[$ci][$t][$m][$a][0] += $min[0]*$cell[$ci][0][0] + $min[1]*$cell[$ci][1][0] + $min[2]*$cell[$ci][2][0];
      $cdata[$ci][$t][$m][$a][1] += $min[0]*$cell[$ci][0][1] + $min[1]*$cell[$ci][1][1] + $min[2]*$cell[$ci][2][1];
      $cdata[$ci][$t][$m][$a][2] += $min[0]*$cell[$ci][0][2] + $min[1]*$cell[$ci][1][2] + $min[2]*$cell[$ci][2][2];
    }
  } elsif($periodic_key[$ci]==6) {
    for($a=1;$a<$mol_numatoms[$fi][$t];$a++) {
      $a_old = $a-1;
      $d[0] = $cdata[$ci][$t][$m][$a][0] - $cdata[$ci][$t][$m][$a_old][0];
      $d[1] = $cdata[$ci][$t][$m][$a][1] - $cdata[$ci][$t][$m][$a_old][1];
      $d[2] = $cdata[$ci][$t][$m][$a][2] - $cdata[$ci][$t][$m][$a_old][2];
      $dminsq = $d[0]*$d[0]+$d[1]*$d[1]+$d[2]*$d[2];
      @min=(0,0);
      $d2[2] = $d[2];
      for($x=-1;$x<2;$x++) {
	for($y=-1;$y<2;$y++) {
	  next if($x==0 and $y==0);
	  $d2[0] = $d[0] + $x*$cell[$ci][0][0] + $y*$cell[$ci][1][0];
	  $d2[1] = $d[1] + $x*$cell[$ci][0][1] + $y*$cell[$ci][1][1];
	  $dsq = $d2[0]*$d2[0]+$d2[1]*$d2[1]+$d2[2]*$d2[2];
	  if($dsq<$dminsq) {
	    $dminsq = $dsq;
	    @min=($x,$y);
	  }
	}
      }
      $cdata[$ci][$t][$m][$a][0] += $min[0]*$cell[$ci][0][0] + $min[1]*$cell[$ci][1][0];
      $cdata[$ci][$t][$m][$a][1] += $min[0]*$cell[$ci][0][1] + $min[1]*$cell[$ci][1][1];
    }
  } else {
    print "**** error: remapping for key $periodic_key[$ci] is not implemented!\n";
    return 1;
  }
  return 0;
} # end subroutine remap_molecule2


sub rotate_molecule {
  # call example:
  # rotate_molecule(\@{$cdata[0][$t][$m]}, \@rotmatrix, \@center);
  my @atoms  = @{$_[0]};
  # format of atoms array must be:
  # index 1: atom index
  # index 2: 0  1  2   3   4   5   6   7   8     9
  #          x  y  z  vx  vy  vz  fx  fy  fz  name
  # (cf. structure of @cdata in subroutine read_config_file)
  my @rotmatrix = @{$_[1]};
  my @center    = @{$_[2]} if($#_>1);
  my @new;
  
  for($a=0;$a<@atoms;$a++) {
    if($#_>1) {
      $atoms[$a][0] -= $center[0];
      $atoms[$a][1] -= $center[1];
      $atoms[$a][2] -= $center[2];
    }
    $new[0]  = $atoms[$a][0]*$rotmatrix[0][0]+$atoms[$a][1]*$rotmatrix[0][1]+$atoms[$a][2]*$rotmatrix[0][2];
    $new[1]  = $atoms[$a][0]*$rotmatrix[1][0]+$atoms[$a][1]*$rotmatrix[1][1]+$atoms[$a][2]*$rotmatrix[1][2];
    $new[2]  = $atoms[$a][0]*$rotmatrix[2][0]+$atoms[$a][1]*$rotmatrix[2][1]+$atoms[$a][2]*$rotmatrix[2][2];
    $atoms[$a][0] = $new[0];
    $atoms[$a][1] = $new[1];
    $atoms[$a][2] = $new[2];
    $new[0]  = $atoms[$a][3]*$rotmatrix[0][0]+$atoms[$a][4]*$rotmatrix[0][1]+$atoms[$a][5]*$rotmatrix[0][2];
    $new[1]  = $atoms[$a][3]*$rotmatrix[1][0]+$atoms[$a][4]*$rotmatrix[1][1]+$atoms[$a][5]*$rotmatrix[1][2];
    $new[2]  = $atoms[$a][3]*$rotmatrix[2][0]+$atoms[$a][4]*$rotmatrix[2][1]+$atoms[$a][5]*$rotmatrix[2][2];
    $atoms[$a][3] = $new[0];
    $atoms[$a][4] = $new[1];
    $atoms[$a][5] = $new[2];
    $new[0]  = $atoms[$a][6]*$rotmatrix[0][0]+$atoms[$a][7]*$rotmatrix[0][1]+$atoms[$a][2]*$rotmatrix[0][8];
    $new[1]  = $atoms[$a][6]*$rotmatrix[1][0]+$atoms[$a][7]*$rotmatrix[1][1]+$atoms[$a][2]*$rotmatrix[1][8];
    $new[2]  = $atoms[$a][6]*$rotmatrix[2][0]+$atoms[$a][7]*$rotmatrix[2][1]+$atoms[$a][2]*$rotmatrix[2][8];
    $atoms[$a][6] = $new[0];
    $atoms[$a][7] = $new[1];
    $atoms[$a][8] = $new[2];
    if($#_>1) {
      $atoms[$a][0] += $center[0];
      $atoms[$a][1] += $center[1];
      $atoms[$a][2] += $center[2];
    }
  }
} # end subroutine rotate_molecule


sub write_config_file {
  my($filehandle,$t,$m,$a,$i,$j);
  my $filename = $_[0];
  my $ci       = $_[1];
  my $title    = $_[2];
  my $mode     = $_[3];
  my $fi       = $_[4];
  $mode = ">" if(not defined($mode));
  $fi   = -1  if(not defined($fi));
  if(not open($filehandle, $mode, $filename)) {
    print "**** error: Can't open output CONFIG file \"$filename\": $!\n";
    return 1;
  }
  $title =~ s/\s+$//;
  print $filehandle "$title\n";
  printf $filehandle "%10u", $config_key[$ci]; # $filehandle file key
  printf $filehandle "%10u\n", $periodic_key[$ci]; # periodic boundary key
  if($periodic_key[$ci] != 0) {
    printf $filehandle " %19.12e %19.12e %19.12e\n", $cell[$ci][0][0], $cell[$ci][0][1], $cell[$ci][0][2];
    printf $filehandle " %19.12e %19.12e %19.12e\n", $cell[$ci][1][0], $cell[$ci][1][1], $cell[$ci][1][2];
    printf $filehandle " %19.12e %19.12e %19.12e\n", $cell[$ci][2][0], $cell[$ci][2][1], $cell[$ci][2][2];
  }
  $i=0;
  if($fi<0) {
    for($t=0;$t<@{$cdata[$ci]};$t++) {
      for($m=0;$m<@{$cdata[$ci][$t]};$m++) {
	for($a=0;$a<@{$cdata[$ci][$t][$m]};$a++) {
	  printf $filehandle "%-8s",       $cdata[$ci][$t][$m][$a][9];
	  printf $filehandle "%10u\n",     $i+1; #atom index
	  printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][0];
	  printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][1];
	  printf $filehandle " %19.12e\n", $cdata[$ci][$t][$m][$a][2];
	  if($config_key[$ci]>0) {
	    printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][3];
	    printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][4];
	    printf $filehandle " %19.12e\n", $cdata[$ci][$t][$m][$a][5];
	    if($config_key[$ci]>1) {
	      printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][6];
	      printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][7];
	      printf $filehandle " %19.12e\n", $cdata[$ci][$t][$m][$a][8];
	    }
	  }
	  $i++;
	} # end for $a
      } # end for $m
    } # end for $t
  } else {
    for($t=0;$t<$field_nummols[$fi];$t++) {
      for($m=0;$m<$mol_numents[$fi][$t];$m++) {
	for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
	  printf $filehandle "%-8s",       $mol_atomdata[$fi][$t][$a][0];
	  printf $filehandle "%10u\n",     $i+1; #atom index
	  printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][0];
	  printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][1];
	  printf $filehandle " %19.12e\n", $cdata[$ci][$t][$m][$a][2];
	  if($config_key[$ci]>0) {
	    printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][3];
	    printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][4];
	    printf $filehandle " %19.12e\n", $cdata[$ci][$t][$m][$a][5];
	    if($config_key[$ci]>1) {
	      printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][6];
	      printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][7];
	      printf $filehandle " %19.12e\n", $cdata[$ci][$t][$m][$a][8];
	    }
	  }
	  $i++;
	} # end for $a
      } # end for $m
    } # end for $t
  } # endif $fi
  close($filehandle);
  return 0;
} # end subroutine write_config_file


sub write_history_timestep {
  my($t,$m,$a,$i,$j);
  my $filehandle = $_[0];
  my $fi         = $_[1]; #index of corresponding FIELD file
  my $ci         = $_[2];
  return 1 if($ci<0);
  printf $filehandle "timestep %9u",$frame_number[$ci]; #current timestep
  printf $filehandle "%10u",$frame_numatoms[$ci]; #number of atoms
  printf $filehandle "%10u", $config_key[$ci];
  printf $filehandle "%10u", $periodic_key[$ci];
  printf $filehandle "%12.6f\n", $frame_timestep[$ci];
  if($periodic_key[$ci] != 0) {
    printf $filehandle " %19.12e %19.12e %19.12e\n", $cell[$ci][0][0], $cell[$ci][0][1], $cell[$ci][0][2];
    printf $filehandle " %19.12e %19.12e %19.12e\n", $cell[$ci][1][0], $cell[$ci][1][1], $cell[$ci][1][2];
    printf $filehandle " %19.12e %19.12e %19.12e\n", $cell[$ci][2][0], $cell[$ci][2][1], $cell[$ci][2][2];
  }
  $i=1;
  for($t=0;$t<@{$cdata[$ci]};$t++) {
    for($m=0;$m<@{$cdata[$ci][$t]};$m++) {
      for($a=0;$a<@{$cdata[$ci][$t][$m]};$a++) {
	printf $filehandle "%-8s",       $cdata[$ci][$t][$m][$a][9];
	printf $filehandle "%10u",       $i; #atom index
	if($fi>=0) {
	  printf $filehandle " %12.6f",   $mol_atomdata[$t][$m][$a][1]; #atomic mass
	  printf $filehandle " %12.6f\n", $mol_atomdata[$t][$m][$a][2]; #atomic charge
	} elsif(defined($cdata[$ci][$t][$m][$a][12])) {
	  printf $filehandle " %12.6f",   $cdata[$ci][$t][$m][$a][12];
	  printf $filehandle " %12.6f\n", $cdata[$ci][$t][$m][$a][13];
	} else {
	  printf $filehandle " %12.6f",   0;
	  printf $filehandle " %12.6f\n", 0;
	}
	printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][0];
	printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][1];
	printf $filehandle " %19.12e\n", $cdata[$ci][$t][$m][$a][2];
	if($config_key[$ci]>0) {
	  printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][3];
	  printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][4];
	  printf $filehandle " %19.12e\n", $cdata[$ci][$t][$m][$a][5];
	  if($config_key[$ci]>1) {
	    printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][6];
	    printf $filehandle " %19.12e",   $cdata[$ci][$t][$m][$a][7];
	    printf $filehandle " %19.12e\n", $cdata[$ci][$t][$m][$a][8];
	  }
	}
	$i++;
      } # end for $a
    } # end for $m
  } # end for $t
  return 0;
} # end subroutine write_history_timestep

sub remove_mol_type {
  my $ci = $_[0];
  my $fi = $_[1];
  my $tr = $_[2];
  my($i,$a,$t,@removetypes);
  return 1 if($fi<0);
  # try to get index from molecule name if necessary
  if(not check_integer($tr)) {
    for($t=0;$t<$field_nummols[$fi];$t++) {
      if($mol_name[$fi][$t] eq $tr) {
	$tr=$t;
      }
    }
    if(not check_integer($tr)) {
      print "**** error in sub remove_mol_type: a molecule named $tr could not be found!\n";
      return 1;
    }
  } elsif(not exists($mol_name[$fi][$tr])) {
    print "**** error in sub remove_mol_type: molecule number $tr does not exist!\n";
    return 1;
  }
  
  # remove atom types if possible
  for($i=0;$i<$field_numatomtypes[$fi];$i++) {
    $removetypes[$i]=$i;
  }
  for($t=0;$t<$field_nummols[$fi];$t++) {
    next if($t==$tr);
    for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
      for($i=0;$i<@removetypes;$i++) {
	if($mol_atomdata[$fi][$t][$a][0] eq $field_atomtypes[$fi][$removetypes[$i]]) {
	  splice(@removetypes,$i,1);
	  $i--;
	}
      }
    }
  }
  for($i=$#removetypes;$i>=0;$i--) {
    # go from highest number to lowest to avoid troubles with indices upon splice
    # remove obsolete vdw data
    for($a=0;$a<$field_numvdw[$fi];$a++) {
      if($field_vdwdata[$fi][$a][0] eq $field_atomtypes[$fi][$removetypes[$i]]
      or $field_vdwdata[$fi][$a][1] eq $field_atomtypes[$fi][$removetypes[$i]]) {
	splice(@{$field_vdwdata[$fi]},$a,1);
	$field_numvdw[$fi]--;
	$a--;
      }
    }
    # remove obsolete tbp data
    for($a=0;$a<$field_numtbp[$fi];$a++) {
      if($field_tbpdata[$fi][$a][0] eq $field_atomtypes[$fi][$removetypes[$i]]
      or $field_tbpdata[$fi][$a][1] eq $field_atomtypes[$fi][$removetypes[$i]]
      or $field_tbpdata[$fi][$a][2] eq $field_atomtypes[$fi][$removetypes[$i]]) {
	splice(@{$field_tbpdata[$fi]},$a,1);
	$field_numtbp[$fi]--;
	$a--;
      }
    }
    # remove obsolete fbp data
    for($a=0;$a<$field_numfbp[$fi];$a++) {
      if($field_fbpdata[$fi][$a][0] eq $field_atomtypes[$fi][$removetypes[$i]]
      or $field_fbpdata[$fi][$a][1] eq $field_atomtypes[$fi][$removetypes[$i]]
      or $field_fbpdata[$fi][$a][2] eq $field_atomtypes[$fi][$removetypes[$i]]
      or $field_fbpdata[$fi][$a][3] eq $field_atomtypes[$fi][$removetypes[$i]]) {
	splice(@{$field_fbpdata[$fi]},$a,1);
	$field_numfbp[$fi]--;
	$a--;
      }
    }
    # remove obsolete metal data
    for($a=0;$a<$field_nummetal[$fi];$a++) {
      if($field_metaldata[$fi][$a][0] eq $field_atomtypes[$fi][$removetypes[$i]]
      or $field_metaldata[$fi][$a][1] eq $field_atomtypes[$fi][$removetypes[$i]]) {
	splice(@{$field_metaldata[$fi]},$a,1);
	$field_nummetal[$fi]--;
	$a--;
      }
    }
    splice(@{$field_atomtypes[$fi]},$removetypes[$i],1);
    $field_numatomtypes[$fi]--;
  }
  $field_numatoms[$fi] -= $mol_numents[$fi][$tr]*$mol_numatoms[$fi][$tr];
  splice(@{$mol_numatoms[$fi]},$tr,1);
  splice(@{$mol_numbonds[$fi]},$tr,1);
  splice(@{$mol_numconstraints[$fi]},$tr,1);
  splice(@{$mol_numangles[$fi]},$tr,1);
  splice(@{$mol_numdihedrals[$fi]},$tr,1);
  splice(@{$mol_numinversions[$fi]},$tr,1);
  splice(@{$mol_numtether[$fi]},$tr,1);
  splice(@{$mol_numents[$fi]},$tr,1);
  splice(@{$mol_atomdata[$fi]},$tr,1);
  splice(@{$mol_bonddata[$fi]},$tr,1);
  splice(@{$mol_bondatoms[$fi]},$tr,1);
  splice(@{$mol_constraintdata[$fi]},$tr,1);
  splice(@{$mol_angledata[$fi]},$tr,1);
  splice(@{$mol_dihedraldata[$fi]},$tr,1);
  splice(@{$mol_inversiondata[$fi]},$tr,1);
  splice(@{$mol_tetherdata[$fi]},$tr,1);
  splice(@{$mol_mass[$fi]},$tr,1);
  splice(@{$mol_name[$fi]},$tr,1);
  splice(@{$mol_charge[$fi]},$tr,1);
  $field_nummols[$fi]--;
  return 0 if(not defined($ci) or $ci<0);
  $frame_numatoms[$ci] -= $mol_numents[$fi][$tr]*$mol_numatoms[$fi][$tr];
  splice(@{$cdata[$ci]},$tr,1);
  return 0
} # end subroutine remove_mol_type

sub remove_mol_entity {
  my $ci = $_[0];
  my $fi = $_[1];
  my $t  = $_[2];
  my $m  = $_[3];
  return -1 if(not @{$cdata[$ci][$t][$m]});
  if($fi<0) { # no FIELD information
    $frame_numatoms[$ci] -= @{$cdata[$ci][$t][$m]};
    splice(@{$cdata[$ci][$t]},$m,1);
  } else {
    splice(@{$cdata[$ci][$t]},$m,1);
    $mol_numents[$fi][$t]--;
    $frame_numatoms[$ci] -= $mol_numatoms[$fi][$t];
    $field_numatoms[$fi] -= $mol_numatoms[$fi][$t];
  }
  return 0
} # end subroutine remove_mol_entity


sub read_statis_timestep {
  my $filehandle = $_[0];
  my $si         = $_[1];
  my($i,$numlines,@linedata);
  if(not $_=<$filehandle>) {return -1} #end of file
  $_ =~ s/^\s+//; $_ =~ s/\s+$//;
  @linedata = split(/\s+/, $_);
  if($#linedata != 2 or not check_real(@linedata)) {
    $_=<$filehandle>; $_=<$filehandle>;
    $_ =~ s/^\s+//; $_ =~ s/\s+$//;
    @linedata = split(/\s+/, $_);
    if($#linedata != 2 or not check_real(@linedata)) {
      print "**** error: expected new frame in STATIS file but found only junk on line $.!\n";
      return 1;
    }
  }
  undef @{$sdata[$si]};
  @{$sdata[$si]} = ();
  $numlines=ceil($linedata[2]/5.0);
  push(@{$sdata[$si]},@linedata[0..1]);
  for($i=0;$i<$numlines;$i++) {
    if(not $_=<$filehandle>) { print "**** error: unexpected end of STATIS file on line $.!"; return 1 } #end of file
    $_ =~ s/^\s+//; $_ =~ s/\s+$//;
    push(@{$sdata[$si]},split(/\s+/, $_));
  }
  return 0;
} # end subroutine read_statis_timestep


sub read_rdfdat_file {
  my $filename = $_[0];
  my $ri       = $_[1];
  my($filehandle,$i,$j,@linedata,$nrdfs,$mxrdf,$a1,$a2);
  if($ri<0) {
    print "**** error: index \$ri must not be smaller than zero in sub read_rdf_file\n";
    return -1;
  }
  if(not open($filehandle, "<", $filename)) {
    print "**** error: Can't open RDFDAT file \"$filename\": $!\n";
    return 1;
  }
  undef @{$rdfdata[$ri]};
  undef @{$rdfnames[$ri]};
  $rdfnames[$ri][0]="radius_[A]";
  if(not $_=<$filehandle>) {close($filehandle); return -1} #config title
  if(not $_=<$filehandle>) {close($filehandle); return -1} # n(RDFs) mxrdf
  $_ =~ s/^\s+//; $_ =~ s/\s+$//;
  ($nrdfs,$mxrdf) = split(/\s+/, $_);
  for($i=1;$i<=$nrdfs;$i++) {
    if(not $_=<$filehandle>) {close($filehandle); return -1} # atom types
    $_ =~ s/^\s+//; $_ =~ s/\s+$//;
    @{$rdfnames[$ri][$i]} = split(/\s+/, $_);
    for($j=0;$j<$mxrdf;$j++) {
      if(not $_=<$filehandle>) {close($filehandle); return -1} # atom types
      $_ =~ s/^\s+//; $_ =~ s/\s+$//;
      @linedata = split(/\s+/, $_);
      $rdfdata[$ri][$j][$i] = $linedata[1];
      $rdfdata[$ri][$j][0]  = $linedata[0] if($i==1);
    }
  }
  close($filehandle);
  return 0;
} # end subroutine read_rdfdat_file


sub print_rdfdat_data {
  my $filename = $_[0];
  my $ri       = $_[1];
  my($filehandle,$i,$j);
  if($ri<0) {
    print "**** error: index \$ri must not be smaller than zero in sub read_rdfdat_file\n";
    return -1;
  }
  if(not open($filehandle, ">", $filename)) {
    print "**** error: Can't open output file \"$filename\": $!\n";
    return 1;
  }
  printf $filehandle "#%16s",$rdfnames[$ri][0];
  for($i=1;$i<@{$rdfnames[$ri]};$i++) {
    printf $filehandle " %17s",join("-",@{$rdfnames[$ri][$i]});
  }
  print $filehandle "\n";
  for($i=0;$i<@{$rdfdata[$ri]};$i++) {
    printf $filehandle "%17s",$rdfdata[$ri][$i][0];
    for($j=1;$j<@{$rdfdata[$ri][$i]};$j++) {
      printf $filehandle " %17.10g",$rdfdata[$ri][$i][$j];
    }
    print $filehandle "\n";
  }
  close($filehandle);
}


sub read_zdndat_file {
  my $filename = $_[0];
  my $zi       = $_[1];
  my($filehandle,$i,$j,@linedata,$nzdns,$mxzdn,$a1,$a2,$iall);
  if($zi<0) {
    print "**** error: index \$zi must not be smaller than zero in sub read_zdndat_file\n";
    return -1;
  }
  if(not open($filehandle, "<", $filename)) {
    print "**** error: Can't open ZDNDAT file \"$filename\": $!\n";
    return 1;
  }
  undef @{$zdndata[$zi]};
  undef @{$zdnnames[$zi]};
  $zdnnames[$zi][0]="zpos_[A]";
  if(not $_=<$filehandle>) {close($filehandle); return -1} #config title
  if(not $_=<$filehandle>) {close($filehandle); return -1} # n(RDFs) mxzdn
  $_ =~ s/^\s+//; $_ =~ s/\s+$//;
  ($nzdns,$mxzdn) = split(/\s+/, $_);
  $iall=$nzdns+1;
  $zdnnames[$zi][$iall]="all_types";
  for($i=1;$i<=$nzdns;$i++) {
    if(not $_=<$filehandle>) {close($filehandle); return -1} # atom types
    $_ =~ s/^\s+//; $_ =~ s/\s+$//;
    $zdnnames[$zi][$i] = $_;
    print $zdnnames[$zi][$i],"\n";
    for($j=0;$j<$mxzdn;$j++) {
      if(not $_=<$filehandle>) {close($filehandle); return -1} # atom types
      $_ =~ s/^\s+//; $_ =~ s/\s+$//;
      @linedata = split(/\s+/, $_);
      $zdndata[$zi][$j][$i]     = $linedata[1];
      $zdndata[$zi][$j][$iall] += $linedata[1];
      $zdndata[$zi][$j][0]      = $linedata[0] if($i==1);
    }
  }
  for($j=0;$j<$mxzdn;$j++) {
    $zdndata[$zi][$j][$iall]/=$nzdns;
  }
  close($filehandle);
  return 0;
} # end subroutine read_zdndat_file


sub print_zdndat_data {
  my $filename = $_[0];
  my $zi       = $_[1];
  my($filehandle,$i,$j);
  if($zi<0) {
    print "**** error: index \$zi must not be smaller than zero in sub read_zdndat_file\n";
    return -1;
  }
  if(not open($filehandle, ">", $filename)) {
    print "**** error: Can't open output file \"$filename\": $!\n";
    return 1;
  }
  printf $filehandle "#%16s",$zdnnames[$zi][0];
  for($i=1;$i<@{$zdnnames[$zi]};$i++) {
    printf $filehandle " %17s",$zdnnames[$zi][$i];
  }
  print $filehandle "\n";
  for($i=0;$i<@{$zdndata[$zi]};$i++) {
    printf $filehandle "%17s",$zdndata[$zi][$i][0];
    for($j=1;$j<@{$zdndata[$zi][$i]};$j++) {
      printf $filehandle " %17.10g",$zdndata[$zi][$i][$j];
    }
    print $filehandle "\n";
  }
  close($filehandle);
}


sub read_xyz_timestep {
  my $filehandle = $_[0];
  my $fi         = $_[1]; #index of corresponding FIELD file
  my $ci         = $_[2];
  my($t,$m,$a,$molindex,@linedata);
  
  return 1 if($ci<0);
  
  # read header
  if(not $_=<$filehandle>) {return -1} #end of file
  $_ =~ s/^\s+//; $_ =~ s/\s+$//;
  @linedata = split(/\s+/, $_);
  if(not check_integer($linedata[0])) {
    print "**** error: expected numatoms in xyz file on line $. but got\n$_\n";
    return 1;
  }
  if(@linedata>1) {
    print "**** warning: expected numatoms but found additional junk on line $.:\n$_\n";
  }
  
  #reset values
  @{$cdata[$ci]} = ();
  @{$cell[$ci]}  = ();
  @{$size[$ci]}  = ();
  $molindex      = 0;
  
  $frame_numatoms[$ci] = $linedata[0];
  $config_key[$ci]     = 0;
  $periodic_key[$ci]   = 0;
  $config_title[$ci]   = <$filehandle>;
  
  # read atom data
  if($fi<0) {
    for($a=0;$a<$frame_numatoms[$ci];$a++) {
      if(not $_=<$filehandle>) {
	print "**** error: unexpected end of xyz file on line $.\n";
	return 1
      } #end of file
      $_ =~ s/^\s+//; $_ =~ s/\s+$//;
      @linedata = split(/\s+/, $_);
      if(/^\s*#/) { $a--; next; }
      if(@linedata<4) {
	print "**** error: expected atom data on line $. but got:\n$_\n";
	return 1;
      }
      ($cdata[$ci][0][0][$a][9],@{$cdata[$ci][0][0][$a]}[0..8]) = @linedata[0..9];
    }
  } else { # read with field information
    if($frame_numatoms[$ci] != $field_numatoms[$fi]) {
      print "**** error: number of atoms in FIELD file ($field_numatoms[$fi]) does not match number of atoms in xyz ($frame_numatoms[$ci]) file\n";
      return 1;
    }
    for($t=0;$t<$field_nummols[$fi];$t++) {
      for($m=0;$m<$mol_numents[$fi][$t];$m++) {
	for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
	  if(not $_=<$filehandle>) {
	    print "**** error: unexpected end of xyz file on line $.\n";
	    return 1
	  } #end of file
	  $_ =~ s/^\s+//; $_ =~ s/\s+$//;
	  @linedata = split(/\s+/, $_);
	  if(/^\s*#/) { $a--; next; }
	  if(@linedata<4) {
	    print "**** error: expected atom data on line $. but got:\n$_\n";
	    return 1;
	  }
	  @{$cdata[$ci][$t][$m][$a]}[0..8] = @linedata[1..9];
	  $cdata[$ci][$t][$m][$a][9] = $mol_atomdata[$fi][$t][$a][0];
	}
      }
    }
  }
  return 0;
} # end subroutine read_xyz_timestep


sub write_xyz_timestep {
  my $filehandle = $_[0];
  my $ci         = $_[1];
  my $title      = $_[2];
  my $fmt        = $_[3]; #output format: 0/undef=element 1=original name
  my($t,$m,$a);
  return 1 if($ci<0);
  $title =~ s/\s+$//;
  print $filehandle "$frame_numatoms[$ci]\n";
  print $filehandle "$title\n";
  if(not defined($fmt) or $fmt==0) {
    for($t=0;$t<@{$cdata[$ci]};$t++) {
      for($m=0;$m<@{$cdata[$ci][$t]};$m++) {
	for($a=0;$a<@{$cdata[$ci][$t][$m]};$a++) {
	  printf $filehandle "%-5s %12.5f %12.5f %12.5f\n", $element_names{uc($cdata[$ci][$t][$m][$a][9])},
	  @{$cdata[$ci][$t][$m][$a]}[0..2];
	}
      }
    }
  } else {
    for($t=0;$t<@{$cdata[$ci]};$t++) {
      for($m=0;$m<@{$cdata[$ci][$t]};$m++) {
	for($a=0;$a<@{$cdata[$ci][$t][$m]};$a++) {
	  printf $filehandle "%-5s %12.5f %12.5f %12.5f\n", $cdata[$ci][$t][$m][$a][9],
	  @{$cdata[$ci][$t][$m][$a]}[0..2];
	}
      }
    }
  }
  return 0;
}

sub read_mol2_file {
  my $filename = $_[0];
  my $mi       = $_[1]; # index for mol2 molecule
  my($filehandle,$keyword,@linedata,$atomcount,$bondcount,$entrynumber,$i,$addtype,%found);
  if(not open($filehandle, "<", $filename)) {
    print "**** error: Can't open MOL2 file $filename: $!\n";
    return 1;
  }
  @{$mol2_atomdata[$mi]}     = ();
  @{$mol2_bonddata[$mi]}     = ();
  @{$mol2_substdata[$mi]}    = ();
  @{$mol2_numatoms[$mi]}     = ();
  @{$mol2_numbonds[$mi]}     = ();
  @{$mol2_numsubst[$mi]}     = ();
  @{$mol2_numfeat[$mi]}      = ();
  @{$mol2_numsets[$mi]}      = ();
  @{$mol2_atomtypes[$mi]}    = ();
  $mol2_name[$mi]            = '';
  $mol2_comment[$mi]         = '';
  $mol2_chargemethod[$mi]    = '';
  $mol2_moltype[$mi]         = '';
  $mol2_statusbits[$mi]      = '';
  $mol2_charge[$mi]          = 0;
  while(<$filehandle>) {
    $_ =~ s/^\s+//; $_ =~ s/\s+$//;
    next if(/^#/ or length($_)==0);
    if(/^\@<TRIPOS>/) {
      ($keyword) = /^\@<TRIPOS>(\S+)/;
      if(exists($found{$keyword})){
	print "**** error in sub read_mol2_file: multiple occurrence of keyword '$keyword'\n";
	return 1;
      }
      $found{$keyword}=1;
      $entrynumber=0;
      next;
    }
    @linedata = split(/\s+/,$_);
    if($keyword=~/^MOLECULE/i) {
      if($entrynumber==0) {
	$mol2_name[$mi]     = $_; }
      elsif($entrynumber==1) {
	$mol2_numatoms[$mi] = $linedata[0];
	$mol2_numbonds[$mi] = $linedata[1];
	$mol2_numsubst[$mi] = $linedata[2];
	$mol2_numfeat[$mi]  = $linedata[3];
	$mol2_numsets[$mi]  = $linedata[4]; }
      elsif($entrynumber==2) {$mol2_moltype[$mi]      = $_;}
      elsif($entrynumber==3) {$mol2_chargemethod[$mi] = $_;}
      elsif($entrynumber==4) {$mol2_statusbits[$mi]   = $_;}
      elsif($entrynumber==5) {$mol2_comment[$mi]      = $_;}
    } elsif($keyword=~/^ATOM/i) {
      $mol2_atomdata[$mi][$entrynumber][0] = $linedata[2];   # x
      $mol2_atomdata[$mi][$entrynumber][1] = $linedata[3];   # y
      $mol2_atomdata[$mi][$entrynumber][2] = $linedata[4];   # z
      $mol2_atomdata[$mi][$entrynumber][3] = $linedata[1];   # name
      $mol2_atomdata[$mi][$entrynumber][4] = $linedata[5];   # type
      $mol2_atomdata[$mi][$entrynumber][5] = $linedata[0];   # atomid
      $mol2_atomdata[$mi][$entrynumber][6] = $linedata[6]-1; # substid
      $mol2_atomdata[$mi][$entrynumber][7] = $linedata[7];   # substname
      $mol2_atomdata[$mi][$entrynumber][8] = $linedata[8];   # charge
      $mol2_atomdata[$mi][$entrynumber][9] = $linedata[9];   # statusbit
      $mol2_charge[$mi] += $linedata[8];
      $addtype = 1;
      for($i=0;$i<@{$mol2_atomtypes[$mi]};$i++) {
	if($linedata[5] eq $mol2_atomtypes[$mi][$i]) {$addtype=0;last}
      }
      push(@{$mol2_atomtypes[$mi]},$linedata[5]) if($addtype);
    } elsif($keyword=~/^BOND/i) {
      $mol2_bonddata[$mi][$entrynumber][0] = $linedata[1]-1; # atom1
      $mol2_bonddata[$mi][$entrynumber][1] = $linedata[2]-1; # atom2
      $mol2_bonddata[$mi][$entrynumber][2] = $linedata[3];   # type
      $mol2_bonddata[$mi][$entrynumber][3] = $linedata[4];   # statusbit
    } elsif($keyword=~/^SUBSTRUCTURE/i) {
      $mol2_substdata[$mi][$entrynumber][0] = $linedata[1];   # name
      $mol2_substdata[$mi][$entrynumber][1] = $linedata[2]-1; # root atom
      $mol2_substdata[$mi][$entrynumber][2] = $linedata[3];   # type
      $mol2_substdata[$mi][$entrynumber][3] = $linedata[4];   # dicttype
      $mol2_substdata[$mi][$entrynumber][4] = $linedata[5];   # chain
      $mol2_substdata[$mi][$entrynumber][5] = $linedata[6];   # subst_type
      $mol2_substdata[$mi][$entrynumber][6] = $linedata[7];   # inter_bonds
      $mol2_substdata[$mi][$entrynumber][7] = $linedata[8];   # status
      $mol2_substdata[$mi][$entrynumber][8] = join(" ",@linedata[9..$#linedata]) if($#linedata>8); # comment
    } else {
      print "**** error: unknown keyword '$keyword'  on line $. of MOL2 file $filename!\n";
      close($filehandle);
      return 1;
    }
    $entrynumber++;
  }
  close($filehandle);
  if(not exists($found{"MOLECULE"})) {
    print "**** error: missing keyword 'MOLECULE' in file $filename\n";
    return 1;
  }
  if(not exists($found{"ATOM"})) {
    print "**** error: missing keyword 'ATOM' in file $filename\n";
    return 1;
  }
  return 0;
} # end subroutine read_mol2_file


sub mol2_to_cdata {
  #example @{$cdata[0][$t][$m]} = mol2_to_cdata(0);
  my $mi=$_[0];
  my ($a,@carr);
  # try to get index from molecule name if necessary
  if(not check_integer($mi)) {
    for($a=0;$a<@mol2_name;$a++) {
      if($mol2_name[$a] eq $mi) {
	$mi=$a;
	last;
      }
    }
    if(not check_integer($mi)) {
      print "**** error in sub mol2_to_cdata: a molecule named $mi could not be found!\n";
      return 1;
    }
  } elsif(not exists($mol2_name[$mi])) {
    print "**** error in sub mol2_to_cdata: molecule number $mi does not exist!\n";
    return 1;
  }
  @carr=();
  for($a=0;$a<$mol2_numatoms[$mi];$a++) {
    $carr[$a][0]  = $mol2_atomdata[$mi][$a][0];
    $carr[$a][1]  = $mol2_atomdata[$mi][$a][1];
    $carr[$a][2]  = $mol2_atomdata[$mi][$a][2];
    $carr[$a][9]  = uc($mol2_atomdata[$mi][$a][4]);
    $carr[$a][10] = 0;  # molid
    $carr[$a][11] = $a; # atomid
    $carr[$a][13] = $mol2_atomdata[$mi][$a][8];
  }
  return @carr;
}

sub write_mol2_file {
# !!!! unfinished !!!!
  my $filename = $_[0];
  my $mi       = $_[1];
  my($filehandle,$a);
  if(not open($filehandle, ">", $filename)) {
    print "**** error: Can't open FIELD file \"$filename\": $!\n";
    return 1;
  }
  # header
  print $filehandle "\@<TRIPOS>MOLECULE\n";
  print $filehandle "$mol2_name[$mi]\n";
  printf $filehandle "%5u %5u %5u %5u %5u\n",$mol2_numatoms[$mi],
         $mol2_numbonds[$mi],$mol2_numsubst[$mi],$mol2_numfeat[$mi],$mol2_numsets[$mi];
  print $filehandle "$mol2_moltype[$mi]\n";
  print $filehandle "$mol2_chargemethod[$mi]\n\n\n";
  # ATOMS section
  print $filehandle "\@<TRIPOS>ATOM";
  for($a=0;$a<$mol2_numatoms[$mi];$a++) {
    printf $filehandle "\n%7u %-8s %9.4f %9.4f %9.4f %-7s %3u %-7s %10.6f",
    $a+1,                         # atomid
    $mol2_atomdata[$mi][$a][3],   # name
    $mol2_atomdata[$mi][$a][0],   # x
    $mol2_atomdata[$mi][$a][1],   # y
    $mol2_atomdata[$mi][$a][2],   # z
    $mol2_atomdata[$mi][$a][4],   # type
    $mol2_atomdata[$mi][$a][6]+1, # substid
    $mol2_atomdata[$mi][$a][7],   # substname
    $mol2_atomdata[$mi][$a][8],   # charge
    $mol2_atomdata[$mi][$a][9];   # status bits
  }
  # BONDS section
  print $filehandle "\n\@<TRIPOS>BOND";
  for($a=0;$a<$mol2_numbonds[$mi];$a++) {
    printf $filehandle "\n%6u %4u %4u %s %s",
    $a+1,
    $mol2_bonddata[$mi][$a][0]+1, # atom 1
    $mol2_bonddata[$mi][$a][1]+1, # atom 2
    $mol2_bonddata[$mi][$a][2],   # bond type
    $mol2_bonddata[$mi][$a][3];   # status bits
  }
  # SUBSTRUCTURE section
  print $filehandle "\n\@<TRIPOS>SUBSTRUCTURE";
  for($a=0;$a<$mol2_numsubst[$mi];$a++) {
    printf $filehandle "\n%6u %-8s %4u %-14s %4u %-5s %-5s %3u %s",
    $a+1,
    $mol2_substdata[$mi][$a][0],   # subst name
    $mol2_substdata[$mi][$a][1]+1, # first atom index
    $mol2_substdata[$mi][$a][2],   # subst type;
    $mol2_substdata[$mi][$a][3],   # dict type;
    $mol2_substdata[$mi][$a][4],   # chain;
    $mol2_substdata[$mi][$a][5],   # sub type;
    $mol2_substdata[$mi][$a][6],   # inter_bonds;
    $mol2_substdata[$mi][$a][7],   # status;
    $mol2_substdata[$mi][$a][8],   # comment;
  }
  # still missing many directives !!!
  close($filehandle);
}

sub connection_depth {
  # !!!! not finished yet !!!!
  # returns the (minimum) number of bonds between two atoms in a molecule
  my $fi = $_[0];
  my $t  = $_[1];
  my $a  = $_[2];
  my $b  = $_[3];
  
}

sub read_gaff_file {
  my $filename = $_[0];
  my($filehandle,$section,$newsection,$i,@linedata);
  if(not open($filehandle, "<", $filename)) {
    print "**** error: Can't open GAFF file $filename:\n$!\n";
    return 1;
  }
  %gaff_atoms     = ();
  %gaff_bonds     = ();
  %gaff_angles    = ();
  %gaff_dihedrals = ();
  %gaff_vdwparam  = ();
  $_=<$filehandle>;
  $section=0;
  while(<$filehandle>) {
    $_ =~ s/^\s+//; $_ =~ s/\s+$//;
    if(length($_) == 0) {
      $newsection=1;
      next;
    }
    next if /^#/;
    last if /^END/i;
    if($newsection) {
      if(/^\S+\s?-\S+/ # e.g. "ha- o"
      or /^\S+\s+[+-]?\d+\.?\d*([eEdD][+-]?\d+)?/ # e.g. "ha -0."
      or /^\S+\s+[+-]?\d*\.?\d+([eEdD][+-]?\d+)?/) { # e.g. "hb -.3"
	$section++;
	$newsection=0;
      } else {
	next;
      }
    }
    if($section == 0) { # atom type data
      @linedata = split(/[\s-]+/,$_);
      # index:   0    1    2
      # format: a1 mass   ??
      @{$gaff_atoms{$linedata[0]}}=@linedata[1..2];
    }
    if($section == 1) { # bond data
      @linedata = split(/[\s-]+/,$_);
      # index:   0   1  2   3    4   5   6
      # format: a1  a2  c  r0  src  ??  ??
      $linedata[0]=uc($linedata[0]);
      $linedata[1]=uc($linedata[1]);
      @linedata[0..1] = sort(@linedata[0..1]);
      @{$gaff_bonds{"$linedata[0]-$linedata[1]"}}=@linedata[2..6];
    } elsif($section == 2) { # angle data
      @linedata = split(/[\s-]+/,$_);
      # index:   0   1   2  3    4    5   6   7
      # format: a1  a2  a3  c  opt  src  ??  ??
      $linedata[0]=uc($linedata[0]);
      $linedata[1]=uc($linedata[1]);
      $linedata[2]=uc($linedata[2]);
      if(($linedata[0] <= $linedata[2])) {
	@{$gaff_bonds{"$linedata[0]-$linedata[1]-$linedata[2]"}}=@linedata[3..7];
      } else {
	@{$gaff_bonds{"$linedata[2]-$linedata[1]-$linedata[0]"}}=@linedata[3..7];
      }
    } elsif($section == 3) { # dihedral data
      @linedata = split(/[\s-]+/,$_);
      # index:   0   1   2   3   4  5    6    7    8
      # format: a1  a2  a3  a4  ??  c  opt  rep  src
      $linedata[0]=uc($linedata[0]);
      $linedata[1]=uc($linedata[1]);
      $linedata[2]=uc($linedata[2]);
      $linedata[3]=uc($linedata[3]);
      if($linedata[0] eq "X" or ($linedata[0] <= $linedata[3])) {
	@{$gaff_bonds{"$linedata[0]-$linedata[1]-$linedata[2]-$linedata[3]"}}=@linedata[4..8];
      } else {
	@{$gaff_bonds{"$linedata[3]-$linedata[2]-$linedata[1]-$linedata[0]"}}=@linedata[4..8];
      }
    } elsif($section == 4) { # improper data
      @linedata = split(/[\s-]+/,$_);
      # index:   0   1   2   3  4    5    6    7
      # format: a1  a2  a3  a4  c  opt  rep  src
      $linedata[0]=uc($linedata[0]);
      $linedata[1]=uc($linedata[1]);
      $linedata[2]=uc($linedata[2]);
      $linedata[3]=uc($linedata[3]);
      if($linedata[0] eq "X" or ($linedata[0] <= $linedata[3])) {
	@{$gaff_bonds{"$linedata[0]-$linedata[1]-$linedata[2]-$linedata[3]"}}=@linedata[4..7];
      } else {
	@{$gaff_bonds{"$linedata[3]-$linedata[2]-$linedata[1]-$linedata[0]"}}=@linedata[4..7];
      }
    } elsif($section == 5) { # vdw data
      @linedata = split(/\s+/,$_);
      # index:   0  1  2
      # format: a1 r0  c
      $linedata[0]=uc($linedata[0]);
      @{$gaff_vdwparam{$linedata[0]}}=@linedata[1..2];
    }
  }
  close($filehandle);
  return 0
} # end subroutine read_gaff_file

sub read_frcmod_file {
  my $filename = $_[0];
  my($filehandle,$section,$newsection,$i,@linedata);
  if(not open($filehandle, "<", $filename)) {
    print "**** error: Can't open frcmod file $filename:\n$!\n";
    return 1;
  }
  $_=<$filehandle>; # title
  while(<$filehandle>) {
    $_ =~ s/^\s+//; $_ =~ s/\s+$//;
    if(length($_) == 0) {
      $newsection=1;
      next;
    }
    next if /^#/;
    last if /^END/i;
    if($newsection) {
      if(/MASS/) {
	$section = 0;
      } elsif(/BOND/) {
	$section = 1;
      } elsif(/ANGLE/) {
	$section = 2;
      } elsif(/DIHE/) {
	$section = 3;
      } elsif(/IMPR/) {
	$section = 4;
      } elsif(/NONBON/) {
	$section = 5;
      } else {
	print "**** error: unknown section \"$_\"  on line $. of frcmod file $filename!\n";
	return 1;
	$section = -1;
      }
      $newsection=0;
    }
    if($section == 0) { # atom type data
      @linedata = split(/[\s-]+/,$_);
      # index:   0    1    2
      # format: a1 mass   ??
      @{$gaff_atoms{$linedata[0]}}=@linedata[1..2];
    }
    if($section == 1) { # bond data
      @linedata = split(/[\s-]+/,$_);
      # index:   0   1  2   3    4   5   6
      # format: a1  a2  c  r0  src  ??  ??
      $linedata[0]=uc($linedata[0]);
      $linedata[1]=uc($linedata[1]);
      @linedata[0..1] = sort(@linedata[0..1]);
      @{$gaff_bonds{"$linedata[0]-$linedata[1]"}}=@linedata[2..6];
    } elsif($section == 2) { # angle data
      @linedata = split(/[\s-]+/,$_);
      # index:   0   1   2  3    4    5   6   7
      # format: a1  a2  a3  c  opt  src  ??  ??
      $linedata[0]=uc($linedata[0]);
      $linedata[1]=uc($linedata[1]);
      $linedata[2]=uc($linedata[2]);
      if(($linedata[0] <= $linedata[2])) {
	@{$gaff_bonds{"$linedata[0]-$linedata[1]-$linedata[2]"}}=@linedata[3..7];
      } else {
	@{$gaff_bonds{"$linedata[2]-$linedata[1]-$linedata[0]"}}=@linedata[3..7];
      }
    } elsif($section == 3) { # dihedral data
      @linedata = split(/[\s-]+/,$_);
      # index:   0   1   2   3   4  5    6    7    8
      # format: a1  a2  a3  a4  ??  c  opt  rep  src
      $linedata[0]=uc($linedata[0]);
      $linedata[1]=uc($linedata[1]);
      $linedata[2]=uc($linedata[2]);
      $linedata[3]=uc($linedata[3]);
      if($linedata[0] eq "X" or ($linedata[0] <= $linedata[3])) {
	@{$gaff_bonds{"$linedata[0]-$linedata[1]-$linedata[2]-$linedata[3]"}}=@linedata[4..8];
      } else {
	@{$gaff_bonds{"$linedata[3]-$linedata[2]-$linedata[1]-$linedata[0]"}}=@linedata[4..8];
      }
    } elsif($section == 4) { # improper data
      @linedata = split(/[\s-]+/,$_);
      # index:   0   1   2   3  4    5    6    7
      # format: a1  a2  a3  a4  c  opt  rep  src
      $linedata[0]=uc($linedata[0]);
      $linedata[1]=uc($linedata[1]);
      $linedata[2]=uc($linedata[2]);
      $linedata[3]=uc($linedata[3]);
      if($linedata[0] eq "X" or ($linedata[0] <= $linedata[3])) {
	@{$gaff_bonds{"$linedata[0]-$linedata[1]-$linedata[2]-$linedata[3]"}}=@linedata[4..7];
      } else {
	@{$gaff_bonds{"$linedata[3]-$linedata[2]-$linedata[1]-$linedata[0]"}}=@linedata[4..7];
      }
    } elsif($section == 5) { # vdw data
      @linedata = split(/\s+/,$_);
      # index:   0  1  2
      # format: a1 r0  c
      $linedata[0]=uc($linedata[0]);
      @{$gaff_vdwparam{$linedata[0]}}=@linedata[1..2];
    }
  }
  close($filehandle);
  return 0
} # end subroutine read_frcmod_file


sub lj_ab_to_sigeps {
  # convert AB form to epsilon/sigma
  my @res;
  $res[0] = $_[1]*$_[1]/(4.0*$_[0]); # epsilon
  $res[1] = ($_[0]/$_[1])**(1.0/6.0); # sigma
  return @res;
}

sub calc_vdw_params {
  my $kcal = 4.184;
  my $err = 0;
  my($sigma,$epsilon);
  if(not exists $gaff_vdwparam{$_[0]}[0]) {
    print "**** error: did not find vdw parameters for atom type $_[0]!\n"; $err=1;
  }
  if(not exists $gaff_vdwparam{$_[1]}[0]) {
    print "**** error: did not find vdw parameters for atom type $_[1]!\n"; $err=1;
  }
  return(undef, undef) if $err;
  $sigma   = ($gaff_vdwparam{$_[0]}[0] + $gaff_vdwparam{$_[1]}[0]) / (2**(1/6));
  $epsilon = sqrt($gaff_vdwparam{$_[0]}[1] * $gaff_vdwparam{$_[1]}[1]) * $kcal;
  return($epsilon,$sigma);
} # end subroutine calc_vdw_params

sub calc_dihedral_angle {
  # calculates the ABSOLUTE of the dihedral angle (without sign)
  # (for dihedral angle with sign use calc_dihedral_angle2)
  my $ci = $_[0];
  my $fi = $_[1];
  my $d  = $_[2]; # index for @mol_dihedraldata
  my $t  = $_[3];
  my $m  = $_[4];
  my ($a1,$a2,$a3,$a4)=@{$mol_dihedraldata[$fi][$t][$d]}[1..4];
  my (@vec,@vecc,@vecn1,@vecn2);
  $vecc[0] = $cdata[$ci][$t][$m][$a3][0]-$cdata[$ci][$t][$m][$a2][0];
  $vecc[1] = $cdata[$ci][$t][$m][$a3][1]-$cdata[$ci][$t][$m][$a2][1];
  $vecc[2] = $cdata[$ci][$t][$m][$a3][2]-$cdata[$ci][$t][$m][$a2][2];
  $vec[0]  = $cdata[$ci][$t][$m][$a1][0]-$cdata[$ci][$t][$m][$a2][0];
  $vec[1]  = $cdata[$ci][$t][$m][$a1][1]-$cdata[$ci][$t][$m][$a2][1];
  $vec[2]  = $cdata[$ci][$t][$m][$a1][2]-$cdata[$ci][$t][$m][$a2][2];
  @vecn1=vector_product(\@vec,\@vecc);
  @vecn1=normalize_vector(\@vecn1);
  $vec[0]  = $cdata[$ci][$t][$m][$a4][0]-$cdata[$ci][$t][$m][$a3][0];
  $vec[1]  = $cdata[$ci][$t][$m][$a4][1]-$cdata[$ci][$t][$m][$a3][1];
  $vec[2]  = $cdata[$ci][$t][$m][$a4][2]-$cdata[$ci][$t][$m][$a3][2];
  @vecn2=vector_product(\@vec,\@vecc);
  @vecn2=normalize_vector(\@vecn2);
  return acos(dot_product(\@vecn1,\@vecn2));
}


sub calc_dihedral_angle2 {
  # calculates the dihedral angle with a sign
  my $ci = $_[0];
  my $fi = $_[1];
  my $d  = $_[2]; # index for @mol_dihedraldata
  my $t  = $_[3];
  my $m  = $_[4];
  my ($a1,$a2,$a3,$a4)=@{$mol_dihedraldata[$fi][$t][$d]}[1..4];
  my (@vec1,@vec2,@vec3,@cross1,@cross2,@vectmp);
  $vec1[0] = $cdata[$ci][$t][$m][$a2][0]-$cdata[$ci][$t][$m][$a1][0];
  $vec1[1] = $cdata[$ci][$t][$m][$a2][1]-$cdata[$ci][$t][$m][$a1][1];
  $vec1[2] = $cdata[$ci][$t][$m][$a2][2]-$cdata[$ci][$t][$m][$a1][2];
  $vec2[0] = $cdata[$ci][$t][$m][$a3][0]-$cdata[$ci][$t][$m][$a2][0];
  $vec2[1] = $cdata[$ci][$t][$m][$a3][1]-$cdata[$ci][$t][$m][$a2][1];
  $vec2[2] = $cdata[$ci][$t][$m][$a3][2]-$cdata[$ci][$t][$m][$a2][2];
  $vec3[0] = $cdata[$ci][$t][$m][$a4][0]-$cdata[$ci][$t][$m][$a3][0];
  $vec3[1] = $cdata[$ci][$t][$m][$a4][1]-$cdata[$ci][$t][$m][$a3][1];
  $vec3[2] = $cdata[$ci][$t][$m][$a4][2]-$cdata[$ci][$t][$m][$a3][2];
  @cross1 = vector_product(\@vec1,\@vec2);
  @cross2 = vector_product(\@vec2,\@vec3);
  @vectmp = vector_product(\@cross1,\@cross2);
  return atan2(dot_product(\@vectmp,\@vec2)/vector_length(@vec2),dot_product(\@cross1,\@cross2));
}


sub calc_angle {
  my $ci = $_[0];
  my $fi = $_[1];
  my $d  = $_[2]; # index for @mol_angledata
  my $t  = $_[3];
  my $m  = $_[4];
  my ($a1,$a2,$a3)=@{$mol_angledata[$fi][$t][$d]}[1..3];
  my (@vec1,@vec2);
  $vec1[0] = $cdata[$ci][$t][$m][$a1][0]-$cdata[$ci][$t][$m][$a2][0];
  $vec1[1] = $cdata[$ci][$t][$m][$a1][1]-$cdata[$ci][$t][$m][$a2][1];
  $vec1[2] = $cdata[$ci][$t][$m][$a1][2]-$cdata[$ci][$t][$m][$a2][2];
  $vec2[0] = $cdata[$ci][$t][$m][$a3][0]-$cdata[$ci][$t][$m][$a2][0];
  $vec2[1] = $cdata[$ci][$t][$m][$a3][1]-$cdata[$ci][$t][$m][$a2][1];
  $vec2[2] = $cdata[$ci][$t][$m][$a3][2]-$cdata[$ci][$t][$m][$a2][2];
  @vec1 = normalize_vector(\@vec1);
  @vec2 = normalize_vector(\@vec2);
  return acos(dot_product(\@vec1,\@vec2));
}


sub calc_angle_orthocell {
  # calc_angle_orthocell($ci, $t1,$m1,$a1, $t2,$m2,$a2, $t3,$m3,$a3);
  my $ci = $_[0];
  my $t1 = $_[1];
  my $m1 = $_[2];
  my $a1 = $_[3];
  my $t2 = $_[4];
  my $m2 = $_[5];
  my $a2 = $_[6];
  my $t3 = $_[7];
  my $m3 = $_[8];
  my $a3 = $_[9];
  my @vec1 = calc_dvec_orthocell($ci, $t1, $m1, $a1, $ci, $t2, $m2, $a2);
  my @vec2 = calc_dvec_orthocell($ci, $t3, $m3, $a3, $ci, $t2, $m2, $a2);
  @vec1 = normalize_vector(\@vec1);
  @vec2 = normalize_vector(\@vec2);
  return acos(dot_product(\@vec1,\@vec2));
}


sub calc_temperature {
  my $fi = $_[0];
  my $ci = $_[1];
  my $t  = $_[2]; # optional (otherwise average is returned)
  my $m  = $_[3]; # optional (otherwise average for mol type $t is returned
  my($a,$vsq,$ekin,$numatoms);
  $ekin     = 0;
  $vsq      = 0;
  $numatoms = 0;
  if(defined($t)) {
    if(defined($m)) {
      for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
	next if($mol_atomdata[0][$t][$a][4]>0); # ommit frozen atoms
	$vsq  = 0;
	$vsq += $cdata[$ci][$t][$m][$a][3]*$cdata[$ci][$t][$m][$a][3];
	$vsq += $cdata[$ci][$t][$m][$a][4]*$cdata[$ci][$t][$m][$a][4];
	$vsq += $cdata[$ci][$t][$m][$a][5]*$cdata[$ci][$t][$m][$a][5];
# 	$vsq = dot_product([@{$cdata[$ci][$t][$m][$a]}[3..5]],[@{$cdata[$ci][$t][$m][$a]}[3..5]]);
	$ekin += $mol_atomdata[$fi][$t][$a][1]*$vsq;
	$numatoms++;
      }
    } else {
      for($m=0;$m<$mol_numents[$fi][$t];$m++) {
	for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
	  next if($mol_atomdata[0][$t][$a][4]>0); # ommit frozen atoms
	  $vsq  = 0;
	  $vsq += $cdata[$ci][$t][$m][$a][3]*$cdata[$ci][$t][$m][$a][3];
	  $vsq += $cdata[$ci][$t][$m][$a][4]*$cdata[$ci][$t][$m][$a][4];
	  $vsq += $cdata[$ci][$t][$m][$a][5]*$cdata[$ci][$t][$m][$a][5];
# 	  $vsq = dot_product([@{$cdata[$ci][$t][$m][$a]}[3..5]],[@{$cdata[$ci][$t][$m][$a]}[3..5]]);
	  $ekin += $mol_atomdata[$fi][$t][$a][1]*$vsq;
	  $numatoms++;
	}
      }
    }
  } else {
    for($t=0;$t<$field_nummols[$fi];$t++) {
      for($m=0;$m<$mol_numents[$fi][$t];$m++) {
	for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
	  next if($mol_atomdata[0][$t][$a][4]>0); # ommit frozen atoms
	  $vsq  = 0;
	  $vsq += $cdata[$ci][$t][$m][$a][3]*$cdata[$ci][$t][$m][$a][3];
	  $vsq += $cdata[$ci][$t][$m][$a][4]*$cdata[$ci][$t][$m][$a][4];
	  $vsq += $cdata[$ci][$t][$m][$a][5]*$cdata[$ci][$t][$m][$a][5];
# 	  $vsq = dot_product([@{$cdata[$ci][$t][$m][$a]}[3..5]],[@{$cdata[$ci][$t][$m][$a]}[3..5]]);
	  $ekin += $mol_atomdata[$fi][$t][$a][1]*$vsq;
	  $numatoms++;
	}
      }
    }
  }
  return $ekin/(3.0*$numatoms-6.0)/0.831451115;
}


sub scale_temperature {
  # scale the temperature of molecule(s) to desired value
  my $target = $_[0]; #target temperature
  my $fi = $_[1];
  my $ci = $_[2];
  my $t  = $_[3]; # optional (otherwise average is returned)
  my $m  = $_[4]; # optional (otherwise average for mol type $t is returned
  my($a,$factor);
  if(defined($t)) {
    if(defined($m)) {
      $factor = sqrt($target/calc_temperature($fi,$ci,$t,$m));
      for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
	next if($mol_atomdata[0][$t][$a][4]>0); # ommit frozen atoms
	@{$cdata[$ci][$t][$m][$a]}[3..5] = scale_vector([@{$cdata[$ci][$t][$m][$a]}[3..5]],$factor);
      }
    } else {
      $factor = sqrt($target/calc_temperature($fi,$ci,$t));
      for($m=0;$m<$mol_numents[$fi][$t];$m++) {
	for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
	  next if($mol_atomdata[0][$t][$a][4]>0); # ommit frozen atoms
	  @{$cdata[$ci][$t][$m][$a]}[3..5] = scale_vector([@{$cdata[$ci][$t][$m][$a]}[3..5]],$factor);
	}
      }
    }
  } else {
    $factor = sqrt($target/calc_temperature($fi,$ci));
    for($t=0;$t<$field_nummols[$fi];$t++) {
      for($m=0;$m<$mol_numents[$fi][$t];$m++) {
	for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
	  next if($mol_atomdata[0][$t][$a][4]>0); # ommit frozen atoms
	  @{$cdata[$ci][$t][$m][$a]}[3..5] = scale_vector([@{$cdata[$ci][$t][$m][$a]}[3..5]],$factor);
	}
      }
    }
  }
}


sub calc_rmsd {
  # calculate the root mean square displacement of a configuration
  # !!!! caution: molecules have to be remapped !!!!
  my $ci1 = $_[0]; # original configuration
  my $ci2 = $_[1]; # new configuration
  my $t   = $_[2]; # optional (otherwise average is returned)
  my $m   = $_[3]; # optional (otherwise average for mol type $t is returned
  my($a,@distvec,$dsq,$numatoms);
  $dsq      = 0;
  $numatoms = 0;
  if(defined($t)) {
    if(defined($m)) {
      for($a=0;$a<@{$cdata[$ci1][$t][$m]};$a++) {
	next if($mol_atomdata[0][$t][$a][4]>0); # ommit frozen atoms
	@distvec = vector_subst([@{$cdata[$ci2][$t][$m][$a]}[0..2]],[@{$cdata[$ci1][$t][$m][$a]}[0..2]]);
	$dsq += dot_product(\@distvec,\@distvec);
	$numatoms++;
      }
    } else {
      for($m=0;$m<@{$cdata[$ci1][$t]};$m++) {
	for($a=0;$a<@{$cdata[$ci1][$t][$m]};$a++) {
	  next if($mol_atomdata[0][$t][$a][4]>0); # ommit frozen atoms
	  @distvec = vector_subst([@{$cdata[$ci2][$t][$m][$a]}[0..2]],[@{$cdata[$ci1][$t][$m][$a]}[0..2]]);
	  $dsq += dot_product(\@distvec,\@distvec);
	  $numatoms++;
	}
      }
    }
  } else {
    for($t=0;$t<@{$cdata[$ci1]};$t++) {
      for($m=0;$m<@{$cdata[$ci1][$t]};$m++) {
	for($a=0;$a<@{$cdata[$ci1][$t][$m]};$a++) {
	  next if($mol_atomdata[0][$t][$a][4]>0); # ommit frozen atoms
	  @distvec = vector_subst([@{$cdata[$ci2][$t][$m][$a]}[0..2]],[@{$cdata[$ci1][$t][$m][$a]}[0..2]]);
	  $dsq += dot_product(\@distvec,\@distvec);
	  $numatoms++;
	}
      }
    }
  }
  return sqrt($dsq/$numatoms);
}

sub calc_minmaxpos {
  my $ci = $_[0];
  my $t  = $_[1]; # optional
  my ($m,$a);
  @{$minpos[$ci]}=( 9e20, 9e20, 9e20);
  @{$maxpos[$ci]}=(-9e20,-9e20,-9e20);
  if(defined($t)) {
    for($m=0;$m<@{$cdata[$ci][$t]};$m++) {
      for($a=0;$a<@{$cdata[$ci][$t][$m]};$a++) {
	$minpos[$ci][0] = min($minpos[$ci][0],$cdata[$ci][$t][$m][$a][0]);
	$minpos[$ci][1] = min($minpos[$ci][1],$cdata[$ci][$t][$m][$a][1]);
	$minpos[$ci][2] = min($minpos[$ci][2],$cdata[$ci][$t][$m][$a][2]);
	$maxpos[$ci][0] = max($maxpos[$ci][0],$cdata[$ci][$t][$m][$a][0]);
	$maxpos[$ci][1] = max($maxpos[$ci][1],$cdata[$ci][$t][$m][$a][1]);
	$maxpos[$ci][2] = max($maxpos[$ci][2],$cdata[$ci][$t][$m][$a][2]);
      }
    }
  } else {
    for($t=0;$t<@{$cdata[$ci]};$t++) {
      for($m=0;$m<@{$cdata[$ci][$t]};$m++) {
	for($a=0;$a<@{$cdata[$ci][$t][$m]};$a++) {
	  $minpos[$ci][0] = min($minpos[$ci][0],$cdata[$ci][$t][$m][$a][0]);
	  $minpos[$ci][1] = min($minpos[$ci][1],$cdata[$ci][$t][$m][$a][1]);
	  $minpos[$ci][2] = min($minpos[$ci][2],$cdata[$ci][$t][$m][$a][2]);
	  $maxpos[$ci][0] = max($maxpos[$ci][0],$cdata[$ci][$t][$m][$a][0]);
	  $maxpos[$ci][1] = max($maxpos[$ci][1],$cdata[$ci][$t][$m][$a][1]);
	  $maxpos[$ci][2] = max($maxpos[$ci][2],$cdata[$ci][$t][$m][$a][2]);
	}
      }
    }
  }
}

sub enlarge_system {
  # example: &enlarge_system(0,0, 1,1, 4,4,1, "rs")
  my $cii=$_[0];
  my $fii=$_[1];
  my $cio=$_[2];
  my $fio=$_[3];
  my @rep=@_[4..6];
  my $flags=$_[7]; # optional flags: m=omit remap, s=omit shift
  my @spacing=@{$_[8]} if($#ARGV>7); # optional spacing
  @spacing=(0,0,0) if (not @spacing);
  my ($x,$y,$z,$t,$m,$mi,$a,@remapaxes,@shift,$omitremap,@tmp);
  if(not @{$cell[$cii]} or $periodic_key[$cii]==0) {
    print "****error: cell information is required for system replication!\n";
    return 1;
  }
  &copy_field($fii,$fio);
  
  # shift for centering
  if($flags=~/s/i) {
    @shift=(0,0,0);
  } else {
    $shift[0] = (-($cell[$cii][0][0]+$spacing[0])*($rep[0]-1)-$cell[$cii][1][0]*($rep[1]-1)-$cell[$cii][2][0]*($rep[2]-1))/2.0;
    $shift[1] = (-$cell[$cii][0][1]*($rep[0]-1)-($cell[$cii][1][1]+$spacing[1])*($rep[1]-1)-$cell[$cii][2][1]*($rep[2]-1))/2.0;
    $shift[2] = (-$cell[$cii][0][2]*($rep[0]-1)-$cell[$cii][1][2]*($rep[1]-1)-($cell[$cii][2][2]+$spacing[2])*($rep[2]-1))/2.0;
  }
  if($flags=~/r/i) {
    $omitremap=1;
  } else {
    $omitremap=0;
    if($periodic_key[$cii]==6) {
      @remapaxes=(0,1);
    } else {
      @remapaxes=(0,1,2);
    }
  }
  for($t=0;$t<$field_nummols[$fio];$t++) {
    $mi=0;
    @{$tmp[$t]}=();
    for($m=0;$m<$mol_numents[$fii][$t];$m++) {
      remap_molecule2($fii,$cii,$t,$m) unless($omitremap);
      for($x=0;$x<$rep[0];$x++) {
	for($y=0;$y<$rep[1];$y++) {
	  for($z=0;$z<$rep[2];$z++) {
	    for($a=0;$a<$mol_numatoms[$fii][$t];$a++) {
	      $tmp[$t][$mi][$a][0] = $cdata[$cii][$t][$m][$a][0] + $shift[0]
	      + $x*($cell[$cii][0][0]+$spacing[0]) + $y*$cell[$cii][1][0] + $z*$cell[$cii][2][0];
	      $tmp[$t][$mi][$a][1] = $cdata[$cii][$t][$m][$a][1] + $shift[1]
	      + $x*$cell[$cii][0][1] + $y*($cell[$cii][1][1]+$spacing[1]) + $z*$cell[$cii][2][1];
	      $tmp[$t][$mi][$a][2] = $cdata[$cii][$t][$m][$a][2] + $shift[2]
	      + $x*$cell[$cii][0][2] + $y*$cell[$cii][1][2] + $z*($cell[$cii][2][2]+$spacing[2]);
	      @{$tmp[$t][$mi][$a]}[3..13]=@{$cdata[$cii][$t][$m][$a]}[3..13];
	    }
	    $mi++;
	  }
	}
      }
    }
    $mol_numents[$fio][$t] *= $rep[0]*$rep[1]*$rep[2];
  }
  @{$cdata[$cio]}=@tmp;
  $cell[$cio][0][0] = ($cell[$cii][0][0]+$spacing[0])*$rep[0];
  $cell[$cio][0][1] = $cell[$cii][0][1]*$rep[0];
  $cell[$cio][0][2] = $cell[$cii][0][2]*$rep[0];
  $cell[$cio][1][0] = $cell[$cii][1][0]*$rep[1];
  $cell[$cio][1][1] = ($cell[$cii][1][1]+$spacing[1])*$rep[1];
  $cell[$cio][1][2] = $cell[$cii][1][2]*$rep[1];
  $cell[$cio][2][0] = $cell[$cii][2][0]*$rep[2];
  $cell[$cio][2][1] = $cell[$cii][2][1]*$rep[2];
  $cell[$cio][2][2] = ($cell[$cii][2][2]+$spacing[2])*$rep[2];
  $periodic_key[$cio] = $periodic_key[$cii];
  $config_key[$cio]   = $config_key[$cii];
  $frame_numatoms[$cio] = $frame_numatoms[$cii]*$rep[0]*$rep[1]*$rep[2];
  $field_numatoms[$fio] = $field_numatoms[$fii]*$rep[0]*$rep[1]*$rep[2];
  $size[$cio][0] = ($cell[$cio][0][0] + $cell[$cio][1][0] + $cell[$cio][2][0])/2;
  $size[$cio][1] = ($cell[$cio][0][1] + $cell[$cio][1][1] + $cell[$cio][2][1])/2;
  $size[$cio][2] = ($cell[$cio][0][2] + $cell[$cio][1][2] + $cell[$cio][2][2])/2;
  &calc_minmaxpos($cio);
  return 0;
}

sub get_atom_list {
  # !!! does not handle ring systems well !!!
  # recursive function that returns a list of atom IDs connected to an atom
  # in the direction of a certain bond
  # call examples: obtain the chain beginning at at1 in direction of at2 (at1 not incl.)
  # @chain = get_atom_list($fi, $t, $at1, $at2);
  
  my $fi    = $_[0];
  my $t     = $_[1];
  my $at1   = $_[2];
  my $at2   = $_[3];
  my @chain = @_[4..$#_];
  my($i);
  push(@chain,$at2);
  foreach $i (@{$mol_bondatoms[$fi][$t][$at2]}) {
    next if($i==$at1);
    next if(contains(@chain,$i));
    @chain = get_atom_list($fi,$t,$at1,$i,@chain);
  }
  return @chain
}

sub rotate_dihedral {
  # rotate_dihedral($fi, $ci, $t, $m, $at1, $at2, $angle);
  # rotates the part of the molecule bonded to at2 around the vector (at1 -> at2)
  my $fi    = $_[0];
  my $ci    = $_[1];
  my $t     = $_[2];
  my $m     = $_[3];
  my $at1   = $_[4];
  my $at2   = $_[5];
  my $angle = $_[6];
  my @chain = get_atom_list($fi,$t,$at1,$at2);
  my(@axis,@rotmatrix,$a,@new,$i);
  $axis[0] = $cdata[$ci][$t][$m][$at2][0]-$cdata[$ci][$t][$m][$at1][0];
  $axis[1] = $cdata[$ci][$t][$m][$at2][1]-$cdata[$ci][$t][$m][$at1][1];
  $axis[2] = $cdata[$ci][$t][$m][$at2][2]-$cdata[$ci][$t][$m][$at1][2];
  @axis = normalize_vector(\@axis);
  @rotmatrix = gen_rot_matrix(\@axis,$angle);
  for($i=1;$i<@chain;$i++) {
    $a=$chain[$i];
    $cdata[$ci][$t][$m][$a][0] -= $cdata[$ci][$t][$m][$at2][0];
    $cdata[$ci][$t][$m][$a][1] -= $cdata[$ci][$t][$m][$at2][1];
    $cdata[$ci][$t][$m][$a][2] -= $cdata[$ci][$t][$m][$at2][2];
    $new[0] = $cdata[$ci][$t][$m][$a][0]*$rotmatrix[0][0]+$cdata[$ci][$t][$m][$a][1]*$rotmatrix[0][1]+$cdata[$ci][$t][$m][$a][2]*$rotmatrix[0][2];
    $new[1] = $cdata[$ci][$t][$m][$a][0]*$rotmatrix[1][0]+$cdata[$ci][$t][$m][$a][1]*$rotmatrix[1][1]+$cdata[$ci][$t][$m][$a][2]*$rotmatrix[1][2];
    $new[2] = $cdata[$ci][$t][$m][$a][0]*$rotmatrix[2][0]+$cdata[$ci][$t][$m][$a][1]*$rotmatrix[2][1]+$cdata[$ci][$t][$m][$a][2]*$rotmatrix[2][2];
    $cdata[$ci][$t][$m][$a][0] = $new[0];
    $cdata[$ci][$t][$m][$a][1] = $new[1];
    $cdata[$ci][$t][$m][$a][2] = $new[2];
    $cdata[$ci][$t][$m][$a][0] += $cdata[$ci][$t][$m][$at2][0];
    $cdata[$ci][$t][$m][$a][1] += $cdata[$ci][$t][$m][$at2][1];
    $cdata[$ci][$t][$m][$a][2] += $cdata[$ci][$t][$m][$at2][2];
  }
}

sub calc_dvec_orthocell {
  # calc_dvec_orthocell($ci1, $t1, $m1, $a1, $ci2, $t2, $m2, $a2)
  my (@d);
  $d[0] = $cdata[$_[0]][$_[1]][$_[2]][$_[3]][0] - $cdata[$_[4]][$_[5]][$_[6]][$_[7]][0];
  $d[1] = $cdata[$_[0]][$_[1]][$_[2]][$_[3]][1] - $cdata[$_[4]][$_[5]][$_[6]][$_[7]][1];
  $d[2] = $cdata[$_[0]][$_[1]][$_[2]][$_[3]][2] - $cdata[$_[4]][$_[5]][$_[6]][$_[7]][2];
  if($d[0] < -$size[$_[0]][0]) {
    $d[0] += 2*$size[$_[0]][0];
  } elsif($d[0] > $size[$_[0]][0]) {
    $d[0] -= 2*$size[$_[0]][0];
  }
  if($d[1] < -$size[$_[0]][1]) {
    $d[1] += 2*$size[$_[0]][1];
  } elsif($d[1] > $size[$_[0]][1]) {
    $d[1] -= 2*$size[$_[0]][1];
  }
  if($periodic_key[$_[0]]!=6) {
    if($d[2] < -$size[$_[0]][2]) {
      $d[2] += 2*$size[$_[0]][2];
    } elsif($d[2] > $size[$_[0]][2]) {
      $d[2] -= 2*$size[$_[0]][2];
    }
  }
  return @d;
}

sub calc_dsq {
  # calc_dsq($ci1, $t1, $m1, $a1, $ci2, $t2, $m2, $a2)
  my ($dx,$dy,$dz);
  $dx = $cdata[$_[0]][$_[1]][$_[2]][$_[3]][0] - $cdata[$_[4]][$_[5]][$_[6]][$_[7]][0];
  $dy = $cdata[$_[0]][$_[1]][$_[2]][$_[3]][1] - $cdata[$_[4]][$_[5]][$_[6]][$_[7]][1];
  $dz = $cdata[$_[0]][$_[1]][$_[2]][$_[3]][2] - $cdata[$_[4]][$_[5]][$_[6]][$_[7]][2];
  return $dx*$dx + $dy*$dy + $dz*$dz;
}

sub calc_dsq_orthocell {
  # calc_dsq_orthocell($ci, $t1, $m1, $a1, $t2, $m2, $a2)
  my ($dx,$dy,$dz);
  $dx = $cdata[$_[0]][$_[1]][$_[2]][$_[3]][0] - $cdata[$_[0]][$_[4]][$_[5]][$_[6]][0];
  $dy = $cdata[$_[0]][$_[1]][$_[2]][$_[3]][1] - $cdata[$_[0]][$_[4]][$_[5]][$_[6]][1];
  $dz = $cdata[$_[0]][$_[1]][$_[2]][$_[3]][2] - $cdata[$_[0]][$_[4]][$_[5]][$_[6]][2];
  if($periodic_key[$_[0]]!=0) {
    if($dx < -$size[$_[0]][0]) {
      $dx += 2*$size[$_[0]][0];
    } elsif($dx > $size[$_[0]][0]) {
      $dx -= 2*$size[$_[0]][0];
    }
    if($dy < -$size[$_[0]][1]) {
      $dy += 2*$size[$_[0]][1];
    } elsif($dy > $size[$_[0]][1]) {
      $dy -= 2*$size[$_[0]][1];
    }
    if($periodic_key[$_[0]]!=6) {
      if($dz < -$size[$_[0]][2]) {
	$dz += 2*$size[$_[0]][2];
      } elsif($dz > $size[$_[0]][2]) {
	$dz -= 2*$size[$_[0]][2];
      }
    }
  }
  return $dx*$dx + $dy*$dy + $dz*$dz;
}

sub calc_dvec_orthocell_vec {
  # calc_dsq_orthocell_vec($ci, \@vec1, \@vec2);
  return undef if(not @{$size[$_[0]]});
  my @d;
  for(my $c=0;$c<@{$_[1]};$c++) {
    $d[$c]=$_[1][$c]-$_[2][$c];
    if($periodic_key[$_[0]]!=0 and ($c!=2 or $periodic_key[$_[0]]!=6)) {
      if($d[$c] < -$size[$_[0]][$c]) {
        $d[$c] += 2*$size[$_[0]][$c];
      } elsif($d[$c] > $size[$_[0]][$c]) {
        $d[$c] -= 2*$size[$_[0]][$c];
      }
    }
  }
  return @d;
}

sub calc_center_of_mass_orthocell {
  # calc_center_of_mass_orthocell($ci, $fi, $t, $m [@atids])
  # !!! for use with remapped molecules only !!!
  my ($dx,$dy,$dz,@com,$t,$m,$a,$ci,$fi,$totmass);
  $ci      = $_[0];
  $fi      = $_[1];
  $t       = $_[2];
  $m       = $_[3];
  @com=(0,0,0);
  if(@_>4) {
    $totmass=0;
    foreach $a (@_[4..$#_]) {
      $com[0] += $cdata[$ci][$t][$m][$a][0]*$mol_atomdata[$fi][$t][$a][1];
      $com[1] += $cdata[$ci][$t][$m][$a][1]*$mol_atomdata[$fi][$t][$a][1];
      $com[2] += $cdata[$ci][$t][$m][$a][2]*$mol_atomdata[$fi][$t][$a][1];
      $totmass += $mol_atomdata[$fi][$t][$a][1];
    }
    $com[0] /= $totmass;
    $com[1] /= $totmass;
    $com[2] /= $totmass;
  } else {
    for($a=0;$a<$mol_numatoms[$ci][$t];$a++) {
      $com[0] += $cdata[$ci][$t][$m][$a][0]*$mol_atomdata[$fi][$t][$a][1];
      $com[1] += $cdata[$ci][$t][$m][$a][1]*$mol_atomdata[$fi][$t][$a][1];
      $com[2] += $cdata[$ci][$t][$m][$a][2]*$mol_atomdata[$fi][$t][$a][1];
    }
    $com[0]/=$mol_mass[$fi][$t];
    $com[1]/=$mol_mass[$fi][$t];
    $com[2]/=$mol_mass[$fi][$t];
  }
  if($periodic_key[$ci]!=0) {
    if($com[0]<-$size[$ci][0]) {
      $com[0]+=2*$size[$ci][0];
    } elsif($com[0]>$size[$ci][0]) {
      $com[0]-=2*$size[$ci][0];
    }
    if($com[1]<-$size[$ci][1]) {
      $com[1]+=2*$size[$ci][1];
    } elsif($com[1]>$size[$ci][1]) {
      $com[1]-=2*$size[$ci][1];
    }
    if($periodic_key[$ci]!=6) {
      if($com[2]<-$size[$ci][2]) {
	$com[2]+=2*$size[$ci][2];
      } elsif($com[2]>$size[$ci][2]) {
	$com[2]-=2*$size[$ci][2];
      }
    }
  }
  return @com;
}

sub calc_dipole_moment {
  # !!!! caution: molecules have to be remapped !!!!
  my ($t,$m,$a,$ci,$fi,@vec);
  $ci      = $_[0];
  $fi      = $_[1];
  $t       = $_[2];
  $m       = $_[3];
  @vec=(0,0,0);
  for($a=0;$a<$mol_numatoms[$ci][$t];$a++) {
    $vec[0]+=$cdata[$ci][$t][$m][$a][0]*$mol_atomdata[$fi][$t][$a][2];
    $vec[1]+=$cdata[$ci][$t][$m][$a][1]*$mol_atomdata[$fi][$t][$a][2];
    $vec[2]+=$cdata[$ci][$t][$m][$a][2]*$mol_atomdata[$fi][$t][$a][2];
  }
  return @vec;
}

sub cut_octahedron {
  my ($t,$m,$a,$i,$x,$y,$z,@vec,@d,$dist,$buffer,$ci,$fi,$lface);
  $ci      = $_[0];
  $fi      = $_[1];
  $lface   = $_[2]; # center-face distance = sqrt(6)/6 a
  $buffer  = $_[3]; # optional buffer zone to avoid overlap of molecules at boundaries
  $buffer  = 0 if(not defined($buffer));
  $x=1.0/sqrt(3);
  @{$vec[0]} = ( $x, $x, $x);
  @{$vec[1]} = ( $x, $x,-$x);
  @{$vec[2]} = ( $x,-$x, $x);
  @{$vec[3]} = ( $x,-$x,-$x);
  @{$vec[4]} = (-$x, $x, $x);
  @{$vec[5]} = (-$x, $x,-$x);
  @{$vec[6]} = (-$x,-$x, $x);
  @{$vec[7]} = (-$x,-$x,-$x);
  @d=($lface,$lface,$lface,$lface,$lface,$lface,$lface,$lface);
  for($t=0;$t<$field_nummols[$ci];$t++) {
    loopmol:for($m=0;$m<$mol_numents[$ci][$t];$m++) {
      for($a=0;$a<$mol_numatoms[$ci][$t];$a++) {
	for($i=0;$i<@vec;$i++) {
	  $dist  = $buffer-$d[$i];
	  $dist += $vec[$i][0]*$cdata[$ci][$t][$m][$a][0];
	  $dist += $vec[$i][1]*$cdata[$ci][$t][$m][$a][1];
	  $dist += $vec[$i][2]*$cdata[$ci][$t][$m][$a][2];
	  if($dist>0) {
	    remove_mol_entity($ci,$fi,$t,$m);
	    $m--;
	    next loopmol;
	  }
	}
      }
    }
  }
}

sub cut_truncated_octahedron {
  my ($t,$m,$a,$i,$x,$y,$z,@vec,@d,$dist,$buffer,$ci,$fi,$lsquare,$lhex);
  $ci      = $_[0];
  $fi      = $_[1];
  $lsquare = $_[2]; # center-square distance = sqrt(2)*/3 a
  $buffer  = $_[3]; # optional buffer zone to avoid overlap of molecules at boundaries
  $buffer  = 0 if(not defined($buffer));
  $lhex    = $lsquare*sqrt(3)/2.0; # center-hexagon distance = sqrt(6)/6 a
  @{$vec[0]}   = (0,0,1);
  @{$vec[1]}   = (0,0,-1);
  @{$vec[2]}   = (0,1,0);
  @{$vec[3]}   = (0,-1,0);
  @{$vec[4]}   = (1,0,0);
  @{$vec[5]}   = (-1,0,0);
  @d = ($lsquare,$lsquare,$lsquare,$lsquare,$lsquare,$lsquare);
  $x=1/sqrt(3);
  @{$vec[6]}   = ( $x, $x, $x);
  @{$vec[7]}   = ( $x, $x,-$x);
  @{$vec[8]}   = ( $x,-$x, $x);
  @{$vec[9]}   = ( $x,-$x,-$x);
  @{$vec[10]}  = (-$x, $x, $x);
  @{$vec[11]}  = (-$x, $x,-$x);
  @{$vec[12]}  = (-$x,-$x, $x);
  @{$vec[13]}  = (-$x,-$x,-$x);
  push(@d,$lhex,$lhex,$lhex,$lhex,$lhex,$lhex,$lhex,$lhex);
  for($t=0;$t<$field_nummols[$ci];$t++) {
    loopmol:for($m=0;$m<$mol_numents[$ci][$t];$m++) {
      for($a=0;$a<$mol_numatoms[$ci][$t];$a++) {
	for($i=0;$i<@vec;$i++) {
	  $dist  = $buffer-$d[$i];
	  $dist += $vec[$i][0]*$cdata[$ci][$t][$m][$a][0];
	  $dist += $vec[$i][1]*$cdata[$ci][$t][$m][$a][1];
	  $dist += $vec[$i][2]*$cdata[$ci][$t][$m][$a][2];
	  if($dist>0) {
	    remove_mol_entity($ci,$fi,$t,$m);
	    $m--;
	    next loopmol;
	  }
	}
      }
    }
  }
  @{$cell[$ci]}=([2.0*$lsquare,0,0],[0,2.0*$lsquare,0],[0,0,2.0*$lsquare]);
  $periodic_key[$ci] = 4;
}

sub calc_cell_abc {
  my $ci=$_[0];
  my(@celldata);
  $celldata[0] = vector_length(@{$cell[$ci][0]}); # a
  $celldata[1] = vector_length(@{$cell[$ci][1]}); # b
  $celldata[2] = vector_length(@{$cell[$ci][2]}); # c
  $celldata[3] = 180.0/pi*acos(($cell[$ci][0][0]*$cell[$ci][2][0]+$cell[$ci][0][1]*$cell[$ci][2][1]+$cell[$ci][0][2]*$cell[$ci][2][2])/($celldata[0]*$celldata[2])); # alpha
  $celldata[4] = 180.0/pi*acos(($cell[$ci][1][0]*$cell[$ci][2][0]+$cell[$ci][1][1]*$cell[$ci][2][1]+$cell[$ci][1][2]*$cell[$ci][2][2])/($celldata[1]*$celldata[2])); # beta
  $celldata[5] = 180.0/pi*acos(($cell[$ci][0][0]*$cell[$ci][1][0]+$cell[$ci][0][1]*$cell[$ci][1][1]+$cell[$ci][0][2]*$cell[$ci][1][2])/($celldata[0]*$celldata[1])); # gamma
  return @celldata
}

sub calc_cell_vecs {
  # input: calc_cell_vecs(a,b,c [alpha, beta, gamma])
  # output: 2d array with vectors
  my @res=();
  $res[0][0]=$_[0];
  $res[0][1]=0;
  $res[0][2]=0;
  $res[1][0]=$_[1]*cos($_[3]);
  $res[1][1]=$_[1]*sin($_[3]);
  $res[1][2]=0;
  $res[2][0]=$_[2]*cos($_[4]);
  $res[2][1]=($_[1]*$_[2]*cos($_[5])-$res[1][0]*$res[2][0])/$res[1][1];
  $res[2][2]=sqrt(abs($_[2]*$_[2]-$res[2][0]*$res[2][0]-$res[2][1]*$res[2][1]));
  $res[1][0]=0 if(abs($res[1][0])<1e-14);
  $res[2][0]=0 if(abs($res[2][0])<1e-14);
  $res[2][1]=0 if(abs($res[2][1])<1e-14);
  return @res;
}

sub print_statis_data_header {
  my $fh = $_[0]; # file handle
  my $fi = $_[1]; # field file
  my @stress=("x","y","z");
  my($i,$j,$type);
  printf $fh "%14s","# 1";
  $j=28+@{$field_atomtypes[$fi]}+9+1;
  for($i=2;$i<=$j;$i++) {
    printf $fh "%14u", $i;
  }
  print  $fh "\n";
  printf $fh "%14s","time";
  
  printf $fh "%14s","E_tot";
  printf $fh "%14s","T_tot";
  printf $fh "%14s","E_cfg";
  printf $fh "%14s","E_vdw,metal";
  printf $fh "%14s","E_coul";

  printf $fh "%14s","E_bnd";
  printf $fh "%14s","E_ang,tbp";
  printf $fh "%14s","E_dih,fbp";
  printf $fh "%14s","E_tether";
  printf $fh "%14s","H=E_tot+P*V";

  printf $fh "%14s","T_rot";
  printf $fh "%14s","vir_tot";
  printf $fh "%14s","vir_vdw,metal";
  printf $fh "%14s","vir_coul";
  printf $fh "%14s","vir_bnd";

  printf $fh "%14s","vir_ang";
  printf $fh "%14s","vir_constr";
  printf $fh "%14s","vir_teth";
  printf $fh "%14s","volume";
  printf $fh "%14s","T_core-shell";

  printf $fh "%14s","E_pot,core-sh";
  printf $fh "%14s","vir_core-shel";
  printf $fh "%14s","cell_alpha";
  printf $fh "%14s","cell_beta";
  printf $fh "%14s","cell_gamma";

  printf $fh "%14s","vir_pmf";
  printf $fh "%14s","pressure";
  foreach $type (@{$field_atomtypes[$fi]}) {
    printf $fh "%14s","msd($type)";
  }
  for($i=0;$i<3;$i++) {
    for($j=0;$j<3;$j++) {
      printf $fh "%14s","stress($stress[$i],$stress[$j])";
    }  
  }

  printf $fh "%14s","E_extfld";
}

sub rotate_cell_vmd {
  my $ci         = $_[0];
  my $lomitcoord = $_[1]; # set to 1 to omit rotating the coordinates
  return 1 if(not @{$cell[0]});
  return 2 if($periodic_key[$ci] == 0);
  my($angle,@axis,@mat1,@mat2,$t,$m);
  # align first vector to x-axis
  $angle=acos($cell[$ci][0][0]/vector_length(@{$cell[$ci][0]}));
  if(abs($angle>1.0e-10)) {
    @axis=(1,0,0);
    @axis=vector_product(\@axis,\@{$cell[$ci][0]});
    @axis=normalize_vector(\@axis);
# 	print "axis (",join(" ",@axis),")\n";
    @mat1 = gen_rot_matrix(\@axis,-$angle);
    @{$cell[$ci][0]}=rotate_vector($cell[$ci][0],\@mat1);
    @{$cell[$ci][1]}=rotate_vector($cell[$ci][1],\@mat1);
    @{$cell[$ci][2]}=rotate_vector($cell[$ci][2],\@mat1);
  }
  # bring second vector into xy plane with y positive
  $angle=acos($cell[$ci][1][1]/vector_length($cell[$ci][1][1],$cell[$ci][1][2]));
  if(abs($angle>1.0e-10)) {
    @axis=(1,0,0);
    @mat2 = gen_rot_matrix(\@axis,-$angle);
    @{$cell[$ci][0]}=rotate_vector($cell[$ci][0],\@mat2);
    @{$cell[$ci][1]}=rotate_vector($cell[$ci][1],\@mat2);
    @{$cell[$ci][2]}=rotate_vector($cell[$ci][2],\@mat2);
  }
  # clear any negligibly small vector components
  for($t=0;$t<3;$t++) {
    for($m=0;$m<3;$m++) {
      $cell[$ci][$t][$m]=0 if(abs($cell[$ci][$t][$m])<1.0e-10);
    }
  }
  # invert last vector if the z-component is negative
  if($cell[$ci][2][2]<0) {
    $cell[$ci][2][0]=-$cell[$ci][2][0];
    $cell[$ci][2][1]=-$cell[$ci][2][1];
    $cell[$ci][2][2]=-$cell[$ci][2][2];
  }
  unless($lomitcoord) {
    @mat1 = matmul(\@mat2,\@mat1);
    for($t=0;$t<@{$cdata[$ci]};$t++) {
      for($m=0;$m<@{$cdata[$ci][$t]};$m++) {
	rotate_molecule(\@{$cdata[$ci][$t][$m]}, \@mat1);
      }
    }
  }
  return 0;
}

sub print_statis_data {
  my $fh       = $_[0]; # file handle
  my $si       = $_[1]; # sdata index
  my $firstrow = $_[2]; # entry for first row, may be time/frame/etc.
  my($i);
  printf $fh "\n%14.7g",$firstrow;
  for($i=2;$i<@{$sdata[$si]};$i++) {
    printf $fh "%14.6e",$sdata[$si][$i];
  }
}

sub read_lammps_thermo {
  my $filename = $_[0];
  my $ti       = $_[1];
  my($filehandle,@linedata,@keys);
  if(not open($filehandle, "<", $filename)) {
    print "**** error: Can't open lammps log file \"$filename\": $!\n";
    return 1;
  }
  undef @{$thermodata[$ti]};
  @{$thermodata[$ti]} = ();
  undef @{$thermokeys[$ti]};
  @{$thermokeys[$ti]} = ();
  my $isthermo=0;
  my $i=0;
  while(<$filehandle>) {
    if($isthermo) {
	$_ =~ s/^\s+//; $_ =~ s/\s+$//;
	@linedata = split(/\s+/, $_);
	if(not check_integer($linedata[0])) {
	  $isthermo=0;
	  $i++;
	  next;
	}
	push(@{$thermodata[$ti][$i]},[@linedata]);
    } elsif(/^Step /) {
	$_ =~ s/^\s+//; $_ =~ s/\s+$//;
	@linedata = split(/\s+/, $_);
	@{$thermodata[$ti][$i]} = ();
	@{$thermokeys[$ti][$i]} = @linedata;
	$isthermo=1;
    }
  }
  close($filehandle);
}

sub print_thermo {
  my $fh = $_[0];
  my $ti = $_[1];
  my($i,$j,$k);
  for(my $i=0;$i<@{$thermokeys[$ti]};$i++) {
    for(my $j=0;$j<@{$thermokeys[$ti][$i]};$j++) {
      printf $fh " %12s",$thermokeys[$ti][$i][$j];
    }
    print "\n";
    for(my $j=0;$j<@{$thermodata[$ti][$i]};$j++) {
      for(my $k=0;$k<@{$thermodata[$ti][$i][$j]};$k++) {
	printf $fh " %12.5g",$thermodata[$ti][$i][$j][$k];
      }
      print "\n";
    }
    print "\n\n";
  }
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
