#!/usr/bin/perl

use hanno_utility;
use dlpoly_utility;
use Switch;
use strict;


if($#ARGV<3) {
  print "input format:\n";
  print "1. name of CONFIG file\n";
  print "2. name of FIELD file\n";
  print "3. name of CONTROL file (may be empty string)\n";
  print "4. name of output lammps files (*.data, *.in)\n";
  print "optional flags:\n";
  print "-d <str> dump as dcd (default), lammpstrj or xyz\n";
  print "-0       ignore pair potentials that are zero\n";
  print "-ii      print only ii pairwise parameters and use mixing\n";
  exit 1;
}

my($i);
my $dump   = "dcd";
my $onlyii = 0;
my $ignorezero = 0;
for($i=4;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-d$/i) {
      $i++;
      $dump=$ARGV[$i];
      if($dump!~/^(custom(\/gz|\/mpiio)?|cfg(\/gz|\/mpiio)?|dcd|xtc|xyz(\/gz|\/mpiio)?)$/) {
        print "**** error: invalid dump style\!\n";
        exit 1;
      }
    } case (/^-ii$/i) {
      $onlyii=1;
    } case (/^-0$/i) {
      $ignorezero=1;
    } else {
      print "**** error: expected optional flag but got \"$ARGV[$i]\"!\n";
      exit 1;
    }
  }
}

exit 1 if(read_field_file($ARGV[1],0)!=0);
exit 1 if(read_config_file($ARGV[0],0,0)!=0);
$ARGV[2]=~s/(^\s+|\s+$)//g;
if(length($ARGV[2]) > 0) {
  exit 1 if(read_control_file($ARGV[2],0)!=0);
  write_lammps_data($ARGV[3],0,0,0,"metal","lammps test",$dump,$onlyii,$ignorezero);
} else {
  write_lammps_data($ARGV[3],0,0,-1,"metal","lammps test",$dump,$onlyii,$ignorezero);
}

sub write_lammps_data {
  my($fh_data,$fh_in,$fh_mols,$t,$m,$k,$l,$i,$j,$indexoffset,$ifix,$fenergy,
  $fvel,$fefield,$group_integrate,$cutglob,$rvdw,$coulstyle,$couladd);
  my $filename  = $_[0];
  my $ci        = $_[1]; # config data
  my $fi        = $_[2]; # field data
  my $oi        = $_[3]; # control data
  my $lunits    = $_[4];
  my $title     = $_[5];
  my $dumpstyle = $_[6];
  my $onlyii    = 0;
  $onlyii=$_[7] if($#_>6);
  my $ignorezero = 0;
  $ignorezero=$_[8] if($#_>6);
  
  my $numatoms=0;
  my $numbonds=0;
  my $numangles=0;
  my $numdihedrals=0;
  my $numimpropers=0;
  my @atomtypes=();     # 1=name      2=mass
  my @bondcoeff=();     # same as mol_bonddata
  my @anglecoeff=();    # same as mol_angledata
  my @dihedralcoeff=(); # same as mol_dihedraldata
  my @paircoeff=();     # same as field_vdwdata except that the atom names are lammps atomtype indices
  my @frozenatoms=();
  my @frozenatomsmol=();
  my @bondstyles=();
  my @anglestyles=();
  my @dihedralstyles=();
  my $specialbonds=undef;
  my @pairstyles=();
  my @impropercoeff=(); 
  my @mol_atomtype=();
  my @mol_bondcoeff=();
  my @mol_anglecoeff=();
  my @mol_dihedralcoeff=();
  my @usedtypes=();
  my @molgroups=();
  my $ev = 1.60217657e-19; # j
  # determine unit conversion factors
  my $ftime   = 1.0e-12;
  my $fmass   = 1.6605402e-27;
  my $flength = 1.0e-10;
  my $fpress  = 101325000.0;
  my $fcharge = 1.60217733e-19;
  my $fdipole = 1.0;
  if($field_units[$fi]=~/^dl/i) {
    $fenergy *= 1.6605402e-23;
  } elsif($field_units[$fi]=~/^ev/i) {
    $fenergy = $fcharge;
  } elsif($field_units[$fi]=~/^kj/i) {
    $fenergy = 1.6605402e-21;
  } elsif($field_units[$fi]=~/^kcal/i) {
    $fenergy = 6.94769484559868e-21;
  }
  if($lunits=~/^si/i) {
  } elsif($lunits=~/^real/i) {
    $fmass      /= 1.6605402e-27; # g/mole;
    $ftime      /= 1.0e-15; # fs;
    $fenergy    /= 6.94769484559868e-21; #kcal/mol
    $fcharge    /= 1.60217733e-19; # e
    $flength    /= 1.0e-10; # A
    $fpress     /= 101325.0; # atm
  } elsif($lunits=~/^metal/i) {
    $fmass   /= 1.6605402e-27; # g/mole;
    $ftime   /= 1.0e-12; # ps;
    $fenergy /= 1.60217733e-19; # eV
    $fcharge /= 1.60217733e-19; # e
    $flength /= 1.0e-10; # A
    $fpress  /= 100000.0; # bar
  } elsif($lunits=~/^electron/i) {
    $fmass   /= 1.6605402e-27; # atomic mass unit
    $ftime   /= 1.0e-15; # fs;
    $fenergy /= 4.35974434e-18; # Ha
    $fcharge /= 1.60217733e-19; # e
    $flength /= 0.52917721092e-10; # bohr
  } elsif($lunits=~/^micro/i) {
    $fmass   /= 1.0e-12; # picogram;
    $ftime   /= 1.0e-6; # microsec;
    $fenergy /= 1.0e-12; # picogram*micrometer**2/microsec**2
    $fcharge /= 1.0e-12; # picocoulomb
    $flength /= 1.0e-6; # micrometer
    $fpress  /= 1.0e3; # picogram/(micrometer*microsec**2)
  } elsif($lunits=~/^nano/i) {
    $fmass   /= 1.0e-18; # attogram;
    $ftime   /= 1.0e-9; # nanosec;
    $fenergy /= 1.0e-18; # attogram*nanometer**2/nanosec**2
    $fcharge /= 1.0e-12; # picocoulomb
    $flength /= 1.0e-9; # nanometer
    $fpress  /= 1.0e9; # attogram/(nanometer*nanosec**2)
  } elsif($lunits=~/^cgs/i) {
    print "**** error: Are you serious? Who the hell uses CGS units?!\n";
    print "            I'm not gonna do this, go fuck yourself!\n";
    return 1;
  } else {
    print "**** error: unrecognized unit system $lunits!\n"; return 1;
  }
  $fefield = 1.03642691902687e+18;
  $fvel    = $flength/$ftime;
  if($lunits=~/^electron/i) {
    $fvel     = $flength/(1.0e-12/1.03275e-15); # bohr/atomic time units
    $fefield /= 1.0e-2; # V/cm
  } else {
    $fefield /= $flength;
  }
  
  if($oi<0 or $control_cutglob[$oi] == undef) {
    $cutglob    = 12;
  } else {
    $cutglob    = $control_cutglob[$oi];
  }
  if($oi<0 or $control_rvdw[$oi] == undef) {
    $rvdw       = $cutglob;
  } else {
    $rvdw       = $control_rvdw[$oi];
  }
  if($oi<0) {
    $coulstyle = "coul/dsf";
    $couladd   = "0 $cutglob";
  } elsif(not defined($control_coulmethod[$oi])) {
    $coulstyle = undef;
  } elsif($control_coulmethod[$oi]=~/^coul/) {
    $coulstyle = "coul/cut";
    $couladd   = "$cutglob";
  } elsif($control_coulmethod[$oi]=~/^shift/) {
    if($control_coulmethod[$oi]=~/^damp/) { # damped shifted force
      $coulstyle = "coul/dsf";
      $couladd   = "$control_coulprm[0][0] $cutglob";
    } elsif($control_coulmethod[$oi]=~/prec/) { # damped shifted force w/ automatic parameter optim
      $i = min(0.5,abs($control_coulprm[0][0])); #eps in define_system_module.f
      $j = sqrt(abs(log($i*$cutglob))); #tol
      $k = sqrt(abs(log($i*$cutglob*$j)))/$cutglob;
      $coulstyle = "coul/dsf";
      $couladd   = sprintf("%.5g %s",$k,$cutglob);
    } else { # shifted force
      $coulstyle = "coul/dsf";
      $couladd   = "0 $cutglob";
    }
  } elsif($control_coulmethod[$oi]=~/^(ewald|spme|hke)/) {
    print "**** error: ewald summation is not implemented (yet)!\n"; return 1;
  } else {
    print "**** warning: unknown coulomb method \"$coulstyle.\"\n";
    print "              using shifted coulomb summation instead\n";
    $coulstyle = "coul/dsf";
    $couladd   = "0 $cutglob";
  }
  # check periodic cell
  if($periodic_key[$ci]>0) {
    if($cell[$ci][0][1] != 0 # 1st vector on x-axis
    or $cell[$ci][0][2] != 0 # 1st vector on x-axis
    or $cell[$ci][1][2] != 0 # 2nd vector in xy-plane
    or $cell[$ci][0][0] <= 0 # 1st vector positive x
    or $cell[$ci][1][1] <= 0 # 2nd vector positive y
    or $cell[$ci][2][2] <= 0) { # 3rd vector positive z
      print "**** error: incompatible periodic cell!\n";
      return 1;
    }
  }
  
  if($fi<0 or not check_integer($fi)) {
    print "**** error: invalid FIELD information \"$filename\": $!\n";
    return 1;
  }
  @{$minpos[$ci]}=( 9e20, 9e20, 9e20);
  @{$maxpos[$ci]}=(-9e20,-9e20,-9e20);
  for($t=0;$t<@{$cdata[$ci]};$t++) {
    for($m=0;$m<@{$cdata[$ci][$t]};$m++) {
      for($k=0;$k<@{$cdata[$ci][$t][$m]};$k++) {
	$numatoms++;
	$minpos[$ci][0] = $cdata[$ci][$t][$m][$k][0] if($cdata[$ci][$t][$m][$k][0]<$minpos[$ci][0]);
	$minpos[$ci][1] = $cdata[$ci][$t][$m][$k][1] if($cdata[$ci][$t][$m][$k][1]<$minpos[$ci][1]);
	$minpos[$ci][2] = $cdata[$ci][$t][$m][$k][2] if($cdata[$ci][$t][$m][$k][2]<$minpos[$ci][2]);
	$maxpos[$ci][0] = $cdata[$ci][$t][$m][$k][0] if($cdata[$ci][$t][$m][$k][0]>$maxpos[$ci][0]);
	$maxpos[$ci][1] = $cdata[$ci][$t][$m][$k][1] if($cdata[$ci][$t][$m][$k][1]>$maxpos[$ci][1]);
	$maxpos[$ci][2] = $cdata[$ci][$t][$m][$k][2] if($cdata[$ci][$t][$m][$k][2]>$maxpos[$ci][2]);
      }
    }
  }
  if($numatoms!=$field_numatoms[$fi]) {
    print "**** error: number of atoms in FIELD file ($field_numatoms[$fi]) does not match number of",
    "atoms in CONFIG file ($numatoms)!\n";
    return 1;
  } elsif($numatoms==0) {
    print "**** error: no atoms found in CONFIG file!\n";
    return 1;
  }
  $indexoffset=1;
  for($t=0;$t<$field_nummols[$fi];$t++) {
    next if($mol_numents[$fi][$t]==0);
    if($mol_numinversions[$fi][$t]>0) {
      print "**** error: inversions are not implemented (yet)!\n";
      return 1;
    }
    if($mol_numconstraints[$fi][$t]>0) {
      print "**** error: constraints are not implemented (yet)!\n";
      return 1;
    }
    $numbonds     += $mol_numents[$fi][$t]*$mol_numbonds[$fi][$t];
    $numangles    += $mol_numents[$fi][$t]*$mol_numangles[$fi][$t];
    $numdihedrals += $mol_numents[$fi][$t]*$mol_numdihedrals[$fi][$t];
    # build atom type list and frozen list
    @frozenatomsmol = ();
    for($k=0;$k<$mol_numatoms[$fi][$t];$k++) {
      $i=find_lammps_atomtype($fi,$t,$k,\@atomtypes);
      if($i<0) {
        push(@usedtypes,$mol_atomdata[$fi][$t][$k][0]);
        push(@atomtypes,[@{$mol_atomdata[$fi][$t][$k]}[0..1]]);
        $i=$#atomtypes;
      }
      if($mol_atomdata[$fi][$t][$k][4] == 1) { # frozen
        push(@frozenatomsmol,$k);
      }
      $mol_atomtype[$t][$k] = $i;
    }
    if(@frozenatomsmol==$mol_numatoms[$fi][$t]) {
      push (@frozenatoms,$indexoffset.":".($indexoffset+$mol_numents[$fi][$t]*$mol_numatoms[$fi][$t]-1))
    } else {
      foreach $k (@frozenatomsmol) {
        push (@frozenatoms,($k+$indexoffset).":".
        ($k+$indexoffset+$mol_numents[$fi][$t]*$mol_numatoms[$fi][$t]-1).":".$mol_numatoms[$fi][$t])
      }
    }
    # build bond type list
    for($k=0;$k<$mol_numbonds[$fi][$t];$k++) {
      $i=find_lammps_bondtype($fi,$t,$k,\@bondcoeff);
      if($i<0) {
	push(@bondcoeff,[@{$mol_bonddata[$fi][$t][$k]}]);
	push(@bondstyles,$mol_bonddata[$fi][$t][$k][0])
          if(not contains(@bondstyles,$mol_bonddata[$fi][$t][$k][0]));
	$i=$#bondcoeff;
      }
      $mol_bondcoeff[$t][$k] = $i;
    }
    for($k=0;$k<$mol_numangles[$fi][$t];$k++) {
      $i=find_lammps_angletype($fi,$t,$k,\@anglecoeff);
      if($i<0) {
	push(@anglecoeff,[@{$mol_angledata[$fi][$t][$k]}]);
	push(@anglestyles,$mol_angledata[$fi][$t][$k][0])
          if(not contains(@anglestyles,$mol_angledata[$fi][$t][$k][0]));
	$i=$#anglecoeff;
      }
      $mol_anglecoeff[$t][$k] = $i;
    }
    for($k=0;$k<$mol_numdihedrals[$fi][$t];$k++) {
      $i=find_lammps_dihedraltype($fi,$t,$k,\@dihedralcoeff);
      if($i<0) {
	push(@dihedralcoeff,[@{$mol_dihedraldata[$fi][$t][$k]}]);
	push(@dihedralstyles,$mol_dihedraldata[$fi][$t][$k][0])
          if(not contains(@dihedralstyles,$mol_dihedraldata[$fi][$t][$k][0]));
	$i=$#dihedralcoeff;
      }
      $mol_dihedralcoeff[$t][$k] = $i;
      # check 1-4 scaling factors:
      if($mol_dihedraldata[$fi][$t][$k][0]=~/^(ryck|rbf)/i) {
        $i=6;
      } elsif($mol_dihedraldata[$fi][$t][$k][0]=~/^(harm|hcos)/i) {
        $i=7;
      } elsif($mol_dihedraldata[$fi][$t][$k][0]=~/^cos/i) {
        $i=8;
      } elsif($mol_dihedraldata[$fi][$t][$k][0]=~/^opls/i) {
        $i=9;
      }
      if(abs($mol_dihedraldata[$fi][$t][$k][$i]-0.833)<0.01 and $mol_dihedraldata[$fi][$t][$k][$i+1]==0.5) {
        $j="amber";
      } elsif($mol_dihedraldata[$fi][$t][$k][$i]==0 and $mol_dihedraldata[$fi][$t][$k][$i+1]==0) {
        $j="charmm";
      } elsif($mol_dihedraldata[$fi][$t][$k][$i]==1 and $mol_dihedraldata[$fi][$t][$k][$i+1]==1) {
        $j="dreiding";
      } else {
        print "**** warning: unknown 1-4 scaling factors. Output for 'special_bonds' must be revised!\n";
        $j="NEEDS REVISION";
      }
      if(not defined($specialbonds)) {
        $specialbonds=$j;
      } elsif($specialbonds!~/^NEEDS/ and $specialbonds ne $j) {
        print "**** warning: mixed 1-4 scaling factors. Output for 'special_bonds' must be revised!\n";
        $specialbonds="NEEDS REVISION";
      }
    }
    $i=$indexoffset+$mol_numatoms[$fi][$t]*$mol_numents[$fi][$t]-1;
    push(@molgroups,"group $mol_name[$fi][$t] id $indexoffset:$i");
    $indexoffset=$i+1;
  }
  
  # check vdw pairs
  for($i=0;$i<@{$field_vdwdata[$fi]};$i++) {
    next if($onlyii and $field_vdwdata[$fi][$i][0] ne $field_vdwdata[$fi][$i][1]);
    next if(not contains(@usedtypes,$field_vdwdata[$fi][$i][0]) or not contains(@usedtypes,$field_vdwdata[$fi][$i][1]));
    if($ignorezero and $field_vdwdata[$fi][$i][2]=~/lj/i and $field_vdwdata[$fi][$i][3]==0) {
      print "ignoring lj-potential $field_vdwdata[$fi][$i][0]-$field_vdwdata[$fi][$i][1] with epsilon=0\n";
      next;
    }
    if(find_dlpoly_pairtype($fi,$i)>1) {
      print "**** error: more than one nonbonded parameter for pair ".
      $field_vdwdata[$fi][$i][0]." and ".$field_vdwdata[$fi][$i][1]."!\n";
      return 1;
    }
    if($field_vdwdata[$fi][$i][2]=~/^lj/i or $i=~/^12-6/) {
      $j="lj/cut";
    } elsif($field_vdwdata[$fi][$i][2]=~/^buck/) {
      $j="buck";
    } elsif($field_vdwdata[$fi][$i][2]=~/^mors/) {
      $j="morse";
    } else {
      print "**** error: unknown VDW type $i\n"; return 1;
    }
    push(@pairstyles,$j) if(not contains(@pairstyles,$j));
    foreach $k (find_lammps_atomtype_pair($field_vdwdata[$fi][$i][0],\@atomtypes)) {
      foreach $l (find_lammps_atomtype_pair($field_vdwdata[$fi][$i][1],\@atomtypes)) {
        $m=$#paircoeff+1;
	@{$paircoeff[$m]} = sort {$a <=> $b} ($k+1,$l+1);
	$paircoeff[$m][2] = $j;
	if($field_vdwdata[$fi][$i][2] =~ /^lj/i) {
          $paircoeff[$m][3] = $field_vdwdata[$fi][$i][3]*$fenergy;
          $paircoeff[$m][4] = $field_vdwdata[$fi][$i][4]*$flength;
        } elsif($field_vdwdata[$fi][$i][2] =~ /^12-6/) {
          $paircoeff[$m][3] = (($field_vdwdata[$fi][$i][3]-$field_vdwdata[$fi][$i][4])/
            (4.0*(($field_vdwdata[$fi][$i][3]/$field_vdwdata[$fi][$i][4])**2
            -($field_vdwdata[$fi][$i][3]/$field_vdwdata[$fi][$i][4]))))*$fenergy; #epsilon
          $paircoeff[$m][4] = (($field_vdwdata[$fi][$i][3]/($field_vdwdata[$fi][$i][4])**(1.0/6.0)))*$flength; #sigma
        } elsif($field_vdwdata[$fi][$i][2] =~ /^buck/) {
          $paircoeff[$m][3] = $paircoeff[$i][3]*$fenergy;
          $paircoeff[$m][4] = $paircoeff[$i][4]*$flength;
          $paircoeff[$m][5] = $paircoeff[$i][5]*$fenergy;
        } elsif($field_vdwdata[$fi][$i][2] =~ /^mors/) {
          $paircoeff[$m][3] = $field_vdwdata[$fi][$i][3]*$fenergy;
          $paircoeff[$m][4] = $field_vdwdata[$fi][$i][5];
          $paircoeff[$m][5] = $field_vdwdata[$fi][$i][4]*$flength;
        }
      }
    }
  }
  @paircoeff=sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @paircoeff;
    
  ############ generate data file #############################################
  
  if(not open($fh_data, ">", "$filename.data")) {
    print "**** error: Can't open output file $filename.data: $!\n";
    return 1;
  }
  
  $title =~ s/^\s+//; $title =~ s/\s+$//;
  print $fh_data "$title\n\n";
  printf $fh_data "%8u atoms\n", $numatoms;
  printf $fh_data "%8u bonds\n", $numbonds if($numbonds>0);
  printf $fh_data "%8u angles\n", $numangles if($numangles>0);
  printf $fh_data "%8u dihedrals\n", $numdihedrals if($numdihedrals>0);
  
  printf $fh_data "\n%8u atom types\n", $#atomtypes+1;
  printf $fh_data "%8u bond types\n", $#bondcoeff+1 if(@bondcoeff>0);
  printf $fh_data "%8u angle types\n", $#anglecoeff+1 if(@anglecoeff>0);
  printf $fh_data "%8u dihedral types\n", $#dihedralcoeff+1 if(@dihedralcoeff>0);
  
  print $fh_data "\n";
  if($periodic_key[$ci]==0) { # non-periodic cell
      printf $fh_data "%12.6f  %12.6f  xlo xhi\n", $minpos[$ci][0]-50, $maxpos[$ci][0]+50;
      printf $fh_data "%12.6f  %12.6f  ylo yhi\n", $minpos[$ci][1]-50, $maxpos[$ci][1]+50;
      printf $fh_data "%12.6f  %12.6f  zlo zhi\n", $minpos[$ci][2]-50, $maxpos[$ci][2]+50;
  } elsif(check_orthogonal($ci)) { # orthogonal cell
    printf $fh_data "\n%12.6f  %12.6f  xlo xhi\n", -$cell[$ci][0][0]/2.0, $cell[$ci][0][0]/2.0;
    printf $fh_data   "%12.6f  %12.6f  ylo yhi\n", -$cell[$ci][1][1]/2.0, $cell[$ci][1][1]/2.0;
    if($periodic_key[$ci]!=6) {
      printf $fh_data "%12.6f  %12.6f  zlo zhi\n", -$cell[$ci][2][2]/2.0, $cell[$ci][2][2]/2.0;
    } else {
      printf $fh_data "%12.6f  %12.6f  zlo zhi\n", $minpos[$ci][2]-20, $maxpos[$ci][2]+20;
    }
  } else { # triclinic cell
      $cell[$ci][1][0],$cell[$ci][2][0],$cell[$ci][2][1];
    my $xlo = -($cell[$ci][0][0]+$cell[$ci][1][0]+$cell[$ci][2][0])/2.0;
    my $ylo = -($cell[$ci][1][1]+$cell[$ci][2][1])/2.0;
    printf $fh_data "%12.6f  %12.6f  xlo xhi\n", $xlo, $xlo + $cell[$ci][0][0];
    printf $fh_data "%12.6f  %12.6f  ylo yhi\n", $ylo, $ylo + $cell[$ci][1][1];
    printf $fh_data "%12.6f  %12.6f  zlo zhi\n", -$cell[$ci][2][2]/2.0, $cell[$ci][2][2]/2.0;
    printf $fh_data "%12.6f  %12.6f  %12.6f xy xz yz\n",
    $cell[$ci][1][0],$cell[$ci][2][0],$cell[$ci][2][1];
  }
  
  print $fh_data "\nMasses\n\n";
  for($i=0;$i<@atomtypes;$i++) {
    printf $fh_data "%8u  %10.4f  # %-8s\n", $i+1,$atomtypes[$i][1]*$fmass,$atomtypes[$i][0];
  }
  if(@bondcoeff) {
    print $fh_data "\nBond Coeffs\n\n";
    for($i=0;$i<@bondcoeff;$i++) {
      if($bondcoeff[$i][0] =~ /^harm/) {
	printf $fh_data "%8u", $i+1;
	printf $fh_data " %-15s", "harmonic" if(@bondstyles>1);
	printf $fh_data " %15.6f %15.6f\n" ,0.5*$fenergy*$bondcoeff[$i][3],$flength*$bondcoeff[$i][4];
      } elsif($bondcoeff[$i][0] =~ /^quar/) {
	printf $fh_data "%8u", $i+1;
	printf $fh_data " %-15s", "class2" if(@bondstyles>1);
	printf $fh_data " %15.6f %15.6f %15.6f %15.6f\n",$flength*$bondcoeff[$i][4],
	0.5*$fenergy*$bondcoeff[$i][3],$bondcoeff[$i][5]*$fenergy/3.0,0.25*$fenergy*$bondcoeff[$i][6];
      } else {
	print "**** error: bond potential \"$bondcoeff[$i][0]\" is not implemented (yet)!\n";
	return 1;
      }
    }
  }
  if(@anglecoeff) {
    print $fh_data "\nAngle Coeffs\n\n";
    for($i=0;$i<@anglecoeff;$i++) {
      if($anglecoeff[$i][0] =~ /^harm/) {
	printf $fh_data "%8u", $i+1;
	printf $fh_data " %-15s", "harmonic" if(@anglestyles>1);
	printf $fh_data " %15.6f %15.6f\n" ,0.5*$fenergy*$anglecoeff[$i][4],$anglecoeff[$i][5];
      } elsif($anglecoeff[$i][0] =~ /^hcos/) {
	printf $fh_data "%8u", $i+1;
	printf $fh_data " %-15s", "cosine/squared" if(@anglestyles>1);
	printf $fh_data " %15.6f %15.6f\n" ,0.5*$fenergy*$anglecoeff[$i][4],$anglecoeff[$i][5];
#       } elsif($anglecoeff[$i][0] =~ /^quar/) {
# 	printf $fh_data "%8u", $i+1;
# 	printf $fh_data " %-15s","quartic" if(@anglestyles>1);
# 	printf $fh_data " %15.6f %15.6f %15.6f %15.6f\n",$anglecoeff[$i][5],
# 	0.5*$fenergy*$anglecoeff[$i][4],0.5*$fenergy*$anglecoeff[$i][6],0.5*$fenergy*$anglecoeff[$i][7];
      } else {
	print "**** error: angle potential \"$anglecoeff[$i][0]\" is not implemented (yet)!\n";
	return 1;
      }
    }
  }
  if(@dihedralcoeff) {
    print $fh_data "\nDihedral Coeffs\n\n";
    for($i=0;$i<@dihedralcoeff;$i++) {
      printf $fh_data "%8u", $i+1;
      if($dihedralcoeff[$i][0] =~ /^harm/) {
	printf $fh_data " %-15s","quadratic" if(@dihedralstyles>1);
	printf $fh_data " %15.6g %15.6g\n",0.5*$fenergy*$dihedralcoeff[$i][5],$dihedralcoeff[$i][6];
      } elsif($dihedralcoeff[$i][0] =~ /^cos3$/) {
	printf $fh_data " %-15s","opls" if(@dihedralstyles>1);
	printf $fh_data " %15.6g %15.6g %15.6g %15.6g\n",
	$dihedralcoeff[$i][5]*$fenergy,
	$dihedralcoeff[$i][6]*$fenergy,
	$dihedralcoeff[$i][7]*$fenergy,0;
      } elsif($dihedralcoeff[$i][0] =~ /^cos$/) {
	printf $fh_data " %-15s","charmm" if(@dihedralstyles>1);
	printf $fh_data " %15.6g %4u %5.0f %7.5f\n",
	$dihedralcoeff[$i][5]*$fenergy,$dihedralcoeff[$i][7],$dihedralcoeff[$i][6],0;
      } else {
	print "**** error: dihedral potential \"$dihedralcoeff[$i][0]\" is not implemented (yet)!\n";
	return 1;
      }
    }
  }
  
  print $fh_data "\nAtoms\n\n";
  $i=0;
  for($t=0;$t<$field_nummols[$fi];$t++) {
    for($m=0;$m<$mol_numents[$fi][$t];$m++) {
      for($k=0;$k<$mol_numatoms[$fi][$t];$k++) {
	$i++;
	printf $fh_data "%8u %4u %4u %15.6f %15.6f %15.6f %15.6f\n",$i,$t+1,$mol_atomtype[$t][$k]+1,
	$mol_atomdata[$fi][$t][$k][2]*$fcharge,$cdata[$fi][$t][$m][$k][0]*$flength,
	$cdata[$fi][$t][$m][$k][1]*$flength,$cdata[$fi][$t][$m][$k][2]*$flength;
      }
    }
  }
  print $fh_data "\n";
  if($numbonds>0) {
    print $fh_data "\Bonds\n\n";
    $i=0;
    $indexoffset=1;
    for($t=0;$t<$field_nummols[$fi];$t++) {
      for($m=0;$m<$mol_numents[$fi][$t];$m++) {
	for($k=0;$k<$mol_numbonds[$fi][$t];$k++) {
	  $i++;
	  printf $fh_data "%8u %4u %8u %8u\n",$i,$mol_bondcoeff[$t][$k]+1,
	  $mol_bonddata[$fi][$t][$k][1]+$indexoffset,
	  $mol_bonddata[$fi][$t][$k][2]+$indexoffset;
	}
	$indexoffset+=$mol_numatoms[$fi][$t];
      }
    }
    print $fh_data "\n";
  }
  if($numangles>0) {
    print $fh_data "\Angles\n\n";
    $i=0;
    $indexoffset=1;
    for($t=0;$t<$field_nummols[$fi];$t++) {
      for($m=0;$m<$mol_numents[$fi][$t];$m++) {
	for($k=0;$k<$mol_numangles[$fi][$t];$k++) {
	  $i++;
	  printf $fh_data "%8u %4u %8u %8u %8u\n",$i,$mol_anglecoeff[$t][$k]+1,
	  $mol_angledata[$fi][$t][$k][1]+$indexoffset,
	  $mol_angledata[$fi][$t][$k][2]+$indexoffset,
	  $mol_angledata[$fi][$t][$k][3]+$indexoffset;
	}
	$indexoffset+=$mol_numatoms[$fi][$t];
      }
    }
    print $fh_data "\n";
  }
  if($numdihedrals>0) {
    print $fh_data "\Dihedrals\n\n";
    $i=0;
    $indexoffset=1;
    for($t=0;$t<$field_nummols[$fi];$t++) {
      for($m=0;$m<$mol_numents[$fi][$t];$m++) {
	for($k=0;$k<$mol_numdihedrals[$fi][$t];$k++) {
	  $i++;
	  printf $fh_data "%8u %4u %8u %8u %8u %8u\n",$i,$mol_dihedralcoeff[$t][$k]+1,
	  $mol_dihedraldata[$fi][$t][$k][1]+$indexoffset,
	  $mol_dihedraldata[$fi][$t][$k][2]+$indexoffset,
	  $mol_dihedraldata[$fi][$t][$k][3]+$indexoffset,
	  $mol_dihedraldata[$fi][$t][$k][4]+$indexoffset;
	}
	$indexoffset+=$mol_numatoms[$fi][$t];
      }
    }
    print $fh_data "\n";
  }
  
  if($config_key[$ci]>0) {
    print $fh_data "\nVelocities\n\n";
    $i=0;
    for($t=0;$t<$field_nummols[$fi];$t++) {
      for($m=0;$m<$mol_numents[$fi][$t];$m++) {
	for($k=0;$k<$mol_numatoms[$fi][$t];$k++) {
	  $i++;
	  printf $fh_data "%8u %15.6f %15.6f %15.6f\n",$i,@{$cdata[$fi][$t][$m][$k]}[3..5];
	}
      }
    }
  }
  
  close($fh_data);
  
  ############ generate input file ############################################
  if(not open($fh_in, ">", "$filename.in")) {
    print "**** error: Can't open output file $filename.in: $!\n";
    return 1;
  }
#   print $fh_in "log $filename.log append\n";
  if($oi<0) {
    print $fh_in "# $config_title[$ci]\n";
  } else {
    print $fh_in "# $control_title[$oi]\n";
  }
  print $fh_in "units      $lunits\n";
  if($periodic_key[$ci]==0) { #0d
    print $fh_in "boundary f f f\n";
  } elsif($periodic_key[$ci]==6) { #2d
    print $fh_in "boundary p p f\n";
  } else { # 3d
    print $fh_in "boundary p p p\n";
  }
  print $fh_in "atom_style full\n";
  
  if($numbonds>0) {
    print $fh_in "bond_style"; 
    print $fh_in " hybrid" if(@bondstyles>1); 
    foreach $i (@bondstyles) {
      if($i=~/^harm/) {
	print $fh_in " harmonic";
      } elsif($i=~/^quar/) {
	print $fh_in " class2";
      }
    }
    print $fh_in "\n";
  }
  print $fh_in "special_bonds $specialbonds\n";
  if($numangles>0) {
    print $fh_in "angle_style"; 
    print $fh_in " hybrid" if(@anglestyles>1); 
    foreach $i (@anglestyles) {
      if($i=~/^harm/) {
	print $fh_in " harmonic";
      } elsif($i=~/^hcos/) {
	print $fh_in " cosine/squared";
      } elsif($i=~/^quar/) {
	print $fh_in " quartic";
      }
    }
    print $fh_in "\n";
  }
  if($numdihedrals>0) {
    print $fh_in "dihedral_style"; 
    print $fh_in " hybrid" if(@dihedralstyles>1); 
    foreach $i (@dihedralstyles) {
      if($i=~/^harm/) {
	print $fh_in " quadratic";
      } elsif($i=~/^cos$/) {
	print $fh_in " charmm";
      } elsif($i=~/^cos3$/) {
	print $fh_in " opls";
      }
    }
    print $fh_in "\n";
  }
  print $fh_in "read_data $filename.data\n";
  if(@pairstyles>1) {
    print $fh_in "pair_style hybrid/overlay";
    #coulomb potential
    print $fh_in " $coulstyle $couladd" if(defined($coulstyle));
    # vdw potentials
    foreach $i (@pairstyles) {
      print $fh_in " $i $rvdw";
    }
    print $fh_in "\n";
    print $fh_in "pair_coeff     *    *   $coulstyle\n" if(defined($coulstyle));
  } elsif(@pairstyles>0) {
    print $fh_in "pair_style $pairstyles[0]";
    if(defined($coulstyle)) {
      print $fh_in "/$coulstyle $couladd";
      print $fh_in " $rvdw" if($rvdw != $cutglob);
    } else {
      print $fh_in " $rvdw"
    }
    print $fh_in "\n";
  }
  if($onlyii) {
    print $fh_in "pair_modify shift yes mix arithmetic\n";
  } else {
    print $fh_in "pair_modify shift yes\n";
  }
  for($i=0;$i<@paircoeff;$i++) {
    print $fh_in "pair_coeff ";
    printf $fh_in " %4u %4u",$paircoeff[$i][0],$paircoeff[$i][1];
    printf $fh_in " %-10s",$paircoeff[$i][2] if(@pairstyles>1);
    for($j=3;$j<@{$paircoeff[$i]};$j++) {
      printf $fh_in "%15.6g",$paircoeff[$i][$j];
    }
    print $fh_in "\n";
  }
  for($i=0;$i<@molgroups;$i++) {
    print $fh_in "# $molgroups[$i]\n";
  }
  if(@frozenatoms) {
    $group_integrate="integrate";
    print $fh_in "group frozen id ",join(" ",@frozenatoms),"\n";
    print $fh_in "group integrate subtract all frozen\n";
    print $fh_in "neigh_modify exclude group frozen frozen\n";
  } else {
    $group_integrate="all";
  }
  
  if($oi<0) {
    close($fh_in);
    return 0;
  }
  
  # data from CONTROL file
  $ifix=1; # index for fixes
  print $fh_in "timestep ",$control_timestep[$oi]*$ftime,"\n";
  if($control_ensemble[$oi]=~/^nvt/) {
    if($control_thermostat[$oi]=~/^ber/) {
      print $fh_in "fix ensemble $group_integrate nve\n"; $ifix++;
      print $fh_in "fix thermostat $group_integrate temp/berendsen $control_temp[$oi] $control_temp[$oi] ",
                   $control_ensembleprm[$oi][0]*$ftime,"\n"; $ifix++;
    } elsif($control_thermostat[$oi]=~/^hoover/) {
      print $fh_in "fix thermostat $group_integrate nvt temp $control_temp[$oi] $control_temp[$oi] ",
                   $control_ensembleprm[$oi][0]*$ftime,"\n"; $ifix++;
    } elsif($control_thermostat[$oi]=~/^evans/) {
      print "**** warning: evans thermostat is not included in Lammps! using hoover thermostat\n";
      print $fh_in "fix thermostat $group_integrate nvt temp $control_temp[$oi] $control_temp[$oi] ",
                   $control_ensembleprm[$oi][0]*$ftime,"\n"; $ifix++;
    }
  } elsif($control_ensemble[$oi]=~/^npt/) {
    if($control_thermostat[$oi]=~/^ber/) {
      print $fh_in "fix ensemble $group_integrate nve\n"; $ifix++;
      print $fh_in "fix thermostat $group_integrate temp/berendsen $control_temp[$oi] $control_temp[$oi] ",
                   $control_ensembleprm[$oi][0]*$ftime,"\n"; $ifix++;
      print $fh_in "fix barostat $group_integrate press/berendsen iso ",
                   $control_pressure[$oi]*$fpress," ",
                   $control_pressure[$oi]*$fpress," ",
                   $control_ensembleprm[$oi][1]*$ftime,"\n"; $ifix++;
    } elsif($control_thermostat[$oi]=~/^hoover/) {
      print $fh_in "fix ensemble $group_integrate npt temp $control_temp[$oi] $control_temp[$oi] ",
                   $control_ensembleprm[$oi][0]*$ftime,
                   "iso",
                   $control_pressure[$oi]*$fpress," ",
                   $control_pressure[$oi]*$fpress," ",
                   $control_ensembleprm[$oi][1]*$ftime,"\n"; $ifix++;
    } elsif($control_thermostat[$oi]=~/^evans/) {
      print "**** warning: evans thermostat is not included in Lammps! using hoover thermostat\n";
      print $fh_in "fix ensemble $group_integrate npt temp $control_temp[$oi] $control_temp[$oi] ",
                   $control_ensembleprm[$oi][0]*$ftime,
                   "iso",
                   $control_pressure[$oi]*$fpress," ",
                   $control_pressure[$oi]*$fpress," ",
                   $control_ensembleprm[$oi][1]*$ftime,"\n"; $ifix++;
    } elsif($control_thermostat[$oi]=~/^piston/) {
      print $fh_in "fix ensemble $group_integrate nve\n"; $ifix++;
      print $fh_in "fix thermostat $group_integrate temp/berendsen $control_temp[$oi] $control_temp[$oi] ",
                   $control_ensembleprm[$oi][0]*$ftime,"\n"; $ifix++;
      print $fh_in "fix barostat $group_integrate wall/pressure/harmonic zhi",
                   " ",($control_ensembleprm[$oi][2]+5.0)*$flength, # 5.0 cutoff
                   " ",$control_ensembleprm[$oi][4]*$fenergy/(2*$flength**2),
                   " 0.0",
                   " ",5.0*$flength, #cutoff
                   " ",$control_pressure[$oi]*$fpress, # pressure
                   " ",$control_ensembleprm[$oi][3]*$fmass/($flength**2), # mass
                   " ",$control_ensembleprm[$oi][1],"\n"; $ifix++;
    } else {
      print "**** error: unknown barostat/thermostat!\n"; return 1;
    }
  }
  if(defined($control_nofic[$oi])) {
    print $fh_in "fix nofic $group_integrate momentum $control_nofic[$oi] linear 1 1 1 rescale\n";
  }
  if(defined($control_mult[$oi])) {
    print "**** warning: directive \"mult\" is not implemented (yet)! neglecting it\n";
  }
  if(defined($control_shake[$oi])) {
    print "**** warning: directive \"shake\" is not implemented (yet)! neglecting it\n";
    #print $fh_in "fix shake $group_integrate shake $control_shake[$oi] 20 0\n"; $ifix++;
  }
  if(@{$control_traj[$oi]}) {
    if(not defined($dumpstyle) or $dumpstyle=~/dcd/) {
      print $fh_in "dump 1 all $dumpstyle $control_traj[$oi][1] $filename.dcd\n";
    } elsif($dumpstyle=~/^cfg/) {
      print $fh_in "dump 1 all $dumpstyle $control_traj[$oi][1] $filename.cfg\n";
    } elsif($dumpstyle=~/^xyz/) {
      print $fh_in "dump 1 all $dumpstyle $control_traj[$oi][1] $filename.xyz\n";
    } elsif($dumpstyle=~/^xtc/) {
      print $fh_in "dump 1 all $dumpstyle $control_traj[$oi][1] $filename.xtc\n";
    } else {
      if($control_traj[$oi][2]==1) {
        print $fh_in "dump 1 all custom $control_traj[$oi][1] $filename.lammpstrj ",
        "element x y z vx vy vz\n";
      } elsif($control_traj[$oi][2]==1) {
        print $fh_in "dump 1 all custom $control_traj[$oi][1] $filename.lammpstrj ",
        "element x y z vx vy vz fx fy fz\n";
      } else {
        print $fh_in "# dump 1 all custom $control_traj[$oi][1] $filename.lammpstrj ",
        "element x y z\n";
      }
      print $fh_in "# dump_modify 1 sort id element";
      for($i=0;$i<@atomtypes;$i++) {
        print $fh_in " $atomtypes[$i][0]";
      }
    }
    print $fh_in "\n";
    print $fh_in "restart 1000 $filename.restart1 $filename.restart2\n";
  }
  if(defined($control_stats[$oi])) {
    print $fh_in "thermo $control_stats[$oi]\n";
    print $fh_in "thermo_style custom step etotal temp pe evdwl ecoul ebond ",
    "eangle edihed eimp enthalpy vol density\n";
    print $fh_in "thermo_modify flush yes\n";
  }
  if($control_equil[$oi]) {
    print "**** warning: directive \"equil\" is not implemented (yet)! neglecting it\n";
  }
  if($control_cap[$oi]) {
    print "**** warning: directive \"cap\" is not implemented (yet)! neglecting it\n";
  }
  if(defined($control_minimmethod[$oi])) {
    print "**** warning: directive \"mimimize\" is not implemented (yet)! neglecting it\n";
  }
  if(defined($control_optimmethod[$oi])) {
    print "**** warning: directive \"optimize\" is not implemented (yet)! neglecting it\n";
  }
  if(defined($control_restartmethod[$oi])) {
    if($control_restartmethod[$oi]=~/^restart\s+scale/i) {
      print $fh_in "velocity $group_integrate scale $control_temp[$oi]\n";
    } elsif($control_restartmethod[$oi]=~/^restart\s+scale/i) {
      print $fh_in "velocity $group_integrate scale $control_temp[$oi]\n";
    } elsif($control_restartmethod[$oi]=~/^restart$/i) {
      print $fh_in "read_restart $filename.restart\n";
    } 
  } elsif($control_temp[$oi]>0) {
    print $fh_in "velocity all create $control_temp[$oi] ",int(rand(1.0e6))," rot yes dist gaussian\n";
  } else {
    print "**** error: don't know what to do for T=0K\n"; return 1;
  }
  
  print $fh_in "run $control_steps[$oi]\n";
  print $fh_in "write_restart $filename.restart\n";
  close($fh_in);
  
}

sub find_lammps_atomtype {
  my $fi = $_[0];
  my $t  = $_[1];
  my $k  = $_[2];
  my @arr  = @{$_[3]};
  for(my $i=0;$i<@arr;$i++) {
    if( $mol_atomdata[$fi][$t][$k][0] eq $arr[$i][0]
    and $mol_atomdata[$fi][$t][$k][1] == $arr[$i][1]) {
      return $i
    }
  }
  return -1;
}
sub find_lammps_bondtype {
  my $fi = $_[0];
  my $t  = $_[1];
  my $k  = $_[2];
  my @arr  = @{$_[3]};
  my ($i,$j);
  looparr:for($i=0;$i<@arr;$i++) {
    next if(@{$arr[$i]}!=@{$mol_bonddata[$fi][$t][$k]});
    next if($arr[$i][0] ne $mol_bonddata[$fi][$t][$k][0]);
    for($j=3;$j<@{$arr[$i]};$j++) {
#       print "$i $j asdf $arr[$i][$j] $mol_bonddata[$fi][$t][$k][$j]\n";
      next looparr if($arr[$i][$j] != $mol_bonddata[$fi][$t][$k][$j]);
    }
    return $i;
  }
  return -1;
}
sub find_lammps_angletype {
  my $fi = $_[0];
  my $t  = $_[1];
  my $k  = $_[2];
  my @arr  = @{$_[3]};
  my ($i,$j);
  looparr:for($i=0;$i<@arr;$i++) {
    next if(@{$arr[$i]}!=@{$mol_angledata[$fi][$t][$k]});
    next if($arr[$i][0] ne $mol_angledata[$fi][$t][$k][0]);
    for($j=4;$j<@{$arr[$i]};$j++) {
      next looparr if($arr[$i][$j] != $mol_angledata[$fi][$t][$k][$j]);
    }
    return $i;
  }
  return -1;
}
sub find_lammps_dihedraltype {
  my $fi = $_[0];
  my $t  = $_[1];
  my $d  = $_[2];
  my @arr  = @{$_[3]};
  my ($i,$j);
  looparr:for($i=0;$i<@arr;$i++) {
    next if(@{$arr[$i]}!=@{$mol_dihedraldata[$fi][$t][$d]});
    next if($arr[$i][0] ne $mol_dihedraldata[$fi][$t][$d][0]);
    for($j=5;$j<@{$arr[$i]};$j++) {
      next looparr if($arr[$i][$j] != $mol_dihedraldata[$fi][$t][$d][$j]);
    }
    return $i;
  }
  return -1;
}
sub find_dlpoly_pairtype {
  # returns the number of vdw data for a certain pair of atoms
  my $fi = $_[0];
  my $k  = $_[1];
  my ($i);
  my $n=0;
  for($i=0;$i<@{$field_vdwdata[$fi]};$i++) {
    if(($field_vdwdata[$fi][$i][0] eq $field_vdwdata[$fi][$k][0] and $field_vdwdata[$fi][$i][1] eq $field_vdwdata[$fi][$k][1]) or
       ($field_vdwdata[$fi][$i][0] eq $field_vdwdata[$fi][$k][1] and $field_vdwdata[$fi][$i][1] eq $field_vdwdata[$fi][$k][0])) {
      $n++;
    }
  }
  return $n;
}
sub find_lammps_atomtype_pair {
  my $name      = $_[0];
  my @atomtypes = @{$_[1]};
  my @arr       = ();
  my ($i);
  for($i=0;$i<@atomtypes;$i++) {
    push(@arr,$i) if($atomtypes[$i][0] eq $name);
  }
  return @arr;
}