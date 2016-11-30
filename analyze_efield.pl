#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;

use constant {
  u  => 1.6605402e-27,  # unit of mass
  e  => 1.60217733e-19, # unit of charge
  t0 => 1e-12,          # unit of time
  l0 => 1e-10,          # unit of length
};
use constant {
  F0 => u*u*l0*l0/(t0*t0),  # unit of force
  E0 => u*u*l0*l0/(t0*t0*e) # unit of electric field
};

if($#ARGV<2) {
  print "input format:\n";
  print "1. name of HISTORY-file with fullerenes\n";
  print "2. name of corresponding FIELD-file\n";
  print "3. name of output-directory\n";
  print "[4. start timestep]\n";
  print "[5. end timestep]\n";
  exit 1;
}

$histfilename  = $ARGV[0];
$fieldfilename = $ARGV[1];
$outdir    = $ARGV[2];
if($#ARGV>2) {
  $starttimestep = $ARGV[3];
} else {
  $starttimestep = 0;
}

if($#ARGV>3) {
  $endtimestep = $ARGV[4];
} else {
  $endtimestep = 99999999;
}

$field_read_success = &read_field_file($fieldfilename,0);
open(HISTORY,"<",$histfilename) or die "**** error: Can't open HISTORY-file:\n$!\n";

@analyze_mols = ();
for($t=0;$t<$field_nummols[0];$t++) {
  if($mol_name[0][$t] =~ /C60/i) {
    push(@analyze_mols,$t);
    for($a=0;$a<$mol_numatoms[0][$t];$a++) {
      if($mol_atomdata[0][$t][$a][0] =~ /CA/i or $mol_atomdata[0][$t][$a][0] =~ /CX/i) {
	push(@{$analyze_atoms[$#analyze_mols]},$a);
      }
    }
  }
}

# indices for hdata array:
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

$f=0;
while(<HISTORY>) {
  if(not /timestep/i) {
    print "**** error reading HISTORY-file: \"timestep\" expected but found:\n",$_;
    exit 1;
  } else {
    ($timestepnumber, $numatoms, $config_key, $periodic_key, $timestep) =
    /timestep\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    if($config_key<2) {
      print "**** error: no force information in timestep $f of HISTORY-file!\n";
    }
    if($timestepnumber > $endtimestep) {
      last;
    } elsif($timestepnumber < $starttimestep) {
      # skip cell vectors
      $_=<HISTORY>;
      $_=<HISTORY>;
      $_=<HISTORY>;
      # skip atom data
      for($t=0;$t<$field_nummols[0];$t++) {
	for($m=0;$m<$mol_numents[0][$t];$m++) {
	  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	    $_=<HISTORY>;
	    $_=<HISTORY>;
	    $_=<HISTORY>;
	  }
	}
      }
      next;
    } else { # read data if timestepnumber is within specified range
      print "analyzing frame $timestepnumber\n";
      $_=<HISTORY>;
      ($cell[0][0],$cell[0][1],$cell[0][2])=/\s*(\S+)\s+(\S+)\s+(\S+)/;
      $_=<HISTORY>;
      ($cell[1][0],$cell[1][1],$cell[1][2])=/\s*(\S+)\s+(\S+)\s+(\S+)/;
      $_=<HISTORY>;
      ($cell[2][0],$cell[2][1],$cell[2][2])=/\s*(\S+)\s+(\S+)\s+(\S+)/;
      $size[0]=($cell[0][0]+$cell[1][0]+$cell[2][0])/2;
      $size[1]=($cell[0][1]+$cell[1][1]+$cell[2][1])/2;
      $size[2]=($cell[0][2]+$cell[1][2]+$cell[2][2])/2;
      for($t=0;$t<$field_nummols[0];$t++) {
	for($m=0;$m<$mol_numents[0][$t];$m++) {
	  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	    $_=<HISTORY>;
	    ($hdata[$t][$m][$a][0]) = /^\s*(\S+)/;
# 	    if($hdata[$t][$m][$a][0] ne $mol_atomdata[0][$t][$a][0]) {
# 	      print "error: atom name ".$hdata[$t][$m][$a][0]." in CONFIG-file $filenameincfg does not match\n";
# 	      print "name ".$mol_atomdata[0][$t][$a][0]." in FIELD file $filenameinfld for molecule ".$mol_name[0][$t]."\n";
# 	      exit 1;
# 	    }
	    $_=<HISTORY>;
	    ($hdata[$t][$m][$a][1],$hdata[$t][$m][$a][2],$hdata[$t][$m][$a][3]) = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
	    $_=<HISTORY>;
	    ($hdata[$t][$m][$a][4],$hdata[$t][$m][$a][5],$hdata[$t][$m][$a][6]) = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
	    $_=<HISTORY>;
	    ($hdata[$t][$m][$a][7],$hdata[$t][$m][$a][8],$hdata[$t][$m][$a][9]) = /^\s*(\S+)\s+(\S+)\s+(\S+)/;
	  }
	}
      }
    }
  }
  ############### analyze timestep ############################
  $time[$f] = $timestepnumber*$timestep;
  foreach $i (keys @analyze_mols) {
    @{$av_all[$f][$i]} = (0,0,0,0);
    foreach $j (keys @{$analyze_atoms[$i]}) {
      @{$av_atom[$f][$i][$j]} = (0,0,0,0);
    }
    @{$av_full[$f][$i]} = ();
    $t = $analyze_mols[$i];
    for($m=0;$m<$mol_numents[0][$t];$m++) {
      @{$av_full[$f][$i][$m]} = (0,0,0,0);
      foreach $j (keys @{$analyze_atoms[$i]}) {
	$a = $analyze_atoms[$i][$j];
	if($mol_atomdata[0][$t][$a][2] != 0) { # if charge != 0 (otherwise, the E-field is 0)
	  $efield[0] = $hdata[$t][$m][$a][7]/$mol_atomdata[0][$t][$a][2];
	  $efield[1] = $hdata[$t][$m][$a][8]/$mol_atomdata[0][$t][$a][2];
	  $efield[2] = $hdata[$t][$m][$a][9]/$mol_atomdata[0][$t][$a][2];
	  $efield[3] = sqrt($efield[0]**2+$efield[1]**2+$efield[2]**2);
	  $av_full[$f][$i][$m][0] += $efield[0];
	  $av_full[$f][$i][$m][1] += $efield[1];
	  $av_full[$f][$i][$m][2] += $efield[2];
	  $av_atom[$f][$i][$j][0] += $efield[0];
	  $av_atom[$f][$i][$j][1] += $efield[1];
	  $av_atom[$f][$i][$j][2] += $efield[2];
	  $av_atom[$f][$i][$j][3] += $efield[3];
	  $av_all[$f][$i][0] += $efield[0];
	  $av_all[$f][$i][1] += $efield[1];
	  $av_all[$f][$i][2] += $efield[2];
	  $av_all[$f][$i][3] += $efield[3];
	}
      }
      $av_full[$f][$i][$m][0] /= @{$analyze_atoms[$i]};
      $av_full[$f][$i][$m][1] /= @{$analyze_atoms[$i]};
      $av_full[$f][$i][$m][2] /= @{$analyze_atoms[$i]};
      $av_full[$f][$i][$m][3] = sqrt($av_full[$f][$i][$m][0]**2+$av_full[$f][$i][$m][1]**2+$av_full[$f][$i][$m][2]**2);
    }
    foreach $j (keys @{$analyze_atoms[$i]}) {
      $av_atom[$f][$i][$j][0] /= $mol_numents[0][$t];
      $av_atom[$f][$i][$j][1] /= $mol_numents[0][$t];
      $av_atom[$f][$i][$j][2] /= $mol_numents[0][$t];
      $av_atom[$f][$i][$j][3] /= $mol_numents[0][$t];
    }
    $av_all[$f][$i][0] /= @{$analyze_atoms[$i]}*$mol_numents[0][$t];
    $av_all[$f][$i][1] /= @{$analyze_atoms[$i]}*$mol_numents[0][$t];
    $av_all[$f][$i][2] /= @{$analyze_atoms[$i]}*$mol_numents[0][$t];
    $av_all[$f][$i][3] /= @{$analyze_atoms[$i]}*$mol_numents[0][$t];
  }
  $f++;
}

print "printing output...\n";

$fmax = $f;

if(not -d $outdir) {
  mkdir $outdir or die "cannot create output-directory: $!";
}

foreach $i (keys @analyze_mols) {
  $t = $analyze_mols[$i];
  open(AV_ALL,">","$outdir/av_all_$mol_name[0][$t]");
  open(AV_FULL,">","$outdir/av_full_$mol_name[0][$t]");
  open(AV_ATOMS,">","$outdir/av_atom_$mol_name[0][$t]");
  for($f=0;$f<$fmax;$f++) {
    printf AV_ALL "%20.10f %20.10f %20.10f %20.10f %20.10f\n", $time[$f], @{$av_all[$f][$i]};
  }
  for($m=0;$m<$mol_numents[0][$t];$m++) {
    for($f=0;$f<$fmax;$f++) {
      printf AV_FULL "%20.10f %20.10f %20.10f %20.10f %20.10f\n", $time[$f], @{$av_full[$f][$i][$m]};
    }
    print AV_FULL "\n\n";
  }
  foreach $j (keys @{$analyze_atoms[$i]}) {
    $a = $analyze_atoms[$i][$j];
    for($f=0;$f<$fmax;$f++) {
      printf AV_ATOMS "%20.10f %20.10f %20.10f %20.10f %20.10f\n", $time[$f], @{$av_atom[$f][$i][$j]};
    }
    print AV_ATOMS "\n\n";
  }
  close(AV_ALL);
  close(AV_FULL);
  close(AV_ATOMS);
}
