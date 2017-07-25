#!/usr/bin/perl

use Switch;
use dlpoly_utility;
use hanno_utility;
use lammps_utility;
use Math::Trig;

if($#ARGV<1) {
  print " 1. name of input CONFIG-file\n";
  print " 2. name of output CONFIG-file\n";
  print " input/output options:\n";
  print " -a            append output to file instead of overwrite\n";
  print " -b xyz <real> crop system (e.g. -b x -200 50 y -80 60)\n";
  print " -c [xyz]+     center the block in x,y and/or z-direction\n";
  print " -e <3*real>   scale system in x,y and z-direction by factor\n";
  print " -f <string>   name of corresponding FIELD file\n";
  print " -fo <string>  name of output FIELD file\n";
  print " -fr <int>     frame number/timestep in xyz/HISTORY file (first=1, default: last)\n";
  print "               (starting from 1, default: last)\n";
  print " -g            rotate system to make cell match VMD reqirements\n";
  print " -i <string>   format of input config file (xyz,hist,cfg,dcd)\n";
  print " -o <string>   format of output config file (xyz,hist,cfg)\n";
  print "operations (executed in the order of appearance):\n";
  print " -k [012]      new config_key\n";
  print " -l <9*real>   new cell vectors\n";
  print " -m [xyz]+     mirror cell in x,y and/or z-direction?\n";
  print " -n <n*str>    list of molecule names to include\n";
  print " -p [0-6]      new periodic key\n";
  print " -r            remap/unwrap molecules\n";
  print " -s <3*real>   shift atoms by vector\n";
  print " -sph <4*real> cut a sphere with given radius and center\n";
  print " -isp <4*real> cut a spherical hole with given radius and center\n";
  print " -to <2*real>  cut a truncated octahedron with specified radius and buffer width\n";
  print " -t <n*str>    list of atom types to include\n";
  print " -uc <3*real>  rotation center\n";
  print " -u xyz <real> rotate system by axis by angle in degree (e.g. -u x 50 y 60)\n";
  print " -w            wrap atoms back into cell\n";
  print " -x <2*int>    extract molecule (int1=type int2=id starting from 0)\n";
  print " -y <3*int>    multiply system in x,y and z-direction by factor\n";
  print " -ys <3*real>  add spacing between entities when scaling\n";
  exit;
}

$incfg  = $ARGV[0];
$outcfg = $ARGV[1];
$writemode   = ">";
$targetframe = -1;
$informat    = 'cfg';
$outformat   = 'cfg';
@axisvecs    = ([1,0,0],[0,1,0],[0,0,1]);
@rotcenter   = (0,0,0);
@spacing     = (0,0,0);
@xyz         = ("x","y","z");

$lincludenames = 0;
$lincludetypes = 0;
$lneedfield    = 0;
$lwarnremap    = 0;
$lremap        = 0;
$linfield      = 0;

for($i=2;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-a$/ or /^--write[-_]mode/) {
      $writemode=">>";
    } case (/^-b$/ or /^--crop/ or /^--cut/) {
      $lneedfield = 1;
      while(not $ARGV[$i+1] =~ /^-/ and $i<$#ARGV) {
	$i++;
	if(not $ARGV[$i] =~ /^[xyz]$/i) {
	  print "**** error: could not interpret \"$ARGV[$i]\" as crop axis!\n";
	  exit 1;
	}
	if(not check_real(@ARGV[$i+1..$i+2])) {
	  print "**** error: could not interpret \"$ARGV[$i+1]\" and/or \"$ARGV[$i+2]\"  as crop values!\n";
	  exit 1;
	}
	$i+=2;
      }
    } case (/^-c$/ or /^--center/) {
      $i++;
      if(not $ARGV[$i] =~ /^[xyz]+$/) {
	print "**** error: input \"$ARGV[$i]\" is not valid for centering!\n";
	exit 1;
      }
    } case (/^-e$/ or /^--scale/) { # scale system
      if(not check_real(@ARGV[$i+1..$i+3])) {
	print "**** error: scale factors must be three real numbers!\n";
	exit 1;
      }
      $i += 3;
    } case (/^-f$/ or /^--field/) { # input FIELD
      $i++;
      $infld = $ARGV[$i];
      $linfield = 1;
    } case (/^-fo$/ or /^--field[-_]out/ or /^--out[-_]field/) { # output FIELD
      $i++;
      $outfld = $ARGV[$i];
      $loutfield = 1;
    } case (/^-fr$/) { # target frame
      $i++;
      $targetframe = $ARGV[$i];
      if(not check_integer($targetframe)) {
	print "**** error: target frame must be integer number!\n";
	exit 1;
      }
    } case (/^-g$/ or /^--rotvmd/) {
    } case (/^-i$/ or /^--config/) { # output format for config file
      $i++;
      $informat = $ARGV[$i];
      if(not $informat=~/^(cfg|xyz|hist|dcd)$/i) {
	print "**** error: output format must be 'xyz', 'cfg', 'hist' or 'dcd'!\n";
	exit 1;
      }
    } case (/^-isp$/ or /^--sphere[-_]hole/) { # spherical hole
      $lwarnremap = 1 if(not $lremap);
      if(not check_real(@ARGV[$i+1..$i+4])) {
	print "**** error: radius and sphere center must be real numbers!\n";
	exit 1;
      }
      $i += 4;
    } case (/^-k$/ or /^--config[-_]key/) { # config_key
      $i++;
      if(not $ARGV[$i] =~ /^[012]$/) {
	print "**** error: config_key must be '0', '1' or '2'!\n";
	exit 1;
      }
    } case (/^-l$/ or /^--cell/) { # set cell vectors
      for($j=0;$j<9;$j++) {
	$i++;
	if(not check_real($ARGV[$i])) {
	  print "**** error: cell vector components must be real numbers!\n";
	  exit 1;
	}
      }
    } case (/^-m$/ or /^--mirror/) { # mirror
      $i++;
      if(not $ARGV[$i] =~ /^[xyz]+$/) {
	print "**** error: input \"$ARGV[$i]\" is not valid for mirroring!\n";
	exit 1;
      }
    } case (/^-n$/ or /^--(include)?.*names/) { # include molecule names
      $lincludenames = 1;
      while(not $ARGV[$i+1] =~ /^-/ and $i<$#ARGV) {
	$i++;
      }
    } case (/^-o$/ or /^--config[-_]out/ or /^--out[-_]config/) { # output format for config file
      $i++;
      $outformat = $ARGV[$i];
      if(not $outformat=~/^(cfg|xyz|hist)$/i) {
	print "**** error: output format must be 'xyz', 'cfg' or 'hist'!\n";
	exit 1;
      }
    } case (/^-p$/ or /^--periodic[-_]key/) { #periodic_key
      $i++;
      if(not $ARGV[$i] =~ /^[0-6]$/) {
	print "**** error: config_key must be an integer between 0 and 6!\n";
	exit 1;
      }
    } case (/^-r$/ or /^--remap/ or /^--unwrap/) {
      $lremap = 1;
    } case (/^-s$/ or /^--shift/) { # shift
      if(not check_real(@ARGV[$i+1..$i+3])) {
	print "**** error: shift vector must consist of three real numbers!\n";
	exit 1;
      }
      $i += 3;
    } case (/^-v$/ or /^--cell[-_]vec/) { #cell vectors
      if(not check_real(@ARGV[$i+1..$i+9])) {
	print "**** error: cell vectors must be real numbers!\n";
	exit 1;
      }
    } case (/^-sph$/ or /^--sphere/) { # cut sphere
      $lwarnremap = 1 if(not $lremap);
      if(not check_real(@ARGV[$i+1..$i+4])) {
	print "**** error: radius and sphere center must be real numbers!\n";
	exit 1;
      }
      $i += 4;
    } case (/^-to$/ or /^--truncated/) { # cut truncated octahedron
      $lwarnremap = 1 if(not $lremap);
      if(not check_real(@ARGV[$i+1..$i+2])) {
	print "**** error: radius and buffer width for truncated octahedron must be real numbers!\n";
	exit 1;
      }
      $i+=2;
    } case (/^-t$/ or /^--types/) { # include types
      $lincludetypes = 1;
      while(not $ARGV[$i+1] =~ /^-/ and $i<$#ARGV) { $i++; }
    } case (/^-uc/ or /^--rot.*center/) { # rotation center
      if(not check_real(@ARGV[$i+1..$i+3])) {
	print "**** error: rotation center must consist of three real numbers!\n";
	exit 1;
      }
      $i += 3;
    } case (/^-u$/ or /^--rotat.*$/) { # rotate system
      while(not $ARGV[$i+1] =~ /^-/ and $i<$#ARGV) {
	$i++;
	if(not $ARGV[$i] =~ /^[xyz]$/i) {
	  print "**** error: could not interpret \"$ARGV[$i]\" as rotation axis!\n";
	  exit 1;
	}
	$i++;
	if(not check_real($ARGV[$i])) {
	  print "**** error: could not interpret \"$ARGV[$i]\" as rotation anlge!\n";
	  exit 1;
	}
      }
    } case (/^-w$/ or /^--wrap/) {
    } case (/^-x$/ or /^--extract/) { # extract single molecule
      $lwarnremap = 1 if(not $lremap);
      if(not check_integer(@ARGV[$i+1..$i+2])) {
	print "**** error: molecule type and ID must be two integer numbers!\n";
	exit 1;
      }
      $i += 2;
    } case (/^-y$/ or /^--multiply/) { # multiply system
      $lwarnremap = 1 if(not $lremap);
      $lneedfield = 1;
      if(not check_integer(@ARGV[$i+1..$i+3])) {
	print "**** error: multiplication factors must be three integer numbers!\n";
	exit 1;
      }
      $i += 3;
    } case (/^-ys$/ or /^--spacing/) { # spacing for system multiplication
      if(not check_real(@ARGV[$i+1..$i+3])) {
	print "**** error: spacing must be three real numbers!\n";
	exit 1;
      }
      $i += 3;
    } else {
      print "**** error: expected optional flag but got \"$ARGV[$i]\"!\n";
      exit 1;
    }
  }
}

##### check dependencies ########################################

$lneedfield=1 if($lremap or $lincludenames or $lincludetypes);

if($lneedfield and not $linfield) {
  print "**** error: FIELD file is required! Please specify with -f option!\n";
  exit 1;
}

if($loutfield and $lincludetypes) {
  print "**** error: cannot output FIELD file with option -t";
  exit 1;
}

if($lwarnremap) {
  print "**** warning: no unwrap/remap command found! Output might be unreliable!\n";
}

#### read files #################################################

if($linfield) {
  $fi =  0;
  exit 1 if(read_field_file($infld,$fi) != 0);
} else {
  $fi = -1;
}
if($informat =~ /^hist$/i) { # read HISTORY file
  if(not open($fhhist,"<",$incfg)) {
    print "**** error: cannot open HISTORY file $incfg: $!\n";
    exit 1;
  }
  if($targetframe>=0) {
    exit 1 if(not find_history_timestep($fhhist,$targetframe));
    exit 1 if(read_history_timestep($fhhist,$fi,0)!=0);
  } else { # read the last frame
    $err=0;
    while($err==0) {
      $err=read_history_timestep($fhhist,$fi,0);
    }
    exit 1 if($err>0);
  }
  close($fhhist);
} elsif($informat =~ /^dcd/i) { # read dcd file
  if(not open($fhhist,"<",$incfg)) {
    print "**** error: cannot open DCD file $incfg: $!\n";
    exit 1;
  }
  if($targetframe>=0) {
    $err=0;
    while($err==0 and $frame_number[$ci]!=$targetframe) {
      $err=read_dcd_timestep($fhhist,$fi,0);
    }
    if($frame_number[$ci]!=$targetframe or $err>0) {
      print "**** error: timestep $targetframe was not found in $incfg!\n";
      exit 1;
    }
  } else { # read the last frame
    exit 1 unless(goto_last_dcd_timestep($fhhist)==0);
    exit 1 unless(read_dcd_timestep($fhhist,$fi,0)==0);
  }
  close($fhhist);
} elsif($informat =~ /^xyz$/i) { # read xyz file
  if(not open($fhxyz,"<",$incfg)) {
    print "**** error: cannot open XYZ file $incfg: $!\n";
    exit 1;
  }
  if($targetframe>=0) {
    for($i=0;$i<$targetframe;$i++) {
      $err=read_xyz_timestep($fhxyz,$fi,0);
      if($err<0) {
	print "**** error: XYZ file ended before frame $targetframe\n";
      }
      exit 1 if($err!=0);
    }
  } else { # read the last frame
    $err=0;
    while($err==0) {
      $err=read_xyz_timestep($fhxyz,$fi,0);
    }
    exit 1 if($err>0);
  }
  close($fhxyz);
} else { # read CONFIG file
  exit 1 if(read_config_file($incfg,$fi,0)!=0);
}

########## perform operations #########################################################

for($i=2;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-b$/ or /^--crop/ or /^--cut/) {
      print "cropping system\n";
      while(not $ARGV[$i+1] =~ /^-/ and $i<$#ARGV) {
	$i++;
	if($ARGV[$i] =~/^x/i) {
	  $c=0;
	} elsif($ARGV[$i] =~/^y/i) {
	  $c=1;
	} elsif($ARGV[$i] =~/^z/i) {
	  $c=2;
	}
	$min=$ARGV[$i+1];
	$max=$ARGV[$i+2];
	$i+=2;
	for($t=0;$t<@{$cdata[0]};$t++) {
	  for($m=0;$m<@{$cdata[0][$t]};$m++) {
	    for($a=0;$a<@{$cdata[0][$t][$m]};$a++) {
	      if($cdata[0][$t][$m][$a][$c]<$min or
		$cdata[0][$t][$m][$a][$c]>$max) {
		remove_mol_entity(0,0,$t,$m);
		$m--;
		last;
	      }
	    }
	  }
	}
      }
    } case (/^-c$/ or /^--center/) {
      print "centering system\n";
      &calc_minmaxpos(0);
      $i++;
      for($c=0;$c<3;$c++) {
	next if($ARGV[$i] !~ /$xyz[$c]/i);
	$center[$c]  =   ($maxpos[0][$c] + $minpos[0][$c])/2;
	for($t=0;$t<@{$cdata[0]};$t++) {
	  for($m=0;$m<@{$cdata[0][$t]};$m++) {
	    for($a=0;$a<@{$cdata[0][$t][$m]};$a++) {
	      $cdata[0][$t][$m][$a][$c] -= $center[$c];
	    }
	  }
	}
      }
    } case (/^-e$/ or /^--scale/) { # scale system
      print "scaling system\n";
      @scalefactors = @ARGV[$i+1..$i+3];
      $i += 3;
      for($t=0;$t<@{$cdata[0]};$t++) {
	for($m=0;$m<@{$cdata[0][$t]};$m++) {
	  for($a=0;$a<@{$cdata[0][$t][$m]};$a++) {
	    $cdata[0][$t][$m][$a][0] *= $scalefactors[0];
	    $cdata[0][$t][$m][$a][1] *= $scalefactors[1];
	    $cdata[0][$t][$m][$a][2] *= $scalefactors[2];
	  }
	}
      }
      for($j=0;$j<3;$j++) {
	for($c=0;$c<3;$c++) {
	  $cell[0][$j][$c] *= $scalefactors[$c];
	}
      }
    } case (/^-g$/ or /^--rotvmd/) {
      $err=rotate_cell_vmd(0);
      if($err==2) {
	print "**** warning: No cell information given, cannot rotate cell!\n";
	next;
      }
    } case (/^-isp$/ or /^--sphere[-_]hole/) { # spherical hole
      print "cutting spherical hole\n";
      $holeradius = $ARGV[$i+1];
      @holecenter = @ARGV[$i+2..$i+4];
      $i += 4;
      $holeradiussq = $holeradius*$holeradius;
      for($t=0;$t<@{$cdata[0]};$t++) {
	for($m=0;$m<@{$cdata[0][$t]};$m++) {
	  for($a=0;$a<@{$cdata[0][$t][$m]};$a++) {
	    $distsq=($cdata[0][$t][$m][$a][0]-$holecenter[0])**2
		  +($cdata[0][$t][$m][$a][1]-$holecenter[1])**2
		  +($cdata[0][$t][$m][$a][2]-$holecenter[2])**2;
# 	    print "$t $m $a $distsq\n";
	    if($distsq<$holeradiussq) {
	      print "removing molecule $t $m\n";
	      remove_mol_entity(0,0,$t,$m);
	      $m--;
	      last;
	    }
	  }
	}
      }
      
    } case (/^-k$/ or /^--config[-_]key/) { # config_key
      print "changing config key\n";
      $i++;
      $config_key[0] = $ARGV[$i];
      
    } case (/^-l$/ or /^--cell/) { # new cell vectors
      print "setting new cell vectors\n";
      for($j=0;$j<3;$j++) {
	for($k=0;$k<3;$k++) {
	  $i++;
	  $cell[0][$j][$k]=$ARGV[$i];
	}
      }
      
    
    } case (/^-m$/ or /^--mirror/) { # mirror
      print "mirroring system\n";
      &calc_minmaxpos(0);
      $dcenter[$c] = 2*($maxpos[$c] - $minpos[$c]);
      $i++;
      for($c=0;$c<3;$c++) {
	next if($ARGV[$i] !~ /$xyz[$c]/i);
	$v=$c+3; $f=$c+6;
	for($t=0;$t<@{$cdata[0]};$t++) {
	  for($m=0;$m<@{$cdata[0][$t]};$m++) {
	    for($a=0;$a<@{$cdata[0][$t][$m]};$a++) {
	      $cdata[0][$t][$m][$a][$c] = - $cdata[0][$t][$m][$a][$c] + $dcenter[$c];
	      $cdata[0][$t][$m][$a][$v] = - $cdata[0][$t][$m][$a][$v];
	      $cdata[0][$t][$m][$a][$f] = - $cdata[0][$t][$m][$a][$f];
	    }
	  }
	}
      }
      
    } case (/^-n$/ or /^--(include)?.*names/) { # include molecule names
      print "removing molecule types\n";
      @include=();
      while(not $ARGV[$i+1] =~ /^-/ and $i<$#ARGV) {
	$i++;
	push(@include,$ARGV[$i]);
      }
      looptype:for($t=0;$t<$field_nummols[0];$t++) {
	for($c=0;$c<@include;$c++) {
	  if($mol_name[0][$t] =~ /$include[$c]/i) {
	    next looptype;
	  }
	}
	remove_mol_type(0,0,$t);
	$t--;
      }
      
    } case (/^-p$/ or /^--periodic[-_]key/) { #periodic_key
      print "changing periodic key\n";
      $i++;
      $periodic_key[0] = $ARGV[$i];
      
    } case (/^-r$/ or /^--remap/ or /^--unwrap/) {
      for($t=0;$t<$field_nummols[0];$t++) {
      print "remapping $mol_name[0][$t]\n";
        if($mol_name[0][$t]=~/PA/) {
          for($m=0;$m<$mol_numents[0][$t];$m++) {
            remap_molecule2(0,0,$t,$m,10);
          }
        } else {
          for($m=0;$m<$mol_numents[0][$t];$m++) {
            remap_molecule2(0,0,$t,$m);
          }
        }
      }
      
    } case (/^-v/ or /^--cell[-_]vec/) { #cell vectors
      print "changing cell vectors\n";
      for($k=0;$k<3;$k++) {
	for($l=0;$l<3;$l++) {
	  $i++;
	  $cell[0][$k][$l] = $ARGV[$i];
	}
      }
      
    } case (/^-s$/ or /^--shift/) { # shift
      print "shifting system\n";
      @shiftvec = @ARGV[$i+1..$i+3];
      $i += 3;
      for($t=0;$t<@{$cdata[0]};$t++) {
	for($m=0;$m<@{$cdata[0][$t]};$m++) {
	  for($a=0;$a<@{$cdata[0][$t][$m]};$a++) {
	    $cdata[0][$t][$m][$a][0] += $shiftvec[0];
	    $cdata[0][$t][$m][$a][1] += $shiftvec[1];
	    $cdata[0][$t][$m][$a][2] += $shiftvec[2];
	  }
	}
      }
      
    } case (/^-sph$/ or /^--sphere/) { # cut sphere
      print "cutting sphere\n";
      $radius  = $ARGV[$i+1];
      @sphcenter = @ARGV[$i+2..$i+4];
      $i += 4;
      $radiussq = $radius*$radius;
      for($t=0;$t<@{$cdata[0]};$t++) {
	for($m=0;$m<@{$cdata[0][$t]};$m++) {
	  for($a=0;$a<@{$cdata[0][$t][$m]};$a++) {
	    $distsq=($cdata[0][$t][$m][$a][0]-$sphcenter[0])**2
		  +($cdata[0][$t][$m][$a][1]-$sphcenter[1])**2
		  +($cdata[0][$t][$m][$a][2]-$sphcenter[2])**2;
	    if($distsq>$radiussq) {
	      remove_mol_entity(0,0,$t,$m);
	      $m--;
	      last;
	    }
	  }
	}
      }
      
    } case (/^-to/ or /^--truncated/) { # cut truncated octahedron
      print "cutting truncated octahedron\n";
      &cut_truncated_octahedron(0,0,$ARGV[$i+1],$ARGV[$i+2]);
      $i+=2;
      
    } case (/^-t/ or /^--types/) { # include types
      print "removing atom types\n";
      @include=();
      while(not $ARGV[$i+1] =~ /^-/ and $i<$#ARGV) {
	$i++;
	push(@include,uc($ARGV[$i]));
      }
      for($t=0;$t<@{$cdata[0]};$t++) {
	for($m=0;$m<@{$cdata[0][$t]};$m++) {
	  loopatom:for($a=0;$a<@{$cdata[0][$t][$m]};$a++) {
	    for($c=0;$c<@include;$c++) {
	      if($cdata[0][$t][$m][$a][9] =~ /^$include[$c]$/i) {
		$remove=0;
		next loopatom
	      }
	    }
	    $frame_numatoms[0]--;
	    splice(@{$cdata[0][$t][$m]},$a,1);
	    $a--;
	  }
	}
      }
      
    } case (/^-u$/ or /^--rotat.*$/) { # rotate system
      print "rotating system\n";
      while(not $ARGV[$i+1] =~ /^-/ and $i<$#ARGV) {
	$i++;
	if($ARGV[$i] =~ /^x$/i) {
	  $c=0;
	} elsif($ARGV[$i] =~ /^y$/i) {
	  $c=1;
	} elsif($ARGV[$i] =~ /^z$/i) {
	  $c=2;
	}
	$i++;
	@rotmatrix = gen_rot_matrix(\@{$axisvecs[$c]},$ARGV[$i]/180.0*pi);
	for($t=0;$t<@{$cdata[0]};$t++) {
	  for($m=0;$m<@{$cdata[0][$t]};$m++) {
	    rotate_molecule(\@{$cdata[0][$t][$m]}, \@rotmatrix, \@rotcenter);
	  }
	}
	if($periodic_key[0]>0) {
	  @{$cell[0][0]} = rotate_vector(\@{$cell[0][0]}, \@rotmatrix);
	  @{$cell[0][1]} = rotate_vector(\@{$cell[0][1]}, \@rotmatrix);
	  @{$cell[0][2]} = rotate_vector(\@{$cell[0][2]}, \@rotmatrix);
	}
      }
    } case (/^-uc$/ or /^--rot.*center/) { # rotation center
      print "setting rotation center\n";
      @rotcenter = @ARGV[$i+1..$i+3];
      $i += 3;
      
    } case (/^-w$/ or /^--wrap/) {
      print "wrapping molecules back into cell\n";
      if($periodic_key[0]==6) {
	@wrap=(0,1);
      } elsif($periodic_key[0]==1 or $periodic_key[0]==2) {
	@wrap=(0,1,2);
      } elsif($periodic_key[0]>2) {
	print "**** error: only implemented for orthogonal cells!\n";exit 1;
      }
      for($t=0;$t<@{$cdata[0]};$t++) {
	for($m=0;$m<@{$cdata[0][$t]};$m++) {
	  for($a=0;$a<@{$cdata[0][$t][$m]};$a++) {
	    foreach $c (@wrap) {
	      $cdata[0][$t][$m][$a][$c] -= int($cdata[0][$t][$m][$a][$c]/$size[0][$c])*$size[0][$c]*2.0;
	    }
	  }
	}
      }
    } case (/^-x$/ or /^--extract/) { # extract single molecule
      $te=$ARGV[$i+1];
      $me=$ARGV[$i+2];
      if($te>=$field_nummols[0]) {
	print "**** error: only $field_nummols[0] molecules in FIELD file!\n";
	exit 1;
      }
      if($me>=$mol_numents[0][$te]) {
	print "**** error: only $mol_numents[0][$te] entities for molecule $mol_name[0][$te]!\n";
	exit 1;
      }
      for($t=0;$t<@{$cdata[$ci]};$t++) {
	for($m=$#{$cdata[$ci][$t]};$m>=0;$m--) {
	  next if($t==$te and $m==$me);
	  remove_mol_entity(0,0,$t,$m);
	}
      }
    } case (/^-y$/ or /^--multiply/) { # multiply system
      print "multiplying system\n";
      @scalefactors = @ARGV[$i+1..$i+3];
      $i+=3;
      &enlarge_system(0,0, 0,0, @scalefactors, "rms", \@spacing);
    
    } case (/^-ys$/ or /^--spacing/) { # spacing for system multiplication
      print "setting spacing\n";
      @spacing = @ARGV[$i+1..$i+3];
      $i += 3;
      
    }
  }
}


######## write output #######################################

write_field_file($outfld,0) if($loutfield);
if($outformat =~ /^hist$/i) { # output HISTORY file
  if(not open($fhhist,$writemode,$outcfg)) {
    print "**** error: cannot open HISTORY file $outcfg: $!\n";
    exit 1;
  }
  exit 1 if (write_history_timestep($fhhist,-1,0) != 0);
  close($fhhist);
} elsif($outformat =~ /^xyz$/i) { # output xyz file
  if(not open($fhxyz,$writemode,$outcfg)) {
    print "**** error: cannot open XYZ file $outcfg: $!\n";
    exit 1;
  }
  exit 1 if (write_xyz_timestep($fhxyz,0,$config_title[0]) != 0);
  close($fhxyz);
} else { #output CONFIG file
  write_config_file($outcfg,0,$config_title[0],$writemode);
}
