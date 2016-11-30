#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;
use Math::Trig;

if($#ARGV<2) {
  print " 1. name of output folder\n";
  print " 2. start frame and step\n";
  print "[3. end frame]\n";
  exit 1;
}

$ang_poh=pi/3;

$outdir     = $ARGV[0];
$startframe = $ARGV[1];
$stepframe  = $ARGV[2];
if($#ARGV>2) {
  $endframe = $ARGV[3];
} else {
  $endframe = 9e20;
}


if(not check_integer($startframe,$stepframe,$endframe)) {
  print "**** error: startframe and step must be integer values\n";
  exit 1;
}
@directories = get_numeric_directories(".",$startframe,$endframe);
push(@directories,'.') if(-e './HISTORY');

if(not -d $outdir) {
  mkdir $outdir or die "cannot create output-directory: $!";
}

foreach $dir (@directories) {
  print "reading directory $dir\n";
  exit 1 if(read_field_file("$dir/FIELD",0)!=0);
  open($fhhist,"<","$dir/HISTORY") or die "**** error: cannot open $dir/HISTORY: $!\n";
  while(read_history_timestep($fhhist,0,0)==0) {
    last if($frame_number[0]>$endframe);
    next if($frame_number[0]<$startframe);
    next if(($frame_number[0]-$startframe)%$stepframe!=0);
    if($periodic_key[0]==6) {
      @remapaxis=(0,1);
    } elsif($periodic_key[0]>0) {
      @remapaxis=(0,1,2);
    } else {
      @remapaxis=();
    }
    for($t=0;$t<$field_nummols[0];$t++) {
      if(not $mol_name[0][$t] =~/PA/) {
	$frame_numatoms[0] -= $mol_numents[0][$t]*$mol_numatoms[0][$t];
	@{$cdata[0][$t]} = ();
# 	$mol_numents[0][$t] = 0;
	next;
      }
      # remap molecules
      if($mol_name[0][$t] =~/C60/) { $full=1; } else { $full=0; }
      for($m=0;$m<$mol_numents[0][$t];$m++) {
	@numatoms=(0,0,0,0); # the four quadrants
	@remapatoms=(undef,undef,undef,undef);
	for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	  next if($full and $mol_atomdata[0][$t][$a][0] ne 'CA');
	  if($cdata[0][$t][$m][$a][0]>0) {
	    if($cdata[0][$t][$m][$a][1]>0) {
	      $remapatoms[0]=$a;
	      $numatoms[0]++
	    } else {
	      $remapatoms[1]=$a;
	      $numatoms[1]++
	    }
	  } else {
	    if($cdata[0][$t][$m][$a][1]>0) {
	      $remapatoms[2]=$a;
	      $numatoms[2]++
	    } else {
	      $remapatoms[3]=$a;
	      $numatoms[3]++
	    }
	  }
	}
	$maxnumatoms=0;
	for($c=0;$c<4;$c++) {
	  if($maxnumatoms<$numatoms[$c]) {
	    $cm          = $c;
	    $maxnumatoms = $numatoms[$c];
	  }
	}
	remap_molecule(\@{$cdata[0][$t][$m]},\@remapaxis,\@{$size[0]},$remapatoms[$cm])
	if($maxnumatoms!=$mol_numatoms[0][$t]);
      }
      # add the hydrogen atom
      if($mol_name[0][$t] =~ /PA-d/) {
	# add a hydrogen atom
	for($m=0;$m<$mol_numents[0][$t];$m++) {
	  $minz = 9e20;
	  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	    if($cdata[0][$t][$m][$a][2]<$minz) {
	      $mina=$a;
	      $minz=$cdata[0][$t][$m][$a][2];
	    }
	  }
	  $dx = $cdata[0][$t][$m][$mina][0]-$cdata[0][$t][$m][$mol_bondatoms[0][$t][$mina][0]][0];
	  $dy = $cdata[0][$t][$m][$mina][1]-$cdata[0][$t][$m][$mol_bondatoms[0][$t][$mina][0]][1];
	  $dz = $cdata[0][$t][$m][$mina][2]-$cdata[0][$t][$m][$mol_bondatoms[0][$t][$mina][0]][2];
	  @bondaxis  = ($dx,$dy,$dz);
	  @bondnorm  = normalize_vector(\@bondaxis);
	  @test      = vector_product(\@bondnorm,[1,0,0]);
	  @rotmatrix = gen_rot_matrix(\@test,$ang_poh);
	  @ohvec     = rotate_vector(\@bondnorm,\@rotmatrix);
	  $angle     = rand(2*pi);
	  @rotmatrix = gen_rot_matrix(\@bondnorm,$angle);
	  @ohvec     = rotate_vector(\@ohvec,\@rotmatrix);
	  @rotmatrix = gen_rot_matrix(\@bondaxis,$angle);
	  $x = $cdata[0][$t][$m][$mina][0]+$ohvec[0];
	  $y = $cdata[0][$t][$m][$mina][1]+$ohvec[1];
	  $z = $cdata[0][$t][$m][$mina][2]+$ohvec[2];
	  push(@{$cdata[0][$t][$m]},[$x,$y,$z,0,0,0,0,0,0,"HG"]);
	  $frame_numatoms[0]++;
	}
      }
    }
    $fn = sprintf("%08u",$frame_number[0]);
    open($fhxyz,">","$outdir/$fn.xyz");
    write_xyz_timestep($fhxyz,0,"");
    close($fhxyz);
  }
  close($fhhist);
}
