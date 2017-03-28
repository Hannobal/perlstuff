#!/usr/bin/perl

use dlpoly_utility;
use lammps_utility;
use hanno_utility;
use POSIX;

# extracts the information from dcd files distributed
# over several folders with numbered names

use Cwd;
use Switch;

$endframe   = 9e20;
$startframe = 0;
$histresdens   = 0.25;
$heatup     = 0;
$splitjob   = 0;
$stepframe  = 1;
$outdir     = "analysis";
$lbulkbound = (0,0,0);
$axis = ('x','y','z');
$timestep   = 0.001;

if($#ARGV<2) {
  print "the input format is:\n";
  print "1. FIELD file\n";
  print "2. name of dcd file in each directory\n";
  print "-t <real>    timestep (default: $timestep)\n";
  print "-f <3*int>   start, end and step frame (default: $startframe, $endframe, 1)\n";
  print "-o <str>     output-directory (default: $outdir)\n";
  print "-s           toggle split job mode\n";
  print "-bx <2*real> minimum and maximum x-value for calculation of bulk-properties\n";
  print "-by <2*real> minimum and maximum y-value for calculation of bulk-properties\n";
  print "-bz <2*real> minimum and maximum z-value for calculation of bulk-properties\n";
#   print "-c60         perform analysis of C60 cages\n";
  print "-h           toggle heatup job mode\n";
  print "-d           analyze bulk density\n";
  print "-dn          include only specified molecule names for density calculation\n";
  print "-dh          create density histogram for each step\n";
  print "-dhav        create time-averaged density histogram\n";
  print "-dr <real>   resolution of z-profile histogram (default: $histresdens)\n";
  print "-v <n*3*int> analyze orientation of vector between atoms in molecule\n";
  print "             format: n*(<mol_name> <atomid 1> <atomid 2>) atom indices starting from 1\n";
  print "-vz <2*real> minimum and maximum z-value for analysis of orientation vectors\n";
  exit 1;
}

$infld=$ARGV[0];
$dcdname=$ARGV[1];

@denstypes=();
for($i=2;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    # general flags
    case ("-o") {
      $i++;
      $outdir = $ARGV[$i];
    } case ("-s") {
      $splitjob = 1;
    } case ("-dn") {
      $lmolnames = 1;
      while(not $ARGV[$i+1] =~ /^-/ and $i<$#ARGV) {
	$i++;
	push(@denstypes,$ARGV[$i]);
      }
    } case ("-h") {
      $splitjob = 1;
      $heatup   = 1;
    } case ("-f") {
      $i++;
      $startframe = $ARGV[$i];
      $i++;
      $endframe   = $ARGV[$i];
      $i++;
      $stepframe  = $ARGV[$i];
      check_integer($startframe) or die "**** error: $startframe is not an integer number!\n";
      check_integer($endframe)   or die "**** error: $endframe is not an integer number!\n";
      die "**** error: end frame must not be smaller than start frame\n" if($endframe<$startframe);
      die "**** error: frame step must be greater than zero\n" if($stepframe<=0);
    } case ("-bx") {
      $lbulkbound[0] = 1;
      $i++;
      $bulkmin[0] = $ARGV[$i];
      $i++;
      $bulkmax[0] = $ARGV[$i];
      check_real($bulkmin[0]) or die "**** error: $bulkmin[0] is not a real number!\n";
      check_real($bulkmax[0]) or die "**** error: $bulkmax[0] is not a real number!\n";
      die "**** error: max x-value must not be smaller than min x-value for bulk-boundaries\n"
        if($bulkmax[0]<$bulkmin[0]);
    } case ("-by") {
      $lbulkbound[1] = 1;
      $i++;
      $bulkmin[1] = $ARGV[$i];
      $i++;
      $bulkmax[1] = $ARGV[$i];
      check_real($bulkmin[1]) or die "**** error: $bulkmin[1] is not a real number!\n";
      check_real($bulkmax[1]) or die "**** error: $bulkmax[1] is not a real number!\n";
      die "**** error: max y-value must not be smaller than min y-value for bulk-boundaries\n"
        if($bulkmax[1]<$bulkmin[1]);
    } case ("-bz") {
      $lbulkbound[2] = 1;
      $i++;
      $bulkmin[2] = $ARGV[$i];
      $i++;
      $bulkmax[2] = $ARGV[$i];
      check_real($bulkmin[2]) or die "**** error: $bulkmin[2] is not a real number!\n";
      check_real($bulkmax[2]) or die "**** error: $bulkmax[2] is not a real number!\n";
      die "**** error: max z-value must not be smaller than min z-value for bulk-boundaries\n"
        if($bulkmax[2]<$bulkmin[2]);
    } case ("-c60") {
      $lfull=1;
    
    # flags for density calculation
    } case ("-d") {
      $ldensbulk = 1;
    } case ("-dh") {
      $ldenshist = 1;
    } case ("-dhav") {
      $ldenshistav = 1;
    } case ("-dr") {
      $i++;
      $histresdens = $ARGV[$i];
      check_real($histresdens) or die "**** error: $histresdens is not a real number!\n";
      die "**** error: histogram step must be greater than zero\n" if($histresdens<=0);
      
    # flag for calculation of orientation vector
    } case ("-vz") {
      $lovecz = 1;
      $i++;
      $oveczmin = $ARGV[$i];
      $i++;
      $oveczmax = $ARGV[$i];
      check_real($oveczmin) or die "**** error: $oveczmin is not a real number!\n";
      check_real($oveczmax) or die "**** error: $oveczmax is not a real number!\n";
      die "**** error: max z-value must not be smaller than min z-value for flag '-vz'\n"
        if($oveczmin>$oveczmax);
    } case ("-v") {
      die "**** error: too few arguments for flag '-v'\n" if($#ARGV<$i+3);
      while($ARGV[$i+1] !~ /^-/ and $#ARGV>=$i+3) {
	($a,$b,$c) = @ARGV[$i+1..$i+3];
	die "**** error: expected molecule name but got flag!\n" if($a =~ /^-/);
	check_integer($b) or die "**** error: expected integer but got '$b'!\n";
	check_integer($c) or die "**** error: expected integer but got '$c'!\n";
	die "**** error: atom ids must be greater than zero\n" if($b<1 or $c<1);
	push(@orientvecs,[-1,$b-1,$c-1,$a]);
	$i+=3;
      }
    } else {
      print "**** error: expected optional flag but got \"$ARGV[$i]\"!\n";
      exit 1;
    }
  } # end switch
} # end for

############### find the relevant directories to read #########################
if(not $splitjob) {
  @directories = (".")
} else {
  @directories = get_numeric_directories(".",$startframe,$endframe);
  push(@directories,'.') if(-e "./$dcdname");
}
die "**** error: no numeric directories and no dcd file in work directory were found" if($#directories<0);

############### read the FIELD file from the first directory ##################
exit 1 if(read_field_file($infld,0)!=0);
for($i=0;$i<@orientvecs;$i++) {
  for($t=0;$t<$field_nummols[0];$t++) {
    $orientvecs[$i][0]=$t if($mol_name[0][$t] =~ $orientvecs[$i][3]);
    last;
  }
  die "**** error: did not find molecule named $orientvecs[$i][3]!\n" if($orientvecs[$i][0]<0);
}
findtypes:for($i=0;$i<@denstypes;$i++) {
  for($t=0;$t<$field_nummols[0];$t++) {
    if($mol_name[0][$t]=~/^$denstypes[$i]$/) {
      $denstypes[$i]=$t;
      next findtypes;
    }
  }
  die "**** error: did not find molecule named $denstypes[$i]!\n";
}
############### generate the output files #####################################
if(not -d $outdir) {
  mkdir $outdir or die "cannot create output-directory: $!";
}

# for density calculation
if($ldensbulk) {
  open(DENS_AV, ">", "$outdir/DENS_AV") or die "Can't create file DENS_AV: $!";
  if(contains(@bulkbound,1)) {
    print DENS_AV "# density between ";
    print DENS_AV "$bulkmin[0] and $bulkmass[0] A on the x-axis " if($lbulkbound[0]);
    print DENS_AV "$bulkmin[1] and $bulkmass[1] A on the y-axis " if($lbulkbound[1]);
    print DENS_AV "$bulkmin[2] and $bulkmass[2] A on the z-axis " if($lbulkbound[2]);
    print DENS_AV "\n";
    print DENS_AV "of molecules",join(" ",@denstypes),"\n" if(@denstypes);
  }
  printf DENS_AV "%14s%14s","time_[ps]","dens_[g/cm^3]";
}
if($ldenshist) {
  open(DENS_HISTZ, ">", "$outdir/DENS_HISTZ") or die "Can't create file DENS_HISTZ: $!";
  print DENS_HISTZ "# density from molecules",join(" ",@denstypes),"\n" if(@denstypes);
}
for($i=0;$i<@orientvecs;$i++) {
  $fn="$outdir/OVEC_$orientvecs[$i][3]_".($orientvecs[$i][1]+1)."-".($orientvecs[$i][2]+1);
  open($ovecfile[$i], ">", $fn) or die "Can't create file $fn: $!";
  $file=$ovecfile[$i];
  print $file "# average vector between atoms ",$orientvecs[$i][1]+1," (",
  $mol_atomdata[0][$orientvecs[$i][0]][$orientvecs[$i][1]][0],") and ",$orientvecs[$i][2]+1," (",
  $mol_atomdata[0][$orientvecs[$i][0]][$orientvecs[$i][2]][0],")\n";
  if($lovecz) {
    print $file "# between z=$oveczmin and $oveczmax\n";
  } else {
    print $file "# for all molecules\n";
  }
  printf $file "\n#%13s%14s%14s%14s%10s","time_[ps]","vec_x","vec_y","vec_z","numvecs for average";
  if($lovecz) {
    print $file "# within z-range $oveczmin A and $oveczmax A\n";
  }
}

############### define some constants in advance ###############################
use constant pi  => 3.14159265358979323846;
use constant eps => 8.854187817620e-12;
use constant ec  => 1.602176565e-19;
use constant l0  => 1e-10;
use constant V0  => 1e-30;
use constant t0  => 1e-12;
use constant u   => 1.660538921e-27;
use constant V1  => 1.0e-24; # to convert nm^-3 to cm^-3
use constant d0  => u*1000/V1; # to convert density

############### read the dcd file from each directory ##########################
$time=0;
$time_old=0;
$starttime=0;
$numframes=0;
foreach $dir (@directories) {
  print "reading directory $dir\n";
  
  open($fhhist, "<", "$dir/$dcdname") or print "**** warning: file $dcdname was not found in $dir\n";
  $starttime=$time if($heatup);
  whilereadhist : while(1) {
    last if(read_dcd_timestep($fhhist,0,0) != 0);
    last if($frame_number[0]>$endframe and not $heatup);
    $time  = $frame_number[0]*$timestep;
    $time += $starttime;
    if($time<$time_old) {
      $starttime = $time_old;
      $time += $starttime;
    }
    $numframes++;
    
    if($periodic_key[0]==6) {
      @remapaxis=(0,1);
    } elsif($periodic_key[0]>0) {
      @remapaxis=(0,1,2);
    } else {
      @remapaxis=();
    }
    
    ############### analyze density ############################################
    if($ldensbulk) {
      # calc bulk volume
      $bulkvol=1;
      for($c=0;$c<3;$c++) {
	if($lbulkbound[$c]) {
	  $bulkvol *= $bulkmax[$c]-$bulkmin[$c];
	} else {
	  $bulkvol *= $size[0][$c]*2;
	}
      }
      # calc bulk mass
      $bulkmass=0;
      for($t=0;$t<$field_nummols[0];$t++) {
	if(@denstypes) {
	  next if(not contains(@denstypes,$t));
	}
	for($m=0;$m<$mol_numents[0][$t];$m++) {
	  densatom : for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	    for($c=0;$c<3;$c++) {
	      if($lbulkbound[$c]) {
		next densatom if($cdata[0][$t][$m][$a][$c]<$bulkmin[$c]
		  or $cdata[0][$t][$m][$a][$c]>$bulkmax[$c])
	      }
	    }
	    $bulkmass += $mol_atomdata[0][$t][$a][1]
	  } # end for $a
	} # end for $m
      } # end for $t
      $densbulk = d0*$bulkmass/$bulkvol;
      printf DENS_AV "\n%14.4f%14.8f",$time,$densbulk;
    } # end if($ldensbulk)
    
    if($ldenshist or $ldenshistav) {
      %denshist = ();
      for($t=0;$t<$field_nummols[0];$t++) {
	if(@denstypes) {
	  next if(not contains(@denstypes,$t));
	}
	for($m=0;$m<$mol_numents[0][$t];$m++) {
	  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
	    $i = floor($cdata[0][$t][$m][$a][2]/$histresdens);
	    $denshist{$i} += $mol_atomdata[0][$t][$a][1];
	  } # end for $a
	} # end for $m
      } # end for $t
      $denshistvol = $histresdens*$size[0][0]*$size[0][1]*4;
      print DENS_HISTZ "\n\n\n# histogram for $time ps";
      printf DENS_HISTZ "\n%14s%14s", "zpos_[A]","dens_[g/cm^3]";
      foreach $i (sort {$a<=> $b} keys %denshist) {
	$denshist{$i} *= d0/$denshistvol;
	printf DENS_HISTZ "\n%14.6f%14.8f",($i*$histresdens),$denshist{$i};
	if($ldenshistav) {
	  $denshistav{$i} += $denshist{$i};
	}
      }
    }
    
    ############### analyze orientation vectors ################################
    
    for($i=0;$i<@orientvecs;$i++) {
      $file=$ovecfile[$i];
      $numvecs=0;
      @resvec=(0,0,0);
      ($t,$a,$b)=@{$orientvecs[$i]};
      for($m=0;$m<$mol_numents[0][$t];$m++) {
	# check boundaries if necessary
	if($lovecz) {
	  next if($cdata[0][$t][$m][$a][2]>$oveczmax);
	  next if($cdata[0][$t][$m][$a][2]<$oveczmin);
	  next if($cdata[0][$t][$m][$b][2]>$oveczmax);
	  next if($cdata[0][$t][$m][$b][2]<$oveczmin);
	}
	$numvecs++;
	remap_molecule(\@{$cdata[0][$t][$m]},\@remapaxes,\@{$size[0]}) if(@remapaxis);
	$vec[0] = $cdata[0][$t][$m][$b][0] - $cdata[0][$t][$m][$a][0];
	$vec[1] = $cdata[0][$t][$m][$b][1] - $cdata[0][$t][$m][$a][1];
	$vec[2] = $cdata[0][$t][$m][$b][2] - $cdata[0][$t][$m][$a][2];
	$veclen = sqrt($vec[0]*$vec[0]+$vec[1]*$vec[1]+$vec[2]*$vec[2]);
	$resvec[0] += $vec[0]/$veclen;
	$resvec[1] += $vec[1]/$veclen;
	$resvec[2] += $vec[2]/$veclen;
      }
      next if($numvecs==0);
      $resvec[0] /= $numvecs;
      $resvec[1] /= $numvecs;
      $resvec[2] /= $numvecs;
      printf $file "\n%14.4f%14.8f%14.8f%14.8f%10u",$time,@resvec,$numvecs;
    }
    
    ############### analyze fullerene positions ################################
    
#     $tot_full_mindist_bottom[$numframes] = 9e20;
#     $b = 0; #buckyball index
#     for($t=0;$t<$field_nummols[0];$t++) {
#       next if(not $mol_name[0][$t] =~ /C60/)
#       for($m=0;$m<$mol_numents[0][$t];$m++) {
# 	$fullbottom[$numframes][$b]=9e20;
# 	for($a=0;$a<$mol_numatoms[0][$t];$a++) {
# 	  next if(not $mol_atomdata[0][$t][$a] =~ /CA/);
# 	  $fullcenter[$numframes][$b][0] += $cdata[0][$t][$m][$a][0];
# 	  $fullcenter[$numframes][$b][1] += $cdata[0][$t][$m][$a][1];
# 	  $fullcenter[$numframes][$b][2] += $cdata[0][$t][$m][$a][2];
# 	}
#       }
#     }
    
    $time_old = $time;
    
    for($i=1;$i<$stepframe;$i++) {
      last whilereadhist if (not read_lammpstrj_timestep($fhhist)==0);
    }
  } # end while
  close_dcd_file($fhhist);
} # end for directories
# write average density histogram
if($ldenshistav) {
print "asdf\n";
  open(DENS_HISTZ_AV, ">", "$outdir/DENS_HISTZ_AV") or die "Can't create file DENS_HISTZ_AV: $!";
  print DENS_HISTZ_AV "# time averaged density histogram";
  print DENS_HISTZ_AV "# $numframes analyzed\n";
  print DENS_HISTZ_AV "# density from molecules ",join(" ",@denstypes),"\n" if(@denstypes);
  foreach $i (sort {$a <=> $b} keys %denshistav) {
    printf DENS_HISTZ_AV "\n%14.6f%14.8f",($i*$histresdens),$denshistav{$i}/$numframes;
  }
  close(DENS_HISTZ_AV)
}

close(DENS_AV) if($ldensbulk);
close(DENS_HISTZ) if($ldenshist);
close(ORIENTVEC) if(@orientvecs);
for($i=0;$i<@orientvecs;$i++) {
  close($ovecfile[$i]);
}
