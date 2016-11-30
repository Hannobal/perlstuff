#!/usr/bin/perl

$|=1;

use hanno_utility;
use dlpoly_utility;
use List::Util 'shuffle';
use Storable qw(dclone);
use Switch;
use Math::Trig;
use IO::Handle;

$pihalf     = pi/2.0;
$cancel     = 1;
$lessoutput = 0;
$histname   = "HISTORY";
$fldname    = "FIELD";
$outdir     = $ARGV[0];
$histresxy  =  2;
$histresdeg = 10;
$startframe =  0;
$endframe   =  9e20;
$tolerance  =  0.1;
$radiuscutoff=9999;
$langular   = 0;
$lspatial   = 0;

if($#ARGV<0) {
  print "input format:\n";
  print " 1. output directory for the output files\n";
  print "-rs <real>   spatial histogram resolution (default: no spatial analysis)\n";
  print "-ra <real>   angular histogram resolution (default: no angular analysis)\n";
  print "-c <real     cutoff around center (default: no cutoff)\n";
  print "-t tolerance for removing inversion symmetric vibrations (default: $tolerance)]\n";
  print "-s <int>     start frame\n";
  print "-e <int>     end frame\n";
  print "-n n*<str>   list of molecules to analyze\n";
  print "-nocancel\n";
  print "-lessoutput\n";
  exit;
}

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    # general flags
    case (/^-rs/i) {
      $i++;
      $histresxy = $ARGV[$i];
      $lspatial=1;
      &error(6) if($histresxy<=0);
    } case (/^-ra/i) {
      $i++;
      $histresangdeg = $ARGV[$i];
      &error(6) if($histresdeg<=0);
      &error(9) if(360%$histresdeg!=0);
      $langular=1;
    } case (/^-t/i) {
      $i++;
      $tolerance = $ARGV[$i];
      &error(6) if($histresxy<=0);
    } case (/^-c/i) {
      $i++;
      $radiuscutoff = $ARGV[$i];
      &error(6) if($radiuscutoff<=0);
    } case(/^-s/i) {
      $i++;
      $startframe = $ARGV[$i];
      &error(1) if not check_integer($startframe);
    } case(/^-e/i) {
      $i++;
      $endframe = $ARGV[$i];
      &error(2) if not check_integer($endframe);
    } case(/^-n/i) {
      $i++;
      until($ARGV[$i]=~/^-/ or $i>$#ARGV) {
	push(@includetypes,$ARGV[$i]);
	$i++;
      }
      $i--;
    } case(/^-lessoutput/i) {
      $lessoutput=1;
    } case(/^-nocancel/i) {
      $cancel=0;
    } else {
      &error(10,$ARGV[$i]);
    }
  }
}

print "xy histogram resolution:      $histresxy\n" if($lspatial);
print "angular histogram resolution: $histresdeg\n" if($langular);
print "start frame: $startframe\n";
print "end   frame: $endframe\n";
print "tolerance:   $tolerance\n";

$histresrad = $histresdeg/180*pi;

@directories=get_numeric_directories(".");
if(@directories) {
  @directories=get_numeric_directories(".",$startframe,$endframe);
}
push(@directories,".") if(-e "./$histname");
&error(7) if(not @directories);

exit 1 if(read_field_file("$directories[0]/$fldname",0)!=0);

if(not @includetypes) {
  @includetypes = @{$mol_name[0]};
  @tlst         = (0..$#{$mol_name[0]});
} else {
  for($i=0;$i<@includetypes;$i++) {
    $tlst[$i] = -1;
    for($t=0;$t<$field_nummols[0];$t++) {
      if($mol_name[0][$t]=~/^$includetypes[$i]$/i) {
	$tlst[$i] = $t;
	last;
      }
    }
    &error(11,$includetypes[$i]) if($tlst[$i]==-1);
  }
}

foreach $t (@tlst) {
  for($a=0;$a<$mol_numatoms[0][$t];$a++) {
    if($mol_atomdata[0][$t][$a][0]=~/C3/) {
      $nh=0;
      @group=();
      $group[0]=$a;
      for($i=0;$i<@{$mol_bondatoms[0][$t][$a]};$i++) {
	$b=$mol_bondatoms[0][$t][$a][$i];
	if($mol_atomdata[0][$t][$b][0]=~/^HC/) {
	  $nh++;
	  $group[$nh]=$b;
	}
      }
      if($nh==2) {
	push(@ch2s,[$group[0],$group[1],$group[2],$t])
      } elsif($nh==3) {
	push(@ch3s,[$group[0],$group[1],$group[2],$group[3],$t])  
      } else {
	&error(3,$nh);
      }
    }
  }
}

mkdir($outdir) if(not -d $outdir);

if($cancel) {
  &error(4) if(not open($fhavch2z,">","$outdir/AV_CH2_Z_CANCELLED"));
  &error(4) if(not open($fhavch3z,">","$outdir/AV_CH3_Z_CANCELLED"));
  &error(4) if(not open($fhch2vecs,">","$outdir/CH2_VECTOR_DIRECTION.xyz"));
  &error(4) if(not open($fhch3vecs,">","$outdir/CH3_VECTOR_DIRECTION.xyz"));
  &error(4) if(not open($fhch2vecsdel,">","$outdir/CH2_VECTOR_DIRECTION_cancelled.xyz"));
  &error(4) if(not open($fhch3vecsdel,">","$outdir/CH3_VECTOR_DIRECTION_cancelled.xyz"));
  &error(4) if(not open($fhch2pos,">","$outdir/CH2pos_cancelled.xyz"));
} else {
  &error(4) if(not open($fhavch2z,">","$outdir/AV_CH2_Z"));
  &error(4) if(not open($fhavch3z,">","$outdir/AV_CH3_Z"));
}
print $fhavch2z "# average z-component of CH2 vectors within $radiuscutoff angstroms from center\n";
printf $fhavch2z "#%9s %17s %17s %7s\n","frame","accum. values [A]","area [nm2]","# values";
print $fhavch3z "# average z-component of CH3 vectors within $radiuscutoff angstroms from center\n";
printf $fhavch3z "#%9s %17s %17s %7s\n","frame","accum. values [A]","area [nm2]","# values";

if($lessoutput) {
  print "begin analysis...\n";
} else {
  print "\n";
}

$firsthist=1;
$f=0;
foreach $dir (@directories) {
  exit 1 if(not open($fhhist,"<","$dir/$histname"));
  if($startframe>0 and $firsthist) {
    &error(5,$startframe) if(not find_history_timestep($fhhist,$startframe));
  }
  $firsthist=0;
  $err=0;
  while($err==0) {
    $err=read_history_timestep($fhhist,0,0);
    last   if($err<0);
    exit 1 if($err>0);
    last if($frame_number[0]>$endframe);
    print "\ranalyzing timestep $frame_number[0]" unless($lessoutput);
    $fn[$f]=$frame_number[0];
    $avch2z[$f]=0;
    $avch3z[$f]=0;
    $jch2=0;
    $jch3=0;
    $numch2cut=0;
    $numch3cut=0;
    # indices :  0  1  2   3   4     5        6         7       8       9        10
    #           dx dy dz dxy phi theta molindex atomindex mindist partner typeindex
    @valch2=();
    @valch3=();
    $minsize = min($size[0][0],$size[0][1],$radiuscutoff);
    $areatot = calc_area($radiuscutoff,0,@{$size[0]})/100;
    for($i=0;$i<@ch2s;$i++) {
      $t=$ch2s[$i][3];
      for($m=0;$m<$mol_numents[0][$t];$m++) {
	@cpos  = @{$cdata[0][$t][$m][$ch2s[$i][0]]}[0..2];
	@h1pos = @{$cdata[0][$t][$m][$ch2s[$i][1]]}[0..2];
	@h2pos = @{$cdata[0][$t][$m][$ch2s[$i][2]]}[0..2];
	@h1pos = unwrap(\@h1pos,\@cpos);
	@h2pos = unwrap(\@h2pos,\@cpos);
	$d[0]  = $h1pos[0] + $h2pos[0] - 2*$cpos[0];
	$d[1]  = $h1pos[1] + $h2pos[1] - 2*$cpos[1];
	$d[2]  = $h1pos[2] + $h2pos[2] - 2*$cpos[2];
	@d   = normalize_vector(\@d);
	$dxy = sqrt($d[0]*$d[0]+$d[1]*$d[1]);
	$phi    = atan2($d[1],$d[0]);
	$phi   += 2*pi if($phi<0);
	$theta  = $pihalf-atan2($d[2],$dxy);
	@{$valch2[$jch2]} = (@d,$dxy,$phi,$theta,$m,$ch2s[$i][0],2,undef,$t);
	$jch2++;
      }
    }
    for($i=0;$i<@ch3s;$i++) {
      $t=$ch3s[$i][4];
      for($m=0;$m<$mol_numents[0][$t];$m++) {
	@cpos  = @{$cdata[0][$t][$m][$ch3s[$i][0]]}[0..2];
	@h1pos = @{$cdata[0][$t][$m][$ch3s[$i][1]]}[0..2];
	@h2pos = @{$cdata[0][$t][$m][$ch3s[$i][2]]}[0..2];
	@h3pos = @{$cdata[0][$t][$m][$ch3s[$i][3]]}[0..2];
	@h1pos = unwrap(\@h1pos,\@cpos);
	@h2pos = unwrap(\@h2pos,\@cpos);
	@h3pos = unwrap(\@h3pos,\@cpos);
	$d[0]  = $h1pos[0] + $h2pos[0] + $h3pos[0] - 3*$cpos[0];
	$d[1]  = $h1pos[1] + $h2pos[1] + $h3pos[1] - 3*$cpos[1];
	$d[2]  = $h1pos[2] + $h2pos[2] + $h3pos[2] - 3*$cpos[2];
	@d   = normalize_vector(\@d);
	$dxy = sqrt($d[0]*$d[0]+$d[1]*$d[1]);
	$phi    = atan2($d[1],$d[0]);
	$phi   += 2*pi if($phi<0);
	$theta  = $pihalf-atan2($dxy,$d[2]);
	@{$valch3[$jch3]} = (@d,$dxy,$phi,$theta,$m,$ch3s[$i][0],2,undef,$t);
	$jch3++;
      }
    }
    if($cancel) {
      # output points before removal
      $numch2=$jch2;
      print $fhch2vecs "$numch2\n\n";
      for($i=0;$i<$numch2;$i++) {
	printf $fhch2vecs "H %10.7f %10.7f %10.7f\n",@{$valch2[$i]}[0..2];
      }
      $numch3=$jch3;
      print $fhch3vecs "$numch3\n\n";
      for($i=0;$i<$numch3;$i++) {
	printf $fhch3vecs "H %10.7f %10.7f %10.7f\n",@{$valch3[$i]}[0..2];
      }
      
      #remove centrosymmetric ch2 vibrations
      @loopch1 = shuffle(0..$#valch2);
      @loopch2 = shuffle(0..$#valch2);
      @usedch2 = ();
      @remch2  = ();
      $numiter=0;
      do {
	$newconn=0;
	foreach $i (@loopch1) {
	  next if(defined($valch2[$i][9]));
	  foreach $j (@loopch2) {
	    next if($i==$j);
	    next if(defined($valch2[$j][9]));
	    $sum[0] = $valch2[$i][0] + $valch2[$j][0];
	    $sum[1] = $valch2[$i][1] + $valch2[$j][1];
	    $sum[2] = $valch2[$i][2] + $valch2[$j][2];
	    $r = vector_length(@sum);
	    next if($r>$tolerance);
	    if($r<$valch2[$i][8] and $r<$valch2[$j][8]) {
	      $newconn++;
	      $valch2[$valch2[$i][9]][9] = undef; # remove previous partner;
	      $valch2[$valch2[$j][9]][9] = undef; # remove previous partner;
	      $valch2[$i][8] = $r; # reset minimum distance
	      $valch2[$j][8] = $r; # reset minimum distance
	      $valch2[$i][9] = $j; # set new partner
	      $valch2[$j][9] = $i; # set new partner
	    }
	  }
	}
	for($i=0;$i<@valch2;$i++) {
	  $valch2[$i][8]=2 if(not defined($valch2[$i][9]));
	}
      } until($newconn==0);
      for($i=0;$i<@valch2;$i++) {
	if(not defined($valch2[$i][9])) { # if vibration has no partner
	  push(@usedch2,$i) if(not contains(@usedch2,$i));
	} else {
	  push(@remch2,$i) if(not contains(@remch2,$i));
	}
      }
      print $fhch2vecsdel $#usedch2+1,"\n\n";
      print $fhch2pos $#usedch2+1,"\n\n";
      foreach $i (@usedch2) {
	printf $fhch2vecsdel "H %10.7f %10.7f %10.7f\n", @{$valch2[$i]}[0..2];
	printf $fhch2pos     "H %10.7f %10.7f %10.7f\n",
	  @{$cdata[0][$valch2[$i][10]][$valch2[$i][6]][$valch2[$i][7]]}[0..2];
      }
      
      #remove centrosymmetric ch3 vibrations
      @loopch1 = shuffle(0..$#valch3);
      @loopch2 = shuffle(0..$#valch3);
      $numiter=0;
      do {
	$newconn=0;
	foreach $i (@loopch1) {
	  next if(defined($valch3[$i][9]));
	foreach $i (@loopch2) {
	    next if($i==$j);
	    next if(defined($valch3[$j][9]));
	    $sum[0] = $valch3[$i][0] + $valch3[$j][0];
	    $sum[1] = $valch3[$i][1] + $valch3[$j][1];
	    $sum[2] = $valch3[$i][2] + $valch3[$j][2];
	    $r = vector_length(@sum);
	    next if($r>$tolerance);
	    if($r<$valch3[$i][8] and $r<$valch3[$j][8]) {
	      $newconn++;
	      $valch3[$valch3[$i][9]][9] = undef; # remove previous partner;
	      $valch3[$valch3[$j][9]][9] = undef; # remove previous partner;
	      $valch3[$i][8] = $r; # reset minimum distance
	      $valch3[$j][8] = $r; # reset minimum distance
	      $valch3[$i][9] = $j; # set new partner
	      $valch3[$j][9] = $i; # set new partner
	    }
	  }
	}
	for($i=0;$i<@valch3;$i++) {
	  $valch3[$i][8]=2 if(not defined($valch3[$i][9]));
	}
      } until($newconn==0);
      for($i=0;$i<@valch3;$i++) {
	if(not defined($valch3[$i][9])) { # if vibration has no partner
	  push(@usedch3,$i) if(not contains(@usedch3,$i));
	} else {
	  push(@remch3,$i) if(not contains(@remch3,$i));
	}
      }
      print $fhch3vecsdel $#usedch3+1,"\n\n";
      foreach $i (@usedch3) {
	printf $fhch3vecsdel "H %10.7f %10.7f %10.7f\n", @{$valch3[$i]}[0..2]
      }
      $fhch2vecs->autoflush;
      $fhch2vecsdel->autoflush;
      $fhch3vecs->autoflush;
      $fhch3vecsdel->autoflush;
    } else {
      @usedch2=(0..$#valch2);
      @usedch3=(0..$#valch3);
    }
    # output ch2
    foreach $i (@usedch2) {
      $dx  = $cdata[0][$valch2[$i][10]][$valch2[$i][6]][$valch2[$i][7]][0];
      $dy  = $cdata[0][$valch2[$i][10]][$valch2[$i][6]][$valch2[$i][7]][1];
      $dxy = sqrt($dx*$dx+$dy*$dy);
      if($dxy<$radiuscutoff) {
	$avch2z[$f] += abs($valch2[$i][2]);
	$numch2cut++;
      };
      histogram_add_value(0,$dxy,$histresxy,abs($valch2[$i][2])) if($lspatial);
      if($langular and $dxy<=$minsize) {
	$angle=atan2($dy,$dx);
	$angle+=2*pi if($angle<0);
	histogram_add_value(2,$angle,$histresrad,abs($valch2[$i][2]));
      }
    }
#     $avch2z[$f] /= $areatot;#@ch2s*$mol_numents[0][0];
    printf $fhavch2z "%10u %17.10g %17.10g %10u\n",$fn[$f],$avch2z[$f],$areatot,$numch2cut;
    
    # output ch3
    foreach $i (@usedch3) {
      $dx = $cdata[0][$valch3[$i][10]][$valch3[$i][6]][$valch3[$i][7]][0];
      $dy = $cdata[0][$valch3[$i][10]][$valch3[$i][6]][$valch3[$i][7]][1];
      $dxy = sqrt($dx*$dx+$dy*$dy);
      if($dxy<$radiuscutoff) {
	$avch3z[$f] += abs($valch3[$i][2]);
	$numch3cut++;
      };
      histogram_add_value(1,$dxy,$histresxy,abs($valch3[$i][2])) if($lspatial);
      if($langular and $dxy<=$minsize) {
	$angle=atan2($dy,$dx);
	$angle+=2*pi if($angle<0);
	histogram_add_value(3,$angle,$histresrad,abs($valch3[$i][2]));
      }
    }
#     $avch3z[$f] /= $areatot;#@ch3s*$mol_numents[0][0];
    printf $fhavch3z "%10u %17.10g %17.10g %10u\n",$fn[$f],$avch3z[$f],$areatot,$numch3cut;
    $f++;
  }
  close($fhhist);
}

close($fhavch2z,$fhavch3z,$fhch2vecs,$fhch3vecs,$fhch2vecsdel,$fhch3vecsdel,$fhch2pos);

if($lessoutput) {
  print "last frame analyzed: $frame_number[0]\n";
} else {
  print "\n";
}

histogram_normalize_value(0,1,2,3,$f);

# write histograms
if($lspatial) {
  if($cancel) {
    open($fhhist,">","$outdir/AV_CH2_Z_RADIAL_DISTRIB_CANCELLED");
  } else {
    open($fhhist,">","$outdir/AV_CH2_Z_RADIAL_DISTRIB");
  }
  print $fhhist "# radial distribution of average z-component of CH2 vector averaged over $f frames.\n";
  printf $fhhist "#%16s %17s %17s %17s %10s\n","radius [A]","value [A]","stdev [A]","area [nm2]","# values";
  for($i=$histminindex[0];$i<=$histmaxindex[0];$i++) {
    ($av,$stdev) = calc_stdev(@{$histvals[0][$i]});
    $area = calc_area(($i+1)*$histresxy,$i*$histresxy,@{$size[0]})/100; # in nm^2
    printf $fhhist "%17.10g %17.10g %17.10g %17.10g %10u\n",$histresxy*$i,$histogram[0][$i],$stdev,$area,$histnum[0][$i];
  }
  close($fhhist);
  if($cancel) {
    open($fhhist,">","$outdir/AV_CH3_Z_RADIAL_DISTRIB_CANCELLED");
  } else {
    open($fhhist,">","$outdir/AV_CH3_Z_RADIAL_DISTRIB");
  }
  print $fhhist "# radial distribution of average z-component of CH3 vector averaged over $f frames.\n";
  printf $fhhist "#%16s %17s %17s %17s %10s\n","radius [A]","value [A]","stdev [A]","area [nm2]","# values";
  for($i=$histminindex[1];$i<=$histmaxindex[1];$i++) {
    ($av,$stdev) = calc_stdev(@{$histvals[1][$i]});
    $area = calc_area(($i+1)*$histresxy,$i*$histresxy,@{$size[0]})/100; # in nm^2
    printf $fhhist "%17.10g %17.10g %17.10g %17.10g %10u\n",$histresxy*$i,$histogram[1][$i],$stdev,$area,$histnum[1][$i];
  }
  close($fhhist);
}

if($langular) {
  $area = pi*$minsize*$minsize*$histresdeg/360.0/100;
  if($cancel) {
    open($fhhist,">","$outdir/AV_CH2_Z_ANGULAR_DISTRIB_CANCELLED");
  } else {
    open($fhhist,">","$outdir/AV_CH2_Z_ANGULAR_DISTRIB");
  }
  print $fhhist "# angular distribution of average z-component of CH2 vector within $minsize A of the center averaged over $f frames.\n";
  printf $fhhist "#%16s %17s %17s %10s\n","angle [deg]","value [A]","stdev [A]","area [nm2]","# values";
  for($i=$histminindex[2];$i<=$histmaxindex[2];$i++) {
    ($av,$stdev) = calc_stdev(@{$histvals[2][$i]});
    printf $fhhist "%17.10g %17.10g %17.10g %17.10g %10u\n",$histresdeg*$i,$histogram[2][$i],$stdev,$area,$histnum[2][$i];
  }
  ($av,$stdev) = calc_stdev(@{$histvals[2][0]});
  printf $fhhist "%17.10g %17.10g %17.10g %17.10g %10u\n",360,$histogram[2][0],$stdev,$area,$histnum[2][0];
  close($fhhist);
  if($cancel) {
    open($fhhist,">","$outdir/AV_CH3_Z_ANGULAR_DISTRIB_CANCELLED");
  } else {
    open($fhhist,">","$outdir/AV_CH3_Z_ANGULAR_DISTRIB");
  }
  print $fhhist "# angular distribution of average z-component of CH3 vector within $minsize A of the center averaged over $f frames.\n";
  printf $fhhist "#%16s %17s %17s %10s\n","angle [deg]","value [A]","stdev [A]","area [nm2]","# values";
  for($i=$histminindex[3];$i<=$histmaxindex[3];$i++) {
    ($av,$stdev) = calc_stdev(@{$histvals[3][$i]});
    printf $fhhist "%17.10g %17.10g %17.10g %17.10g %10u\n",$histresdeg*$i,$histogram[3][$i],$stdev,$area,$histnum[3][$i];
  }
  ($av,$stdev) = calc_stdev(@{$histvals[3][0]});
  printf $fhhist "%17.10g %17.10g %17.10g %17.10g %10u\n",360,$histogram[3][0],$stdev,$area,$histnum[3][0];
  close($fhhist);
}

sub unwrap {
  my @pos1 = @{$_[0]};
  my @pos2 = @{$_[1]};
  my($d,$c);
  for($c=0;$c<2;$c++) {
    $d = $pos1[$c] - $pos2[$c];
    if($d > $size[0][$c]) {
      $pos1[$c] -= 2 * $size[0][$c];
    } elsif($d <-$size[0][$c]) {
      $pos1[$c] += 2 * $size[0][$c];
    }
  }
  return @pos1;
}

sub error {
  if($_[0]==1) {
    print "**** error: start frame must be an integer number!\n";
  } elsif($_[0]==2) {
    print "**** error: end frame must be an integer number!\n";
  } elsif($_[0]==3) {
    print "**** error: carbon atom has $_[1] hydrogen atoms bound!\n";
  } elsif($_[0]==4) {
    print "**** error: could not open output file!\n";
  } elsif($_[0]==5) {
    print "**** error: did not find timestep $_[1] in $dir/HISTORY!\n";
  } elsif($_[0]==6) {
    print "**** error: histogram resolution must be greater than zero!\n";
  } elsif($_[0]==7) {
    print "**** error: no suitable directories found!\n";
  } elsif($_[0]==8) {
    print "**** error: tolerance must be greater than zero!\n";
  } elsif($_[0]==9) {
    print "**** error: angular resolution must be divisor of 360!\n";
  } elsif($_[0]==10) {
    print "**** error: unknown flag $_[1]!\n";
  } elsif($_[0]==11) {
    print "**** error: molecule $_[1] was not found in FIELD data!\n";
  } else {
    print "**** unknown error!\n";
    exit 99;
  }
  exit $_[0]
}

sub rem_centrosymmetric {
  my $data   = @{$_[0]};
  my $cutoff = $_[1];
  my($i,$j,@sum,$r);
  loopouter : for($i=0;$i<@{$data};$i++) {
    for($j=$i+1;$j<@{$data};$j++) {
      $sum[0] = ${$data}[$i][0] + ${$data}[$j][0];
      $sum[1] = ${$data}[$i][1] + ${$data}[$j][1];
      $sum[2] = ${$data}[$i][2] + ${$data}[$j][2];
      $r = vector_length(@sum);
      if($r<$cutoff) {
	splice(@{$data},$j,1);
	splice(@{$data},$i,1);
	$i--;
	next loopouter;
      }
    }
  }
}

sub calc_area {
  my(@syssize,$router,$rinner,$area,$halfdiag,$h);
  $router = $_[0];
  $rinner = $_[1];
  @syssize = @_[2..3];
  $halfdiag = sqrt($syssize[0]*$syssize[0]+$syssize[1]*$syssize[1]);
  return 0 if($rinner>$halfdiag);
  $router=$halfdiag if($router>$halfdiag);
  $area = pi*($router*$router-$rinner*$rinner);
  # check x-axis with outer radius
  $h = $router-$syssize[0];
  $area -= 2*segment_area($router,$h) if($h>0);
  # check y-axis with outer radius
  $h = $router-$syssize[1];
  $area -= 2*segment_area($router,$h) if($h>0);
  # check x-axis with inner radius
  $h = $rinner-$syssize[0];
  $area += 2*segment_area($rinner,$h) if($h>0);
  # check y-axis with inner radius
  $h = $rinner-$syssize[1];
  $area += 2*segment_area($rinner,$h) if($h>0);
  return $area;
}

sub segment_area {
  my $r = $_[0];
  my $h = $_[1];
  return $r*$r*acos(1-$h/$r)-sqrt(2*$r*$h-$h*$h)*($r-$h);
}