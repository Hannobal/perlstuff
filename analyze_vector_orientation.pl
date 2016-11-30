#!/usr/bin/perl

use hanno_utility;
use dlpoly_utility;
use Switch;
use Storable qw(dclone);
use Math::Trig;

$lessoutput = 0;
$histresthetadeg = 1;
$histresphideg   = 1;
$startframe      = 0;
$outdir          = "analysis";
$endframe        = 9e20;
@lboundaries = (0,0,0);
@boundaries  = ([-9.0e20,9.0e20],[-9.0e20,9.0e20],[-9.0e20,9.0e20]);
@vectors = ();
$pihalf  = pi/2.0;
$twopi   = pi*2.0;

if($#ARGV<0) {
  print "required input:\n";
  print " -a n*<3*int> IDs (starting from 1): 1. molecule 2. from-atom 3. to-atom\n";
  print "optional flags:\n";
  print " -o output directory for the output files\n";
  print "-rt <real>   histogram resolution for theta (default: $histresthetadeg deg)\n";
  print "-rp <real>   histogram resolution for phi   (default: $histresphideg deg)\n";
  print "-x <2*real>  cutoff for x position\n";
  print "-y <2*real>  cutoff for y position\n";
  print "-z <2*real>  cutoff for z position\n";
  print "-s <int>     start frame\n";
  print "-e <int>     end frame\n";
  print "-lessoutput  do not print current frame\n";
  exit;
}

for($i=0;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-a/i) {
       $i++;
       until($ARGV[$i]=~/^-/ or $i>$#ARGV) {
	&error(4) if(not check_integer($ARGV[$i]));
	push(@vectors,[$ARGV[$i]-1,$ARGV[$i+1]-1,$ARGV[$i+2]-1]);
	$i+=3;
      };
      $i--;
    } case (/^-rt/i) {
      $i++;
      $histresthetadeg = $ARGV[$i];
      &error(1) if($histresthetadeg<=0);
    } case (/^-rp/i) {
      $i++;
      $histresphideg = $ARGV[$i];
      &error(1) if($histresphideg<=0);
    } case (/^-x/i) {
      $i++;
      $boundaries[0][0] = $ARGV[$i];
      $i++;
      $boundaries[0][1] = $ARGV[$i];
      &error(2) if(not check_real(@{$boundaries[0]}));
      $lboundaries[0]=1;
    } case (/^-y/i) {
      $i++;
      $boundaries[1][0] = $ARGV[$i];
      $i++;
      $boundaries[1][1] = $ARGV[$i];
      &error(2) if(not check_real(@{$boundaries[1]}));
      $lboundaries[1]=1;
    } case (/^-z/i) {
      $i++;
      $boundaries[2][0] = $ARGV[$i];
      $i++;
      $boundaries[2][1] = $ARGV[$i];
      &error(2) if(not check_real(@{$boundaries[2]}));
      $lboundaries[2]=1;
    } case(/^-s/i) {
      $i++;
      $startframe = $ARGV[$i];
      &error(3) if not check_integer($startframe);
    } case(/^-e/i) {
      $i++;
      $endframe = $ARGV[$i];
      &error(3) if not check_integer($endframe);
    } case(/^-lessoutput/i) {
      $lessoutput=1;
    } else {
      &error(10,$ARGV[$i]);
    }
  }
}

$histresphirad   = $histresphideg/180.0*pi;
$histresthetarad = $histresthetadeg/180.0*pi;

@directories=get_numeric_directories(".");
if(@directories) {
  @directories=get_numeric_directories(".",$startframe,$endframe);
}
push(@directories,".") if(-e "./HISTORY");
&error(7) if(not @directories);

exit 1 if(read_field_file("$directories[0]/FIELD",0)!=0);

# check validity of molecule/atom indices
for($i=0;$i<@vectors;$i++) {
  &error(5,$vectors[$i][0]+1) if($vectors[$i][0]>=$field_nummols[0]);
  &error(6,$vectors[$i][1]+1) if($vectors[$i][1]>=$mol_numatoms[0][$vectors[$i][0]]);
  &error(6,$vectors[$i][2]+1) if($vectors[$i][2]>=$mol_numatoms[0][$vectors[$i][0]]);
}

mkdir($outdir) if(not -d $outdir);
print "vectors for analysis:\n";
for($i=0;$i<@vectors;$i++) {
  $t  = $vectors[$i][0];
  $a1 = $vectors[$i][1];
  $a2 = $vectors[$i][2];
  # histogram indices
  $hit[$i] = $i*2;
  $hip[$i] = $i*2+1;
  
  &error(8) if(not open($fhavtheta[$i],">","$outdir/AV_THETA_".($t+1)."_".($a1+1)."-".($a2+1)));
  print {$fhavtheta[$i]} "# vectors from ",$t+1," (",$mol_name[0][$t],
    ") from atom ",$a1+1," (",$mol_atomdata[0][$t][$a1][0],
    ") to atom ",$a2+1," (",$mol_atomdata[0][$t][$a2][0],")\n";
    
  print  {$fhavtheta[$i]} "# " if(contains(@lboundaries,1));
  print  {$fhavtheta[$i]} "from x=$boundaries[0][0] to x=$boundaries[0][1]" if($lboundaries[0]);
  print  {$fhavtheta[$i]} "from y=$boundaries[1][0] to y=$boundaries[1][1]" if($lboundaries[1]);
  print  {$fhavtheta[$i]} "from z=$boundaries[2][0] to z=$boundaries[2][1]" if($lboundaries[2]);
  print  {$fhavtheta[$i]} "\n" if(contains(@lboundaries,1));
  printf {$fhavtheta[$i]} "# %8s %17s %17s\n","frame","angle [deg]","stdev [deg]";
  
  &error(8) if(not open($fhavphi[$i],">","$outdir/AV_PHI_".($t+1)."_".($a1+1)."-".($a2+1)));
  print {$fhavphi[$i]} "# vectors from ",$t+1," (",$mol_name[0][$t],
    ") from atom ",$a1+1," (",$mol_atomdata[0][$t][$a1][0],
    ") to atom ",$a2+1," (",$mol_atomdata[0][$t][$a2][0],")\n";
    
  print  {$fhavphi[$i]} "# " if(contains(@lboundaries,1));
  print  {$fhavphi[$i]} "from x=$boundaries[0][0] to x=$boundaries[0][1]" if($lboundaries[0]);
  print  {$fhavphi[$i]} "from y=$boundaries[1][0] to y=$boundaries[1][1]" if($lboundaries[1]);
  print  {$fhavphi[$i]} "from z=$boundaries[2][0] to z=$boundaries[2][1]" if($lboundaries[2]);
  print  {$fhavphi[$i]} "\n" if(contains(@lboundaries,1));
  printf {$fhavphi[$i]} "# %8s %17s %17s\n","frame","angle [deg]","stdev [deg]";
  
  print $t+1," (",$mol_name[0][$t],
    ") from atom ",$a1+1," (",$mol_atomdata[0][$t][$a1][0],
    ") to atom ",$a2+1," (",$mol_atomdata[0][$t][$a2][0],")\n";
}


$firsthist=1;
print "begin analysis...\n" if($lessoutput);

foreach $dir (@directories) {
  &error(11,"$dir/HISTORY") if(not open($fhhist,"<","$dir/HISTORY"));
  if($startframe>0 and $firsthist) {
    &error(9,$startframe) if(not find_history_timestep($fhhist,$startframe));
  }
  $firsthist=0;
  $err=0;
  while($err==0) {
    $err=read_history_timestep($fhhist,0,0);
    last if($err!=0);
    last if($frame_number[0]>$endframe);
    print "\ranalyzing timestep $frame_number[0]" unless($lessoutput);
    for($i=0;$i<@vectors;$i++) {
      undef @{$arrphi[$i]};
      undef @{$arrtheta[$i]};
      $t  = $vectors[$i][0];
      $a1 = $vectors[$i][1];
      $a2 = $vectors[$i][2];
      loopmols:for($m=0;$m<$mol_numents[0][$t];$m++) {
	# check boundaries
	for($c=0;$c<3;$c++) {
	  if($lboundaries[$c]) {
	    next loopmols if($cdata[0][$t][$m][$a1][$c]<$boundaries[$c][0]);
	    next loopmols if($cdata[0][$t][$m][$a1][$c]>$boundaries[$c][1]);
	    next loopmols if($cdata[0][$t][$m][$a2][$c]<$boundaries[$c][0]);
	    next loopmols if($cdata[0][$t][$m][$a2][$c]>$boundaries[$c][1]);
	  }
	}
	# calc theta, phi
	@d = calc_dvec_orthocell(0,$t,$m,$a1,0,$t,$m,$a2);
	$dxy = sqrt($d[0]*$d[0]+$d[1]*$d[1]);
	$theta = $pihalf-atan2($d[2],$dxy);
	$phi   = atan2($d[2],$d[1]);
	$phi  += $twopi if($phi<0);
	histogram_add_one($hit[$i],$theta,$histresthetarad);
	histogram_add_one($hip[$i],$phi,$histresphirad);
	push(@arrtheta,$theta);
	push(@arrphi,$phi);
      }
      ($av,$stdev) = calc_stdev(@arrtheta);
      printf {$fhavtheta[$i]} "%10u %17.10g %17.10g\n",$frame_number[0],$av/pi*180.0,$stdev/pi*180.0;
      ($av,$stdev) = calc_stdev(@arrphi);
      printf {$fhavphi[$i]} "%10u %17.10g %17.10g\n",$frame_number[0],$av/pi*180.0,$stdev/pi*180.0;
    }
  }
  close($fhhist);
}


for($i=0;$i<@vectors;$i++) {
  $t  = $vectors[$i][0];
  $a1 = $vectors[$i][1];
  $a2 = $vectors[$i][2];
  &histogram_normalize_integral($hit[$i],$histresthetarad);
  &histogram_normalize_integral($hip[$i],$histresphirad);
  &error(8) if(not open($fhhist,">","$outdir/HIST_THETA_".($t+1)."_".($a1+1)."-".($a2+1)));
  print $fhhist "# histogram for theta values of vectors from ",$t+1," (",$mol_name[0][$t],
    ") from atom ",$a1+1," (",$mol_atomdata[0][$t][$a1][0],
    ") to atom ",$a2+1," (",$mol_atomdata[0][$t][$a2][0],")\n";
    
  print  $fhhist "# " if(contains(@lboundaries,1));
  print  $fhhist "from x=$boundaries[0][0] to x=$boundaries[0][1]" if($lboundaries[0]);
  print  $fhhist "from y=$boundaries[1][0] to y=$boundaries[1][1]" if($lboundaries[1]);
  print  $fhhist "from z=$boundaries[2][0] to z=$boundaries[2][1]" if($lboundaries[2]);
  print  $fhhist "\n" if(contains(@lboundaries,1));
  printf $fhhist "#%16s %17s %10s\n","theta [deg]","normalized count","raw count";
  for($j=$histminindex[$hit[$i]];$j<=$histmaxindex[$hit[$i]];$j++) {
    printf $fhhist "%17.10g %17.10g %10u\n",$histresthetadeg*$j,$histogram[$hit[$i]][$j],$histnum[$hit[$i]][$j];
  }
  close($fhhist);
  &error(8) if(not open($fhhist,">","$outdir/HIST_PHI_".($t+1)."_".($a1+1)."-".($a2+1)));
  print $fhhist "# histogram for phi values of vectors from ",$t+1," (",$mol_name[0][$t],
    ") from atom ",$a1+1," (",$mol_atomdata[0][$t][$a1][0],
    ") to atom ",$a2+1," (",$mol_atomdata[0][$t][$a2][0],")\n";
    
  print  $fhhist "# " if(contains(@lboundaries,1));
  print  $fhhist "from x=$boundaries[0][0] to x=$boundaries[0][1]" if($lboundaries[0]);
  print  $fhhist "from y=$boundaries[1][0] to y=$boundaries[1][1]" if($lboundaries[1]);
  print  $fhhist "from z=$boundaries[2][0] to z=$boundaries[2][1]" if($lboundaries[2]);
  print  $fhhist "\n" if(contains(@lboundaries,1));
  printf $fhhist "#%16s %17s %10s\n","phi [deg]","normalized count","raw count";
  for($j=$histminindex[$hip[$i]];$j<=$histmaxindex[$hip[$i]];$j++) {
    printf $fhhist "%17.10g %17.10g %10u\n",$histresphideg*$j,$histogram[$hip[$i]][$j],$histnum[$hip[$i]][$j];
  }
  close($fhhist,$fhavphi[$i],$fhavtheta[$i]);
}
print "\n" unless($lessoutput);


sub error {
  if($_[0]==1) {
    print "**** error: histogram resolution must be greater than zero!\n";
  } elsif($_[0]==2) {
    print "**** error: boundaries must be real numbers!\n";
  } elsif($_[0]==3) {
    print "**** error: start and end frame must be integer numbers!\n";
  } elsif($_[0]==4) {
    print "**** error: molecule and atom indices must be integer numbers!\n";
  } elsif($_[0]==5) {
    print "**** error: molecule index $_[1] is greater than the number of molecules in FIELD file!\n";
  } elsif($_[0]==6) {
    print "**** error: atom index $_[1] is greater than the number of molecules in FIELD file!\n";
  } elsif($_[0]==7) {
    print "**** error: no suitable directories found!\n";
  } elsif($_[0]==8) {
    print "**** error: could not open output file!\n";
  } elsif($_[0]==9) {
    print "**** error: did not find timestep $_[1]!\n";
  } elsif($_[0]==10) {
    print "**** error: unknown flag $_[1]!\n";
  } elsif($_[0]==11) {
    print "**** error: could not open file $_[1]!\n";
  } else {
    print "**** unknown error!\n";
    exit 99;
  }
  exit $_[0]
}
