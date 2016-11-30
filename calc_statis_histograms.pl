#!/usr/bin/perl

# calculates a histogram from the energies in the STATIS file

use dlpoly_utility;
use hanno_utility;
use Switch;
use POSIX qw(ceil floor);

$startframe = 0;
$endframe   = 9e20;

if($#ARGV<0) {
  print "the input format is:\n",
  " 1. output directory\n",
  "-s <int>   start frame\n",
  "-e <int>   end frame\n",
  "possible Histograms (each followed by resoultion):\n",
  "-etot      total energy\n",
  "-ecfg      potential energy\n",
  "-evdw      van-der-Waals energy\n",
  "-ecoul     Coulomb energy\n",
  "-ebnd      bond energy\n",
  "-eang      angle energy\n",
  "-edih      torsion energy\n",
  "-etet      tether energy\n",
  "-htot      enthalpy\n",
  "-ttot      total temperature\n",
  "-trot      rotational temperature\n",
  "-vol       volume\n",
  "-pres      pressure\n";
  exit 1;
}

$outdir  = $ARGV[0];
@indices = ();
@names = ("step","time","ETOT","TTOT","ECFG","EVDW","ECOUL","EBND","EANG","EDIH","ETET","HTOT",
"TROT","VIRTOT","VIRVDW","VIRCOUL","VIRBND","VIRANG","VIRCONST","VIRTET","VOL","TCORESHELL",
"EPOTCORESHELL","VIRCORESHELL","CELL_ALPHA","CELL_BETA","CELL_GAMMA","VIRPMF","PRESS");

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-s$/) { # start frame
      $i++;
      $startframe = $ARGV[$i];
    } case (/^-e$/) { # end frame
      $i++;
      $endframe = $ARGV[$i];
    } case (/^-etot/) {
      push(@indices,2);
      $i++;
      $histres[2] = $ARGV[$i];
      die "**** error: histogram resoultion must be a real number!\n" if(not check_real($histres[2]));
    } case (/^-ecfg/) {
      push(@indices,4);
      $i++;
      $histres[4] = $ARGV[$i];
      die "**** error: histogram resoultion must be a real number!\n" if(not check_real($histres[4]));
    } case (/^-evdw/) {
      push(@indices,5);
      $i++;
      $histres[5] = $ARGV[$i];
      die "**** error: histogram resoultion must be a real number!\n" if(not check_real($histres[5]));
    } case (/^-ecoul/) {
      push(@indices,6);
      $i++;
      $histres[6] = $ARGV[$i];
      die "**** error: histogram resoultion must be a real number!\n" if(not check_real($histres[6]));
    } case (/^-ebnd/) {
      push(@indices,7);
      $i++;
      $histres[7] = $ARGV[$i];
      die "**** error: histogram resoultion must be a real number!\n" if(not check_real($histres[7]));
    } case (/^-eang/) {
      push(@indices,8);
      $i++;
      $histres[8] = $ARGV[$i];
      die "**** error: histogram resoultion must be a real number!\n" if(not check_real($histres[8]));
    } case (/^-edih/) {
      push(@indices,9);
      $i++;
      $histres[9] = $ARGV[$i];
      die "**** error: histogram resoultion must be a real number!\n" if(not check_real($histres[9]));
    } case (/^-etet/) {
      push(@indices,10);
      $i++;
      $histres[10] = $ARGV[$i];
      die "**** error: histogram resoultion must be a real number!\n" if(not check_real($histres[10]));
    } case (/^-htot/) {
      push(@indices,11);
      $i++;
      $histres[11] = $ARGV[$i];
      die "**** error: histogram resoultion must be a real number!\n" if(not check_real($histres[11]));
    } case (/^-ttot/) {
      push(@indices,3);
      $i++;
      $histres[3] = $ARGV[$i];
      die "**** error: histogram resoultion must be a real number!\n" if(not check_real($histres[3]));
    } case (/^-trot/) {
      push(@indices,12);
      $i++;
      $histres[12] = $ARGV[$i];
      die "**** error: histogram resoultion must be a real number!\n" if(not check_real($histres[12]));
    } case (/^-vol/) {
      push(@indices,20);
      $i++;
      $histres[20] = $ARGV[$i];
      die "**** error: histogram resoultion must be a real number!\n" if(not check_real($histres[20]));
    } case (/^-pres/) {
      push(@indices,28);
      $i++;
      $histres[28] = $ARGV[$i];
      die "**** error: histogram resoultion must be a real number!\n" if(not check_real($histres[28]));
    } else {
      die "**** error: expected optional flag but got \"$ARGV[$i]\"!\n";
    }
  }
}
check_integer($startframe) or die "error: $startframe is not an integer number!\n";
check_integer($endframe)   or die "error: $endframe is not an integer number!\n";

if($endframe<$startframe) {
  die "**** error: end frame must not be smaller than start frame\n";
}
if(not @indices) {
  die "**** error: nothing to do!\n";
}

@directories = get_numeric_directories(".",$startframe,$endframe);
push(@directories,'.') if(-e './STATIS');

if(not @directories) {
  die "**** error: no suitable directories found!\n";
}
# find minimum and maximum values
$numframes=0;
foreach $i (@indices) {
  $minval[$i] =  9e20;
  $maxval[$i] = -9e20;
}
foreach $dir (@directories) {
  open($fhstatis, "<", "$dir/STATIS") or print "**** error: file STATIS was not found in $dir\n";
  while(1) {
    last if(read_statis_timestep($fhstatis,0) != 0);
    next if($sdata[0][0]<$startframe);
    last if($sdata[0][0]>$endframe);
    foreach $i (@indices) {
      $minval[$i] = min($sdata[0][$i],$minval[$i]);
      $maxval[$i] = max($sdata[0][$i],$maxval[$i]);
    }
    $numframes++;
  } # end while
  close($fhstatis);
} # end for directories
foreach $i (@indices) {
  $minval[$i] =  $histres[$i]*(floor($minval[$i]/$histres[$i])-1);
  $maxval[$i] =  $histres[$i]*ceil($maxval[$i]/$histres[$i]);
}

if($numframes==0) {
  die "**** error: no valid frames found!\n";
}

foreach $dir (@directories) {
  open($fhstatis, "<", "$dir/STATIS") or print "**** error: file STATIS was not found in $dir\n";
  print "reading $dir/STATIS\n";
  while(1) {
    last if(read_statis_timestep($fhstatis,0) != 0);
    next if($sdata[0][0]<$startframe);
    last if($sdata[0][0]>$endframe);
    foreach $i (@indices) {
      histogram_add_one($i,$sdata[0][$i]-$minval[$i],$histres[$i]);
    }
  } # end while
  close($fhstatis);
} # end for directories

print "estimates for average and standard deviation:\n";
foreach $i (@indices) {
  histogram_normalize_integral($i,$histres[$i]);
  ($average,$stdev) = histogram_estimate_stdev($i,$histres[$i]);
  print "$names[$i] ",$average+$minval[$i]," +/- $stdev\n";
}

########## write output files ###################

if(not -d $outdir) {
  mkdir $outdir or die "**** error: cannot create output-directory: $!";
}

foreach $i (@indices) {
  open($fhhist,">","$outdir/HISTOGRAM_$names[$i]") or die "**** error: cannot open output file!\n";
  print $fhhist "# histogram from $numframes frames in between $startframe and $endframe with resolution $histres[$i]\n";
  $histminindex[$i]--;
  $histmaxindex[$i]++;
  for($j=$histminindex[$i];$j<=$histmaxindex[$i];$j++) {
    printf $fhhist "%20.8g %10.7f %10u\n", $histres[$i]*$j+$minval[$i], $histogram[$i][$j], $histnum[$i][$j];
  }
  close($fhhist);
}

# sub error {
#   print $_[0];
#   if(not defined $_[1]) {
#     exit 1 
#   } else {
#     exit $_[1];
#   }
# }