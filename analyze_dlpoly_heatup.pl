#!/usr/bin/perl

# extracts the information from DL_POLY OUTPUT files distributed
# over several folders with numbered names

if($#ARGV==-1) {
  print "please enter an output directory for the output files\n";
  print "optional: start and end temperature and temperature step\n";
  print "optional argument: format-flag (frame/step/time)\n";
  exit;
}

$outdir = $ARGV[0];

if($ARGV[1] eq '') {
  $T = 1;
} else {
  $T = $ARGV[1];
}
if($ARGV[2] eq '') {
  $end = 5000;
} else {
  $end = $ARGV[2];
}
if($ARGV[3] eq '') {
  $step = 1;
} else {
  $step = $ARGV[3];
}

if($ARGV[4] eq "frame") {
  $outfmt=0;
} elsif($ARGV[4] eq "step") {
  $outfmt=1;
} elsif($ARGV[4] eq "time") {
  $outfmt=2;
# } elsif($ARGV[4] eq "trajframe") {
#   $outfmt=3;
} else {
  print "using frames as output for ordinate\n";
  $outfmt=0;
}

if(not chdir($outdir)) {
  mkdir $outdir or die "cannot create output-directory: $!";
  chdir $outdir or die "cannot access output-directory: $!";
}
open(E_TOT, ">", "E_TOT") or die "Can't create file E_TOT: $!";
open(T_TOT, ">", "T_TOT") or die "Can't create file T_TOT: $!";
open(E_CFG, ">", "E_CFG") or die "Can't create file E_CFG: $!";
open(E_VDW, ">", "E_VDW") or die "Can't create file E_VDW: $!";
open(E_COU, ">", "E_COU") or die "Can't create file E_COU: $!";
open(E_BND, ">", "E_BND") or die "Can't create file E_BND: $!";
open(E_ANG, ">", "E_ANG") or die "Can't create file E_ANG: $!";
open(E_DIH, ">", "E_DIH") or die "Can't create file E_DIH: $!";
open(E_TET, ">", "E_TET") or die "Can't create file E_TET: $!";
open(E_PV, ">", "E_PV") or die "Can't create file E_PV: $!";
open(T_ROT, ">", "T_ROT") or die "Can't create file T_ROT: $!";
open(VOLUME, ">", "VOLUME") or die "Can't create file VOLUME: $!";
open(ALPHA, ">", "ALPHA") or die "Can't create file ALPHA: $!";
open(BETA, ">", "BETA") or die "Can't create file BETA: $!";
open(GAMMA, ">", "GAMMA") or die "Can't create file GAMMA: $!";
open(PRESS, ">", "PRESS") or die "Can't create file PRESS: $!";
chdir ".." or die "cannot access root directory: $!";
$error = 0;
$framestot = 0;
$stepstot = 0;
$timetot = 0;
$trajframestot = 0;
while($T<=$end and $error==0) {
  if(-e "./$T/OUTPUT") {
#     print "analyzing $T K\n";
    open(OUTPUT, "<", "./$T/OUTPUT");
    open(CONTROL, "<", "./$T/CONTROL") or die "Can't open file CONTROL: $!";
    while(<CONTROL>) {
      if(/\s*timestep/) {
	($timestep) = /\s*timestep\s+(\S*)/;
      } elsif(/traj/) {
	($trajstart,$trajinterval) = /\s*traj\S*\s+(\S+)\s+(\S+)/;
      }
    }
    close(CONTROL);
    $i=0;
    $nexttrajframe=$trajstart;
    while(<OUTPUT>) {
      if(/---------------/) {
        $_ = <OUTPUT>;
        ($step[$i], $etot[$i], $ttot[$i], $ecfg[$i], $evdw[$i], $ecou[$i], $ebnd[$i], $eang[$i], $edih[$i], $etet[$i]) =
        /\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)/;
        if($step[$i]!=0) {
          $_ = <OUTPUT>;
          ($time[$i], $epv[$i], $trot[$i], $vcfg[$i], $vvdw[$i], $vcou[$i], $vbnd[$i], $vang[$i], $vcou[$i], $vtet[$i]) =
          /\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)/;
	  $time[$i]=$step[$i]*$timestep;
          $_ = <OUTPUT>;
          ($cpu[$i], $vol[$i], $tshl[$i], $eshl[$i], $vshl[$i], $alpha[$i], $beta[$i], $gamma[$i], $vpmf[$i], $press[$i]) =
          /\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)/;
	  $ordinate[$i][0]=$i+$framestot;
	  $ordinate[$i][1]=$step[$i]+$stepstot;
	  $ordinate[$i][2]=$time[$i]+$timetot;
# 	  $ordinate[$i][3]=($step[$i]-$trajstart)/$trajinterval+$trajframestot;
          $i++;
        }
      } elsif(/run terminating/) {
        last;
      }
    }
    $imax = $i;
    for($i=0;$i<$imax;$i++) {
      printf E_TOT  "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$etot[$i];
      printf T_TOT  "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$ttot[$i];
      printf E_CFG  "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$ecfg[$i];
      printf E_VDW  "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$evdw[$i];
      printf E_COU  "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$ecou[$i];
      printf E_BND  "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$ebnd[$i];
      printf E_ANG  "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$eang[$i];
      printf E_DIH  "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$edih[$i];
      printf E_TET  "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$etet[$i];
      printf E_PV   "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$epv[$i];
      printf T_ROT  "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$trot[$i];
      printf VOLUME "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$vol[$i];
      printf ALPHA  "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$alpha[$i];
      printf BETA   "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$beta[$i];
      printf GAMMA  "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$gamma[$i];
      printf PRESS  "%15.5f %10.4e\n",$ordinate[$i][$outfmt],$press[$i];
    }
    $framestot += $imax;
    $timetot += $time[$imax-1];
    $stepstot += $step[$imax-1];
  } else {
#     print "file $T/OUTPUT was not found\n";
  }
  $T += $step;
}
close(E_TOT, T_TOT, E_CFG, E_VDW, E_COU, E_BND, E_ANG, E_DIH, E_TET, E_PV, T_ROT, VOLUME, ALPHA, BETA, GAMMA, PRESS);