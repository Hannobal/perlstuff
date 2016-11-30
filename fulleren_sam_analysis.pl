#!/usr/bin/perl
if($#ARGV<1) {
  print "1. HISTORY file to analyze\n";
  print "2. corresponding FIELD file\n";
  print "3. output folder (optional)\n";
  exit;
}

open(HISTORY, "<", $ARGV[0]) or die "cannot open HISTORY file: $!";
open(FIELD, "<", $ARGV[1]) or die "cannot open FIELD file: $!";

$outdir = $ARGV[2];
if(not ($outdir eq "")) {
  if(not chdir($outdir)) {
    mkdir $outdir or die "cannot create output-directory: $!";
    chdir $outdir or die "cannot access output-directory: $!";
  }
}

#open(FULLERENES, ">", "FUL_DISTANCE") or die "cannot open file FUL_DISTANCE: $!";
open(FULLERENES_AV_DIST_MIN, ">", "FULLERENES_AV_DIST_MIN") or die "cannot open file FULLERENES_AV_DIST_MIN: $!";
open(FULLERENES_AV_DIST_CENTER, ">", "FULLERENES_AV_DIST_CENTER") or die "cannot open file FULLERENES_AV_DIST_CENTER: $!";
open(FULLERENES_MIN_DIST, ">", "FULLERENES_MIN_DIST") or die "cannot open file FULLERENES_MIN_DIST: $!";
open(TEST, ">", "fullerenes.xyz") or die "cannot open fullerenes.xyz: $!";

#read out molecules from FIELD file
while(<FIELD>) {
  if(/molecules/ or /MOLECULES/) { last; }
}
($field_nummols) = /\s*\S+\s+(\S+)/;
if($field_nummols==0) {
  print "no molecules found in FIELD file\n";
  exit;
}
for($i=0;$i<$field_nummols;$i++) {
  $mol_name[$i]=<FIELD>;
  chomp($mol_name[$i]);
  $_=<FIELD>;
  ($test1)=/\s*\S+\s+(\S+)/;
  $mol_nummols[$i]=$test1;
  $_=<FIELD>;
  ($test2)=/\s*\S+\s+(\S+)/;
  $mol_numatoms[$i]=$test2;
  $field_numatoms += $mol_nummols[$i]*$mol_numatoms[$i];
  until(/FINISH/ or /finish/ or /CLOSE/ or /close/) {
    $_=<FIELD>;
  }
}
close(FIELD);
# read HISTORY file
$frame=-1;
while(<HISTORY>) {
# reading frame out of history file
  if(/timestep/) {
    $frame++;
    print "evaluating frame ".($frame+1)."\n";
    ($timestep[$frame], $numatoms[$frame], $trajkey[$frame], $periodkey[$frame], $integration[$frame]) =
    /timestep\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    if($numatoms[$frame] != $field_numatoms) {
      print "error in HISOTRY file: $field_numatoms expected, $numatoms[$frame] given in HISTORY!\n";
      exit;
    }
    $surfzmax[$frame]=-99999;
    $_=<HISTORY>;
    ($cell[$frame][0][0],$cell[$frame][0][1],$cell[$frame][0][2])=/\s*(\S+)\s+(\S+)\s+(\S+)/;
    $_=<HISTORY>;
    ($cell[$frame][1][0],$cell[$frame][1][1],$cell[$frame][1][2])=/\s*(\S+)\s+(\S+)\s+(\S+)/;
    $_=<HISTORY>;
    ($cell[$frame][2][0],$cell[$frame][2][1],$cell[$frame][2][2])=/\s*(\S+)\s+(\S+)\s+(\S+)/;
    for($i=0;$i<$field_nummols;$i++) {
      for($j=0;$j<$mol_nummols[$i];$j++) {
        for($k=0;$k<$mol_numatoms[$i];$k++) {
          $_=<HISTORY>;
          ($at_name[$i][$frame][$j][$k], $index[$i][$frame][$j][$k], $mass[$i][$frame][$j][$k],$charge[$i][$frame][$j][$k]) = 
          /\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
          $_=<HISTORY>;
          ($coords[$i][$frame][$j][$k][0],$coords[$i][$frame][$j][$k][1],$coords[$i][$frame][$j][$k][2]) =
          /\s*(\S+)\s+(\S+)\s+(\S+)/;
          if($coords[$i][$frame][$j][$k][2]<0) {
            $coords[$i][$frame][$j][$k][2]+=$cell[$frame][2][2];
          }
          if(uc($at_name[$i][$frame][$j][$k]) eq "AL" or uc($at_name[$i][$frame][$j][$k]) eq "OA") {
            if($surfzmax[$frame]<$coords[$i][$frame][$j][$k][2]) {
              $surfzmax[$frame]=$coords[$i][$frame][$j][$k][2];
            }
          }
          if($trajkey[$frame] > 0) {
            $_=<HISTORY>;
            ($vel[$i][$frame][$j][$k][0],$vel[$i][$frame][$j][$k][1],$vel[$i][$frame][$j][$k][2]) =
            /\s*(\S+)\s+(\S+)\s+(\S+)/;
            if($trajkey[$frame] > 0) {
              $_=<HISTORY>;
              ($force[$i][$frame][$j][$k][0],$force[$i][$frame][$j][$k][1],$force[$i][$frame][$j][$k][2]) =
              /\s*(\S+)\s+(\S+)\s+(\S+)/;
            }
          }
        }
      }
    }
  } else {
    print "error while reading HISTORY file!\n";
    close(HISTORY);
  }
}

$totframes = $frame+1;

#evaluate averages and write output
for($frame=0;$frame<$totframes;$frame++) {
  $numfullatoms=0;
  $numfulls=0;
  $fullerenezmintot[$frame]=99999;
  for($i=0;$i<$field_nummols;$i++) {
    if($mol_name[$i] =~ "C60") {
      for($j=0;$j<$mol_nummols[$i];$j++) {
      $fullerenezmin[$frame][$numfulls]=99999;
        for($k=0;$k<$mol_numatoms[$i];$k++) {
          if(uc($at_name[$i][$frame][$j][$k]) eq "CA" or 
             uc($at_name[$i][$frame][$j][$k]) eq "CX") {
            $numfullatoms++;
            $coords[$i][$frame][$j][$k][2] -= $surfzmax[$frame];
            $fullerenecenter[$frame][$numfulls][0] += $coords[$i][$frame][$j][$k][0];
            $fullerenecenter[$frame][$numfulls][1] += $coords[$i][$frame][$j][$k][1];
            $fullerenecenter[$frame][$numfulls][2] += $coords[$i][$frame][$j][$k][2];
            if($coords[$i][$frame][$j][$k][2]<$fullerenezmin[$frame][$numfulls]) {
              $fullerenezmin[$frame][$numfulls] = $coords[$i][$frame][$j][$k][2];
            }
            $num_ful_atoms[$frame][$numfulls] += 1;
          }
        }
        $fullerenecenter[$frame][$numfulls][0] /= $num_ful_atoms[$frame][$numfulls];
        $fullerenecenter[$frame][$numfulls][1] /= $num_ful_atoms[$frame][$numfulls];
        $fullerenecenter[$frame][$numfulls][2] /= $num_ful_atoms[$frame][$numfulls];
        $avfullerenecenter[$frame][0] += $fullerenecenter[$frame][$numfulls][0];
        $avfullerenecenter[$frame][1] += $fullerenecenter[$frame][$numfulls][1];
        $avfullerenecenter[$frame][2] += $fullerenecenter[$frame][$numfulls][2];
        $avfullerenezmin[$frame] += $fullerenezmin[$frame][$numfulls];
        if($fullerenezmintot[$frame]>$fullerenezmin[$frame][$numfulls]) {
          $fullerenezmintot[$frame] = $fullerenezmin[$frame][$numfulls];
        }
        $numfulls++;
      }
      $avfullerenecenter[$frame][0] /= $numfulls;
      $avfullerenecenter[$frame][1] /= $numfulls;
      $avfullerenecenter[$frame][2] /= $numfulls;
      $avfullerenezmin[$frame] /= $numfulls;
    }
  }
# write an xyz file with the fullerene coordinates + centers
  print TEST ($numfulls+$numfullatoms)."\n";
  print TEST "\n";
  for($i=0;$i<$numfulls;$i++) {
    print TEST "N $fullerenecenter[$frame][$i][0] $fullerenecenter[$frame][$i][1] $fullerenecenter[$frame][$i][2]\n"
  }
   for($i=0;$i<$field_nummols;$i++) {
     if($mol_name[$i] =~ "C60") {
       for($j=0;$j<$mol_nummols[$i];$j++) {
       $fullerenezmin[$frame][$numfulls]=99999;
         for($k=0;$k<$mol_numatoms[$i];$k++) {
           if(uc($at_name[$i][$frame][$j][$k]) eq "CA" or 
              uc($at_name[$i][$frame][$j][$k]) eq "CX") {
              print TEST "C $coords[$i][$frame][$j][$k][0] $coords[$i][$frame][$j][$k][1] $coords[$i][$frame][$j][$k][2]\n";
           }
         }
       }
     }
   }

  if($frame==0) {
    $timetot[$frame]=0;
  } elsif($timestep[$frame]>$timestep[$frame-1]) {
    $timetot[$frame] = $timetot[$frame-1]+($timestep[$frame]-$timestep[$frame-1])*$integration[$frame];
  } else {
    $timetot[$frame] = $timetot[$frame-1]+$timestep[$frame]*$integration[$frame];
}
  print FULLERENES_AV_DIST_CENTER "$frame $avfullerenecenter[$frame][2]\n";
  print FULLERENES_AV_DIST_MIN "$frame $avfullerenezmin[$frame]\n";
  print FULLERENES_MIN_DIST "$frame $fullerenezmintot[$frame]\n";
}