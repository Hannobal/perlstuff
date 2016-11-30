#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;

if($#ARGV<1) {
  print "1. name of CONFIG-file\n";
  print "2. name of HISTORY-file\n";
  print "[3. append? (y/n)]\n";
  print "[4. frame number (default: 0)]\n";
  print "[5. timestep (default: 0.001)]\n";
  print "[6. list of atom types to inclue (default: all)]\n";
  exit;
}
$anzges = 0; #number of atoms counted in CONFIG

if($#ARGV>1 and not $ARGV[2] =~ /^[yn]/i) {
  print "**** error: could not interpret input $ARGV[2] for append key!\n";
  exit 1;
}

$configfilename = $ARGV[0];
$histfilename   = $ARGV[1];

$starti = 3;
if(check_integer($ARGV[3])) {
  $framenumber = $ARGV[3];
  $starti++;
} else {
  $framenumber = 0;
}

if(check_real($ARGV[4])) {
  $timestep = $ARGV[4];
  $starti++;
} else {
  $timestep = 0.001;
}

$numtypes = 0;
for($i=$starti;$i<@ARGV;$i++) {
  $includetype[$numtypes] = $ARGV[$i];
  $numtypes++;
}

############# read input file #############################################

if($configfilename =~ "CONFIG") {
  $fieldfilename = $configfilename;
  $fieldfilename =~ s/CONFIG/FIELD/i;
} elsif($configfilename =~ "REVCON") {
  $fieldfilename = $configfilename;
  $fieldfilename =~ s/REVCON/FIELD/i;
} elsif($configfilename =~ /.cfg/i) {
  $fieldfilename = $configfilename;
  $fieldfilename =~ s/.cfg/.fld/i;
} else {
  $field_read_err = 1;
}

if(defined($fieldfilename)) {
  if(-e $fieldfilename) {
    $field_read_err = read_field_file($fieldfilename,0);
  } else {
    $field_read_err = 1;
  }
}

if($field_read_err==0) { $fi=0; } else { $fi=-1; }
exit 1 if(read_config_file($configfilename,$fi,0) > 0);

if($field_read_err == 0) {
  for($t=0;$t<@{$cdata[0]};$t++) {
    for($m=0;$m<@{$cdata[0][$t]};$m++) {
      for($a=0;$a<@{$cdata[0][$t][$m]};$a++) {
	$cdata[0][$t][$m][$a][12] = $mol_atomdata[0][$t][$a][1];
	$cdata[0][$t][$m][$a][13] = $mol_atomdata[0][$t][$a][2];
      }
    }
  }
}

if($numtypes>0) {
  for($t=0;$t<@{$cdata[0]};$t++) {
    for($m=0;$m<@{$cdata[0][$t]};$m++) {
      fora : for($a=0;$a<@{$cdata[0][$t][$m]};$a++) {
	for($i=0;$i<$numtypes;$i++) {
	  next fora if($includetype[$i] eq $cdata[0][$t][$m][$a][9])
	}
	splice(@{$cdata[0][$t][$m]},$a,1);
	$frame_numatoms[0]--;
	$a--;
      }
    }
  }
}

############# write output file ###########################################

if($ARGV[2] =~ /^y/i) {
  open($fhhist, ">>", $histfilename) or die "**** error: Can't create HISTORY-File: $!";
} else {
  open($fhhist, ">", $histfilename) or die "**** error: Can't create HISTORY-File: $!";
}
exit 1 if (write_history_timestep($fhhist,-1,0) > 0);
close($fhhist)
