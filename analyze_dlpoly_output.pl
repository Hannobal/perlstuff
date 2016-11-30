#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;

# extracts the information from DL_POLY OUTPUT files distributed
# over several folders with numbered names

use Cwd;
if($#ARGV<1) {
  print "the input format is:\n";
  print "<output-directory> <split job? (y/n/h)> [begin frame] [end frame]\n";
  exit 1;
}

$outdir = $ARGV[0];

if(not $ARGV[1] =~ /^[ynh]/i) {
  print "**** error: did not recognize input for \"split job?\"!";
  exit 1;
} elsif($ARGV[1] =~ /^y/i) {
  $splitjob = 1;
  $heatup   = 0;
} elsif($ARGV[1] =~ /^h/i) {
  $heatup   = 1;
  $splitjob = 1;
} else {
  $splitjob = 0;
}

if ($#ARGV==1) {
  $startframe = 0;
  $endframe = 999999999;
} elsif ($#ARGV==2) {
  $startframe = $ARGV[2];
  $endframe = 999999999;
} else {
  $startframe = $ARGV[2];
  $endframe = $ARGV[3];
}

check_integer($startframe) or die "error: $startframe is not an integer number!\n";
check_integer($endframe)   or die "error: $endframe is not an integer number!\n";

if($endframe<$startframe) {
  print "**** error: end frame must not be smaller than start frame\n";
  exit;
}

############### find the relevant directories to read #########################
if(not $splitjob) {
  @directories = (".")
} else {
  @directories = get_numeric_directories(".",$startframe,$endframe);
  push(@directories,'.') if(-e './OUTPUT');
}

############### generate the output files #####################################
if(not -d $outdir) {
  mkdir $outdir or die "cannot create output-directory: $!";
}
open(OUTDATA,    ">", "$outdir/OUTPUT_DATA")    or die "Can't create file OUTPUT_DATA: $!";
open(OUTDATA_AV, ">", "$outdir/OUTPUT_DATA_AV") or die "Can't create file OUTPUT_DATA_AV: $!";
printf OUTDATA    "%14s","# 1";
printf OUTDATA_AV "%14s","# 1";
for($i=2;$i<=28;$i++) {
  printf OUTDATA    " %11u", $i;
  printf OUTDATA_AV " %11u", $i;
}
printf OUTDATA "\n%14s","time";
printf OUTDATA " %11s","E_tot";
printf OUTDATA " %11s","T_tot";
printf OUTDATA " %11s","E_cfg";
printf OUTDATA " %11s","E_vdw,metal";
printf OUTDATA " %11s","E_coul";

printf OUTDATA " %11s","E_bnd";
printf OUTDATA " %11s","E_ang,tbp";
printf OUTDATA " %11s","E_dih,fbp";
printf OUTDATA " %11s","E_tether";
printf OUTDATA " %11s","H=E_tot+P*V";

printf OUTDATA " %11s","T_rot";
printf OUTDATA " %11s","vir_cfg";
printf OUTDATA " %11s","vir_vdw,met";
printf OUTDATA " %11s","vir_coul";
printf OUTDATA " %11s","vir_bnd";

printf OUTDATA " %11s","vir_ang";
printf OUTDATA " %11s","vir_constr";
printf OUTDATA " %11s","vir_teth";
printf OUTDATA " %11s","volume";
printf OUTDATA " %11s","T_core-sh";

printf OUTDATA " %11s","E_core-sh";
printf OUTDATA " %11s","vir_core-sh";
printf OUTDATA " %11s","cell_alpha";
printf OUTDATA " %11s","cell_beta";
printf OUTDATA " %11s","cell_gamma";

printf OUTDATA " %11s","vir_pmf";
printf OUTDATA " %11s","pressure";


printf OUTDATA_AV "\n%14s","time";
printf OUTDATA_AV " %11s","E_tot";
printf OUTDATA_AV " %11s","T_tot";
printf OUTDATA_AV " %11s","E_cfg";
printf OUTDATA_AV " %11s","E_vdw,metal";
printf OUTDATA_AV " %11s","E_coul";

printf OUTDATA_AV " %11s","E_bnd";
printf OUTDATA_AV " %11s","E_ang,tbp";
printf OUTDATA_AV " %11s","E_dih,fbp";
printf OUTDATA_AV " %11s","E_tether";
printf OUTDATA_AV " %11s","H=E_tot+P*V";

printf OUTDATA_AV " %11s","T_rot";
printf OUTDATA_AV " %11s","vir_cfg";
printf OUTDATA_AV " %11s","vir_vdw,met";
printf OUTDATA_AV " %11s","vir_coul";
printf OUTDATA_AV " %11s","vir_bnd";

printf OUTDATA_AV " %11s","vir_ang";
printf OUTDATA_AV " %11s","vir_constr";
printf OUTDATA_AV " %11s","vir_teth";
printf OUTDATA_AV " %11s","volume";
printf OUTDATA_AV " %11s","T_core-sh";

printf OUTDATA_AV " %11s","E_core-sh";
printf OUTDATA_AV " %11s","vir_core-sh";
printf OUTDATA_AV " %11s","cell_alpha";
printf OUTDATA_AV " %11s","cell_beta";
printf OUTDATA_AV " %11s","cell_gamma";

printf OUTDATA_AV " %11s","vir_pmf";
printf OUTDATA_AV " %11s","pressure";

foreach $type (@{$field_atomtypes[0]}) {
  printf OUTDATA "%14s","msd($type)";
}

############### read the OUTPUT file from each directory ######################
$time=0;
$time_old=0;
foreach $dir (@directories) {
  print "reading $dir/OUTPUT\n";
  open(OUTPUT, "<", "$dir/OUTPUT") or print "file OUTPUT was not found in $dir\n";
  open(CONTROL, "<", "$dir/CONTROL") or print "file CONTROL was not found in $dir\n";
  $i=0;
  $starttime=$time if($heatup);
  # read the time step interval out of the CONTROL file
  while(<CONTROL>) {
    if(/\s*timestep/) {
      ($timestep) = /\s*timestep\s+(\S*)/;
    } elsif(/traj/) {
      ($trajinterval) = /\s*traj\S*\s+\S+\s+(\S*)/;
    }
  }
  close(CONTROL);
  while(<OUTPUT>) {
    if(/---------------/) {
      @data    = ();
      @data_av = ();
      $_ = <OUTPUT>;
      $_ =~ s/^\s+//; $_ =~ s/\s+$//;
      push(@data,split(/\s+/,$_));
      if($data[0]!=0 and $data[0]>=$startframe and $data[0]<=$endframe) {
	$time = $data[0] * $timestep;
	$_ = <OUTPUT>;
	$_ =~ s/^\s+//; $_ =~ s/\s+$//;
	push(@data,split(/\s+/,$_));
	$_ = <OUTPUT>;
	$_ =~ s/^\s+//; $_ =~ s/\s+$//;
	push(@data,split(/\s+/,$_));
	splice(@data,20,1);
	splice(@data,10,1);
	
	# read rolling averages
	$_ = <OUTPUT>; # empty line
	$_ = <OUTPUT>;
	last if(/r.m.s./);
	$_ =~ s/^\s*(rolling)*\s+//; $_ =~ s/\s+$//;
	push(@data_av,split(/\s+/,$_));
	$_ = <OUTPUT>;
	$_ =~ s/^\s*(averages)*\s+//; $_ =~ s/\s+$//;
	push(@data_av,split(/\s+/,$_));
	$_ = <OUTPUT>;
	$_ =~ s/^\s+//; $_ =~ s/\s+$//;
	push(@data_av,split(/\s+/,$_));
	if($heatup or $data[0]>=$startframe) {
	  printf OUTDATA    "\n%14.4f",$time;
	  printf OUTDATA_AV "\n%14.4f",$time;
	  for($i=1;$i<@data;$i++) {
	    printf OUTDATA    " %11.4e",$data[$i];
	  }
	  for($i=0;$i<@data_av;$i++) {
	    printf OUTDATA_AV " %11.4e",$data_av[$i];
	  }
	}
	$time_old = $time;
      }
    }
  } # end while
  close(OUTPUT);
} # end for directories
close(OUTDATA);
