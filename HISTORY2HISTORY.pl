#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;
use Switch;

# join restarted DLPOLY HISTORY files

if($#ARGV==-1) {
  print "the input format is:\n";
  print " 1. <output-filename>\n";
  print " optional flags:\n";
  print "-a           append output to file instead of overwrite\n";
  print "-ca <3*real> reset cell vector a\n";
  print "-cb <3*real> reset cell vector b\n";
  print "-cc <3*real> reset cell vector c\n";
  print "-fixyz       use xyz input format (does not work yet)\n";
  print "-foxyz       use xyz output format (does not work yet)\n";
  print "-fovel       use vel output format\n";
  print "-fld <str>   specify field file\n";
  print "-f <2*int>   specify start and end frame\n";
  print "-i  <str>    input filename (default: HISTORY)\n";
  print "-k [012]     new config_key\n";
  print "-n <n*str>   list of molecule names to include\n";
  print "-r1          remap molecules according to cell in first frame\n";
  print "-r [xyz]+    remap molecules in x,y and/or z-direction (all frames)?\n";
  print "-t <n*str>   list of atom types to include\n";
  print "-g           rotate system to make cell match VMD reqirements\n";
  exit 1;
}

$outhist  = $ARGV[0];

$histname = "HISTORY";
$fldname  = "FIELD";
@remap       = ();
$writemode   = ">";
$startframe  = 0;
$endframe    = 9e20;
%newcell     = ();
$xyzin       = 0;
$xyzout      = 0;
$velout      = 0;
$lincludetypes = 0;
$lincludenames = 0;
$lrotatevmd    = 0;
$fldreaderror  = 0;
$lremapfirst   = 0;

for($i=1;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-a/) {
      $writemode=">>";
    } case (/^-fixyz/) {
      $xyzin  = 1;
    } case (/^-foxyz/) {
      $xyzout = 1;
    } case (/^-fovel/) {
      $velout = 1;
    } case (/^-fld/) {
      $i++;
      $fldname = $ARGV[$i];
    } case (/^-f$/) { # frame range
      $i++;
      $startframe = $ARGV[$i];
      &error(1) if(not check_integer($startframe));
      $i++;
      $endframe = $ARGV[$i];
      &error(2) if(not check_integer($endframe));
    } case (/^-ca/) { # reset cell vectors
      for($j=0;$j<3;$j++) {
	$i++;
	&error(9) if(not check_real($ARGV[$i]));
	$newcell{0}[$j] = $ARGV[$i];
      }
    } case (/^-cb/) { # reset cell vectors
      for($j=0;$j<3;$j++) {
	$i++;
	&error(9) if(not check_real($ARGV[$i]));
	$newcell{1}[$j] = $ARGV[$i];
      }
    } case (/^-cc/) { # reset cell vectors
      for($j=0;$j<3;$j++) {
	$i++;
	&error(9) if(not check_real($ARGV[$i]));
	$newcell{2}[$j] = $ARGV[$i];
      }
    } case (/^-i/) { # input file
      $i++;
      $histname = $ARGV[$i];
    } case (/^-k/) { # config_key
      $i++;
      $newconfkey = $ARGV[$i];
      error(8) if(not $newconfkey =~ /^[012]$/);
    } case (/^-n/) { # include molecule names
      $lincludenames = 1;
      while(not $ARGV[$i+1] =~ /^-/ and $i<$#ARGV) {
	$i++;
	push(@includemols,$ARGV[$i]);
      }
    } case (/^-r1/ or /^--remapfirst/) {
      $lremapfirst=1;
    } case (/^-r/ or /^--remap/) {
      $i++;
      & error(3,$ARGV[$i]) if(not $ARGV[$i] =~ /^[xyz]+$/);
      if($ARGV[$i] =~ /x/i) {
	push(@remap,0);
      }
      if($ARGV[$i] =~ /y/i) {
	push(@remap,1);
      }
      if($ARGV[$i] =~ /z/i) {
	push(@remap,2);
      }
    } case (/^-t/) { # include types
      $lincludetypes = 1;
      while(not $ARGV[$i+1] =~ /^-/ and $i<$#ARGV) {
	$i++;
	push(@includetypes,uc($ARGV[$i]));
      }
    } case (/^-g/) { # rotate cell to match VMD requirements
      $lrotatevmd = 1;
    } 
  }
}

if($histname=~/^(.*)\/(.*)$/) {
  @directories=($1);
  $histname=$2;
  print $directories[0],"   ",$histname,"\n";
} else {
  @directories = get_numeric_directories(".",$startframe,$endframe);
  push(@directories,'.') if(-e "./$histname");
  &error(7) if(not @directories);

}

open($fhhistout, $writemode, $outhist) or die "Cant open $outhist: $!";

if($fldname=~/^(.*)\/(.*)$/) {
  $fldreaderror=read_field_file($fldname,0);
} else {
  $fldreaderror=read_field_file("$directories[0]/$fldname",0);
}

if($fldreaderror and ($lincludenames or @remap or $lremapfirst)) {
    &error(10);
}

$error = 0;
$firsthist  = 1;
$firstframe = 1;
foreach $dir (@directories) {
  exit 1 if(not open($fhhistin,"<","$dir/$histname"));
  if($startframe>0 and $firsthist) {
    &error(4,$startframe) if(not find_history_timestep($fhhistin,$startframe));
  }
  $firsthist=0;

  $err=0;
  while($err==0) {
    if($fldreaderror) {
      $err=read_history_timestep($fhhistin,-1,0);
    } else {
      $err=read_history_timestep($fhhistin,0,0);
    }
    last   if($err<0);
    exit 1 if($err>0);
    last if($frame_number[0]>$endframe);
    print "$frame_number[0]\n";
    #remap first frame
    if($lremapfirst and $firstframe) {
      print "remapping first frame\n";
      if($periodic_key[0] != 0) {
	for($t=0;$t<$field_nummols[0];$t++) {
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
      }
      $firstframe=0;
    }
    $config_key[0]=$newconfkey if(defined($newconfkey));
    foreach $key (keys %newcell) {
      @{$cell[0][$key]} = @{$newcell{$key}};
    }
    #remap if desired------------------------------------------------
    if(@remap>0) {
      if($periodic_key[0] != 1 and $periodic_key[0] != 2 and $periodic_key[0] != 6) {
	for($t=0;$t<$field_nummols[0];$t++) {
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
      }
      elsif($periodic_key[0]==6 and contains(@remap,2)) {
	&error(6)
      } else {
	for($t=0;$t<$field_nummols[0];$t++) {
	  if($mol_name[0][$t]=~/PA/) {
	    for($m=0;$m<$mol_numents[0][$t];$m++) {
	      remap_molecule(\@{$cdata[0][$t][$m]},\@remap,\@{$size[0]},10);
	    }
	  } else {
	    for($m=0;$m<$mol_numents[0][$t];$m++) {
	      remap_molecule(\@{$cdata[0][$t][$m]},\@remap,\@{$size[0]});
	    }
	  }
	}
      }
    }
    if($lincludetypes) {
      for($t=0;$t<$field_nummols[0];$t++) {
	next if(contains(@includemols,$mol_name[0][$t]));
	remove_mol_type(0,0,$t);
	$t--;
      }
    }
    if($lincludeatoms) {
      for($t=0;$t<@{$cdata[$ci]};$t++) {
	for($m=0;$m<@{$cdata[$ci][$t]};$m++) {
	  for($a=0;$a<@{$cdata[$ci][$t][$m]};$a++) {
	    next if(contains(@includetypes,$cdata[0][$t][$m][$a][9]));
	    splice(@{$cdata[0][$t][$m]},$a,1);
	    $a--;
	  }
	}
      }
    }
    if($lrotatevmd) {
      rotate_cell_vmd(0);
    }
    if($xyzout) {
      write_xyz_timestep($fhhistout,0,0);
    } elsif($velout) {
      print $fhhistout $frame_number[0]," ",$frame_number[0]*$frame_timestep[0],"\n";
      for($t=0;$t<@{$cdata[$ci]};$t++) {
	for($m=0;$m<@{$cdata[$ci][$t]};$m++) {
	  for($a=0;$a<@{$cdata[$ci][$t][$m]};$a++) {
	  printf $fhhistout "%15u %17.10g %17.10g %17.10g\n",@{$cdata[0][$t][$m][$a]}[3..5];
	  } # end for $a
	} # end for $m
      } # end for $t
    } else {
      write_history_timestep($fhhistout,0,0);
    }
  }
  close($fhhistin);
}
close($fhhistout);

sub error {
  if($_[0]==1) {
    print "**** error: start frame must be an integer number!\n";
  } elsif($_[0]==2) {
    print "**** error: end frame must be an integer number!\n";
  } elsif($_[0]==3) {
    print "**** error: input $_[1] is not valid for remapping!\n";
  } elsif($_[0]==4) {
    print "**** error: did not find timestep $_[1] in $dir/HISTORY!\n";
  } elsif($_[0]==5) {
    print "**** error: remapping is only possible for periodic systems with orthogonal cells!\n";
  } elsif($_[0]==6) {
    print "**** error: remapping in z-direction is only possible for 3d-periodic systems!\n";
  } elsif($_[0]==7) {
    print "**** error: no suitable directories found!\n";
  } elsif($_[0]==8) {
    print "**** error: config_key must be '0', '1' or '2'!\n";
  } elsif($_[0]==9) {
    print "**** error: cell vectors must be specified as three real numbers!\n";
  } elsif($_[0]==10) {
    print "**** error: unable to read required FIELD file!\n";
  } else {
    print "**** unknown error!\n";
    exit 99;
  }
  exit $_[0]
}