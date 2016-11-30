#!/usr/bin/perl

if($#ARGV<1) {
  print "1. name of HISTORY-file\n";
  print "2. name of XYZ-file\n";
  print "[4. start timestep]\n";
  print "[5. end timestep]\n";
  print "[6. step timestep]\n";
  exit;
}

if($#ARGV>1) {
  $startframe = $ARGV[2];
} else {
  $startframe = 0;
}
if($#ARGV>2) {
  $endframe = $ARGV[3];
} else {
  $endframe = 999999999;
}
if($#ARGV>3) {
  $step = $ARGV[4];
} else {
  $step = 1;
}

if(not (check_integer($startframe) and check_integer($endframe) and check_integer($step))) {
  print "**** error: start-, endframe and frame step must be integer values!\n";
  exit 1;
}

open(HISTORY, "<", $ARGV[0]) or die "Can't open HISTORY-File: $!";
open(XYZ,">",$ARGV[1]);

# necessary step because some HISTORY files have a title line and some don't:
$_=<HISTORY>;
if(not /^\s*timestep/) {
  $title=$_;
  $title =~ s/\s+$//;
  $_=<HISTORY>;
} else {
  $title="title";
  seek(HISTORY,0,0); #rewind file;
}

while(<HISTORY>) {
  if( not /timestep/) {
    print "error reading HISTORY-file: \"timestep\" expected on line:\n$_";
    exit 1;
  }
  ($framenumber, $numatoms, $config_key, $periodic_key, $timestep) =
  /timestep\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
  print "converting frame number $framenumber\n";
  if($framenumber>$startframe and $framenumber<$endframe and ($framenumber-$startframe)%$step==0) {
    if($periodic_key>0) {
      for($j=0;$j<=2;$j++) {
	$_=<HISTORY>;
	($vec[$j][0], $vec[$j][1], $vec[$j][2]) = /\s+(\S+)\s+(\S+)\s+(\S+)/;
      }
    }
    print XYZ "$numatoms\n";
    print XYZ "$title frame $framenumber\n";
    for($i=0;$i<$numatoms;$i++) {
      &read_atom($i);
      printf XYZ "%-3s", name_atom($cdata[$i][0]);
      printf XYZ  " %19.12e %19.12e %19.12e\n", @{$cdata[$i]}[1..3];
    }
  } else {
    $imax = $numatoms*($config_key+2); # all lines for atoms
    $imax +=3 if($periodic_key>0);
    for($i=0;$i<$imax;$i++) {
      $_ = <HISTORY>;
    }
  }
}
close(HISTORY, XYZ);

# for($i=0;$i<$anzges;$i++) {
#   printf XYZ "\n%-3s", $type{lc($cdata[$i][0])};
#   printf XYZ  " %19.12e %19.12e %19.12e", @{$cdata[$i]}[1..3];
#   if($print_key==2) {
#     if($divide) {
#       if($cdata[$i][11] == 0) {
# 	printf XYZ " %19.12e", 0;
# 	printf XYZ " %19.12e", 0;
# 	printf XYZ " %19.12e", 0;
#       } else {
# 	printf XYZ " %19.12e", $cdata[$i][7]/$cdata[$i][11]*$scalefactor;
# 	printf XYZ " %19.12e", $cdata[$i][8]/$cdata[$i][11]*$scalefactor;
# 	printf XYZ " %19.12e", $cdata[$i][9]/$cdata[$i][11]*$scalefactor;
#       }
#     } else {
#       printf XYZ " %19.12e", $cdata[$i][7]*$scalefactor;
#       printf XYZ " %19.12e", $cdata[$i][8]*$scalefactor;
#       printf XYZ " %19.12e", $cdata[$i][9]*$scalefactor;
#     }
#   } elsif($print_key==1) {
#     printf XYZ " %19.12e", $cdata[$i][4]*$scalefactor;
#     printf XYZ " %19.12e", $cdata[$i][5]*$scalefactor;
#     printf XYZ " %19.12e", $cdata[$i][6]*$scalefactor;
#   }
# }


sub read_atom {
  my $i = $_[0];
  $_=<HISTORY>;
  ($cdata[$i][0],$cdata[$i][10],$cdata[$i][11]) = /\s*(\S+)\s+\S+\s+(\S+)\s+(\S+)/;
  $_=<HISTORY>;
  ($cdata[$i][1],$cdata[$i][2],$cdata[$i][3]) = /\s*(\S+)\s+(\S+)\s+(\S+)/;
  if($config_key>0) {
    $_ = <HISTORY>;
    ($cdata[$i][4],$cdata[$i][5],$cdata[$i][6]) = /\s*(\S+)\s+(\S+)\s+(\S+)/;
    if($config_key>1) {
      $_ = <HISTORY>;
      ($cdata[$i][7],$cdata[$i][8],$cdata[$i][9]) = /\s*(\S+)\s+(\S+)\s+(\S+)/;
    } else {
      ($cdata[$i][7],$cdata[$i][8],$cdata[$i][9]) = (0,0,0);
    }
  } else {
    ($cdata[$i][4],$cdata[$i][5],$cdata[$i][6]) = (0,0,0);
    ($cdata[$i][7],$cdata[$i][8],$cdata[$i][9]) = (0,0,0);
  }
}

sub check_integer {
  my $check= 1;
  for(my $i=0;$i<@_;$i++) {
    if (not $_[$i] =~ /^[+-]?\d+$/ ) {
      $check = 0;
      break;
    }
  }
  return $check;
}

sub name_atom {
  my $namein = $_[0];
  my $nameout = "X";
  if($namein =~ /^Cl/) {
    $nameout = "Cl";
  } elsif($namein =~ /^Cl/i) {
    $nameout = "Cl";
  } elsif($namein =~ /^Na/i) {
    $nameout = "Na";
  } elsif($namein =~ /^C/i) {
    $nameout = "C";
  } elsif($namein =~ /^N/i) {
    $nameout = "N";
  } elsif($namein =~ /^H/i) {
    $nameout = "H";
  } elsif($namein =~ /^O/i) {
    $nameout = "O";
  } elsif($namein =~ /^S/i) {
    $nameout = "S";
  } elsif($namein =~ /^P/i) {
    $nameout = "P";
  } elsif($namein =~ /^F/i) {
    $nameout = "F";
  } else {
    print "atom name $namein not known!\n";
    exit 1;
  }
}