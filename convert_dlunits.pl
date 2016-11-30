#!/usr/bin/perl

# convert the units to/from DLPOLY FIELD file with "kJ" as units

use hanno_utility;

$pi  = 3.14159265358979323846;
$eps = 8.854187817620e-12;
$ec  = 1.602176565e-19;
$l0  = 1e-10;
$t0  = 1e-12;
$u   = 1.660538921e-27;
$factor{"E"} = ($t0/$u)*($t0/$l0)*$ec/100;
$factor{"F"} = ($t0/$u)*($t0/$l0);
$factor{"v"} = ($t0/$l0);
$factor{"a"} = ($t0/$l0)*$t0/100;
$factor{"e"} = 1000*($t0/$l0)^2/$u;

$siunit{"E"} = "V/m";
$siunit{"F"} = "N";
$siunit{"v"} = "m/s";
$siunit{"a"} = "m/s^2";
$siunit{"e"} = "kJ";

$dlunit{"E"} = "u*l0/t0^2/e0";
$dlunit{"F"} = "u*l0/t0^2";
$dlunit{"v"} = "l0/t0";
$dlunit{"a"} = "l0/t0^2";
$dlunit{"e"} = "u*l0^2/t0^2";

if($#ARGV<2) {
  print "The input format is:\n";
  print "1. either of the following:\n";
  print "   E for electric field (V/m   or u*l0/t0^2/e0)\n";
  print "   F for force          (N     or u*l0/t0^2)\n";
  print "   v for velocity       (m/s   or l0/t0)\n";
  print "   a for acceleration   (m/s^2 or l0/t0^2)\n";
  print "   e for energy         (kJ    or u*l0^2/t0^2)\n";
  print "2. 'si' for conversion to SI units OR 'dl' for conv. to dl_poly units\n";
  print "3. the actual number to convert\n";
  print "!!! dlpoly_units using kj as energy base in FIELD file!!!\n";
  exit 1;
}

if(not defined $factor{$ARGV[0]}) {
  print "**** error: unknown conversion keyword '$ARGV[0]'\n";
  exit 1;
} elsif($ARGV[1] =~ /^si$/i) {
  $result = $ARGV[2]/$factor{$ARGV[0]};
  print "$ARGV[2] $dlunit{$ARGV[0]} equals $result $siunit{$ARGV[0]}\n"
} elsif($ARGV[1] =~ /^dl$/i) {
  $result = $ARGV[2]*$factor{$ARGV[0]};
  print "$ARGV[2] $siunit{$ARGV[0]} equals $result $dlunit{$ARGV[0]}\n"
} else {
  print "**** error: unknown conversion keyword '$ARGV[1]'\n";
}
