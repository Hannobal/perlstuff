#!/usr/bin/perl

if($#ARGV!=2 and $#ARGV!=4) {
  print "input format:\n";
  print "1. number of molecules\n";
  print "2. chemical formula\n";
  print "3. volume or box dimensions\n";
  exit 1;
}

$N       = $ARGV[0];
$formula = $ARGV[1];
if($#ARGV==2) {
  $V=$ARGV[2];
} else {
  $V=$ARGV[2]*$ARGV[3]*$ARGV[4];
}

$pi  = 3.14159265358979323846;
$eps = 8.854187817620e-12;
$ec  = 1.602176565e-19;
$l0  = 1e-10;
$t0  = 1e-12;
$u   = 1.660538921e-27;

%mass = (
  "H" => 1.00797,
  "He" => 4.00260,
  "Li" => 6.941,
  "Be" => 9.01218,
  "B" => 10.81,
  "C" => 12.011,
  "N" => 14.0067,
  "O" => 15.9994,
  "F" => 18.998403,
  "Ne" => 20.179,
  "Na" => 22.98977,
  "Mg" => 24.305,
  "Al" => 26.98154,
  "Si" => 28.0855,
  "P" => 30.97376,
  "S" => 32.06,
  "Cl" => 35.453,
  "K" => 39.0983,
  "Ar" => 39.948,
  "Ca" => 40.08,
  "Sc" => 44.9559,
  "Ti" => 47.90,
  "V" => 50.9415,
  "Cr" => 51.996,
  "Mn" => 54.9380,
  "Fe" => 55.847,
  "Ni" => 58.70,
  "Co" => 58.9332,
  "Cu" => 63.546,
  "Zn" => 65.38,
  "Ga" => 69.72,
  "Ge" => 72.59,
  "As" => 74.9216,
  "Se" => 78.96,
  "Br" => 79.904,
  "Kr" => 83.80,
  "Rb" => 85.4678,
  "Sr" => 87.62,
  "Y" => 88.9059,
  "Zr" => 91.22,
  "Nb" => 92.9064,
  "Mo" => 95.94,
  "Tc" => 98,
  "Ru" => 101.07,
  "Rh" => 102.9055,
  "Pd" => 106.4,
  "Ag" => 107.868,
  "Cd" => 112.41,
  "In" => 114.82,
  "Sn" => 118.69,
  "Sb" => 121.75,
  "I" => 126.9045,
  "Te" => 127.60,
  "Xe" => 131.30,
  "Cs" => 132.9054,
  "Ba" => 137.33,
  "La" => 138.9055,
  "Ce" => 140.12,
  "Pr" => 140.9077,
  "Nd" => 144.24,
  "Pm" => 145,
  "Sm" => 150.4,
  "Eu" => 151.96,
  "Gd" => 157.25,
  "Tb" => 158.9254,
  "Dy" => 162.50,
  "Ho" => 164.9304,
  "Er" => 167.26,
  "Tm" => 168.9342,
  "Yb" => 173.04,
  "Lu" => 174.967,
  "Hf" => 178.49,
  "Ta" => 180.9479,
  "W" => 183.85,
  "Re" => 186.207,
  "Os" => 190.2,
  "Ir" => 192.22,
  "Pt" => 195.09,
  "Au" => 196.9665,
  "Hg" => 200.59,
  "Tl" => 204.37,
  "Pb" => 207.2,
  "Bi" => 208.9804,
  "Po" => 209,
  "At" => 210,
  "Rn" => 222,
  "Fr" => 223,
  "Ra" => 226.0254,
  "Ac" => 227.0278,
  "Pa" => 231.0359,
  "Th" => 232.0381,
  "Np" => 237.0482,
  "U" => 238.029,
);

$formulacp=$formula;
$m=0;
while(length($formulacp)>0) {
  if(not $formulacp =~ /^([A-Z][a-z]*)/) {
    print "**** error in chemical formula!\n      (element names must begin with capital letter!)\n";
    exit 1;
  }
  $element=$1;
  $len=length($element);
  if($formulacp =~ /^$element(\d+)/) {
    $anz = $1;
    $len += length($1);
  } else {
    $anz = 1;
  }
  if(not defined($mass{$element})) {
    print "**** error: found unknown element '$element' in formula!\n";
    exit 1;
  }
  $m += $anz * $mass{$element};
  $formulacp=substr($formulacp,$len);
}

$dens=$N*$m*$u/$V*1e27;
print "$dens g/cm^3\n"

# $rho=1e-27;
# $mi=(15.9994+2*1.0079)*$u;
# $A=28.978140323126*28.978140323126;
# $h=$N*$mi/($rho*$A);
# print "$h\n";
