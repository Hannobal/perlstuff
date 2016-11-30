#!/usr/bin/perl

if($#ARGV==-1) {
  print "the input format is:\n";
  print "name of .car- and .mdf-file without ending (must be equivalent)\n";
  exit;
}
#open(CAR, "<", "$ARGV[0].car")  or die "Can't open .car-File: $!";
#open(MDF, "<", "$ARGV[0].mdf")  or die "Can't open .car-File: $!";
open(CAR, "<", "$ARGV[0].car")  or die "Can't open .car-File: $!";
open(MDF, "<", "$ARGV[0].mdf")  or die "Can't open .mdf-File: $!";
open(CONFIG, ">", "CONFIG") or die "Can't create CONFIG-File: $!";
open(FIELD, ">", "FIELD") or die "Can't create FIELD-File: $!";

#reading out the material studio xsd-file------------------------------
for($i=1; $i<=5; $i++) {
  $_ = <CAR>;
}
($length_a, $length_b, $length_c, $alpha, $beta, $gamma, $symmetry) =
/(\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)\s*\((\w*)\)/;
print "$length_a, $length_b, $length_c, $alpha, $beta, $gamma\n";
$moleculeid=1;
$i=0;
$j=0;
while(<CAR>) {
  if(/end/) {
    $moleculeid++;
  } else {
    
  }
}
print "$moleculeid\n";
#print "$j\n";
#for($i=90; $i<$#atomid; $i++) {
#  printf "%-4s ", $type[$i];
#  printf "%-4s ", $name[$i];
#  printf "%12.7f ", $x[$i];
#  printf "%12.7f ", $y[$i];
#  printf "%12.7f ", $z[$i];
#  printf "%12.7f ", $charge[$i];
#  printf "%-10s\n", "$connections[$i]";
#}

#write CONFIG file-----------------------------------------------------
printf CONFIG "%-80s\n", "title";
printf CONFIG "%10u", "0"; # CONFIG file key
printf CONFIG "%10u\n", "6"; # periodic boundary key
printf CONFIG "%20.12e%20.12e%20.12e xyz of Vector1\n", 0,0,0;
printf CONFIG "%20.12e%20.12e%20.12e xyz of Vector2\n", 0,0,0;
printf CONFIG "%20.12e%20.12e%20.12e xyz of Vector3\n", 0,0,0;
for($i=0;$i<=$#oldatomid;$i++) {
  if($name[$i] eq "") {
    printf CONFIG "%-8s%10u\n" ,$type[$i],$i;
  } else {
    printf CONFIG "%-8s%10u\n" ,$name[$i],$i;
  }
  printf CONFIG "%20.12g%20.12g%20.12g\n",$x[$i],$y[$i],$z[$i];
  printf CONFIG "%20.12g%20.12g%20.12g\n",0,0,0;
  printf CONFIG "%20.12g%20.12g%20.12g\n",0,0,0;
}
#write FIELD file------------------------------------------------------
printf FIELD "%-80s\n", "title";
