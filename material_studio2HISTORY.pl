#!/usr/bin/perl

if($#ARGV<1) {
  print "the input format is:\n";
  print "1. name of .xsd-file with ending\n";
  print "2. name of output HISTORY-file\n";
  exit;
}
open(XSD, "<", $ARGV[0])  or die "Can't open .xsd-File: $!";
open(CONFIG, ">", $ARGV[1]) or die "Can't create HISTORY-File: $!";
#open(FIELD, ">", "FIELD") or die "Can't create FIELD-File: $!";

#reading out the material studio xsd-file------------------------------
<XSD>;
until(/<IdentityMapping/) {
  $_ = <XSD>;
}
$i=0;
$j=0;
$period_key = 999;
%mass = (
  "Al", 26.9815,
  "O",  15.999
);
while(<XSD>) {
  if(/<\/IdentityMapping>/) {
    last;
  } elsif(/Atom3d/) {
    $i++;
    ($nid) = /ID=\"(\d*)\"/;
    $newatomid[$nid]=$i;
    $oldatomid[$i]=$nid;
    ($name[$i]) = /Name=\"(\S*)\"/;
    ($x[$i], $y[$i], $z[$i]) =  /XYZ=\"(\S*),(\S*),(\S*)\"/;
    ($connections[$i]) = /Connections=\"(\S*)\"/;
    ($type[$i], $charge[$i]) = /Components=\"(\S*),(\S*)\S*\"/;
  }elsif(/Bond/) {
    $j++;
    ($nid) = /ID=\"(\d*)\"/;
    $bondid[$nid]=$j;
    ($at1[$j], $at2[$j]) = /Connects=\"(\S*),(\S*)\"/;
  } elsif(/<SpaceGroup/) {
   $period_key = 2;
   ($a[0], $a[1], $a[2]) = /AVector=\"(\S*),(\S*),(\S*)\"/;
   ($b[0], $b[1], $b[2]) = /BVector=\"(\S*),(\S*),(\S*)\"/;
   ($c[0], $c[1], $c[2]) = /CVector=\"(\S*),(\S*),(\S*)\"/;
  } elsif(/<PlaneGroup/) {
   $period_key = 2;
   ($a[0], $a[1], $a[2]) = /AVector=\"(\S*),(\S*),(\S*)\"/;
   ($b[0], $b[1], $b[2]) = /BVector=\"(\S*),(\S*),(\S*)\"/;
   ($c[0], $c[1], $c[2]) = /CVector=\"(\S*),(\S*),(\S*)\"/;
  }
}
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
printf CONFIG "%10u", $period_key; # periodic boundary key
printf CONFIG "%10u\n", $i; # number of atoms
printf CONFIG "%-8s","timestep";
printf CONFIG "%10u",0; #current timestep
printf CONFIG "%10u",$i; #number of atoms
printf CONFIG "%10u", "0"; # CONFIG file key
printf CONFIG "%10u", "6"; # periodic boundary key
printf CONFIG "%12.6f\n", "0.001"; # integration timestep
printf CONFIG "%12.4g%12.4g%12.4g\n", $a[0], $a[1], $a[2];
printf CONFIG "%12.4g%12.4g%12.4g\n", $b[0], $b[1], $b[2];
printf CONFIG "%12.4g%12.4g%12.4g\n", $c[0], $c[1], $c[2];
for($i=1;$i<=$#oldatomid;$i++) {
  $new_pos[0] = $x[$i]*$a[0] + $y[$i]*$b[0] + $z[$i]*$c[0];
  $new_pos[1] = $x[$i]*$a[1] + $y[$i]*$b[1] + $z[$i]*$c[1];
  $new_pos[2] = $x[$i]*$a[2] + $y[$i]*$b[2] + $z[$i]*$c[2];
  if($name[$i] eq "") {
    printf CONFIG "%-8s" ,$type[$i];
  } else {
    printf CONFIG "%-8s" ,$name[$i];
  }
  printf CONFIG "%10u",$i; #atom index
  printf CONFIG "%12.6f",$mass{$type[$i]}; #atomic mass
  printf CONFIG "%12.6f\n",$charge[$i]; #atomic charge
#  printf CONFIG " %11.4e %11.4e %11.4e\n",$x[$i],$y[$i],$z[$i];
  printf CONFIG " %11.4e %11.4e %11.4e\n",$new_pos[0],$new_pos[1],$new_pos[2];
#  printf CONFIG "%20.12g%20.12g%20.12g\n",0,0,0; # xyz of velocity
#  printf CONFIG "%20.12g%20.12g%20.12g\n",0,0,0; # xyz of force
}
#write FIELD file------------------------------------------------------
#printf FIELD "%-80s\n", "title"