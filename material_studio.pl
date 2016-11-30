#!/usr/bin/perl

use dlpoly_utility;

if($#ARGV<1) {
  print "the input format is:\n";
  print "1. name of .xsd-file with ending\n";
  print "2. name of output CONFIG file\n";
  exit;
}
open(XSD, "<", $ARGV[0])  or die "Can't open .xsd-File: $!";

#reading out the material studio xsd-file------------------------------
while(<XSD>) {
  if(/<IdentityMapping/) { last; }
}

$i=0;
$j=0;
$perioic_key[0] = 999;
$zmin=999999;
while(<XSD>) {
  if(/Atom3d/ and /Name/) {
#     ($nid) = /ID=\"\s*(\d+)\s*\"/;
    ($cdata[0][0][0][$i][9]) = /Name=\"\s*(\S+)\s*\"/;
    ($cdata[0][0][0][$i][0],$cdata[0][0][0][$i][1],$cdata[0][0][0][$i][2])
      =  /XYZ=\"\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*\"/;
    ($data[$i][6]) = /Connections=\"(\S*)\"/;
    $i++;
  } elsif(/<SpaceGroup/) {
   $periodic_key[0] = 6;
   ($cell[0][0][0], $cell[0][0][1], $cell[0][0][2]) = /AVector=\"\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*\"/;
   ($cell[0][1][0], $cell[0][1][1], $cell[0][1][2]) = /BVector=\"\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*\"/;
   ($cell[0][2][0], $cell[0][2][1], $cell[0][2][2]) = /CVector=\"\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*\"/;
  } elsif(/<PlaneGroup/) {
   $periodic_key[0] = 6;
   ($cell[0][0][0], $cell[0][0][1], $cell[0][0][2]) = /AVector=\"\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*\"/;
   ($cell[0][1][0], $cell[0][1][1], $cell[0][1][2]) = /BVector=\"\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*\"/;
   ($cell[0][2][0], $cell[0][2][1], $cell[0][2][2]) = /CVector=\"\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*\"/;
  }
}
close(XSD);
print $i,"\n";

for($i=0;$i<@{$cdata[0][0][0]};$i++) {
  $newpos[0] = $cdata[0][0][0][$i][0]*$cell[0][0][0] + $cdata[0][0][0][$i][1]*$cell[0][1][0] + $cdata[0][0][0][$i][2]*$cell[0][2][0];
  $newpos[1] = $cdata[0][0][0][$i][0]*$cell[0][0][1] + $cdata[0][0][0][$i][1]*$cell[0][1][1] + $cdata[0][0][0][$i][2]*$cell[0][2][1];
  $newpos[2] = $cdata[0][0][0][$i][0]*$cell[0][0][2] + $cdata[0][0][0][$i][1]*$cell[0][1][2] + $cdata[0][0][0][$i][2]*$cell[0][2][2];
  @{$cdata[0][0][0][$i]}[0..2] = @newpos;
}

# sort data by atom name and z-coordinate in ascending order
@{$cdata[0][0][0]} = sort { $a->[2] <=> $b->[2] || $a->[9] cmp $b->[9] } @{$cdata[0][0][0]};

write_config_file($ARGV[1],0,'title');

#write FIELD file------------------------------------------------------
#printf FIELD "%-80s\n", "title"
