#!/usr/bin/perl

use dlpoly_utility;
use hanno_utility;

$dist=3;

if($#ARGV<4) {
  print "input format:\n";
  print "1. name of CONFIG/HISTORY file\n";
  print "2. name of corresponding FIELD file\n";
  print "3. name of GAFF file\n";
  print "4. name of output CONFIG file\n";
  print "5. name of output FIELD file\n";
  exit;
}

$incfg = $ARGV[0];
$infld = $ARGV[1];
$ingaff = $ARGV[2];
$outcfg = $ARGV[3];
$outfld = $ARGV[4];

exit 1 if(read_field_file($infld,0) != 0);
exit 1 if(read_config_file($incfg,0,0) != 0);
exit 1 if(read_gaff_file($ingaff) != 0);

print "@{$size[0]}\n";

# remap molecules in z-direction
for($t=0;$t<$field_nummols[0];$t++) {
  for($m=0;$m<$mol_numents[0][$t];$m++) {
    remap_molecule(\@{$cdata[0][$t][$m]},[2],\@{$size[0]});
  }
}

$minpos=9e20;
for($t=0;$t<$field_nummols[0];$t++) {
  for($m=0;$m<$mol_numents[0][$t];$m++) {
    for($a=0;$a<$mol_numatoms[0][$t];$a++) {
#       for($c=0;$c<3;$c++) {
	$minpos = $cdata[0][$t][$m][$a][2] if($cdata[0][$t][$m][$a][2] < $minpos);
#       } # end for c
    } # end for a
  } # end for m
} # end for t

$numx  = int(2*$size[0][0]/$dist);
$numy  = 2*sprintf("%.0f",$size[0][1]/($dist*sqrt(3)/2));
$distx = 2*$size[0][0]/$numx;
$disty = 2*$size[0][1]/$numy;

print "$minpos\n";

push(@{$mol_name[0]},'barrier');
push(@{$mol_numents[0]},$numx*$numy);
push(@{$mol_numatoms[0]},1);
push(@{$mol_atomdata[0]},[['XX',99.9999,0,1,1,0]]);
push(@{$mol_numbonds[0]},0);
push(@{$mol_numangles[0]},0);
push(@{$mol_numdihedrals[0]},0);
push(@{$mol_numinversions[0]},0);
push(@{$mol_numconstraints[0]},0);
push(@{$mol_numatoms[0]},0);
$field_nummols[0]++;
$t    = $#{$mol_name[0]};
$posz = $minpos-5;
$periodic_key[0] = 6;

# update vdw-parameters if necessary
if(not contains(@{$field_atomtypes[0]},'XX')) {
  push(@{$field_atomtypes[0]},'XX');
  $field_numatomtypes[0]++;
  foreach $i ( 0..$#{$field_atomtypes[0]} ) {
    if(not exists $gaff_vdwparam{$field_atomtypes[0][$i]}[0]) {
      print "**** error: did not find entry for ".$field_atomtypes[0][$i]."\n";
      exit;
    }
  }
  for($i=0;$i<@{$field_atomtypes[0]};$i++) {
    loop1:for($j=$i;$j<@{$field_atomtypes[0]};$j++) {
      next if(($field_atomtypes[0][$i] eq 'AL' and $field_atomtypes[0][$j] eq 'AL') or 
              ($field_atomtypes[0][$i] eq 'XX' and $field_atomtypes[0][$j] eq 'XX'));
      for($k=0;$k<@{$field_vdwdata[0]};$k++) {
	if(($field_vdwdata[0][$k][0] eq $field_atomtypes[0][$i] and 
	    $field_vdwdata[0][$k][1] eq $field_atomtypes[0][$j]) or
	   ($field_vdwdata[0][$k][1] eq $field_atomtypes[0][$i] and 
	    $field_vdwdata[0][$k][0] eq $field_atomtypes[0][$j])) {
	  next loop1;
	}
      }
      ($epsilon, $sigma) = calc_vdw_params($field_atomtypes[0][$i], $field_atomtypes[0][$j]);
      push(@{$field_vdwdata[0]},[$field_atomtypes[0][$i],$field_atomtypes[0][$j],'lj',$epsilon,$sigma]);
      $field_numvdw[0]++;
    }
  }
}
# print "$size[0][0] $size[0][1] $size[0][2] $posz\n";
# print "$numx $distx\n";
# print "$numy $disty\n";

for($x=0;$x<$numx;$x++) {
  for($y=0;$y<$numy;$y++) {
    $posx=($x+0.5*($y % 2))*$distx-$size[0][0];
    $posy=$y*$disty-$size[0][1];
    push(@{$cdata[0][$t]},[[$posx,$posy,$posz,0,0,0,0,0,0,'XX']])
  }
}

write_config_file($outcfg,0,$config_title[0]);
write_field_file($outfld,0);
