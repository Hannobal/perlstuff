#!/usr/bin/perl

use hanno_utility;
use dlpoly_utility;
use Math::Trig;

if($#ARGV<0) {
  print "input format\n",
  " 1. name of Gaussian output file\n",
  "[2. name of xyz file]\n";
  exit 1;
}
$infile = $ARGV[0];
if($#ARGV>0) {
  $outfile = $ARGV[1];
} else {
  $outfile = "$infile.xyz";
  $outfile =~ s/.out//;
}

open($fhin, "<", $infile) or die "**** error: can't open input file:\n$!";

our $foundnames = 0;
our @atomnames  = ();
our @coords     = ();

$frame=0;
while($_=<$fhin>) {
  if(/Multiplicity =/) {
    $_=<$fhin>;
    while(&readnames) {
      $_=<$fhin>;
    }
    $foundnames=1;
    if(not @atomnames) {
      print "**** error reading atom names!\n"; exit 1;
    }
  } elsif(/Coordinates/) {
    if(not @atomnames) {
      print "**** error: no atom names defined yet at first coordinates section!\n"; exit 1;
    }
    $_=<$fhin>; $_=<$fhin>; $_=<$fhin>;
    while(&readcoords($frame)) {
      $_=<$fhin>;
    }
    if($#{$coords[$frame]} != $#atomnames) {
      print "**** error: number of atoms does not fit in frame $frame!\n"; exit 1;
    }
    $frame++;
  }
}

# alignment
$c  = 0; # central atom
$o  = 1; # reference atom
$o2 = 2; # second reference atom
# center first atom in first frame
@shiftvec = (-$coords[0][$c][0],-$coords[0][$c][1],-$coords[0][$c][2]);
move_mol(\@{$coords[0]},\@shiftvec);
for($f=1;$f<$frame;$f++) {
  # center first atom
  $g=$f-1;
  @shiftvec = (-$coords[$f][$c][0],-$coords[$f][$c][1],-$coords[$f][$c][2]);
  move_mol(\@{$coords[$f]},\@shiftvec);
  $angle = dot_product(\@{$coords[$g][$o]},\@{$coords[$f][$o]});
  $angle /= vector_length(@{$coords[$g][$o]})*vector_length(@{$coords[$f][$o]});
  $angle = acos($angle);
#   print "$angle\n";
  next if ($angle==0);
  @axis = vector_product(\@{$coords[$g][$o]},\@{$coords[$f][$o]});
#   print join(" ",@axis),"\n";
  @axis = normalize_vector(\@axis);
  @rotmatrix = gen_rot_matrix(\@axis,-$angle);
  rotate_molecule(\@{$coords[$f]},\@rotmatrix,\@{$coords[$f][$c]});
  #rotate around other axis
  @axis = @{$coords[$f][$o]};
  @axis = normalize_vector(\@axis);
#   printf "axis %8.5f %8.5f %8.5f\n",@axis;
  @tmp = vector_projection(\@{$coords[$g][$o2]},\@axis);
#   printf "tmp  %8.5f %8.5f %8.5f\n",@tmp;
  @vec1 = vector_subst(\@{$coords[$g][$o2]},\@tmp);
  @tmp = vector_projection(\@{$coords[$f][$o2]},\@axis);
  @vec2 = vector_subst(\@{$coords[$f][$o2]},\@tmp);
  $angle = dot_product(\@vec1,\@vec2);
  $angle /= vector_length(@vec1)*vector_length(@vec2);
  $angle = acos($angle);
  print $angle/pi*180,"\n";
  @rotmatrix = gen_rot_matrix(\@axis,$angle);
  rotate_molecule(\@{$coords[$f]},\@rotmatrix,\@{$coords[$f][$c]});
}

open($fhout,">",$outfile) or die "**** error: can't open output file:\n$!";
for($f=($frame-1)%2;$f<$frame;$f+=2) {
  print $fhout $#atomnames+1,"\nframe $f\n";
  for($i=0;$i<@atomnames;$i++) {
    printf $fhout "%5s %10.5f %10.5f %10.5f\n",$atomnames[$i],@{$coords[$f][$i]};
  }
}
close($fhout);

sub readnames {
  my(@linedata);
  $_=~s/^\s+//; $_=~s/\s+$//;
  if(/,/) {
    @linedata=split(",");
  } else {
    @linedata=split(/\s+/);
  }
  return 0 if(not check_real(@linedata[1..3]));
  push(@atomnames,$linedata[0]) if(not $foundnames);
#   push(@{$coords[0]},[@linedata[1..3]]);
  return 1;
}

for($f=0;$f<@frame;$f++) {
  
}

sub readcoords {
  my(@linedata);
  $_=~s/^\s+//; $_=~s/\s+$//;
  @linedata=split(/\s+/);
#   print "$_[0]: $_\n";
  return 0 if(not check_real(@linedata));
  push(@{$coords[$_[0]]},[@linedata[3..5]]);
  return 1;
}