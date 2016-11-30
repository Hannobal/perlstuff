#!/usr/bin/perl

use hanno_utility;
use Storable qw(dclone);

if($#ARGV<2) {
  print "1. name of input data-file\n";
  print "2. name of output output-file\n";
  print "3. normal to plane (x/y/z)\n";
  print "4. height value (nearest will be chosen)";
  print "[filters, e.g. x>-15 x<15]";
  print "  keywords for filters: x,y,z,vx,vy,vz,vabs";
  exit;
}

$infilename  = $ARGV[0];
$outfilename = $ARGV[1];
$targetpos   = $ARGV[3];
@filters     = @ARGV[4..$#ARGV];
if($ARGV[2] =~ /^x/i) {
  $surfnorm=0;
} elsif($ARGV[2] =~ /^y/i) {
  $surfnorm=1;
} elsif($ARGV[2] =~ /^z/i) {
  $surfnorm=2;
} else {
  print "**** error: input for surface normal could not be interpreted!\n";
  exit 1;
}
for($i=0;$i<@filters;$i++) {
  $filters[$i] =~ s/vx/\$line[3]/i;
  $filters[$i] =~ s/vy/\$line[4]/i;
  $filters[$i] =~ s/vz/\$line[5]/i;
  $filters[$i] =~ s/vabs/\$line[6]/i;
  $filters[$i] =~ s/x/\$line[0]/i;
  $filters[$i] =~ s/y/\$line[1]/i;
  $filters[$i] =~ s/z/\$line[2]/i;
  print "$filters[$i]\n";
}

# read data, apply filters and determine nearest plane
open(INFILE,"<",$infilename) or die "**** error: could not open file $infilename: $!\n";
$ndist=9e20;
$n=0;
while(<INFILE>) {
  $_ =~ s/^\s+//; $_ =~ s/\s+$//;
  @line=split(/\s+/,$_);
  if(check_real($line[0])) {
    $dist = abs($line[$surfnorm]-$targetpos);
    if($dist<$ndist) {
      $ndist  = $dist;
      $npos   = $line[$surfnorm];
    }
    $add=1;
    for($i=0;$i<@filters;$i++) {
      $add=0 if(not eval($filter[$i]));
    }
    @{$data[@data]} = @line if($add);
  }
  $n++;
}
close(INFILE);

@planenames = ('x','y','z');
print "taking $planenames[$surfnorm]-plane with value $npos\n";

# extract and sort data in plane
@data2=();
for($i=0;$i<@data;$i++) {
  push(@data2,\@{dclone \@{$data[$i]}}) if($data[$i][$surfnorm]==$npos);
}
if($surfnorm==0) {
  @data2 = sort { $a->[2] <=> $b->[2] || $a->[1] <=> $b->[1] } @data2;
  $breakindex=2;
}elsif($surfnorm==1) {
  @data2 = sort { $a->[2] <=> $b->[2] || $a->[0] <=> $b->[0] } @data2;
  $breakindex=2;
}else {
  @data2 = sort { $a->[1] <=> $b->[1] ||  $a->[0] <=> $b->[0] } @data2;
  $breakindex=1;
}

#print output

open(OUTFILE,">",$outfilename) or die "**** error: could not open file $outfilename: $!\n";
$old=$data2[0][$breakindex];
for($i=0;$i<@data2;$i++) {
  if($old != $data2[$i][$breakindex]) {
    print OUTFILE "\n";
    $old=$data2[$i][$breakindex];
  }
  printf OUTFILE "%10.5f %10.5f %10.5f", @{$data2[$i]};
  for($j=3;$j<=$#{$data2[$i]};$j++) {
    printf OUTFILE " %20.13e", $data2[$i][$j];
  }
  print OUTFILE "\n";
}
close(OUTFILE);
