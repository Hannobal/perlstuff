#!/usr/bin/perl

#analyzes how long the execution of a command takes

if($#ARGV==-1) {
  print "enter commands to execute!\n";
  exit 1;
}

$starttime = time();
foreach $command (values @ARGV) {
  print "$command\n";
  system($command);
  exit 1 if($?!=0);
}
$timeelapsed = time2human(time()-$starttime);
print "time elapsed: $timeelapsed\n";

sub time2human {
  my $time = $_[0];
  my $s = $time%60;
  my $m = int(($time/60)%60);
  my $h = int(($time/3600)%24);
  my $d = int($time/86400);
  if($d>1) {
    return sprintf "%s days %02u:%02u:%02u",$d,$h,$m,$s;
  }elsif($d>0) {
    return sprintf "%s day %02u:%02u:%02u",$d,$h,$m,$s;
  }elsif($h>0) {
    return sprintf "%02u:%02u:%02u",$h,$m,$s;
  }elsif($m>0) {
    return sprintf "%02u:%02u",$m,$s;
  }else{
    return "$s sec";
  }
}