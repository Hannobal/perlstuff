#!/usr/bin/perl

# changes the ifrz integer in DL_POLY FIELD files to 1 or 0 respectively
# except for molecules with "frozen" or "fix" in the name and user defined exceptions

if($#ARGV<2) {
  print "input format:\n";
  print "1. name of input FIELD file\n";
  print "2. name of output FIELD file\n";
  print "3. (f/u) for freeze/unfreeze\n";
  print "4. list of molecule names to exclude\n   (molecules with \"frozen\" in the name stay frozen)\n";
  exit 1;
}

$inputfilename  = $ARGV[0];
$outputfilename = $ARGV[1];
if($ARGV[2] =~ /^[fu]/i) {
  if($ARGV[2] =~ /^[f]/i) {
    $newifrz = 1;
  } else {
    $newifrz = 0;
  }
} else {
  print "error: could not interpret input \"$ARGV[2]\"!";
  exit 1;
}

if($#ARGV>2) {
  print "exceptions:\n";
  for($i=3;$i<@ARGV;$i++) {
    $exceptions[@exceptions]=$ARGV[$i];
    $foundexception[@foundexception]=0;
    print "$exceptions[$#exceptions]\n";
  }
}

# if($inputfilename eq $outputfilename) {
#   print "error: input FIELD and output FIELD file must not be the same";
#   exit 1;
# }

open(INFIELD, "<", $inputfilename) or die "Can't open input FIELD file $inputfilename:\n$!\n";

while(<INFIELD>) {
  $line[@line] = $_;
  if(/MOLECULES/i) { last; }
}
($field_nummols) = /^\s*\S+\s+(\S+)/;

for($t=0;$t<$field_nummols;$t++) {
  $_=<INFIELD>;
  $line[@line] = $_;
  $mol_name[$t]=$_;
  $mol_name[$t] =~ s/^\s+//;
  $mol_name[$t] =~ s/\s+$//;
  #check whether molecule has "frozen/fix" in the name or is exception
  $isexception = 0;
  if(/frozen/i or /fix/i) {
    $isexception = 1;
  } else {
    for($i=0;$i<@exceptions;$i++) {
      if(/$exceptions[$i]/i) {
	$foundexception[$i] = 1;
	$isexception        = 1;
	last;
      }
    }
  }
  while(<INFIELD>) {
    $line[@line] = $_;
    if(/ATOMS/i and not $isexception) {
      while (<INFIELD>) {
	if(/BONDS/i or /ANGLES/i or /DIHEDRALS/i or /FINISH/i) {
	  $line[@line] = $_;
	  last;
	}
	$_ =~ s/(\s*\S+\s+\S+\s+\S+\s+\S+\s+)\S+/$1$newifrz/;
	$line[@line] = $_;
      }
    }
    if(/FINISH/i) { last; }
  }
}

while(<INFIELD>) {
  $line[@line] = $_;
}

close(INFIELD);

$foundallexceptions=1;
for($i=0;$i<@foundexception;$i++) {
  if(not $foundexception[$i]) {
    $foundallexceptions = 0;
    print "**** warning: molecule \"$exceptions[$i]\" (declared as exception) was not found in file $inputfilename\n";
  }
}

open(OUTFIELD,">",$outputfilename);
print OUTFIELD @line;