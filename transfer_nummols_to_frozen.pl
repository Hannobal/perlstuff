#!/usr/bin/perl

# changes the number of molecules in a DL_POLY FIELD file so that molecules
# are either transferred to their frozen or unfrozen counterparts respectively

if($#ARGV<3) {
  print "input format:\n";
  print "1. name of input FIELD file\n";
  print "2. name of output FIELD file\n";
  print "list of molecule names and number of molecules to transfer to frozen\n";
  print "      (put negative numbers to unfreeze)\n";
  print "example: ... input.fld output.fld \"my molecue\" 2 \"another molecule\" -5\n";
  exit 1;
}

$inputfilename  = $ARGV[0];
$outputfilename = $ARGV[1];

if(@argv % 2 != 0) {
  print "**** error: uneven number of command arguments!";
  exit 1;
}

$j=0;
for($i=2;$i<@ARGV;$i+=2) {
  $transfer_mol[$j] = $ARGV[$i];
  $transfer_num[$j] = $ARGV[$i+1];
  $found_free[$j]   = 0;
  $found_frozen[$j] = 0;
  $j++;
}

# if($inputfilename eq $outputfilename) {
#   print "**** error: input FIELD and output FIELD file must not be the same";
#   exit 1;
# }

open(INFIELD, "<", $inputfilename) or die "Can't open input FIELD file $inputfilename:\n$!\n";

$error=0;

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
  $mol_addnuments[$t] = 0; # will be the index for the @transfer_mol array
  for($i=0;$i<@transfer_mol;$i++) {
    if(/$transfer_mol[$i]/i) {
      if(/frozen/i or /fix/i) {
	$found_frozen[$i] = 1;
	$sign             = 1; # the sign whether to add or subtract the number
      } else {
	$found_free[$i]  =  1;
	$sign            = -1;
      }
      $mol_addnuments[$t] = $sign*$transfer_num[$i];
      last;
    }
  }
  $_=<INFIELD>;
  ($mol_numents[$t])=/^\s*\S+\s+(\S+)/;
  $mol_numents[$t] += $mol_addnuments[$t];
  if($mol_numents[$t]<0) {
    print "**** error: number of molecules for molecule \"$mol_name[$t]\" ($mol_numents[$t]) is smaller than zero!\n";
    $error=1;
  }
  $line[@line] = "NUMMOLS $mol_numents[$t]\n";
  while(<INFIELD>) {
    $line[@line] = $_;
    if(/FINISH/i) { last; }
  }
}

while(<INFIELD>) {
  $line[@line] = $_;
}

close(INFIELD);

for($i=0;$i<@found_frozen;$i++) {
  if(not @found_frozen[$i]) {
    print "**** error: frozen version of molecule \"$transfer_mol[$i]\" was not found in file $inputfilename\n";
    $error = 1;
  }
  if(not @found_free[$i]) {
    print "**** error: free   version of molecule \"$transfer_mol[$i]\" was not found in file $inputfilename\n";
    $error = 1;
  }
}

exit 1 if($error);

open(OUTFIELD,">",$outputfilename);
print OUTFIELD @line;