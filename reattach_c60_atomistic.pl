#!/usr/bin/perl

use Switch;
use dlpoly_utility;
use hanno_utility;

use Math::Trig;

if($#ARGV<4) {
  print " 1. input CONFIG-file\n";
  print " 2. input FIELD-file with C60-dummies\n";
  print " 3. input FIELD file with atomistic C60\n";
  print " 4. mol2 file with atomistic C60\n";
  print " 5. output CONFIG-file\n";
  exit 1;
}

$incfg  = $ARGV[0];
$dumfld = $ARGV[1];
$atfld  = $ARGV[2];
$mol2   = $ARGV[3];
$outcfg = $ARGV[4];

&error(1) if (read_field_file($dumfld,0)!=0);
&error(1) if (read_field_file($atfld, 1)!=0);
&error(1) if (read_mol2_file($mol2, 0)!=0);
&error(1) if (read_config_file($incfg, 0, 0)!=0);

$tc60=-1;
for($t=0;$t<$field_nummols[0];$t++) {
  &error(2,$t) if($mol_name[0][$t] != $mol_name[1][$t]);
  if($mol_name[0][$t] eq $mol2_name[0]) {
    $tc60=$t;
  } else {
    for($a=0;$a<$mol_numatoms[0][$t];$a++) {
      &error(3,$t,$a) if($mol_atomdata[0][$t][$a][0] ne $mol_atomdata[1][$t][$a][0]);
    }
  }
}

&error(4) if($tc60<0);

# match the cx of molecule with dummy to cx of atomistic molecule
$dumcx  = -1;
$dumful = -1;
$dumref = -1;
for($a=0;$a<$mol_numatoms[0][$tc60];$a++) {
  if($mol_atomdata[0][$tc60][$a][0] =~ /^CX$/i) {
    foreach $b (@{$mol_bondatoms[0][$tc60][$a]}) {
      if($mol_atomdata[0][$tc60][$b][0] =~ /^ful$/i) {
	$dumcx=$a;
	$dumful=$b;
      } else {
	$dumref=$b;
      }
    }
  }
}
&error(5) if($dumcx<0 or $dumful<0);

$atcx  = -1;
$atref = -1;
@c60   = ();
@nonc60 = ();
findatcx:for($a=0;$a<$mol_numatoms[1][$tc60];$a++) {
  if($mol_atomdata[1][$tc60][$a][0] =~ /^CX$/i) {
    foreach $b (@{$mol_bondatoms[1][$tc60][$a]}) {
      if($mol_atomdata[1][$tc60][$b][0] =~ /^CA$/i) {
	push(@c60,$a);
	next findatcx;
      }
    }
    foreach $b (@{$mol_bondatoms[1][$tc60][$a]}) {
      if($mol_atomdata[1][$tc60][$b][0] !~ /^CX$/i) {
	$atref=$b;
      }
    }
    $atcx=$a;
    push(@nonc60,$a);
  } elsif($mol_atomdata[1][$tc60][$a][0] =~ /^CA$/i) {
    push(@c60,$a);
  } else {
    push(@nonc60,$a);
  }
}
&error(6) if($atcx<0);
&error(7,int(@c60)) if(int(@c60)!=60);

# match non-c60 atoms
@match = ();
$n=0;
for($a=0;$a<$mol_numatoms[1][$tc60];$a++) {
  next if($mol_atomdata[1][$tc60][$a][0] =~ /^CA$/i);
  next if($mol_atomdata[1][$tc60][$a][0] =~ /^CX$/i and $a!=$atcx);
  $match[$a]=-1;
  while($match[$a]<0 and $n<$mol_numatoms[0][$tc60]) {
    $match[$a]=$n if($mol_atomdata[0][$tc60][$n][0] eq $mol_atomdata[1][$tc60][$a][0]);
    $n++;
  }
  &error(8,$a) if($match[$a]<0);
}

# calc center postion of atomistic fullerene in mol2
@atpful=(0,0,0);
foreach $a (@c60) {
  $atpful[0] += $mol2_atomdata[0][$a][0];
  $atpful[1] += $mol2_atomdata[0][$a][1];
  $atpful[2] += $mol2_atomdata[0][$a][2];
}
$atpful[0] /= 60;
$atpful[1] /= 60;
$atpful[2] /= 60;
# calc distance of atomistic fullerene to CX in mol2
$atvful[0] = $atpful[0]-$mol2_atomdata[0][$atcx][0];
$atvful[1] = $atpful[1]-$mol2_atomdata[0][$atcx][1];
$atvful[2] = $atpful[2]-$mol2_atomdata[0][$atcx][2];
# $lat = vector_length(@atvful);
# @atvfuln = normalize_vector(@atvful);

# calc reference vector for atomistic molecule in mol2
$atvref[0] = $mol2_atomdata[0][$atref][0] - $mol2_atomdata[0][$atcx][0];
$atvref[1] = $mol2_atomdata[0][$atref][1] - $mol2_atomdata[0][$atcx][1];
$atvref[2] = $mol2_atomdata[0][$atref][2] - $mol2_atomdata[0][$atcx][2];

# open($fhtest,">","test.xyz");
# print $fhtest $mol_numents[0][$tc60]*3,"\n\n";

# do the fitting
copy_config(0,1);

if($periodic_key[0]==0) {
  @remapaxes=();
} elsif($periodic_key[0]==6) {
  @remapaxes=(0,1);
} elsif($periodic_key[0]==1 or $periodic_key[0]==2) {
  @remapaxes=(0,1,2);
} else {
  &error(9);}

for($m=0;$m<$mol_numents[0][$tc60];$m++) {
  remap_molecule(\@{$cdata[0][$tc60][$m]},\@remapaxes,\@{$size[0]});
  $dumvful[0] = $cdata[0][$tc60][$m][$dumful][0] - $cdata[0][$tc60][$m][$dumcx][0];
  $dumvful[1] = $cdata[0][$tc60][$m][$dumful][1] - $cdata[0][$tc60][$m][$dumcx][1];
  $dumvful[2] = $cdata[0][$tc60][$m][$dumful][2] - $cdata[0][$tc60][$m][$dumcx][2];
  @axis=vector_product(\@dumvful, \@atvful);
  @axis=normalize_vector(\@axis);
  $angle1=acos(dot_product(\@dumvful,\@atvful)/(vector_length(@dumvful) * vector_length(@atvful)));
  @rm1=gen_rot_matrix(\@axis,-$angle1);
  
  $dumvref[0]= $cdata[0][$tc60][$m][$dumref][0] - $cdata[0][$tc60][$m][$dumcx][0];
  $dumvref[1]= $cdata[0][$tc60][$m][$dumref][1] - $cdata[0][$tc60][$m][$dumcx][1];
  $dumvref[2]= $cdata[0][$tc60][$m][$dumref][2] - $cdata[0][$tc60][$m][$dumcx][2];
  @dumvref=vector_projection_plane(\@dumvref,\@dumvful);
  @dumvref=normalize_vector(\@dumvref);
  
  @atvref2=rotate_vector(\@atvref,\@rm1);
  @atvref2=vector_projection_plane(\@atvref2,\@dumvful);
  @dumvful=normalize_vector(\@dumvful);
  @atvref2=normalize_vector(\@atvref2);
  @dumvref=normalize_vector(\@dumvref);
  @axis=vector_product(\@dumvref, \@atvref2);
  @axis=normalize_vector(\@axis);
#   $test1[0] = $cdata[0][$tc60][$m][$dumcx][0] + $atvref2[0];
#   $test1[1] = $cdata[0][$tc60][$m][$dumcx][1] + $atvref2[1];
#   $test1[2] = $cdata[0][$tc60][$m][$dumcx][2] + $atvref2[2];
#   $test2[0] = $cdata[0][$tc60][$m][$dumcx][0] + $dumvref[0];
#   $test2[1] = $cdata[0][$tc60][$m][$dumcx][1] + $dumvref[1];
#   $test2[2] = $cdata[0][$tc60][$m][$dumcx][2] + $dumvref[2];
#   print $fhtest "O  ",join(" ",@{$cdata[0][$tc60][$m][$dumcx]}[0..2]),"\n";
#   print $fhtest "H  ",join(" ",@test1),"\n";
#   print $fhtest "Li ",join(" ",@test2),"\n";
  $angle2=acos(dot_product(\@dumvref,\@atvref2)/(vector_length(@dumvref) * vector_length(@atvref2)));
#   print $angle2/pi*180,"\n";
  @rm2=gen_rot_matrix(\@axis,-$angle2);
  
  @{$cdata[1][$tc60][$m]} = mol2_to_cdata(0);
  foreach $a (@c60) {
    $pos[0]=$mol2_atomdata[0][$a][0]-$atpful[0];
    $pos[1]=$mol2_atomdata[0][$a][1]-$atpful[1];
    $pos[2]=$mol2_atomdata[0][$a][2]-$atpful[2];
    @pos=rotate_vector(\@pos,\@rm1);
    @pos=rotate_vector(\@pos,\@rm2);
    $cdata[1][$tc60][$m][$a][0]=$pos[0]+$cdata[0][$tc60][$m][$dumful][0];
    $cdata[1][$tc60][$m][$a][1]=$pos[1]+$cdata[0][$tc60][$m][$dumful][1];
    $cdata[1][$tc60][$m][$a][2]=$pos[2]+$cdata[0][$tc60][$m][$dumful][2];
  }
  foreach $a (@nonc60) {
    $cdata[1][$tc60][$m][$a][0]=$cdata[0][$tc60][$m][$match[$a]][0];
    $cdata[1][$tc60][$m][$a][1]=$cdata[0][$tc60][$m][$match[$a]][1];
    $cdata[1][$tc60][$m][$a][2]=$cdata[0][$tc60][$m][$match[$a]][2];
  }
}

$config_key[1]=0;
write_config_file($outcfg,1,$config_title[0],">",1);

sub error {
  if($_[0]==1) {
    # no message
  }elsif($_[0]==2) {
    print "**** error: molecule name $mol_name[0][$_[1]] does not fit $mol_name[0][$_[2]] !\n";
  }elsif($_[0]==3) {
    print "**** error: atom name $mol_atomdata[0][$_[1]][$_[2]][0] does not fit ",
    "$mol_atomdata[1][$_[1]][$_[2]][0] for molecule $mol_name[0][$_[1]]!\n";
  }elsif($_[0]==4) {
    print "**** error: did not find molecule named $mol2_name[0] (from mol2) in FIELD files!\n";
  }elsif($_[0]==5) {
    print "**** error: did not find dummy C60 or CX atom!\n";
  }elsif($_[0]==6) {
    print "**** error: did not find attachment point for fullerene in atomistic FIELD file!\n";
  }elsif($_[0]==7) {
    print "**** error: C60 has $_[1] atoms instead of 60!\n";
  }elsif($_[0]==8) {
    print "**** error: could not find a match for atom ",($_[1]+1)," ($mol_atomdata[1][$tc60][$_[1]][0])\n";
   }elsif($_[0]==8) {
    print "**** error: script only implemented for orthogonal cells!\n";
  } else {
    print "**** unknown error!\n";
    exit 99;
  }
  exit $_[0]
}