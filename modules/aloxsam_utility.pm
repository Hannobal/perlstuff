package aloxsam_utility;

use 5.014002;
use Math::Trig;
use strict;
use warnings;
use dlpoly_utility;
use hanno_utility;

our $doo=2.755333333;
our (@ohperiodicvec);
$ohperiodicvec[0][0]=$doo;
$ohperiodicvec[0][1]=0;
$ohperiodicvec[1][0]=$doo*cos(pi/3);
$ohperiodicvec[1][1]=$doo*sin(pi/3);

require Exporter;
our @ISA = qw(Exporter);
our %EXPORT_TAGS = ( 'all' => [ qw() ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw(calc_zsurf find_hydroxide_candidates %refcindex %refp5index
    find_pa_phosphorous find_longest_chain find_longest_alkyl_chain
    find_c60_atoms find_acidic_protons find_acidic_oxygens);

our $VERSION = '0.01';

our %refp5index = (
  # neutral PAs
  "C10-PA"     => 10,
  "C14-PA"     => 10,
  "C17-PA"     => 17,
  "C18-PA"     => 10,
  "C1-PA"      =>  1,
  "C60-C18-PA" => 10,
  #single deprotonated PAs
  "BTBT-C12-PA-d"    => 10,
  "C10-PA-d"         => 10,
  "C12-C60-C18-PA-d" => 10,
  "C14-PA-d"         => 10,
  "C17-PA-d"         => 17,
  "C18-PA-d"         => 10,
  "C1-PA-d"          =>  1,
  "C2-4T-C12-PA-d"   => 10,
  "C60-C18-PA-d"     => 10,
  "F15-C18-PA-d"     => 10,
  "PHDA-d"           => 10,
  "C1-gly2-C2-PA-d"  => 10,
  #double deprotonated PAs
  "C10-PA-2d" => 10
);

our %refcindex = (
  # neutral PAs
  "C10-PA"     => 0,
  "C14-PA"     => 38,
  "C17-PA"     => 15,
  "C18-PA"     => 41,
  "C1-PA"      => 0,
  "C60-C18-PA" => 52,
  #single deprotonated PAs
  "BTBT-C12-PA-d"    => 36,
  "C10-PA-d"         => 0,
  "C12-C60-C18-PA-d" => 52,
  "C14-PA-d"         => 38,
  "C17-PA-d"         => 15,
  "C18-PA-d"         => 42,
  "C1-PA-d"          => 0,
  "C2-4T-C12-PA-d"   => 36,
  "C60-C18-PA-d"     => 52,
  "F15-C18-PA-d"     => 52,
  "PHDA-d"           => 48,
  "C1-gly2-C2-PA-d"  =>  2,
  #double deprotonated PAs
  "C10-PA-2d" => 0
);

sub calc_zsurf {
  my($ci,$fi,$ltop,$t,$m,$a,$surftop,$surfbot,$toh);
  $ci   = $_[0];
  $fi   = $_[1];
  $ltop = $_[2]; # must be -1 for bottom, 1 for top or 0 for both
  $toh=-1;
  for($t=0;$t<$field_nummols[$fi];$t++) {
    next if(not $mol_name[$fi][$t] =~ /hydroxide/i);
    $toh=$t;
    $surftop = -9e20;
    $surfbot =  9e20;
    # determine boundaries for z-position
    for($m=0;$m<$mol_numents[$fi][$t];$m++) {
      $surftop=$cdata[$ci][$t][$m][0][2] if($surftop<$cdata[$ci][$t][$m][0][2]);
      $surfbot=$cdata[$ci][$t][$m][0][2] if($surfbot<$cdata[$ci][$t][$m][0][2]);
    }
    if($ltop>0) {
      return $surftop;
    } elsif($ltop<0) {
      return $surfbot;
    } else {
      return $surftop, $surfbot;
    }
  }
  print "**** error: did not find entry for hydroxide ions in FIELD file $fi!\n";
  return undef, undef;
}

sub find_hydroxide_candidates {
  # @candidates = find_hydroxide_candidates($ci,$fi,$ltop,[$mindist],[$tol],[@gridvec],[@gridcenter])
  my($ci,$fi,$ltop,$mindistsq,@candidates,$t,$m,$a,$i,$surftop, $surfbot,
     $toh,$zmintol,$zmaxtol,$tol,$x,$y,$dx,$dy,$dz,$distsq,$shiftx,$shifty,
     $a1,$a2,$b1,$b2,$c1,$c2,$ca,$cb,@pos);
  $ci   = $_[0];
  $fi   = $_[1];
  $ltop = $_[2]; # must be -1 for bottom, 1 for top or 0 for both
  @candidates=();
  if($#_>2) {
    # min distance to neighboring PAs (off if ==0);
    # or radial tolerance for grid positions
    $mindistsq = $_[3]*$_[3];
  } else {
    $mindistsq = 0;
  }
  if($#_>3) {
    $tol = abs($_[4]); # tolerance in angstroms
  } else {
    $tol = 0.5; 
  }
  if($#_>4) {
    if($#_<8) {
      print "**** error in call of find_hydroxide_candidates: too few arguments for grid vectors!\n";
      return @candidates;
    }
    $a1 = $_[5]*$ohperiodicvec[0][0]+$_[6]*$ohperiodicvec[1][0];
    $a2 = $_[5]*$ohperiodicvec[0][1]+$_[6]*$ohperiodicvec[1][1];
    $b1 = $_[7]*$ohperiodicvec[0][0]+$_[8]*$ohperiodicvec[1][0];
    $b2 = $_[7]*$ohperiodicvec[0][1]+$_[8]*$ohperiodicvec[1][1];
  }
  my @center=(0,0);
  if($#_>8) {
    if($#_<10) {
      print "**** error in call of find_hydroxide_candidates: too few arguments for grid center!\n";
      return @candidates;
    }
    $center[0] = $_[9];
    $center[1] = $_[10];
  }
  if(not defined($field_nummols[$fi])) {
    print "**** error: FIELD data $fi is not defined!\n";
    return @candidates;
  }
  if(not @{$cdata[$ci]}) {
    print "**** error: \@{cdata[$ci]} is empty!\n";
    return @candidates;
  }
  
  $surftop = -9e20;
  $surfbot =  9e20;
  $toh=-1;
  for($t=0;$t<$field_nummols[$fi];$t++) {
    next if(not $mol_name[$fi][$t] =~ /hydroxide/i);
    $toh=$t; last;
  }
  if($toh<0) {
    print "**** error: did not find entry for hydroxide ions in FIELD file $fi!\n";
    return @candidates;
  }
  # determine boundaries for z-position
  for($m=0;$m<$mol_numents[$fi][$toh];$m++) {
    $surftop=$cdata[$ci][$toh][$m][0][2] if($surftop<$cdata[$ci][$toh][$m][0][2]);
    $surfbot=$cdata[$ci][$toh][$m][0][2] if($surfbot>$cdata[$ci][$toh][$m][0][2]);
  }
  if($ltop>0) {
    $zmintol=$surftop-$tol;
    $zmaxtol=$surftop+$tol;
  } elsif($ltop<0) {
    $zmintol=$surfbot-$tol;
    $zmaxtol=$surfbot+$tol;
  } else {
    $zmintol=$surfbot-$tol;
    $zmaxtol=$surftop+$tol;
  }
  # build @candidates list
  for($m=0;$m<$mol_numents[$fi][$toh];$m++) {
    if($cdata[$ci][$toh][$m][0][2]>$zmintol and $cdata[$ci][$toh][$m][0][2]<$zmaxtol) {
      push(@candidates,$m);
    }
  }
  if(defined($a1)) {
    checkgrid : for($i=0;$i<@candidates;$i++) {
      $c1 = $cdata[$ci][$toh][$candidates[$i]][0][0]-$center[0];
      $c2 = $cdata[$ci][$toh][$candidates[$i]][0][1]-$center[1];
      if(abs($a1)>1e-14) {
	$cb = ($c2-$c1*$a2/$a1)/($b2-$b1*$a2/$a1);
	$ca = ($c1-$cb*$b1)/$a1;
      } else {
	$ca = ($c2-$c1*$b2/$b1)/($a2-$a1*$b2/$b1);
	$cb = ($c1-$ca*$a1)/$b1;
      }
      $ca = sprintf("%.0f",$ca);
      $cb = sprintf("%.0f",$cb);
      $dx = $c1-($ca*$a1+$cb*$b1);
      $dy = $c2-($ca*$a2+$cb*$b2);
      if($dx*$dx+$dy*$dy>$mindistsq) {
	splice(@candidates,$i,1);
	$i--;
      }
    }
  } elsif($mindistsq>0) { # exclude hydroxides neighboring to PAs
    candcheck:for($i=0;$i<=$#candidates;$i++) {
      for($x=-1;$x<2;$x++) {
        for($y=-1;$y<2;$y++) {
          $shiftx=$x*$cell[$ci][0][0]+$y*$cell[$ci][0][1];
          $shifty=$x*$cell[$ci][1][0]+$y*$cell[$ci][1][1];
          for($t=0;$t<$field_nummols[$ci];$t++) {
            next if(not $mol_name[$ci][$t] =~ /-PA-/);
            for($a=0;$a<$mol_numatoms[$ci][$t];$a++) {
              next if($mol_atomdata[$ci][$t][$a][0] ne 'O');
              for($m=0;$m<$mol_numents[$ci][$t];$m++) {
                $dx=$cdata[$ci][$toh][$candidates[$i]][0][0]-$cdata[$ci][$t][$m][$a][0]-$shiftx;
                $dy=$cdata[$ci][$toh][$candidates[$i]][0][1]-$cdata[$ci][$t][$m][$a][1]-$shifty;
                $dz=$cdata[$ci][$toh][$candidates[$i]][0][2]-$cdata[$ci][$t][$m][$a][2];
                $distsq = $dx*$dx + $dy*$dy + $dz*$dz;
                if($distsq<$mindistsq) {
                  splice(@candidates,$i,1);
                  $i--;
                  next candcheck;
                }
              }
            }
          }
        } # end for $y
      } # end for $x
    } # end candcheck
  } # end if $mindistsq>0
  return @candidates;
}

sub find_pa_phosphorous {
  # searches a PA-molecule for the phosphorous atom
  my $fi = $_[0];
  my $t  = $_[1];
  my($a,$p5);
  $p5=-1;
  for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
    if($mol_atomdata[$fi][$t][$a][0]=~/^P5$/i) {
      $p5=$a;
    }
  }
  return $p5;
}

sub find_longest_chain {
  # !!! does not handle ring systems well !!!
  # recursive function that searches a PA-molecule for the terminal
  # carbon atom in the alkyl chain
  # call examples: 
  # @chain = find_longest_chain($fi, $t, $p5); to obtain the entire chain
  
  my $fi    = $_[0];
  my $t     = $_[1];
  my @chain = @_[2..$#_];
  my($i,$cc,@longest,@other);
  @longest=@chain;
#   print "$chain[$#chain]    ",join(" ",@longest),"\n";
  foreach $i (@{$mol_bondatoms[$fi][$t][$chain[$#chain]]}) {
    next if(contains(@chain,$i));
    next if($mol_atomdata[$fi][$t][$i][0] !~ /^(C|CX|C3|OS)$/);
#     print "$chain[$#chain] doing $i\n";
    @other = find_longest_chain($fi,$t,@chain,$i);
    next if(@other<=@longest);
    undef @longest;
    @longest=@other;
  }
  return @longest;
}

sub find_longest_alkyl_chain {
  # !!! does not handle ring systems well !!!
  # recursive function that searches a PA-molecule for the terminal
  # carbon atom in the alkyl chain
  # call examples: 
  # @chain = find_longest_chain($fi, $t, $p5); to obtain the entire chain
  
  my $fi    = $_[0];
  my $t     = $_[1];
  my @chain = @_[2..$#_];
  my($i,$cc,@longest,@other);
  @longest=@chain;
#   print "$chain[$#chain]    ",join(" ",@longest),"\n";
  foreach $i (@{$mol_bondatoms[$fi][$t][$chain[$#chain]]}) {
    next if(contains(@chain,$i));
    next if($mol_atomdata[$fi][$t][$i][0] !~ /^(C|C3|OS)$/);
#     print "$chain[$#chain] doing $i\n";
    @other = find_longest_alkyl_chain($fi,$t,@chain,$i);
    next if(@other<=@longest);
    undef @longest;
    @longest=@other;
  }
  return @longest;
}

sub find_c60_atoms {
  # returns an array with atomids of the c60 atoms excluding the bing
  my $fi    = $_[0];
  my $t     = $_[1];
  my(@c60,$b,$a);
  loopatoms:for($a=0;$a<$mol_numatoms[$fi][$t];$a++) {
    if($mol_atomdata[$fi][$t][$a][0] =~/^CA$/i) {
      push(@c60,$a);
    } elsif($mol_atomdata[$fi][$t][$a][0] =~/^CX$/i) {
      foreach $b (@{$mol_bondatoms[$fi][$t][$a]}) {
	next loopatoms if($mol_atomdata[$fi][$t][$b][0] !~/^(CA|CX)$/i);
      }
      push(@c60,$a);
    }
  }
  return @c60;
}


sub find_acidic_protons {
  my $fi=$_[0];
  my $t=$_[1];
  my @lst=();
  for(my $a=0;$a<$mol_numatoms[$fi][$t];$a++) {
    if($mol_atomdata[$fi][$t][$a][$fi]=~/^P5$/i) {
      for my $b (@{$mol_bondatoms[$fi][$t][$a]}) {
	for my $c (@{$mol_bondatoms[$fi][$t][$b]}) {
	  if($mol_atomdata[$fi][$t][$c][0]=~/^(HO|HG)$/i) {
	    push(@lst,$c);
	  }
	}
      }
    }
  }
  return @lst;
}

sub find_acidic_oxygens {
  # the third argument is as follows:
  # 0/undef  return only doubly bound oxygen
  # 1        return doubly bound oxygen and hydroxyl oxygen
  # else     return only doubly bound oxygen
  my $fi=$_[0];
  my $t=$_[1];
  my $addoh=0;
  if($#_>1) {
    if($_[2]==0) {
      $addoh=0;
    } elsif($_[2]==1) {
      $addoh=1;
    } else {
      $addoh=2;
    }
  }
  my @lst=();
  for(my $a=0;$a<$mol_numatoms[$fi][$t];$a++) {
    if($mol_atomdata[$fi][$t][$a][$fi]=~/^P5$/i) {
      for my $b (@{$mol_bondatoms[$fi][$t][$a]}) {
        if($addoh<2 and $mol_atomdata[$fi][$t][$b][0]=~/^(O)$/i) {
          push(@lst,$b);
	} elsif($addoh>0 and $mol_atomdata[$fi][$t][$b][0]=~/^(OH)$/i) {
          push(@lst,$b);
	}
      }
    }
  }
  return @lst;
}

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

aloxsam_utility - Perl extension for blah blah blah

=head1 SYNOPSIS

  use aloxsam_utility;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for aloxsam_utility, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Hanno Dietrich, E<lt>dietrich@E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Hanno Dietrich

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.2 or,
at your option, any later version of Perl 5 you may have available.


=cut
