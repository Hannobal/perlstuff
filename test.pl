#!/usr/bin/perl
use feature state;
use hanno_utility;
use dlpoly_utility;
use lammps_utility;
use Math::Trig;
use Cwd;

read_field_file("FIELD",0);
print $mol_dihedraldata[0][0][0][8],"\n";