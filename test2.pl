#!/usr/bin/perl
# use strict vars;
use dlpoly_utility;
use hanno_utility;
use Math::Trig;
# use Math::Matrix;
# use PDL;

# $a= Math::MatrixReal->new_from_rows([[1,0],[cos(pi/3.0),sin(pi/3.0)]]);
# $v= Math::MatrixReal->new_from_rows([[0],[2.0*sin(pi/3.0)]]);
# # 
# # $a = $a->transpose()->invert();
# # $b = $a->invert();
# # $w = $b*$v;
# print $v;


# $a= new Math::Matrix([1,0],[cos(pi/3.0),sin(pi/3.0)]);
# $v= new Math::Matrix([0],[2.0*sin(pi/3.0)]);
# 
# $a = $a->transpose();
# $b = $a->invert();
# $w = $b*$v;
# $w->print();

# $angle=pi/2.0;
# @axis=(0,0,1);
# @mat1=gen_rot_matrix(\@axis,$angle);
# $angle=pi/2.0;
# @axis=(1,0,0);
# @mat2=gen_rot_matrix(\@axis,$angle);
# 
# @vec=(1,0,0);
# @vec2=rotate_vector(\@vec,\@mat1);
# @vec2=rotate_vector(\@vec2,\@mat2);
# print join(" ",@vec2),"\n";
# 
# @mat3=matmul(\@mat2,\@mat1);
# @vec2=rotate_vector(\@vec,\@mat3);
# print join(" ",@vec2),"\n";

$mxblock=5;
@res=();
$ncurves=();
for($b=0;$b<$mxblock;$b++) {
  @{$data[$b]}=read_data_file($ARGV[0],"#",$b);
  for($l=0;$l<@{$data[$b]};$l++) {
    $res[$l]+=$data[$b][$l][1];
    $ncurves[$l]++;
  }
}
$av=$res[0]/$ncurves[0];
for($l=0;$l<@res;$l++) {
  print $data[0][$l][0]," ",$res[$l]," ",$res[$l]/$ncurves[$l]/$av,"\n";
}