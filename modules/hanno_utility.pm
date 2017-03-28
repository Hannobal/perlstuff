package hanno_utility;

use 5.014002;
use strict;
use warnings;
use Math::Trig;
require Exporter;
our @ISA = qw(Exporter);
our %EXPORT_TAGS = ( 'all' => [ qw() ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw(@histogram @histminindex @histmaxindex @histnumtotal @histnum
    @histdata @histvals

    histogram_add_value histogram_add_one histogram_add_value2 histogram_add_one2 
    histogram_normalize_integral histogram_normalize_maxpeak histogram_normalize_individual
    histogram_clear check_integer check_real check_number contains
    sgn time2human find_directoryindex_timestep get_numeric_directories gen_rot_matrix 
    rotate_vector vector_product vector_length normalize_vector scale_vector
    dot_product histogram_normalize_value sum_array calc_average calc_stdev
    subarray min min_nonrecursive max vector_projection vector_projection_plane vector_add
    vector_subst parse_intlist randomgauss histogram_normalize_numdata
    histogram_estimate_stdev read_data_file matmul round generate_octahedron
    generate_prism generate_tetrahedron matrix_det matrix_invert matrix_multiply
    matrix_transpose);

our $VERSION = '0.01';

our(@histdata,@histvals,@histogram,@histminindex,@histmaxindex,@histnumtotal,@histnum);

sub histogram_add_one {
  my $hi  = $_[0]; # index of histogram
  my $val = $_[1]; # value for index
  my $res = $_[2]; # resolution
  my $i;
  $i=int($val/$res+0.5);
  $histogram[$hi][$i]++;
  $histnum[$hi][$i]++;
  $histmaxindex[$hi]=$i if(not defined($histmaxindex[$hi]) or $histmaxindex[$hi]<$i);
  $histminindex[$hi]=$i if(not defined($histminindex[$hi]) or $histminindex[$hi]>$i);
  $histnumtotal[$hi]++;
}

sub histogram_add_one2 {
  # same as above but saves original data into histdata array (!!! may consume huge amounts of RAM !!!)
  my $hi  = $_[0]; # index of histogram
  my $val = $_[1]; # value for index
  my $res = $_[2]; # resolution
  my $i;
  $i=int($val/$res+0.5);
  $histogram[$hi][$i]++;
  $histnum[$hi][$i]++;
  push(@{$histdata[$hi][$i]},$val);
  push(@{$histvals[$hi][$i]},1);
  $histmaxindex[$hi]=$i if(not defined($histmaxindex[$hi]) or $histmaxindex[$hi]<$i);
  $histminindex[$hi]=$i if(not defined($histminindex[$hi]) or $histminindex[$hi]>$i);
  $histnumtotal[$hi]++;
}

sub histogram_add_value {
  # add a value to the histogram
  my $hi   = $_[0]; # index of histogram
  my $val  = $_[1]; # value for index
  my $res  = $_[2]; # resolution
  my $val2 = $_[3]; # value to add
  my $i;
  $i=int($val/$res+0.5);
  $histogram[$hi][$i]+=$val2;
#   push(@{$histdata[$hi][$i]},$val);  # histogram value
#   push(@{$histvals[$hi][$i]},$val2); # added value
  $histmaxindex[$hi]=$i if(not defined($histmaxindex[$hi]) or $histmaxindex[$hi]<$i);
  $histminindex[$hi]=$i if(not defined($histminindex[$hi]) or $histminindex[$hi]>$i);
  $histnum[$hi][$i]++;
  $histnumtotal[$hi]++;
}

sub histogram_add_value2 {
  # same as above but saves original data into histdata array (!!! may consume huge amounts of RAM !!!)
  my $hi   = $_[0]; # index of histogram
  my $val  = $_[1]; # value for index
  my $res  = $_[2]; # resolution
  my $val2 = $_[3]; # value to add
  my $i;
  $i=int($val/$res+0.5);
  $histogram[$hi][$i]+=$val2;
  push(@{$histdata[$hi][$i]},$val);  # histogram value
  push(@{$histvals[$hi][$i]},$val2); # added value
  $histmaxindex[$hi]=$i if(not defined($histmaxindex[$hi]) or $histmaxindex[$hi]<$i);
  $histminindex[$hi]=$i if(not defined($histminindex[$hi]) or $histminindex[$hi]>$i);
  $histnum[$hi][$i]++;
  $histnumtotal[$hi]++;
}

sub histogram_normalize_value {
  my($hi,$i);
  my $val=splice(@_,-1,1);
  foreach $hi (@_) {
    next if(not defined $histmaxindex[$hi]);
    for($i=$histminindex[$hi];$i<=$histmaxindex[$hi];$i++) {
      next if(not defined($histogram[$hi][$i]));
      $histogram[$hi][$i]/=$val;
    }
  }
}

sub histogram_normalize_integral {
  # make the area under the plot 1
  my $hi  = $_[0];
  my $res = $_[1];
  my($i,$factor);
  return if(not defined $histmaxindex[$hi]);
  $factor = 0;
  for($i=$histminindex[$hi];$i<=$histmaxindex[$hi];$i++) {
    next if(not defined($histogram[$hi][$i]));
    $factor+=$histogram[$hi][$i]*$res;
  }
  for($i=$histminindex[$hi];$i<=$histmaxindex[$hi];$i++) {
    next if(not defined($histogram[$hi][$i]));
    $histogram[$hi][$i]/=$factor;
  }
}

sub histogram_normalize_numdata {
  # make the area under the plot 1
  my $hi  = $_[0];
  my($i);
  return if(not defined $histmaxindex[$hi]);
  for($i=$histminindex[$hi];$i<=$histmaxindex[$hi];$i++) {
    next if(not defined($histogram[$hi][$i]));
    $histogram[$hi][$i]/=$histnumtotal[$hi];
  }
}

sub histogram_normalize_maxpeak {
  my($hi,$i,$max);
  foreach $hi (@_) {
    next if(not defined $histmaxindex[$hi]);
    $max=-9e20;
    for($i=$histminindex[$hi];$i<=$histmaxindex[$hi];$i++) {
      next if(not defined($histogram[$hi][$i]));
      $max=$histogram[$hi][$i] if($histogram[$hi][$i]>$max);
    }
    for($i=$histminindex[$hi];$i<=$histmaxindex[$hi];$i++) {
      next if(not defined($histogram[$hi][$i]));
      $histogram[$hi][$i]/=$max;
    }
  }
}

sub histogram_normalize_individual {
  my $hi;
  foreach $hi (@_) {
    next if(not defined $histmaxindex[$hi]);
    for(my $i=$histminindex[$hi];$i<=$histmaxindex[$hi];$i++) {
      next if(not defined($histogram[$hi][$i]));
      next if(not defined($histnum[$hi][$i]));
      $histogram[$hi][$i]/=$histnum[$hi][$i];
    }
  }
}

sub histogram_clear {
  my $hi;
  foreach $hi (@_) {
    next if(not defined $histmaxindex[$hi]);
    undef @{$histogram[$hi]};
    undef @{$histdata[$hi]};
    undef @{$histvals[$hi]};
    undef @{$histnum[$hi]};
    delete $histnumtotal[$hi];
    delete $histmaxindex[$hi];
    delete $histminindex[$hi];
  }
}

sub histogram_estimate_stdev {
  # !!! only estimates for gaussian distribution, no fit !!!
  my $hi  = $_[0];
  my $res = $_[1];
  my ($i,$average,$stdev);
  return undef,undef if(not(defined $histnumtotal[$hi]));
  return undef,undef if($histnumtotal[$hi]==0);
  for($i=$histminindex[$hi];$i<=$histmaxindex[$hi];$i++) {
    next if(not defined($histnum[$hi][$i]));
    $average += $histnum[$hi][$i]*$i*$res;
  }
  $average /= $histnumtotal[$hi];
  $stdev = ($histmaxindex[$hi]-$histminindex[$hi])/8.0*$res; # factor 1/8 is just empirical
  return $average,$stdev;
}

sub sum_array {
  # sums all data in a (multidimensional) array
  my $sum     = 0;
  my $numvals = 0;
  my ($a ,$b);
  foreach my $val (@_) {
    next if(not defined($val));
    if(ref($val) eq "ARRAY") {
      ($a,$b) = sum_array(@{$val});
      $numvals += $a;
      $sum += $b;
    } else {
      $sum += $val;
      $numvals++;
    }
  }
  return $numvals,$sum;
}

sub calc_average {
  # calculates the average of a set of data
  my ($average, $numvals);
  ($numvals,$average) = sum_array(@_);
  return 0,0 if($numvals==0);
  return $numvals,$average/$numvals;
}

sub sum_array_stdev {
  # sums for standard deviation (calc_stdev)
  my $av      = splice(@_,$#_,1);
  my $sum     = 0;
  my $numvals = 0;
  my $a;
  foreach my $val (@_) {
    next if(not defined($val));
    if(ref($val) eq "ARRAY") {
      $a = sum_array_stdev(@{$val},$av);
      $sum += $a;
    } else {
      $a = $val-$av;
      $sum += $a*$a;
    }
  }
  return $sum;
}

sub calc_stdev {
  # calculates the standard deviation of a (multidimensional) array of data
  my($stdev, $av, $val, $j, $numvals, $a, $b);
  ($numvals,$av) = calc_average(@_);
  return 0,0 if($numvals==0);
  return $av,0 if($numvals==1);
  $stdev   = sum_array_stdev(@_,$av);
  return $av,sqrt((1/($numvals-1))*$stdev);
}

sub check_integer {
  my $check = 1;
  for(my $i=0;$i<@_;$i++) {
    return 0 if (not defined($_[$i]));
    return 0 if (not $_[$i] =~ /^[+-]?\d+(e[+]?\d+)?$/ );
  }
  return 1;
}

sub check_real {
  for(my $i=0;$i<@_;$i++) {
    return 0 if (not defined($_[$i]));
    if (not (($_[$i] =~ /^[+-]?\d+\.?\d*([eEdD][+-]?\d+)?$/ 
        or    $_[$i] =~ /^[+-]?\d*\.?\d+([eEdD][+-]?\d+)?$/ )
        and length($_[$i])>0)) {
      return 0;
    }
  }
  return 1;
}

sub check_number {
  for(my $i=0; $i<@_; $i++) {
    return 0 if (not defined($_[$i]));
    if($_[$i]*1 ne $_[$i]) {
      return 0;
    }
  }
  return 1;
}

sub contains {
  # (recursivly) check whether an array contains a certain element
  my $value = $_[$#_];
  for(my $i=0;$i<$#_;$i++) {
    return 0 if (not defined($_[$i]));
    if(ref($_[$i]) eq 'ARRAY'){
      return 1 if(contains(@{$_[$i]},$value));
    } elsif($_[$i] eq $value) {
      return 1;
    }
  }
  return 0;
} # end subroutine contains

sub sgn {
  return ($_[0] <=> 0);
}

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

sub find_directoryindex_timestep {
  my @files       = @{$_[0]};
  my $targetframe = $_[1];
  my($i);
  for($i=0;$i<@files;$i++) {
    if($files[$i]>=$targetframe) {
      return $i;
    }
  }
  return undef;
}

sub get_numeric_directories {
  my $inputdir   = $_[0];
  my $startframe = $_[1]; # optional
  my $endframe   = $_[2]; # optional
  my $prefix     = $_[3]; # optional
  my $suffix     = $_[4]; # optional
  my($i,$dir,@files,$startdirindex,$enddirindex,$tmp);
  if(not opendir($dir, $inputdir)) {
    print "**** error: could not open directory $inputdir\n";
    @files=();
    return @files;
  }
  @files = readdir($dir);
  closedir($dir);
# find the relevant folders
  for($i=0;$i<@files;$i++) {
    $tmp=$files[$i];
    $tmp=~s/^$prefix// if(defined($prefix));
    $tmp=~s/$suffix$// if(defined($suffix));
    if (not check_real($tmp)) {
      splice(@files,$i,1);
      $i--;
    }
  }
  return @files if(not @files);
  @files = sort {$a <=> $b} @files;
  if(defined($startframe)) {
    $startdirindex = find_directoryindex_timestep(\@files,$startframe);
    if(not defined($startdirindex)) {
      print "**** error: specified start frame is larger than last frame in simulation!\n";
      @files=();
      return @files;
    }
  } else {
    $startdirindex = 0;
  }
  if(defined($endframe)) {
    $enddirindex = find_directoryindex_timestep(\@files,$endframe);
    $enddirindex = $#files if(not defined($enddirindex));
  } else {
    $enddirindex = $#files;
  }
  return @files[$startdirindex..$enddirindex];
} # end subroutine get_numeric_directories

sub gen_rot_matrix {
  # args: axis (3 dim vector), angle (radian)
  my @axis  = @{$_[0]};
  my $angle =   $_[1];
  my $cosa  = cos($angle);
  my $sina  = sin($angle);
  my @rotmatrix;
  if($#axis == 2) {
  # 3 dimensional rotation matrix
  #1st index=line   2nd index=row
    $rotmatrix[0][0] = $cosa + ($axis[0]*$axis[0])*(1-$cosa);
    $rotmatrix[1][0] = $axis[1]*$axis[0]*(1-$cosa)+$axis[2]*$sina;
    $rotmatrix[2][0] = $axis[2]*$axis[0]*(1-$cosa)-$axis[1]*$sina;
    $rotmatrix[0][1] = $axis[0]*$axis[1]*(1-$cosa)-$axis[2]*$sina;
    $rotmatrix[1][1] = $cosa + ($axis[1]*$axis[1])*(1-$cosa);
    $rotmatrix[2][1] = $axis[2]*$axis[1]*(1-$cosa)+$axis[0]*$sina;
    $rotmatrix[0][2] = $axis[0]*$axis[2]*(1-$cosa)+$axis[1]*$sina;
    $rotmatrix[1][2] = $axis[1]*$axis[2]*(1-$cosa)-$axis[0]*$sina;
    $rotmatrix[2][2] = $cosa + ($axis[2]*$axis[2])*(1-$cosa);
  } elsif ($#axis == 1) {
    # 2 dimensional rotation matrix
    $rotmatrix[0][0] =   $cosa;
    $rotmatrix[0][1] = - $sina;
    $rotmatrix[1][0] =   $sina;
    $rotmatrix[1][1] =   $cosa;
  } else {
    print "error in subroutine gen_rot_matrix: axis vector has @axis dimensions instead of 2 or 3!\n";
    return undef;
  }
  return @rotmatrix;
} # end subroutine gen_rot_matrix


sub rotate_vector {
  # call example:
  # @vec = rotate_vector(\@vec, \@rotmatrix, \@center);
  my @vec       = @{$_[0]};
  my @rotmatrix = @{$_[1]};
  my @center    = @{$_[2]} if($#_>1);
  my (@new,@tmp);
  if(@center) {
    $tmp[0] = $vec[0]-$center[0];
    $tmp[1] = $vec[1]-$center[1];
    $tmp[2] = $vec[2]-$center[2];
    $new[0] = $tmp[0]*$rotmatrix[0][0]+$tmp[1]*$rotmatrix[0][1]+$tmp[2]*$rotmatrix[0][2]+$center[0];
    $new[1] = $tmp[0]*$rotmatrix[1][0]+$tmp[1]*$rotmatrix[1][1]+$tmp[2]*$rotmatrix[1][2]+$center[1];
    $new[2] = $tmp[0]*$rotmatrix[2][0]+$tmp[1]*$rotmatrix[2][1]+$tmp[2]*$rotmatrix[2][2]+$center[2];
  } else {
    $new[0] = $vec[0]*$rotmatrix[0][0]+$vec[1]*$rotmatrix[0][1]+$vec[2]*$rotmatrix[0][2];
    $new[1] = $vec[0]*$rotmatrix[1][0]+$vec[1]*$rotmatrix[1][1]+$vec[2]*$rotmatrix[1][2];
    $new[2] = $vec[0]*$rotmatrix[2][0]+$vec[1]*$rotmatrix[2][1]+$vec[2]*$rotmatrix[2][2];
  }
  return @new;
} # end subroutine rotate_vector

sub dot_product {
  my @vec1 = @{$_[0]};
  my @vec2 = @{$_[1]};
  my ($res,$i);
  $res=0;
  for($i=0;$i<@vec1;$i++) {
    $res += $vec1[$i]*$vec2[$i];
  }
  return $res;
}

sub vector_product {
  my @vec1 = @{$_[0]};
  my @vec2 = @{$_[1]};
  my @new;
  $new[0] = $vec1[1]*$vec2[2]-$vec1[2]*$vec2[1];
  $new[1] = $vec1[2]*$vec2[0]-$vec1[0]*$vec2[2];
  $new[2] = $vec1[0]*$vec2[1]-$vec1[1]*$vec2[0];
  return @new;
}

sub vector_length {
  my $i;
  my $res = 0;
  foreach my $val (@_) {
    $res += $val*$val;
  }
  return sqrt($res);
}

sub normalize_vector {
  my @vec = @{$_[0]};
  my $veclen=vector_length(@vec);
  my(@res,$i);
  for($i=0;$i<@vec;$i++) {
    $res[$i] = $vec[$i]/$veclen;
  }
  return @res;
}

sub vector_projection_plane {
  # project vec1 onto plane defined by normal vector n (and point p)
  my @v = @{$_[0]};
  my @n = @{$_[1]};
  my (@tmp,$fac);
  if(@_>2) {
    @tmp = vector_subst($_[0],$_[2]);
    $fac = dot_product(\@tmp,$_[1])/dot_product($_[1],$_[1])
  } else {
    $fac = dot_product($_[0],$_[1])/dot_product($_[1],$_[1])
  }
  $tmp[0] = $v[0]-$fac*$n[0];
  $tmp[1] = $v[1]-$fac*$n[1];
  $tmp[2] = $v[2]-$fac*$n[2];
  return @tmp;
}

sub vector_projection {
  # project vec1 onto vec2
  my @vec1 = @{$_[0]};
  my @vec2 = @{$_[1]};
  my $dp = dot_product(\@vec2,\@vec2);
  return undef if($dp==0);
  return scale_vector(\@vec2,dot_product(\@vec1,\@vec2)/$dp);
}

sub scale_vector {
  my @vec    = @{$_[0]};
  my $factor = $_[1];
  my(@res,$i);
  for($i=0;$i<@vec;$i++) {
    $res[$i] = $factor*$vec[$i];
  }
  return @res;
}

sub vector_add {
  my @vec1 = @{$_[0]};
  my @vec2 = @{$_[1]};
  my (@res,$i);
  for($i=0;$i<@vec1;$i++) {
    $res[$i] = $vec1[$i] + $vec2[$i];
  }
  return @res;
}

sub vector_subst {
  my @vec1 = @{$_[0]};
  my @vec2 = @{$_[1]};
  my(@res,$i);
  for($i=0;$i<@vec1;$i++) {
    $res[$i] = $vec1[$i] - $vec2[$i];
  }
  return @res;
}

sub subarray {
  my @arr = @{$_[0]};
  my $start = $_[1];
  my $end   = $_[2];
  my @res=();
  return @res if($start>$#arr);
  $end=$#arr if($end>$#arr);
  return @arr[$start..$end];
}

sub max {
  # find the maximum value in a multidimensional array
  # !!! values must all be numeric !!!
  my ($max, $next, @vars) = @_;
  return $max if not $next;
  return max( $max > $next ? $max : max($max,@{$next}) ) if(ref($next) eq "ARRAY");
  return max( $max > $next ? $max : $next, @vars );
}

sub min {
  # find the minimum value in a multidimensional array
  # !!! values must all be numeric !!!
  my ($min, $next, @vars) = @_;
  return $min if not defined($next);
  return min( $min > $next ? $min : min($min,@{$next}) ) if(ref($next) eq "ARRAY");
  return min( $min < $next ? $min : $next, @vars );
}

sub min_nonrecursive {
  # find the minimum value in a one-dimensional array
  # !!! values must all be numeric !!!
  my $min=$_[0];
  for(my $i=1;$i<@_;$i++) {
    $min = $_[$i] if($_[$i]<$min);
  }
  return $min;
}
sub parse_intlist {
  # parses a list of integers e.g. "1 5-7 13-22" and returns array
  # first argument is a shift
  my $shift = splice(@_,0,1);
  my($arg,@linedata,@arr,$i);
  foreach $arg (@_) {
    @linedata = split("-",$arg);
    if(not check_integer(@linedata)) {
      print "**** error: could not parse integers from expression $arg\n";
      return undef;
    }
    if(@linedata==1) {
      push(@arr,$linedata[0]+$shift);
    }elsif(@linedata==2) {
      @linedata = sort {$a <=> $b} @linedata;
      for($i=$linedata[0];$i<=$linedata[1];$i++) {
	push(@arr,$i+$shift);
      }
    }else{
      print "**** error: could not parse integers from expression $arg\n";
      return undef;
    }
  }
  return @arr;
}

sub randomgauss {
  # returns two random numbers in normal distribution standard deviation $stdev
  # the first one is positive, the second negative
  my $stdev = $_[0];
  my ($x1,$x2,$y1,$y2,$w);
  do {
    $x1 = rand();
    $x2 = rand();
    $w = $x1*$x1 + $x2*$x2;
  } while ($w >=1 || $w==0);
  $w = sqrt( (-2.0*log($w)) / $w);
  $y2 = $x1 *$w;
  $y1 = $x2 *$w;
  if(not defined($stdev)) {
    return ($y1,-$y2);
  } else {
    return ($y1*$stdev,-$y2*$stdev);
  }
}

sub read_data_file {
  my ($filename,$comment,$targetblock,$start,$end,$skiptext)=@_;
  $comment = '#'   if(not defined($comment));
  $start   = -9e20 if(not defined($start));
  $end     =  9e20 if(not defined($end));
  $targetblock = -1 if(not defined($targetblock));
  $skiptext    =  1 if(not defined($skiptext));
  my($filehandle,@data,@linedata,$block,$emptyline);
  
  if(not open($filehandle, "<", $filename)) {
    print "**** error: Can't open data file \"$filename\": $!\n";
    return ();
  }
  
  @data=();
  $block=0;
  $emptyline=0;
  readlines:while($_=<$filehandle>) {
    $_ =~ s/^\s+//;  $_ =~ s/\s+$//;
    if(length($_)==0) {
      $emptyline=1;
      next;
    } elsif($emptyline) {
      $emptyline=0;
      $block++;
    }
    next if($_=~/^$comment/);
    @linedata = split(/\s+/);
    if(check_real($linedata[0])) {
      next if($linedata[0]<$start);
      next if($linedata[0]>$end);
    } elsif($skiptext) {
      next;
    }
    push(@data,[@linedata]) if($targetblock<0 or $block==$targetblock);
  }
  close($filehandle);
  return @data;
}

sub matmul {
  my @mat1 = @{$_[0]};
  my @mat2 = @{$_[1]};
  my $l = @mat1;
  my $m = @{$mat1[0]};
  return undef if(@mat2 != $m);
  my $n = @{$mat2[0]};
  my(@res,$i,$j,$k);
  for($i=0;$i<$l;$i++) {
#     print "i=$i\n";
    for($j=0;$j<$n;$j++) {
#       print "  j=$j\n";
      $res[$i][$j]=0;
      for($k=0;$k<$m;$k++) {
	$res[$i][$j]+=$mat1[$i][$k]*$mat2[$k][$j];
      }
    }
  }
  return @res;
}

sub round {
  if($_[0]>=0) {
    return int($_[0]+0.5);
  } else {
    return int($_[0]-0.5);
  }
}

sub generate_prism {
  my($l,$h,$r,$n,$i,@points,$dang,$staggered,$sx,$sy);
  $n = $_[0]; # number of edges
  $l = $_[1]; # length of polygon edges
  $h = $_[2]; # height
  if($#_>2) { # staggered or eclipsed (0=eclipsed/1=staggered)?
    $staggered=$_[3];
  } else {
    $staggered=0;
  }
  if($#_>3) { # x-shear
    $sx=$_[4];
  } else {
    $sx=0;
  }
  if($#_>4) { # y-shear
    $sy=$_[5];
  } else {
    $sy=0;
  }
  $dang=2.0*pi/$n;
  $r=$l*0.5/sin(pi/$n);
  @points=();
  for($i=0;$i<$n;$i++) {
    push(@points,[$r*cos($i*$dang)-$sx*0.5,$r*sin($i*$dang)-$sy*0.5,$h*0.5]);
    if($staggered) {
      push(@points,[$r*cos(($i+0.5)*$dang)+$sx*0.5,$r*sin(($i+0.5)*$dang)+$sy*0.5,-$h*0.5]);
    } else {
      push(@points,[$r*cos($i*$dang)+$sx*0.5,$r*sin($i*$dang)+$sy*0.5,-$h*0.5]);
    }
  }
  return @points;
}

sub generate_tetrahedron {
  my($l,$r,$x,@points);
  $l = $_[0]; # length of square edges
  if($#_>0) {
    $r = $_[1]; # factor for height (1=no distortion)
  }else {
    $r=1;
  }
  $x=0.5*$l/sqrt(2)*$r;
  @points=();
  push(@points,[-$l*0.5,0,$x]);
  push(@points,[$l*0.5,0,$x]);
  push(@points,[0,-$l*0.5,-$x]);
  push(@points,[0,$l*0.5,-$x]);
  return @points;
}

sub generate_octahedron {
  my($l,$r,$x,@points);
  $l = $_[0]; # length of square edges
  if($#_>0) {
    $r = $_[1]; # factor for height (1=no distortion)
  }else {
    $r=1;
  }
  $x=$l/sqrt(2);
  @points=();
  push(@points,[$x,0,0]);
  push(@points,[-$x,0,0]);
  push(@points,[0,$x,0]);
  push(@points,[0,-$x,0]);
  push(@points,[0,0,$r*$x]);
  push(@points,[0,0,-$r*$x]);
  return @points;
}

sub matrix_multiply {
  my @m1=@{$_[0]};
  my @m2=@{$_[1]};
  return undef if($#{$m1[0]}!=$#m2);
  my @res=();
  if(ref($m2[0]) eq "ARRAY") {
    for(my $i=0;$i<@m1;$i++) {
      for(my $k=0;$k<@{$m2[0]};$k++) {
        $res[$i][$k]=0;
        for(my $j=0;$j<@m2;$j++) {
          $res[$i][$k] += $m1[$i][$j]*$m2[$j][$k];
        }
      }
    }
  } else { # we have a vector
    for(my $i=0;$i<@m1;$i++) {
      $res[$i]=0;
      for(my $j=0;$j<@m2;$j++) {
        $res[$i] += $m1[$i][$j]*$m2[$j];
      }
    }
  }
  return @res;
}

sub matrix_det {
  if($#_==1) { # assume 2x2-matrix
    return $_[0][0]*$_[1][1]-$_[0][1]*$_[1][0];
  } if($#_==2) { # assume 3x3-matrix
    return $_[0][0]*$_[1][1]*$_[2][2]
          +$_[0][1]*$_[1][2]*$_[2][0]
          +$_[0][2]*$_[1][0]*$_[2][1]
          -$_[2][0]*$_[1][1]*$_[0][2]
          -$_[2][1]*$_[1][2]*$_[0][0]
          -$_[2][2]*$_[1][0]*$_[0][1];
  } else {
    return undef;
  }
}

sub matrix_invert {
  my $f=matrix_det(@_);
  return undef if(not defined($f) or $f==0);
  $f=1.0/$f;
  my @inv=();
  if($#_==1) { # assume 2x2-matrix
    $inv[0][0] =  $f*$_[1][1];
    $inv[0][1] = -$f*$_[0][1];
    $inv[1][0] = -$f*$_[1][0];
    $inv[1][1] =  $f*$_[0][0];
  } elsif($#_==2) { # assume 3x3-matrix
    $inv[0][0] = -1*($_[1][1]*$_[2][2]-$_[1][2]*$_[2][1]);
    $inv[0][1] = -1*($_[0][2]*$_[2][1]-$_[0][1]*$_[2][2]);
    $inv[0][2] = -1*($_[0][1]*$_[1][2]-$_[0][2]*$_[1][1]);
    $inv[1][0] = -1*($_[1][2]*$_[2][0]-$_[1][0]*$_[2][2]);
    $inv[1][1] = -1*($_[0][0]*$_[2][2]-$_[0][2]*$_[2][0]);
    $inv[1][2] = -1*($_[0][2]*$_[1][0]-$_[0][0]*$_[1][2]);
    $inv[2][0] = -1*($_[1][0]*$_[2][1]-$_[1][1]*$_[2][0]);
    $inv[2][1] = -1*($_[0][1]*$_[2][0]-$_[0][0]*$_[2][1]);
    $inv[2][2] = -1*($_[0][0]*$_[1][1]-$_[0][1]*$_[1][0]);
  } else {
    return undef;
  }
  return @inv
}

sub matrix_transpose {
  my @res=();
  if(ref($_[0]) eq "ARRAY") {
    for(my $i=0;$i<@_;$i++) {
      for(my $j=0;$j<@{$_[$i]};$j++) {
        $res[$j][$i]=$_[$i][$j];
      }
    }
  } else {
    for(my $i=0;$i<@_;$i++) {
      $res[0][$i]=$_[$i];
    }
  }
  return @res;
}

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

hanno_utility - Perl extension for blah blah blah

=head1 SYNOPSIS

  use hanno_utility;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for hanno_utility, created by h2xs. It looks like the
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
