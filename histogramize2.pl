#!/usr/bin/perl

use hanno_utility;
use Switch;

$startframe = -9e20;
$endframe   = 9e20;

if($#ARGV<0) {
  print "the input format is:\n",
  " 1. input file\n",
  " 2. output directory\n",
  "-n         include numeric subdirectories\n",
  "-s <int>   start value in row 1\n",
  "-e <int>   end value in row 1\n",
  "-m set the point at min value of range instead of centered\n",
  "-o <real>  offset: shifts the values by this number\n",
  "-c <expr>  define a condition (e.g. \$2+$4<3)\n",
  "-<expr>    <resolution> <output filename>\n\n",
  "expr can be a row number starting from 1 (e.g. 2) or an expression with row e.g. '\$2+\$1/2'\n",
  "you may have to put the argument in quotes (e.g. '-(\$2+\$1)/2)'\n";
  exit 1;
}

$infilename = $ARGV[0];
$outdir     = $ARGV[1];
@indices    = ();
@rowflags   = ();
@histres    = ();
@outfiles   = ();
@conditions = ();
@condflags  = ();
$lminval    = 0;
$offset     = 0;
$incndirs   = 0;

for($i=2;$i<@ARGV;$i++) {
  switch($ARGV[$i]) {
    case (/^-n$/) {
      $incndirs=1;
    }case (/^-s$/) { # start frame
      $i++;
      $startframe = $ARGV[$i];
    } case (/^-e$/) { # end frame
      $i++;
      $endframe = $ARGV[$i];
    } case (/^-m$/) { # center point at start of bin instead of center
      $lminval = 1;
    } case (/^-o$/) { # offset
      $i++;
      $offset = $ARGV[$i];
    } case (/^-c$/) {
      $i++;
      $expr = $ARGV[$i];
      push(@condflags,$expr);
      $expr =~ s/\$(\d+)/\$linedata[\$$1]/g;
      $expr =~ s/\$(\d+)/($1-1)/eg;
      push(@conditions,$expr);
    } else {
      push(@rowflags,substr($ARGV[$i],1));
      if($ARGV[$i]=~/^-(\d+)$/) {
	push(@indices,"\$linedata[".($1-1)."]");
      } else {
	$expr = substr($ARGV[$i],1); # cut the "-" from the flag
	$expr =~ s/\$(\d+)/\$linedata[\$$1]/g;
	$expr =~ s/\$(\d+)/($1-1)/eg;
	push(@indices,$expr);
      }
      push(@minval,9e20);
      $i++;
      die "**** error: resolution must be a real number!\n" if(not check_real($ARGV[$i]));
      push(@histres,$ARGV[$i]);
      $i++;
      push(@outfiles,$ARGV[$i]);
    }
  }
}

print "parsed expressions for histograms:\n";
print join("\n",@indices),"\n";
if(@conditions) {
  print "parsed expressions for conditions:\n";
  print join("\n",@conditions),"\n";
}

if($incndirs) {
  @directories=get_numeric_directories(".",$startframe,$endframe);
  if(not @directories) {
    print "**** warning: no suitable directories found!\n";
  }
} else {
  @directories=(".");
}

$numframes=0;
foreach $dir (@directories) {
  if(not -e "$dir/$infilename") {
    print "**** warning: $dir/$infilename does not exist\n";
    next;
  }
  open(INPUT,"<","$dir/$infilename") or die "**** error: Can't open file $dir/$infilename!";

  $minframe=9e90;
  $maxframe=-9e20;
  readfile:while(<INPUT>) {
    $_ =~ s/^\s+//;
    $_ =~ s/\s+$//;
    next if(length($_)==0);
    next if(/^#/);
    @linedata = split(/\s+/,$_);
    next if(not check_real($linedata[0]));
    next if($linedata[0]<$startframe);
    last if($linedata[0]>$endframe);
    for($i=0;$i<@conditions;$i++) {
      next readfile unless(eval($conditions[$i]));
    }
    $minframe=min($minframe,$linedata[0]);
    $maxframe=max($maxframe,$linedata[0]);
    @{$data[$numframes]} = @linedata;
    for($i=0;$i<@indices;$i++) {
      $minval[$i] = min($minval[$i],eval($indices[$i]));
    }
    $numframes++;
  }
  close(INPUT);
}

for($i=0;$i<@indices;$i++) {
  $minval[$i] = $histres[$i]*(int($minval[$i]/$histres[$i])-1);
}

if(not -d $outdir) {
  mkdir $outdir or die "**** error: cannot create output-directory: $!";
}

if($startframe==-9e20) {$startframe=$minframe};
if($endframe  == 9e20) {$endframe  =$maxframe};

for($line=0;$line<@data;$line++) {
  @linedata = @{$data[$line]};
  for($i=0;$i<@indices;$i++) {
    histogram_add_one($i,eval($indices[$i])-$minval[$i]+$offset,$histres[$i]);
  }
}

for($i=0;$i<@indices;$i++) {
  histogram_normalize_integral($i,$histres[$i]);
  open($fhhist,">","$outdir/$outfiles[$i]") or die "**** error: cannot open output file $outfiles[$i]:\n$!";
  print $fhhist "# histogram from $numframes frames in between $startframe and $endframe with resolution $histres[$i]\n",
                "# of row $rowflags[$i] in file $infilename\n";
  printf $fhhist "# %18s %10s %10s\n", "value", "norm. cnt", "count";
  if(@conditions) {
    print $fhhist "# conditions:\n";
    print $fhhist "# ",join("\n# ",@condflags),"\n";
  }
  for($j=$histminindex[$i];$j<=$histmaxindex[$i];$j++) {
    if($lminval) {
      printf $fhhist "%20.8g %10.7f %10u\n", $histres[$i]*($j-0.5)+$minval[$i], $histogram[$i][$j], $histnum[$i][$j];
    } else {
      printf $fhhist "%20.8g %10.7f %10u\n", $histres[$i]*$j+$minval[$i], $histogram[$i][$j], $histnum[$i][$j];
    }
  }
}