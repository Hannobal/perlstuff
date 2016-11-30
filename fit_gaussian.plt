set print "-"

f0(x)=A/sqrt(2*s0**2*pi)*exp(-(x-m0)**2/(2*s0**2))
A=1

set fit logfile '/dev/null'
set fit results
# set fit quiet
if(lfitmagnitude==1) {
#   A=magguess
  fit [fitmin:fitmax] f0(x) filename via A
  fit [fitmin:fitmax] f0(x) filename via m0
  fit [fitmin:fitmax] f0(x) filename via s0
  fit [fitmin:fitmax] f0(x) filename via A,m0,s0
  print "after fit:             ", m0," +/- ", s0, "   magnitude:", A
} else {
  fit [fitmin:fitmax] f0(x) filename via m0
  fit [fitmin:fitmax] f0(x) filename via s0
  fit [fitmin:fitmax] f0(x) filename via m0,s0
  print "after fit:             ", m0," +/- ", s0
}

if(lshowplot==1) plot filename w l t filename, f0(x)