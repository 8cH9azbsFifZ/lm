#!/bin/zsh
NMAX=`ls Hr__pic* | sort -t "c" -n -k 2 | tail -n 1 | sed 's/Hr__pic//g;s/.dat//g' `
RMAX=6
QMAX=10

cat << eof | gnuplot
set term post enhanced color
set encoding iso_8859_1
set xlabel "r/{/Symbol s}_{12}"
set key bottom
set ytics nomirror
set xtics nomirror

set out "Cr.eps"
set ylabel "c(r)"
plot [:$RMAX] "Cr__pic$NMAX.dat" u 1:2 w l t "c_{11}(r)", \
   "" u 1:3 w l t "c_{12}(r)", \
   "" u 1:5 w l t "c_{22}(r)"

set out "Gr.eps"
set ylabel "g(r)"
plot [:$RMAX] "Hr__pic$NMAX.dat" u 1:(\$2+1) w l t "g_{11}(r)", \
   "" u 1:(\$3+1) w l t "g_{12}(r)", \
   "" u 1:(\$5+1) w l t "g_{22}(r)"

set out "Gammar.eps"   
set ylabel "{/Symbol g}(r)"
plot [:$RMAX] "Gammar__pic$NMAX.dat" u 1:2 w l t"{/Symbol g}_{11}", \
   "" u 1:3 w l t"{/Symbol g}_{12}", \
   "" u 1:5 w l t"{/Symbol g}_{22}"

set out "Cq.eps"
set ylabel "C(q)"
set xlabel "q{/Symbol s}"
plot [:$QMAX] "Cq__pic$NMAX.dat" u 1:2 w l t "C_{11}(q)", \
   "" u 1:3 w l t "C_{12}(q)", \
   "" u 1:5 w l t "C_{22}(q)"

set out "Hq.eps"
plot [:$QMAX] "Hq__pic$NMAX.dat" u 1:2 w l t "H_{11}(q)", \
   "" u 1:3 w l t "H_{12}(q)", \
   "" u 1:5 w l t "H_{22}(q)"

set out "Gammaq.eps"
plot [:$QMAX] "Cq__pic$NMAX.dat" u 1:2 w l t "C_{11}(q)", \
   "" u 1:3 w l t "C_{12}(q)", \
   "" u 1:5 w l t "C_{22}(q)"

eof

cat Cr.eps Gr.eps Gammar.eps Cq.eps Hq.eps Gammaq.eps | mpage -k -2 -ba4 > Correlators.eps
