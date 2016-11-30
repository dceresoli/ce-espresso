 reset
 set term post eps
 set output "wannier_bands.eps"
 unset xtics
set yrange [ -0.643: 15.654]
 set style line 1 lt 1 lc rgb "black" lw 2
 set style line 2 lt 2 lc rgb "red" lw 2
 set style line 3 lt 1 lc rgb "green" lw 1
 set ylabel "Energy (eV)"
 plot \
 "original_bands.dat" title "LDA bands" with lines linestyle 1,\
 "wannier_bands.dat" title "Wannier bands" with lines linestyle 2,\
 12.807 title "Fermi energy" with lines linestyle 3
