set terminal png enhanced
set output 'hist.png'
set xlabel 'E'
set ylabel 'P(E)' 
p 'cubic_energy-T019.dat.hist' u 1:2 w l t 'T=38.2500', 'cubic_energy-T021.dat.hist' u 1:2 w l t 'T=38.5000', 'cubic_energy-T023.dat.hist' u 1:2 w l t 'T=38.7500'
