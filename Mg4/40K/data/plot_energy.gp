set terminal pngcairo size 1800,1000 enhanced font 'sans,14'
set output 'plots/fondo_rege_energy.png'
set title 'Espectro de fondo REGe calibrado en energia (E = 0.522 + 0.23633 {/Symbol \264} ch)'
set xlabel 'Energia (keV)'
set ylabel 'Cuentas'
set grid lw 0.3
set key top right
set logscale y

E(ch) = 0.5220 + 0.23633 * ch

set arrow 1 from 1460.82,graph 0.02 to 1460.82,graph 0.4 nohead dt 2 lw 2 lc rgb '#cc0000'
set label 1 '^{40}K: 1460.82 keV' at 1300,0.03 tc rgb '#cc0000' font 'sans,12'
set arrow 2 from 511.0,graph 0.02 to 511.0,graph 0.35 nohead dt 3 lw 1 lc rgb '#006600'
set label 2 'Aniquilacion: 511 keV' at 520,0.03 tc rgb '#006600' font 'sans,10'
set arrow 3 from 186.21,graph 0.02 to 186.21,graph 0.15 nohead dt 3 lw 1 lc rgb '#006600'
set arrow 4 from 609.31,graph 0.02 to 609.31,graph 0.25 nohead dt 3 lw 1 lc rgb '#006600'
set arrow 5 from 1764.52,graph 0.02 to 1764.52,graph 0.10 nohead dt 3 lw 1 lc rgb '#006600'

plot 'data/fondo_rege.txt' using (E($1)):2 with lines lw 1 lc rgb '#003366' title 'Espectro REGe'
