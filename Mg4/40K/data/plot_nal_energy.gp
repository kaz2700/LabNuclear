set terminal pngcairo size 1800,1000 enhanced font 'sans,14'
set output '/home/comoving/LabNuclear/Mg4/40K/plots/nal_energy.png'
set title 'Espectro NaI(Tl) calibrado en energia'
set xlabel 'Energia (keV)'
set ylabel 'Cuentas'
set grid lw 0.5
set logscale y
E(ch) = -9.130820 + 3.141230 * ch
set arrow from 1460.8,graph 0 to 1460.8,graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'
set label 1 '^{40}K: 1460.8 keV' at 1490.8,graph 0.85 tc rgb '#cc0000'
plot 'data/fondo_nal.txt' using (E($1)):2 with lines lw 1 lc rgb '#003366' title 'Fondo NaI(Tl)'
