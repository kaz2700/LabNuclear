set terminal pngcairo size 1400,850 enhanced font 'sans,13'
set output '/home/comoving/LabNuclear/Mg4/40K/plots/k40_rege_fit.png'
set title 'Fondo REGe: pico 40K (ajuste gaussiano + fondo lineal)'
set xlabel 'Canal'
set ylabel 'Cuentas'
set grid lw 0.5
set key top left
plot 'data/k40_fit_curve.dat' using 1:2 with lines lw 2 lc rgb '#d62728' title 'Ajuste (gauss + fondo)', \
     'data/k40_fit_curve.dat' using 1:3 with lines lw 1.5 lc rgb '#2ca02c' title 'Fondo lineal', \
'data/pico.dat' using ($1 + 6165):2 with points pt 7 ps 1.2 lc rgb '#1f77b4' title 'Datos'
