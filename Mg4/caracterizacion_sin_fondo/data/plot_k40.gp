set terminal pngcairo size 1400,850 enhanced font 'sans,13'
set output '/home/comoving/LabNuclear/Mg4/caracterizacion_sin_fondo/plots/k40_hpge_fit.png'
set title 'HPGe Mg4: pico candidato a 40K (datos y ajuste Fortran)'
set xlabel 'Canal'
set ylabel 'Cuentas'
set grid lw 0.5
set key top left
plot 'data/k40_fit_curve.dat' using 1:2 with points pt 7 ps 0.8 lc rgb '#1f77b4' title 'Datos', \
     'data/k40_fit_curve.dat' using 1:3 with lines lw 2 lc rgb '#d62728' title 'Ajuste (gauss + fondo lineal)'
