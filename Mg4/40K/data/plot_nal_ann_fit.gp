set terminal pngcairo size 1200,700 enhanced font 'sans,13'
set output '/home/comoving/LabNuclear/Mg4/40K/plots/nal_ann_fit.png'
set title 'Pico de aniquilacion (511 keV) en NaI(Tl)'
set xlabel 'Canal'
set ylabel 'Cuentas'
set grid lw 0.5
set key top left
set arrow from 165.6,graph 0 to 165.6,graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'
plot 'data/nal_ann_fit_curve.dat' using 1:2 with lines lw 2 lc rgb '#d62728' title 'Ajuste (gauss+fondo)', \
     'data/nal_ann_fit_curve.dat' using 1:3 with lines lw 1.5 lc rgb '#2ca02c' title 'Fondo lineal', \
     'data/fondo_nal.txt' using 1:2 every ::139::194 with points pt 7 ps 1.2 lc rgb '#1f77b4' title 'Datos'
