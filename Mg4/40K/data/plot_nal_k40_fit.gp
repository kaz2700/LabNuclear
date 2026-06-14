set terminal pngcairo size 1400,850 enhanced font 'sans,13'
set output '/home/comoving/LabNuclear/Mg4/40K/plots/nal_k40_fit.png'
set title 'Fondo NaI(Tl): pico K-40 (ajuste gaussiano + fondo lineal)'
set xlabel 'Canal'
set ylabel 'Cuentas'
set grid lw 0.5
set key top left
set arrow from 468.0,graph 0 to 468.0,graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'
plot 'data/nal_fit_curve.dat' using 1:2 with lines lw 2 lc rgb '#d62728' title 'Ajuste (gauss+fondo)', \
     'data/nal_fit_curve.dat' using 1:3 with lines lw 1.5 lc rgb '#2ca02c' title 'Fondo lineal', \
     'data/fondo_nal.txt' using 1:2 every ::429::501 with points pt 7 ps 1.2 lc rgb '#1f77b4' title 'Datos'
