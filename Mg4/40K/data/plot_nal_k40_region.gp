set terminal pngcairo size 1400,700 enhanced font 'sans,14'
set output '/home/comoving/LabNuclear/Mg4/40K/plots/nal_k40_region.png'
set title 'Region del K-40 (canales 430-502)'
set xlabel 'Canal'
set ylabel 'Cuentas'
set grid lw 0.5
set key top left
set arrow from 468.0,graph 0 to 468.0,graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'
plot [430:502] 'data/fondo_nal.txt' using 1:2 with lines lw 1.5 lc rgb '#003366' title 'Datos', \
     'data/nal_fit_curve.dat' using 1:2 with lines lw 2 lc rgb '#d62728' title 'Ajuste'
