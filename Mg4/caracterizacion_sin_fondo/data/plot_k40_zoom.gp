set terminal pngcairo size 1400,850 enhanced font 'sans,13'
set output '/home/comoving/LabNuclear/Mg4/caracterizacion_sin_fondo/plots/k40_hpge_zoom.png'
set title 'Pico de 40K: canal vs energia (calib. con lineas del fondo)'
set xlabel 'Canal'
set ylabel 'Cuentas'
set label 1 'E(keV) = -7.1119 + 0.28289 * canal' at graph 0.03, graph 0.96 tc rgb '#333333' font 'sans,10'
set label 2 'Centro: ch 5189.09 = 1460.84 keV' at graph 0.03, graph 0.91 tc rgb '#cc0000' font 'sans,10'
set grid lw 0.5
set key top left
set xrange [5175:5205]
set xtics add ('{/=8 609.3 keV (Bi-214)\ncanal 2179}' 2179, '{/=8 1460.8 keV (40K)\ncanal 5189}' 5189, '{/=8 2614.5 keV (Tl-208)\ncanal 9205}' 9205)
plot 'data/k40_fit_curve.dat' using 1:2 with points pt 7 ps 1.0 lc rgb '#1f77b4' title 'Datos', \
     'data/k40_fit_curve.dat' using 1:3 with lines lw 2 lc rgb '#d62728' title 'Ajuste gaussiano + fondo lineal'
