set terminal pngcairo size 1500,850 enhanced font 'sans,13'
set output '/home/comoving/LabNuclear/caracterizacion/plots/k40_nai_k40_fuera_rango.png'
set title 'NaI Mg4: imposibilidad de caracterizar 40K (1460.82 keV)'
set xlabel 'Energia (keV)'
set ylabel 'Cuentas'
set xrange [0:1600]
set yrange [1:*]
set grid lw 0.5
set key top right
set logscale y
set object 1 rect from 751.19, graph 0 to 1600, graph 1 fc rgb '#efefef' fs solid 1.0 noborder behind
set arrow 1 from 1460.82, graph 0 to 1460.82, graph 1 nohead dt 2 lw 2 lc rgb '#cc0000'
set label 1 '40K: 1460.82 keV' at 1210, graph 0.88 tc rgb '#cc0000'
set label 2 'Sin datos NaI medidos en esta adquisicion' at 860, graph 0.16 tc rgb '#555555'
plot '< tail -n +6 ../Mg4/datos/Mg4_NaI.txt' using 3:2 with lines lw 1.5 lc rgb '#7f0000' title 'Datos NaI disponibles'
