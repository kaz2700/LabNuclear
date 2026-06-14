set terminal pngcairo size 1600,900 enhanced font 'sans,13'
set output '/home/comoving/LabNuclear/Mg4/caracterizacion_sin_fondo/plots/cal_hpge_verificacion.png'
set title 'Verificacion en energia (calibracion 5 puntos con lineas de fondo)'
set xlabel 'Energia (keV)'
set ylabel 'Cuentas'
set grid lw 0.5
set key top right
set logscale y
set xrange [0:2000]
E(ch) = -17.357252 + 0.285675 * ch
set arrow 11 from 46.539,graph 0 to 46.539,graph 1 nohead dt 2 lw 1.5 lc rgb '#888888'
set arrow 12 from 53.225,graph 0 to 53.225,graph 1 nohead dt 2 lw 1.5 lc rgb '#006600'
set arrow 13 from 66.2,graph 0 to 66.2,graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'
set arrow 14 from 74.816,graph 0 to 74.816,graph 1 nohead dt 2 lw 1.5 lc rgb '#0066cc'
set arrow 15 from 77.105,graph 0 to 77.105,graph 1 nohead dt 2 lw 1.5 lc rgb '#cc6600'
set arrow 16 from 84.343,graph 0 to 84.343,graph 1 nohead dt 2 lw 1.5 lc rgb '#990099'
set arrow 17 from 114.673,graph 0 to 114.673,graph 1 nohead dt 2 lw 1.5 lc rgb '#669900'
set arrow 18 from 186.211,graph 0 to 186.211,graph 1 nohead dt 2 lw 1.5 lc rgb '#888888'
set arrow 19 from 238.632,graph 0 to 238.632,graph 1 nohead dt 2 lw 1.5 lc rgb '#006600'
set arrow 20 from 295.224,graph 0 to 295.224,graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'
set arrow 21 from 351.932,graph 0 to 351.932,graph 1 nohead dt 2 lw 1.5 lc rgb '#0066cc'
set arrow 22 from 609.312,graph 0 to 609.312,graph 1 nohead dt 2 lw 1.5 lc rgb '#cc6600'
set arrow 23 from 1460.82,graph 0 to 1460.82,graph 1 nohead dt 2 lw 1.5 lc rgb '#990099'
set arrow 24 from 1764.52,graph 0 to 1764.52,graph 1 nohead dt 2 lw 1.5 lc rgb '#669900'

set label 11 '46.5 keV (210Pb)' at 51.539,graph 0.85 tc rgb '#888888' font 'sans,9' left
set label 12 '53.2 keV (234U)' at 58.225,graph 0.8 tc rgb '#006600' font 'sans,9' left
set label 13 '66.2 keV (n)' at 71.2,graph 0.75 tc rgb '#cc0000' font 'sans,9' left
set label 14 '74.8 keV (212Pb)' at 79.816,graph 0.7 tc rgb '#0066cc' font 'sans,9' left
set label 15 '77.1 keV (214Bi)' at 82.105,graph 0.65 tc rgb '#cc6600' font 'sans,9' left
set label 16 '84.3 keV (231Th)' at 89.343,graph 0.6 tc rgb '#990099' font 'sans,9' left
set label 17 '114.7 keV (7Be)' at 119.673,graph 0.55 tc rgb '#669900' font 'sans,9' left
set label 18 '186.2 keV (226Ra)' at 191.211,graph 0.5 tc rgb '#888888' font 'sans,9' left
set label 19 '238.6 keV (212Pb)' at 243.632,graph 0.45 tc rgb '#006600' font 'sans,9' left
set label 20 '295.2 keV (214Pb)' at 300.224,graph 0.4 tc rgb '#cc0000' font 'sans,9' left
set label 21 '351.9 keV (214Pb)' at 356.932,graph 0.35 tc rgb '#0066cc' font 'sans,9' left
set label 22 '609.3 keV (214Bi)' at 614.312,graph 0.29999999999999993 tc rgb '#cc6600' font 'sans,9' left
set label 23 '1460.8 keV (40K)' at 1465.82,graph 0.25 tc rgb '#990099' font 'sans,9' left
set label 24 '1764.5 keV (214Bi)' at 1769.52,graph 0.19999999999999996 tc rgb '#669900' font 'sans,9' left

plot '/home/comoving/LabNuclear/Mg4/datos/Mg4_GeHP.txt' every ::4 using (E($1)):2 with lines lw 1 lc rgb '#003366' title 'Espectro HPGe Mg4'
