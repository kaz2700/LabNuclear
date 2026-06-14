set terminal pngcairo size 1800,1000 enhanced font 'sans,14'
set output 'plots/fondo_rege_full.png'
set title 'Espectro de fondo del detector REGe (8192 canales, t_l = 4.86 h)'
set xlabel 'Canal'
set ylabel 'Cuentas'
set grid lw 0.3
set key top right
set logscale y

set arrow 1 from 6179,graph 0.02 to 6179,graph 0.4 nohead dt 2 lw 2 lc rgb '#cc0000'
set label 1 '^{40}K: 1460.82 keV (canal 6179)' at 5500,0.03 tc rgb '#cc0000' font 'sans,12'

set arrow 2 from 2160,graph 0.02 to 2160,graph 0.4 nohead dt 3 lw 1 lc rgb '#006600'
set label 2 'Aniquilacion: 511 keV' at 1300,0.03 tc rgb '#006600' font 'sans,10'

set arrow 3 from 785.9,graph 0.02 to 785.9,graph 0.15 nohead dt 3 lw 1 lc rgb '#006600'
set label 3 '^{226}Ra: 186 keV' at 500,0.03 tc rgb '#006600' font 'sans,10'

plot 'data/fondo_rege.txt' using 1:2 with lines lw 1 lc rgb '#003366' title 'Espectro REGe'
