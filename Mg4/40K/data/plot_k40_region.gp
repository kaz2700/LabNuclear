set terminal pngcairo size 1400,700 enhanced font 'sans,14'
set output 'plots/fondo_rege_k40_region.png'
set title 'Region del ^{40}K (canales 6150-6210, sigma ~3-4 canales)'
set xlabel 'Canal'
set ylabel 'Cuentas'
set grid lw 0.3
set key top left

set arrow 1 from 6179,graph 0 to 6179,graph 1 nohead dt 2 lw 2 lc rgb '#cc0000'
set label 1 '^{40}K: 1460.82 keV' at 6150,65 tc rgb '#cc0000' font 'sans,11'

plot [6150:6210] 'data/fondo_rege.txt' using 1:2 with lines lw 1.5 lc rgb '#003366' title 'Datos'
