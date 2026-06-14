set terminal pngcairo size 1600,900 enhanced font 'sans,13'
set output '/home/comoving/LabNuclear/Mg4/caracterizacion_sin_fondo/plots/k40_verificacion_fondo.png'
set title 'Verificacion de lineas de fondo GeHP BEGe en espectro Mg4'
set xlabel 'Canal'
set ylabel 'Cuentas'
set grid lw 0.5
set key top right
set logscale y
E(ch) = -8.870720 + 0.283520 * ch
set arrow from 5189.1765,graph 0 to 5189.1765,graph 1 nohead dt 2 lw 2 lc rgb '#006600'
set label '40K: 1462.37 keV' at 5309.176544,graph 0.80 tc rgb '#006600' font 'sans,12'
set arrow 1 from 5189,graph 0 to 5189,graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'
set arrow 2 from 6251,graph 0 to 6251,graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'
set label 1 '1460.8 keV (40K)' at 5269,graph 0.90-0.0 tc rgb '#cc0000' font 'sans,11'
set label 2 '1764.5 keV (214Bi)' at 6331,graph 0.90-0.05 tc rgb '#cc0000' font 'sans,11'
plot '< tail -n +5 /home/comoving/LabNuclear/Mg4/datos/Mg4_GeHP.txt' using 1:2 with lines lw 1 lc rgb '#003366' title 'Espectro Mg4 HPGe', \
     '-' using 1:2 with points pt 7 ps 1.5 lc rgb '#006600' title 'Pico 40K'
5189.1765 1717
e
