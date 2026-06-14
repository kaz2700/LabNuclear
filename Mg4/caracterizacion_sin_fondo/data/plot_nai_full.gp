set terminal pngcairo size 1500,850 enhanced font 'sans,13'
set output '/home/comoving/LabNuclear/caracterizacion/plots/k40_nai_full.png'
set title 'NaI Mg4: espectro medido'
set xlabel 'Energia (keV)'
set ylabel 'Cuentas'
set grid lw 0.5
set key top right
set logscale y
plot '< tail -n +6 ../Mg4/datos/Mg4_NaI.txt' using 3:2 with lines lw 1.5 lc rgb '#7f0000' title 'NaI (0-751 keV)'
