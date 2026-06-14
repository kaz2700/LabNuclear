set terminal pngcairo size 1800,1000 enhanced font 'sans,14'
set output '/home/comoving/LabNuclear/Mg4/40K/plots/nal_full.png'
set title 'Espectro de fondo detector NaI(Tl) (512 canales, t_l = 6375 s)'
set xlabel 'Canal'
set ylabel 'Cuentas'
set grid lw 0.5
set logscale y
set arrow from 468.0,graph 0 to 468.0,graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'
set label 1 '^{40}K: ch 468.0' at 478.0,graph 0.85 tc rgb '#cc0000'
plot 'data/fondo_nal.txt' using 1:2 with lines lw 1 lc rgb '#003366' title 'Fondo NaI(Tl)'
