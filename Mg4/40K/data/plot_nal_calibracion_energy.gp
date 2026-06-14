set terminal pngcairo size 2000,1200 enhanced font 'sans,13'
set output '/home/comoving/LabNuclear/Mg4/40K/plots/nal_calibracion_energy.png'
set title 'Verificacion calibracion NaI(Tl): lineas de fondo en energia'
set xlabel 'Energia (keV)'
set ylabel 'Cuentas'
set grid lw 0.5
set logscale y
E(ch) = -9.130820 + 3.141230 * ch
plot 'data/fondo_nal.txt' using (E($1)):2 with lines lw 1 lc rgb '#003366' title 'Fondo NaI(Tl)'
