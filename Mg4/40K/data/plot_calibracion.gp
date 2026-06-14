set terminal pngcairo size 2000,1200 enhanced font "sans,13"
set output "plots/calibracion_verificacion.png"
set title "Verificación de la calibración: líneas de fondo conocidas sobre el espectro"
set xlabel "Canal"
set ylabel "Cuentas"
set grid lw 0.3
set key top right
set logscale y
set xrange [0:8192]
set yrange [0.5:*]

set arrow 1 from 194.7,graph 0 to 194.7,graph 0.6 nohead dt 2 lw 1 lc rgb "#cc0000"
set label 1 "Pb-210 (46.5 keV)" at 194.7,graph 0.62 tc rgb "#cc0000" font "sans,9" rotate by 90
set arrow 2 from 223.0,graph 0 to 223.0,graph 0.6 nohead dt 2 lw 1 lc rgb "#cc0000"
set label 2 "U-234 (53.2 keV)" at 223.0,graph 0.62 tc rgb "#cc0000" font "sans,9" rotate by 90
set arrow 3 from 785.9,graph 0 to 785.9,graph 0.6 nohead dt 2 lw 1 lc rgb "#cc0000"
set label 3 "Ra-226 (186.2 keV)" at 785.9,graph 0.62 tc rgb "#cc0000" font "sans,9" rotate by 90
set arrow 4 from 1007.6,graph 0 to 1007.6,graph 0.6 nohead dt 2 lw 1 lc rgb "#cc0000"
set label 4 "Pb-212 (238.6 keV)" at 1007.6,graph 0.62 tc rgb "#cc0000" font "sans,9" rotate by 90
set arrow 5 from 1247.0,graph 0 to 1247.0,graph 0.6 nohead dt 2 lw 1 lc rgb "#cc0000"
set label 5 "Pb-214 (295.2 keV)" at 1247.0,graph 0.62 tc rgb "#cc0000" font "sans,9" rotate by 90
set arrow 6 from 1487.2,graph 0 to 1487.2,graph 0.6 nohead dt 2 lw 1 lc rgb "#cc0000"
set label 6 "Pb-214 (351.9 keV)" at 1487.2,graph 0.62 tc rgb "#cc0000" font "sans,9" rotate by 90
set arrow 7 from 2160.0,graph 0 to 2160.0,graph 0.6 nohead dt 2 lw 1 lc rgb "#cc0000"
set label 7 "e+e- (511.0 keV)" at 2160.0,graph 0.62 tc rgb "#cc0000" font "sans,9" rotate by 90
set arrow 8 from 2465.8,graph 0 to 2465.8,graph 0.6 nohead dt 2 lw 1 lc rgb "#cc0000"
set label 8 "Tl-208 (583.2 keV)" at 2465.8,graph 0.62 tc rgb "#cc0000" font "sans,9" rotate by 90
set arrow 9 from 2576.2,graph 0 to 2576.2,graph 0.6 nohead dt 2 lw 1 lc rgb "#cc0000"
set label 9 "Bi-214 (609.3 keV)" at 2576.2,graph 0.62 tc rgb "#cc0000" font "sans,9" rotate by 90
set arrow 10 from 3854.5,graph 0 to 3854.5,graph 0.6 nohead dt 2 lw 1 lc rgb "#cc0000"
set label 10 "Ac-228 (911.2 keV)" at 3854.5,graph 0.62 tc rgb "#cc0000" font "sans,9" rotate by 90
set arrow 11 from 4099.0,graph 0 to 4099.0,graph 0.6 nohead dt 2 lw 1 lc rgb "#cc0000"
set label 11 "Ac-228 (969.0 keV)" at 4099.0,graph 0.62 tc rgb "#cc0000" font "sans,9" rotate by 90
set arrow 12 from 4738.7,graph 0 to 4738.7,graph 0.6 nohead dt 2 lw 1 lc rgb "#cc0000"
set label 12 "Bi-214 (1120.3 keV)" at 4738.7,graph 0.62 tc rgb "#cc0000" font "sans,9" rotate by 90
set arrow 13 from 6179.0,graph 0 to 6179.0,graph 0.6 nohead dt 2 lw 1 lc rgb "#cc0000"
set label 13 "K-40 (1460.8 keV)" at 6179.0,graph 0.62 tc rgb "#cc0000" font "sans,9" rotate by 90
set arrow 14 from 7465.6,graph 0 to 7465.6,graph 0.6 nohead dt 2 lw 1 lc rgb "#cc0000"
set label 14 "Bi-214 (1764.5 keV)" at 7465.6,graph 0.62 tc rgb "#cc0000" font "sans,9" rotate by 90

plot "data/fondo_rege.txt" using 1:2 with lines lw 1 lc rgb "#003366" title "Espectro REGe"
