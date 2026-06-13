set terminal pngcairo size 1600,900 enhanced font "sans,12"

# --- Mg4 GeHP ---
set output "../plots/Mg4_GeHP_raw.png"
set multiplot layout 2,1 title "Mg4 - Germanio HPGe (Raw)"
set logscale y
set xlabel "Canal"
set ylabel "Cuentas"
set title "Espectro completo (16384 canales)"
plot "< tail -n +5 ../datos/Mg4_GeHP.txt" using 1:2 with lines lc "#000080" notitle
unset logscale y
set xlabel "Canal"
set ylabel "Cuentas"
set title "Zoom (canales 0-500)"
set xrange [0:500]
plot "< tail -n +5 ../datos/Mg4_GeHP.txt" using 1:2 with lines lc "#000080" notitle
set xrange [*:*]
unset multiplot

# --- Mg4 NaI ---
set output "../plots/Mg4_NaI_raw.png"
set multiplot layout 2,1 title "Mg4 - NaI(Tl) (Raw)"
set logscale y
set xlabel "Energia (keV)"
set ylabel "Cuentas"
set title "Espectro completo (0-750 keV)"
plot "< tail -n +6 ../datos/Mg4_NaI.txt" using 3:2 with lines lc "#8B0000" notitle
unset logscale y
set xlabel "Energia (keV)"
set ylabel "Cuentas"
set title "Zoom lineal (0-100 keV)"
set xrange [0:100]
plot "< tail -n +6 ../datos/Mg4_NaI.txt" using 3:2 with lines lc "#8B0000" notitle
set xrange [*:*]
unset multiplot
