import math, os, glob, subprocess

def linear_ls(x, y):
    n = len(x); sx = sum(x); sy = sum(y); sxx = sum(xi*xi for xi in x); sxy = sum(xi*yi for xi, yi in zip(x, y))
    denom = n*sxx - sx*sx
    b = (n*sxy - sx*sy) / denom if abs(denom) > 1e-30 else 0
    a = (sy - b*sx) / n if n > 0 else 0
    return (a, b)

plots_dir = "plots"
os.makedirs(plots_dir, exist_ok=True)

for itx_path in sorted(glob.glob(os.path.join("datos", "60*.itx"))):
    base = os.path.splitext(os.path.basename(itx_path))[0]
    with open(itx_path, "r", errors="replace") as f:
        lines = f.readlines()
    in_ch0 = False; counts = []
    for line in lines:
        s = line.strip()
        if s.startswith("WAVES") and "MCAch0" in s: in_ch0 = True; continue
        if in_ch0 and s == "BEGIN": continue
        if in_ch0 and s == "END": break
        if in_ch0: counts.append(int(s))

    sub = counts[197:214]
    pk = 197 + sub.index(max(sub))

    best_fit = None
    for hw in [6, 7, 8, 9, 10, 12]:
        for bg_offset in [10, 12, 14, 16, 18, 20]:
            bg_left = list(range(max(0, pk-bg_offset), pk-bg_offset+6))
            bg_right = list(range(pk+8, min(len(counts), pk+8+bg_offset-2)))
            bg_x = bg_left + bg_right
            if len(bg_x) < 6: continue
            bg_y = [counts[i] for i in bg_x]
            a_bg, b_bg = linear_ls(bg_x, bg_y)

            flo = pk - hw; fhi = pk + hw
            x_net, y_net = [], []
            for i in range(flo, fhi+1):
                bg = a_bg + b_bg*i
                net = counts[i] - bg
                if net > 0: x_net.append(float(i)); y_net.append(net)
            if len(x_net) < 5: continue

            y_log = [math.log(v) for v in y_net]
            sig_log = [1.0/math.sqrt(max(v, 0.01)) for v in y_net]

            n = len(x_net)
            s0 = sum(1/(s*s) for s in sig_log); s1 = sum(x/(s*s) for x,s in zip(x_net, sig_log))
            s2 = sum(x*x/(s*s) for x,s in zip(x_net, sig_log)); s3 = sum(x*x*x/(s*s) for x,s in zip(x_net, sig_log))
            s4 = sum(x*x*x*x/(s*s) for x,s in zip(x_net, sig_log))
            t0 = sum(y/(s*s) for y,s in zip(y_log, sig_log)); t1 = sum(x*y/(s*s) for x,y,s in zip(x_net, y_log, sig_log))
            t2 = sum(x*x*y/(s*s) for x,y,s in zip(x_net, y_log, sig_log))

            det = s0*(s2*s4 - s3*s3) - s1*(s1*s4 - s2*s3) + s2*(s1*s3 - s2*s2)
            if abs(det) < 1e-40: continue

            a00 = (s2*s4 - s3*s3)/det; a01 = -(s1*s4 - s2*s3)/det; a02 = (s1*s3 - s2*s2)/det
            a11 = (s0*s4 - s2*s2)/det; a12 = -(s0*s3 - s1*s2)/det; a22 = (s0*s2 - s1*s1)/det
            a0 = a00*t0 + a01*t1 + a02*t2; a1 = a01*t0 + a11*t1 + a12*t2; a2 = a02*t0 + a12*t1 + a22*t2
            if a2 >= 0: continue

            sigma = math.sqrt(-1.0/(2.0*a2))
            x0 = a1 * sigma * sigma
            ampl = math.exp(a0 + x0*x0/(2.0*sigma*sigma))
            if sigma < 1 or sigma > 12 or abs(x0 - pk) > 6: continue

            sum_sq = 0; n_pts = 0
            for i in range(flo, fhi+1):
                bg = a_bg + b_bg*i
                gauss = ampl * math.exp(-(i-x0)**2/(2.0*sigma*sigma))
                ymod = bg + gauss; res = counts[i] - ymod; sum_sq += res*res; n_pts += 1
            rms_v = math.sqrt(sum_sq/n_pts) if n_pts > 0 else 0

            chi2 = sum(((yi-(a0+a1*xi+a2*xi*xi))/si)**2 for xi,yi,si in zip(x_net, y_log, sig_log))
            ndf = len(x_net) - 3

            score = abs(x0 - pk) * 3 + chi2/ndf + (rms_v/ampl)*10
            if best_fit is None or score < best_fit[0]:
                best_fit = (score, hw, bg_offset, x0, sigma, ampl, rms_v, chi2, ndf, a_bg, b_bg, pk, flo, fhi, list(x_net))

    if best_fit:
        sc, hw, bg_off, x0, sigma, ampl, rms_v, chi2, ndf, a_bg, b_bg, pk, flo, fhi, x_net = best_fit
        bg_points = set(range(max(0, pk-bg_off), pk-bg_off+6))
        bg_points |= set(range(pk+8, min(len(counts), pk+8+bg_off-2)))

        dat_path = os.path.join(plots_dir, "_tmp_" + base + ".dat")
        with open(dat_path, "w") as f:
            for i, c in enumerate(counts):
                used = 1 if i in x_net else 0
                is_bg = 1 if i in bg_points else 0
                in_win = 1 if flo <= i <= fhi else 0
                if in_win or is_bg:
                    f.write("%d %d %d %d\n" % (i, c, used, is_bg))

        f_lo = max(0, pk - 25)
        f_hi = min(len(counts)-1, pk + 25)
        fwhm = 2.3548 * sigma
        c2n = chi2 / ndf

        label = "x0 = %.2f\\nsigma = %.2f\\nFWHM = %.2f\\nAmpl = %.0f\\nchi2/ndf = %.1f\\nRMS = %.1f" % (x0, sigma, fwhm, ampl, c2n, rms_v)

        gp = """set terminal pngcairo size 1200,800 enhanced font 'Verdana,10'
set output '%s'
set xlabel 'Channel'
set ylabel 'Counts'
set grid
set title '%s - 1332 keV Gaussian fit'
set xrange [%d:%d]
set style line 1 lc rgb 'blue' pt 7 ps 0.8
set style line 2 lc rgb 'red' pt 5 ps 0.8
set style line 3 lc rgb 'green' pt 9 ps 1.2
bg(x) = %.4f + %.4f*x
g(x) = %.4f * exp(-(x-%.4f)**2 / (2*%.4f**2))
set label 1 at graph 0.03, graph 0.72 "%s" front boxed font 'Verdana,11'
plot '%s' using 1:(column(3)==1?column(2):1/0) with points ls 1 title 'Used in fit', \
     '%s' using 1:(column(3)==0&&column(4)==1?column(2):1/0) with points ls 3 title 'Background only', \
     '%s' using 1:(column(3)==0&&column(4)==0&&column(1)>=%d&&column(1)<=%d?column(2):1/0) with points ls 2 title 'Excluded (net<=0)', \
     bg(x) with lines lw 2 lt 2 title 'Background', \
     bg(x) + g(x) with lines lw 2 lc rgb 'red' title 'Gaussian fit'
""" % (plots_dir + "/" + base + "_gauss.png", base, f_lo, f_hi,
       a_bg, b_bg, ampl, x0, sigma, label,
       dat_path, dat_path, dat_path, flo, fhi)

        gp_path = os.path.join(plots_dir, "_g_" + base + ".gp")
        with open(gp_path, "w") as f:
            f.write(gp)
        subprocess.run(["gnuplot", gp_path], check=True)
        os.remove(gp_path); os.remove(dat_path)
        print("%s: x0=%.2f sigma=%.2f ampl=%.0f chi2/ndf=%.1f" % (base, x0, sigma, ampl, c2n))
    else:
        print("%s: NO FIT" % base)
