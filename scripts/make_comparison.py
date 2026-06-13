import math, os, glob, subprocess

def linear_ls(x, y):
    n = len(x); sx = sum(x); sy = sum(y); sxx = sum(xi*xi for xi in x); sxy = sum(xi*yi for xi, yi in zip(x, y))
    denom = n*sxx - sx*sx; b = (n*sxy - sx*sy) / denom if abs(denom) > 1e-30 else 0; a = (sy - b*sx) / n if n > 0 else 0
    return (a, b)

def fit_gaussian(x_net, y_log, sig_log):
    n = len(x_net)
    s0 = sum(1/(s*s) for s in sig_log); s1 = sum(x/(s*s) for x,s in zip(x_net, sig_log))
    s2 = sum(x*x/(s*s) for x,s in zip(x_net, sig_log)); s3 = sum(x*x*x/(s*s) for x,s in zip(x_net, sig_log))
    s4 = sum(x*x*x*x/(s*s) for x,s in zip(x_net, sig_log))
    t0 = sum(y/(s*s) for y,s in zip(y_log, sig_log)); t1 = sum(x*y/(s*s) for x,y,s in zip(x_net, y_log, sig_log))
    t2 = sum(x*x*y/(s*s) for x,y,s in zip(x_net, y_log, sig_log))
    det = s0*(s2*s4 - s3*s3) - s1*(s1*s4 - s2*s3) + s2*(s1*s3 - s2*s2)
    if abs(det) < 1e-40: return None
    a00 = (s2*s4 - s3*s3)/det; a01 = -(s1*s4 - s2*s3)/det; a02 = (s1*s3 - s2*s2)/det
    a11 = (s0*s4 - s2*s2)/det; a12 = -(s0*s3 - s1*s2)/det; a22 = (s0*s2 - s1*s1)/det
    a0 = a00*t0 + a01*t1 + a02*t2; a1 = a01*t0 + a11*t1 + a12*t2; a2 = a02*t0 + a12*t1 + a22*t2
    if a2 >= 0: return None
    sigma = math.sqrt(-1.0/(2.0*a2)); x0 = a1*sigma*sigma; ampl = math.exp(a0 + x0*x0/(2.0*sigma*sigma))
    chi2 = sum(((yi-(a0+a1*xi+a2*xi*xi))/si)**2 for xi,yi,si in zip(x_net, y_log, sig_log))
    return (x0, sigma, ampl, chi2, n-3)

# Load caracterizacion
itx_files = glob.glob("datos/*.itx")
matches = [f for f in itx_files if "caracterizacion" in f.lower().replace(",","").replace("ó","o").encode("ascii","replace").decode()]
with open(matches[0], "r", errors="replace") as f:
    lines = f.readlines()
in_ch0 = False; counts = []
for line in lines:
    s = line.strip()
    if s.startswith("WAVES") and "MCAch0" in s: in_ch0 = True; continue
    if in_ch0 and s == "BEGIN": continue
    if in_ch0 and s == "END": break
    if in_ch0: counts.append(int(s))

pk = 197 + counts[197:214].index(max(counts[197:214]))
bg_off = 14
bg_left = list(range(max(0, pk-bg_off), pk-bg_off+6))
bg_right = list(range(pk+7, min(len(counts), pk+7+bg_off-2)))
bg_x = bg_left + bg_right
bg_y = [counts[i] for i in bg_x]
a_bg, b_bg = linear_ls(bg_x, bg_y)

hw_values = [4, 6, 10, 15, 20]
fits = []
for hw in hw_values:
    flo = pk - hw; fhi = pk + hw
    x_net, y_net = [], []
    for i in range(flo, fhi+1):
        bg = a_bg + b_bg*i
        net = counts[i] - bg
        if net > 0: x_net.append(float(i)); y_net.append(net)
    if len(x_net) < 5: continue
    y_log = [math.log(v) for v in y_net]
    sig_log = [1.0/math.sqrt(max(v, 0.01)) for v in y_net]
    result = fit_gaussian(x_net, y_log, sig_log)
    if result and 1 <= result[1] <= 15 and abs(result[0]-pk) <= 10:
        fits.append((hw, result[0], result[1], result[2]))

# Write data file
with open("comparison_data.dat", "w") as f:
    for i in range(182, 232):
        f.write("%d %d\n" % (i, counts[i]))

# Build gnuplot script
l = []
l.append('set terminal pngcairo size 1400,900 enhanced font "Verdana,12"')
l.append('set output "comparison_gaussians.png"')
l.append('set xlabel "Channel"')
l.append('set ylabel "Counts"')
l.append('set grid')
l.append('set title "Gaussian fit for 60-Co-caracterizacion - 1332 keV peak (different window widths)"')
l.append('set xrange [182:232]')
l.append('set key outside right top')
l.append('')

bg_fn = "bg(x) = %.4f + %.4f*x" % (a_bg, b_bg)
l.append(bg_fn)
l.append('')

colors = ["red", "blue", "green", "orange", "purple"]
dashes = ["1", "2", "3", "4", "5"]

for idx, (hw, x0, sigma, ampl) in enumerate(fits):
    c = colors[idx]
    fn_name = "g_%d" % hw
    expr = "%s(x) = bg(x) + %.4f * exp(-(x-%.4f)**2 / (2*%.4f**2))" % (fn_name, ampl, x0, sigma)
    l.append(expr)

l.append('')
plot = '"comparison_data.dat" using 1:2 with points pt 7 ps 1.0 lc rgb "black" title "Data"'

for idx, (hw, x0, sigma, ampl) in enumerate(fits):
    c = colors[idx]
    fwhm = 2.3548*sigma
    fn_name = "g_%d" % hw
    label = "hw=%d (x0=%.1f, s=%.2f, FWHM=%.1f)" % (hw, x0, sigma, fwhm)
    plot += ', \\\n    g_%d(x) with lines lw 2 lc rgb "%s" title "%s"' % (hw, c, label)

l.append("plot " + plot)

with open("comparison_gaussians.gp", "w") as f:
    f.write("\n".join(l))

subprocess.run(["gnuplot", "comparison_gaussians.gp"], check=True)
os.remove("comparison_data.dat")
os.remove("comparison_gaussians.gp")
print("Saved comparison_gaussians.png")
