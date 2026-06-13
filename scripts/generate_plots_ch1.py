#!/usr/bin/env python3
"""
Generate all plots for MCAch1 (mobile detector, green, better resolution).
Creates raw spectrum plots and Gaussian fits for Co-60 spectra.
Output: plots/MCAch1/
"""
import os, sys, math, glob
import numpy as np

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots', 'MCAch1')
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'datos')
os.makedirs(OUT_DIR, exist_ok=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10, 'figure.dpi': 150})

def parse_itx_ch(filepath, ch='MCAch1'):
    meta = {}
    counts = []
    with open(filepath, 'r', errors='replace') as f:
        lines = f.readlines()
    in_ch = False
    for line in lines:
        s = line.strip()
        if s.startswith('WAVES') and ch in s:
            in_ch = True
            continue
        if in_ch:
            if s == 'BEGIN': continue
            if s == 'END': break
            try:
                counts.append(int(s))
            except ValueError:
                pass
    return meta, np.array(counts, dtype=float)

def linear_ls(x, y):
    n = len(x); sx = sum(x); sy = sum(y)
    sxx = sum(xi*xi for xi in x); sxy = sum(xi*yi for xi, yi in zip(x, y))
    denom = n*sxx - sx*sx
    b = (n*sxy - sx*sy) / denom if abs(denom) > 1e-30 else 0
    a = (sy - b*sx) / n if n > 0 else 0
    return (a, b)

def fit_gaussian_log(x_net, y_net):
    y_log = [math.log(v) for v in y_net]
    sig_log = [1.0/math.sqrt(max(v, 0.01)) for v in y_net]
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
    sigma = math.sqrt(-1.0/(2.0*a2))
    x0 = a1 * sigma * sigma
    ampl = math.exp(a0 + x0*x0/(2.0*sigma*sigma))
    if sigma < 1 or sigma > 12: return None
    chi2 = sum(((yi-(a0+a1*xi+a2*xi*xi))/si)**2 for xi,yi,si in zip(x_net, y_log, sig_log))
    ndf = len(x_net) - 3
    return (x0, sigma, ampl, chi2, ndf)

def make_raw_spectrum_plot(counts, base, label):
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.step(np.arange(len(counts)), counts, where='mid', color='black', linewidth=0.6)
    nonzero = np.where(counts > 0)[0]
    if len(nonzero) > 0:
        ax.set_xlim(max(0, nonzero[0] - 10), min(len(counts), nonzero[-1] + 10))
    ax.set_xlabel('Canal (raw ADC)', fontsize=11)
    ax.set_ylabel('Cuentas', fontsize=11)
    ax.set_title(f'{label} — MCAch1 (detector móvil)', fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, f'{base}_MCAch1.png'), dpi=150)
    plt.close()

def make_gaussian_fit_plot(counts, base, best_fit, label):
    if best_fit is None: return
    sc, hw, bg_off, x0, sigma, ampl, rms_v, chi2, ndf, a_bg, b_bg, pk, flo, fhi = best_fit
    fwhm = 2.3548 * sigma
    c2n = chi2 / ndf if ndf > 0 else 0

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.step(np.arange(len(counts)), counts, where='mid', color='black', linewidth=0.7, label='Datos')

    x_min = max(0, pk-25)
    x_max = min(len(counts)-1, pk+25)
    full_x = np.linspace(x_min, x_max, 400)
    bg_full = a_bg + b_bg * full_x
    g_y = ampl * np.exp(-(full_x - x0)**2 / (2 * sigma**2))
    ax.plot(full_x, bg_full + g_y, 'r-', lw=2, label='Ajuste: Gaussiana + fondo lineal')
    ax.plot(full_x, bg_full, 'g--', lw=1.5, label='Fondo lineal')

    x_plot = np.arange(flo, fhi+1)
    bg_vals = a_bg + b_bg * x_plot
    net_vals = counts[flo:fhi+1] - bg_vals
    used_mask = net_vals > 0

    ax.plot(x_plot[used_mask], counts[flo:fhi+1][used_mask], 'o', color='blue', ms=4, label='Usado en ajuste')
    ax.plot(x_plot[~used_mask], counts[flo:fhi+1][~used_mask], 'o', color='purple', ms=4, label='Excluido (net<=0)')

    txt = (f'$x_0 = {x0:.2f}$\n'
           f'$\\sigma = {sigma:.2f}$\n'
           f'FWHM = {fwhm:.2f}\n'
           f'Ampl = {ampl:.0f}\n'
           f'$\\chi^2/\\nu = {c2n:.1f}$')
    ax.text(0.03, 0.72, txt, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))

    ax.set_xlim(max(0, pk-25), min(len(counts)-1, pk+25))
    ax.set_xlabel('Canal (raw ADC)', fontsize=11)
    ax.set_ylabel('Cuentas', fontsize=11)
    ax.set_title(f'{label} — Pico 1332 keV (MCAch1, detector móvil)', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, f'{base}_gauss.png'), dpi=150)
    plt.close()

# ─── Main ──────────────────────────────────────────────────────────────────

itx_files = sorted(glob.glob(os.path.join(DATA_DIR, '60*.itx')))

for itx_path in itx_files:
    base = os.path.splitext(os.path.basename(itx_path))[0]
    print(f'Processing {base} (MCAch1)...')

    meta, counts = parse_itx_ch(itx_path, 'MCAch1')

    # Always generate raw spectrum plot
    make_raw_spectrum_plot(counts, base, base)

    sub = counts[197:214]
    pk = 197 + int(np.argmax(sub))

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

            result = fit_gaussian_log(x_net, y_net)
            if result is None: continue
            x0, sigma, ampl, chi2, ndf = result
            if x0 - pk > 6 or pk - x0 > 6: continue
            if sigma < 1 or sigma > 12: continue

            sum_sq = 0; n_pts = 0
            for i in range(flo, fhi+1):
                bg = a_bg + b_bg*i
                gauss = ampl * math.exp(-(i-x0)**2/(2.0*sigma*sigma))
                res = counts[i] - (bg + gauss); sum_sq += res*res; n_pts += 1
            rms_v = math.sqrt(sum_sq/n_pts) if n_pts > 0 else 0

            score = abs(x0 - pk) * 3 + chi2/ndf + (rms_v/ampl)*10
            if best_fit is None or score < best_fit[0]:
                best_fit = (score, hw, bg_offset, x0, sigma, ampl, rms_v, chi2, ndf, a_bg, b_bg, pk, flo, fhi)

    if best_fit:
        make_gaussian_fit_plot(counts, base, best_fit, base)
        x0 = best_fit[3]; sigma = best_fit[4]; ampl = best_fit[5]
        c2n = best_fit[7]/best_fit[8] if best_fit[8] > 0 else 0
        print(f'  -> Fit: x0={x0:.2f}, σ={sigma:.2f}, A={ampl:.0f}, χ²/ν={c2n:.1f}')
    else:
        print('  -> No fit found')

# Also process remaining non-60 files for raw spectra only
for extra in ['22Na-157-coin.itx', '22Na-180-coin.itx', 'F-F-coin.itx']:
    path = os.path.join(DATA_DIR, extra)
    if not os.path.exists(path): continue
    base = os.path.splitext(extra)[0]
    if os.path.exists(os.path.join(OUT_DIR, f'{base}_MCAch1.png')):
        continue
    print(f'Processing {base} (MCAch1)...')
    meta, counts = parse_itx_ch(path, 'MCAch1')
    make_raw_spectrum_plot(counts, base, base)

print(f'\nAll MCAch1 plots saved to: {OUT_DIR}')
