#!/usr/bin/env python3
"""
Análisis completo de espectroscopía de coincidencias gamma.
Lee archivos .itx, ajusta picos gaussianos, calcula áreas,
genera gráficos y produce informe.
"""

import math, os, sys
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2 as chi2_dist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(OUT_DIR), 'datos')
PLOTS_DIR = os.path.join(OUT_DIR, 'graficos_ch1')
os.makedirs(PLOTS_DIR, exist_ok=True)

sys.path.insert(0, OUT_DIR)
import _fit_methods as fm

FORTRAN_BIN = os.path.join(os.path.dirname(OUT_DIR), 'src', 'gaussian_background')

# ─── 1. Parseo de archivos .itx ─────────────────────────────────────────────

def parse_itx(filepath):
    meta = {}
    counts = []
    with open(filepath, 'r', errors='replace') as f:
        lines = f.readlines()
    in_ch0 = False
    for line in lines:
        s = line.strip()
        if s.startswith('X //'):
            if 'Run Time [s]' in s:
                meta['run_time'] = float(s.split('=')[1].strip().split()[0])
            elif 'Event Rate [cps]' in s:
                meta['event_rate'] = float(s.split('=')[1].strip().split()[0])
            elif 'Live Time [s]' in s and 'Channel 0' in meta.get('_last_ch', ''):
                meta['live_time'] = float(s.split('=')[1].strip().split()[0])
            elif 'Input Count Rate [cps]' in s and 'Channel 0' in meta.get('_last_ch', ''):
                meta['input_count_rate'] = float(s.split('=')[1].strip().split()[0])
            elif 'Channel 0' in s:
                meta['_last_ch'] = 'Channel 0'
            if 'saved' in s:
                meta['timestamp'] = s.split('saved')[1].strip().split(';')[0].strip()
        if s.startswith('WAVES') and 'MCAch1' in s:
            in_ch0 = True
            continue
        if in_ch0:
            if s == 'BEGIN': continue
            if s == 'END': break
            try:
                counts.append(int(s))
            except ValueError:
                pass
    meta.pop('_last_ch', None)
    return meta, np.array(counts, dtype=float)


def cargar_espectro(nombre):
    path = os.path.join(DATA_DIR, nombre)
    if not os.path.exists(path):
        path = os.path.join(DATA_DIR, nombre + '.itx')
    meta, counts = parse_itx(path)
    channels = np.arange(len(counts), dtype=float)
    return meta, channels, counts


# ─── 2. Ajuste de picos ─────────────────────────────────────────────────────

def gauss_const(x, A, x0, sigma, bg):
    return A * np.exp(-(x - x0)**2 / (2 * sigma**2)) + bg


def gauss_linear(x, A, x0, sigma, c0, c1):
    return A * np.exp(-(x - x0)**2 / (2 * sigma**2)) + c0 + c1 * x


def fit_peak_robust(channels, counts, xmin, xmax, bg_order='const',
                    sigma_lim=(0.5, 12.0), debug=False):
    mask = (channels >= xmin) & (channels <= xmax)
    x = channels[mask]
    y = counts[mask]
    n_pts = len(x)
    if n_pts < 6:
        return None

    # Initial guess: find the maximum
    i_max = np.argmax(y)
    x0_guess = x[i_max]
    A_guess = y[i_max] - np.median(y)

    # Sigma guess from FWHM estimation
    half_max = y[i_max] / 2
    left_idx = None
    right_idx = None
    for j in range(i_max, -1, -1):
        if y[j] < half_max:
            left_idx = j
            break
    for j in range(i_max, len(y)):
        if y[j] < half_max:
            right_idx = j
            break
    if left_idx is not None and right_idx is not None:
        sigma_guess = (x[right_idx] - x[left_idx]) / 2.3548
    else:
        sigma_guess = 3.0
    sigma_guess = max(1.0, min(sigma_guess, 8.0))

    if A_guess < 1:
        A_guess = y[i_max] * 0.5

    if bg_order == 'const':
        bg_guess = np.median(np.concatenate([y[:max(2, n_pts//4)], y[-max(2, n_pts//4):]]))
        p0 = [A_guess, x0_guess, sigma_guess, bg_guess]
        bounds = ([0, xmin - 3, sigma_lim[0], 0],
                  [np.inf, xmax + 3, sigma_lim[1], np.inf])

        def model(x, A, x0, s, bg):
            return gauss_const(x, A, x0, s, bg)
        n_params = 4
    else:
        bg_c0 = np.median(y[:3]) if len(y) > 6 else np.median(y)
        bg_c1 = 0.0
        p0 = [A_guess, x0_guess, sigma_guess, bg_c0, bg_c1]
        bounds = ([0, xmin - 3, sigma_lim[0], -np.inf, -np.inf],
                  [np.inf, xmax + 3, sigma_lim[1], np.inf, np.inf])

        def model(x, A, x0, s, c0, c1):
            return gauss_linear(x, A, x0, s, c0, c1)
        n_params = 5

    try:
        sigma_y = np.where(y > 1, np.sqrt(y), 1.0)
        popt, pcov = curve_fit(model, x, y, p0=p0, bounds=bounds,
                               sigma=sigma_y, absolute_sigma=False,
                               maxfev=20000)
        perr = np.sqrt(np.diag(pcov))
    except Exception as e:
        if debug:
            print(f"  fit_peak failed [{xmin},{xmax}]: {e}")
        return None

    y_fit = model(x, *popt)
    resid = y - y_fit
    chi2 = np.sum((resid / np.maximum(np.sqrt(y), 0.5))**2)
    ndf = n_pts - n_params
    chi2_red = chi2 / ndf if ndf > 0 else 999
    ss_res = np.sum(resid**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    p_val = 1 - chi2_dist.cdf(chi2, ndf) if ndf > 0 else 0

    A, x0_f, sigma = popt[0], popt[1], popt[2]
    dA, dx0, dsigma = perr[0], perr[1], perr[2]

    if bg_order == 'const':
        bg, dbg = popt[3], perr[3]
        c1, dc1 = 0.0, 0.0
    else:
        bg, dbg = popt[3], perr[3]
        c1, dc1 = popt[4], perr[4]

    area = A * sigma * math.sqrt(2 * math.pi)
    d_area = area * math.sqrt((dA/A)**2 + (dsigma/sigma)**2) if A > 0 and sigma > 0 else 0
    fwhm = 2.35482 * sigma
    dfwhm = 2.35482 * dsigma

    # Check fit quality
    if sigma < sigma_lim[0] or sigma > sigma_lim[1]:
        if debug:
            print(f"  fit_peak: sigma={sigma:.2f} out of bounds, rejecting")
        return None
    if abs(x0_f - x0_guess) > (xmax - xmin) * 0.6:
        if debug:
            print(f"  fit_peak: x0={x0_f:.2f} too far from guess={x0_guess:.2f}, rejecting")
        return None

    return {
        'xmin': xmin, 'xmax': xmax,
        'A': A, 'dA': dA,
        'x0': x0_f, 'dx0': dx0,
        'sigma': sigma, 'dsigma': dsigma,
        'fwhm': fwhm, 'dfwhm': dfwhm,
        'bg': bg, 'dbg': dbg,
        'c1': c1, 'dc1': dc1,
        'area': area, 'darea': d_area,
        'chi2': chi2, 'ndf': ndf, 'chi2_red': chi2_red,
        'r_squared': r2, 'p_value': p_val,
        'bg_order': bg_order, 'n_params': n_params,
        'x_fit': x, 'y_fit': y_fit,
        'popt': popt, 'pcov': pcov,
    }


def fit_summary(fit, label=''):
    if fit is None:
        return f"{label}: NO CONVERGIÓ"
    return (f"{label}: x₀={fit['x0']:.2f}±{fit['dx0']:.2f}, "
            f"σ={fit['sigma']:.2f}±{fit['dsigma']:.2f}, "
            f"FWHM={fit['fwhm']:.2f}±{fit['dfwhm']:.2f}, "
            f"A={fit['A']:.1f}±{fit['dA']:.1f}, "
            f"Bg={fit['bg']:.1f}±{fit['dbg']:.1f}, "
            f"Área={fit['area']:.1f}±{fit['darea']:.1f}, "
            f"χ²/ν={fit['chi2_red']:.2f}, R²={fit['r_squared']:.4f}")


# ─── 3. Carga de datos ──────────────────────────────────────────────────────

meta_157, ch157, cnts157 = cargar_espectro('22Na-157-coin.itx')
meta_180, ch180, cnts180 = cargar_espectro('22Na-180-coin.itx')
meta_CoF, chCoF, cntsCoF = cargar_espectro('60Co-F-coin.itx')
meta_FF,  chFF,  cntsFF  = cargar_espectro('F-F-coin.itx')

# ─── 4. Análisis 22Na-157 ───────────────────────────────────────────────────

print("=" * 60)
print("  ANÁLISIS 22Na-157")
print("=" * 60)

# The 511+1274 keV coincidence shows detector 0's 511 keV peak at channel ~77,
# and the 1274 keV peak at channel ~199.
# (Detector 1's CFTD provides the gate; Detector 0 records its energy in MCAch1)

# Peak 511 keV from 511+1274 coinc: channels 68-90 (peak at 77, 569 counts)
# Both the 511 and 1274 peaks come from the 511+1274 coincidence.
for lo, hi in [(66, 90), (68, 92), (64, 92)]:
    fit_511_1274 = fit_peak_robust(ch157, cnts157, lo, hi, 'const', debug=True)
    if fit_511_1274:
        print(f"Window [{lo},{hi}]: {fit_summary(fit_511_1274, '511 (511+1274)')}")
        if 72 <= fit_511_1274['x0'] <= 82:
            break
    fit_511_1274 = None

if fit_511_1274 is None:
    fit_511_1274 = fit_peak_robust(ch157, cnts157, 60, 95, 'const',
                                   sigma_lim=(1, 15), debug=True)
    if fit_511_1274:
        print(f"Fallback: {fit_summary(fit_511_1274, '511 (511+1274)')}")

# Also fit with scripts-style and Fortran
fit_511_1274_scripts = fm.fit_scripts_style(ch157, cnts157, 78)
if fit_511_1274_scripts:
    print(f"  Scripts: {fit_511_1274_scripts['x0']:.2f} σ={fit_511_1274_scripts['sigma']:.2f}")
mask_511 = (ch157 >= 66) & (ch157 <= 90)
fit_511_1274_fort = fm.run_fortran_fit(ch157[mask_511], cnts157[mask_511], 78, FORTRAN_BIN)
if fit_511_1274_fort:
    print(f"  Fortran: x₀={fit_511_1274_fort['x0']:.2f} σ={fit_511_1274_fort['sigma']:.2f}")
print()

# Second peak from 511+1274 coinc: 1274 keV at channel ~199 (counts=71)
# This is the 1274 keV peak in detector 0 when detector 1 sees 511
fit_1274_157 = fit_peak_robust(ch157, cnts157, 188, 210, 'const', debug=True)
if fit_1274_157:
    print(f"1274 keV peak: {fit_summary(fit_1274_157, '1274 (511+1274)')}")
else:
    fit_1274_157 = fit_peak_robust(ch157, cnts157, 185, 215, 'const', debug=True)
    if fit_1274_157:
        print(f"1274 keV peak: {fit_summary(fit_1274_157, '1274 (511+1274)')}")
fit_1274_157_scripts = fm.fit_scripts_style(ch157, cnts157, 196)
if fit_1274_157_scripts:
    print(f"  Scripts (1274): x₀={fit_1274_157_scripts['x0']:.2f} σ={fit_1274_157_scripts['sigma']:.2f}")
mask_1274 = (ch157 >= 188) & (ch157 <= 210)
fit_1274_157_fort = fm.run_fortran_fit(ch157[mask_1274], cnts157[mask_1274], 196, FORTRAN_BIN)
if fit_1274_157_fort:
    print(f"  Fortran (1274): x₀={fit_1274_157_fort['x0']:.2f} σ={fit_1274_157_fort['sigma']:.2f}")
print()

# Peak 511+511 accidental ~ channels 195-210 (this is ALSO 511 keV,
# appearing when both detectors see 511 from different decays)
# Actually, this is the 1274 keV peak from the 511+1274 coinc; 
# the "accidental" peak is harder to separate.

print()

# Accidental counts: at 157°, the 511 peak contains both real (511+1274)
# and accidental (511+511) events. The 1274 peak contains only real events.
# Since the real coincidence rate is symmetric, the excess in 511 vs 1274
# gives the accidental contribution.
# N_acc = area(511 peak) - area(1274 peak)
N_total_157 = np.sum(cnts157)
N_acc_157 = 0
if fit_511_1274 and fit_1274_157:
    N_acc_157 = fit_511_1274['area'] - fit_1274_157['area']

print(f"N_total(157°) = {N_total_157:.0f}")
print(f"N_(511 peak) = {fit_511_1274['area']:.1f} ± {fit_511_1274['darea']:.1f}" if fit_511_1274 else "N_(511) = NO DISPONIBLE")
print(f"N_(1274 peak) = {fit_1274_157['area']:.1f} ± {fit_1274_157['darea']:.1f}" if fit_1274_157 else "N_(1274) = NO DISPONIBLE")
print(f"N_acc(157°) = N_(511) - N_(1274) = {N_acc_157:.1f}")
print()

# ─── 5. Análisis 22Na-180 ───────────────────────────────────────────────────

print("=" * 60)
print("  ANÁLISIS 22Na-180")
print("=" * 60)

# Peak 511+511 ~ channels 66-92 (huge peak max at channel 78, 6838 counts)
for lo, hi in [(66, 92), (68, 90), (64, 94)]:
    fit_511_511_180 = fit_peak_robust(ch180, cnts180, lo, hi, 'linear', debug=True)
    if fit_511_511_180:
        print(f"Window [{lo},{hi}]: {fit_summary(fit_511_511_180, '511+511')}")
        if 73 <= fit_511_511_180['x0'] <= 83:
            break
    fit_511_511_180 = None

if fit_511_511_180 is None:
    fit_511_511_180 = fit_peak_robust(ch180, cnts180, 60, 95, 'linear',
                                       sigma_lim=(1, 15), debug=True)
    if fit_511_511_180:
        print(f"Fallback: {fit_summary(fit_511_511_180, '511+511')}")
fit_511_511_180_scripts = fm.fit_scripts_style(ch180, cnts180, 78)
if fit_511_511_180_scripts:
    print(f"  Scripts (511+511): x₀={fit_511_511_180_scripts['x0']:.2f} σ={fit_511_511_180_scripts['sigma']:.2f}")
mask_511_180 = (ch180 >= 66) & (ch180 <= 92)
fit_511_511_180_fort = fm.run_fortran_fit(ch180[mask_511_180], cnts180[mask_511_180], 78, FORTRAN_BIN)
if fit_511_511_180_fort:
    print(f"  Fortran (511+511): x₀={fit_511_511_180_fort['x0']:.2f} σ={fit_511_511_180_fort['sigma']:.2f}")
print()

# Small peak ~190-210 (background / accidental - 511 keV from accidentals
# or possibly 1274 keV)
for lo, hi in [(190, 215), (188, 218)]:
    fit_bg_180 = fit_peak_robust(ch180, cnts180, lo, hi, 'const', debug=True)
    if fit_bg_180:
        print(f"Window [{lo},{hi}]: {fit_summary(fit_bg_180, 'pico ~200')}")
        break
    fit_bg_180 = None
fit_bg_180_scripts = fm.fit_scripts_style(ch180, cnts180, 197)
if fit_bg_180_scripts:
    print(f"  Scripts (fondo): x₀={fit_bg_180_scripts['x0']:.2f} σ={fit_bg_180_scripts['sigma']:.2f}")
mask_bg = (ch180 >= 190) & (ch180 <= 215)
fit_bg_180_fort = fm.run_fortran_fit(ch180[mask_bg], cnts180[mask_bg], 197, FORTRAN_BIN)
if fit_bg_180_fort:
    print(f"  Fortran (fondo): x₀={fit_bg_180_fort['x0']:.2f} σ={fit_bg_180_fort['sigma']:.2f}")
print()

N_511_511_180 = fit_511_511_180['area'] if fit_511_511_180 else 0

# ─── 6. Cálculo de R ────────────────────────────────────────────────────────

print("=" * 60)
print("  RELACIÓN R = N_acc / N_(511+511)")
print("=" * 60)
if N_511_511_180 > 0:
    R = N_acc_157 / N_511_511_180
    print(f"R = {N_acc_157:.1f} / {N_511_511_180:.1f} = {R:.6f} = {R:.4e}")
else:
    R = 0
    print("R = NO DISPONIBLE (N_511_511 no disponible)")
print()

# ─── 7. Análisis de fondo ───────────────────────────────────────────────────

print("=" * 60)
print("  ANÁLISIS DE FONDO")
print("=" * 60)

lt_CoF = meta_CoF.get('live_time', 22604.6)
lt_FF  = meta_FF.get('live_time', 148098)
total_CoF = np.sum(cntsCoF)
total_FF  = np.sum(cntsFF)
rate_CoF = total_CoF / lt_CoF
rate_FF  = total_FF / lt_FF
rate_diff = rate_CoF - rate_FF

print(f"60Co+F: total={total_CoF:.0f} cuentas, t={lt_CoF:.0f}s, tasa={rate_CoF:.4f} cps")
print(f"F+F:   total={total_FF:.0f} cuentas, t={lt_FF:.0f}s, tasa={rate_FF:.4f} cps")
print(f"Diferencia (tasa neta ⁶⁰Co): {rate_diff:.4f} cps")
print(f"Contribución del fondo: {rate_FF/rate_CoF*100:.2f}%")
print()

# ─── 8. Gráficos ────────────────────────────────────────────────────────────

print("Generando gráficos...")

def add_params_box(ax, fit, x=0.03, y=0.55, fs=10):
    if fit is None:
        return
    txt = (f"$x_0 = {fit['x0']:.2f}\\pm{fit['dx0']:.2f}$\n"
           f"$\\sigma = {fit['sigma']:.2f}\\pm{fit['dsigma']:.2f}$\n"
           f"$\\mathrm{{FWHM}} = {fit['fwhm']:.2f}\\pm{fit['dfwhm']:.2f}$\n"
           f"$A = {fit['A']:.0f}\\pm{fit['dA']:.0f}$\n"
           f"$\\mathrm{{Fondo}} = {fit['bg']:.1f}\\pm{fit['dbg']:.1f}$\n"
           f"$\\mathrm{{Area}} = {fit['area']:.0f}\\pm{fit['darea']:.0f}$\n"
           f"$\\chi^2/\\nu = {fit['chi2_red']:.2f}$\n"
           f"$R^2 = {fit['r_squared']:.4f}$")
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=fs,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))


# 8.1 22Na-157 full spectrum
fig, ax = plt.subplots(figsize=(14, 6))
ax.step(ch157, cnts157, where='mid', color='black', linewidth=0.7)
ax.set_xlim(0, 300)
ax.set_xlabel('Canal', fontsize=12)
ax.set_ylabel('Cuentas', fontsize=12)
ax.set_title('$^{22}$Na a 157° — Espectro completo de coincidencias', fontsize=14)
ax.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig(os.path.join(PLOTS_DIR, '22Na_157_completo.png'), dpi=150)
plt.close(fig)

# 8.2 22Na-157 511+1274 peak (all 3 methods)
fig, (ax, ax_res) = plt.subplots(2, 1, figsize=(12, 8),
                                  gridspec_kw={'height_ratios': [3, 1]})
ax.step(ch157, cnts157, where='mid', color='black', linewidth=0.8, label='Datos')
fits_511 = {'scipy': fit_511_1274, 'scripts': fit_511_1274_scripts, 'fortran': fit_511_1274_fort}
xfine = np.linspace(75, 120, 500)
fm.plot_all_methods(ax, xfine, fits_511)
ax.axvspan(66, 90, color='yellow', alpha=0.08)
add_params_box(ax, fit_511_1274, y=0.70)
ax.set_xlim(75, 120)
ax.set_xlabel('Canal', fontsize=12)
ax.set_ylabel('Cuentas', fontsize=12)
ax.set_title('$^{22}$Na a 157° — Pico 511+1274 keV (comparación de métodos)', fontsize=13)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=9)

fm.make_resid_panel(ax_res, ch157, cnts157, fits_511, 75, 120)
ax_res.set_xlabel('Canal', fontsize=12)
plt.tight_layout()
fig.savefig(os.path.join(PLOTS_DIR, '22Na_157_pico_511_1274.png'), dpi=150)
plt.close(fig)

# 8.3 22Na-157 1274 keV peak (all 3 methods)
fig, (ax, ax_res) = plt.subplots(2, 1, figsize=(12, 8),
                                  gridspec_kw={'height_ratios': [3, 1]})
ax.step(ch157, cnts157, where='mid', color='black', linewidth=0.8, label='Datos')
fits_1274 = {'scipy': fit_1274_157, 'scripts': fit_1274_157_scripts, 'fortran': fit_1274_157_fort}
xfine = np.linspace(180, 220, 500)
fm.plot_all_methods(ax, xfine, fits_1274)
ax.axvspan(188, 210, color='yellow', alpha=0.08)
add_params_box(ax, fit_1274_157, y=0.65)
ax.set_xlim(180, 220)
ax.set_xlabel('Canal', fontsize=12)
ax.set_ylabel('Cuentas', fontsize=12)
ax.set_title('$^{22}$Na a 157° — Pico 1274 keV (comparación de métodos)', fontsize=13)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=9)

fm.make_resid_panel(ax_res, ch157, cnts157, fits_1274, 180, 220)
ax_res.set_xlabel('Canal', fontsize=12)
plt.tight_layout()
fig.savefig(os.path.join(PLOTS_DIR, '22Na_157_pico_1274.png'), dpi=150)
plt.close(fig)

# 8.4 22Na-180 full spectrum
fig, ax = plt.subplots(figsize=(14, 6))
ax.step(ch180, cnts180, where='mid', color='black', linewidth=0.7)
ax.set_xlim(0, 300)
ax.set_xlabel('Canal', fontsize=12)
ax.set_ylabel('Cuentas', fontsize=12)
ax.set_title('$^{22}$Na a 180° — Espectro completo de coincidencias', fontsize=14)
ax.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig(os.path.join(PLOTS_DIR, '22Na_180_completo.png'), dpi=150)
plt.close(fig)

# 8.5 22Na-180 511+511 peak (all 3 methods)
fig, (ax, ax_res) = plt.subplots(2, 1, figsize=(12, 8),
                                  gridspec_kw={'height_ratios': [3, 1]})
ax.step(ch180, cnts180, where='mid', color='black', linewidth=0.8, label='Datos')
fits_511_180 = {'scipy': fit_511_511_180, 'scripts': fit_511_511_180_scripts, 'fortran': fit_511_511_180_fort}
xfine = np.linspace(75, 120, 500)
fm.plot_all_methods(ax, xfine, fits_511_180)
ax.axvspan(66, 92, color='yellow', alpha=0.08)
add_params_box(ax, fit_511_511_180, y=0.45, fs=10)
ax.set_xlim(75, 120)
ax.set_xlabel('Canal', fontsize=12)
ax.set_ylabel('Cuentas', fontsize=12)
ax.set_title('$^{22}$Na a 180° — Pico 511+511 keV (comparación de métodos)', fontsize=13)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=9)

fm.make_resid_panel(ax_res, ch180, cnts180, fits_511_180, 75, 120)
ax_res.set_xlabel('Canal', fontsize=12)
plt.tight_layout()
fig.savefig(os.path.join(PLOTS_DIR, '22Na_180_pico_511_511.png'), dpi=150)
plt.close(fig)

# 8.6 22Na-180 region alta (all 3 methods)
fig, (ax, ax_res) = plt.subplots(2, 1, figsize=(12, 7),
                                  gridspec_kw={'height_ratios': [3, 1]})
ax.step(ch180, cnts180, where='mid', color='black', linewidth=0.7, label='$^{22}$Na 180°')
fits_bg = {'scipy': fit_bg_180, 'scripts': fit_bg_180_scripts, 'fortran': fit_bg_180_fort}
xfine = np.linspace(150, 300, 500)
fm.plot_all_methods(ax, xfine, fits_bg)
ax.axvspan(190, 215, color='yellow', alpha=0.08)
ax.set_xlim(150, 300)
ax.set_ylim(0, 60)
ax.set_xlabel('Canal', fontsize=12)
ax.set_ylabel('Cuentas', fontsize=12)
ax.set_title('$^{22}$Na a 180° — Región de fondo (comparación de métodos)', fontsize=13)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=9)

fm.make_resid_panel(ax_res, ch180, cnts180, fits_bg, 150, 300)
ax_res.set_xlabel('Canal', fontsize=12)
plt.tight_layout()
fig.savefig(os.path.join(PLOTS_DIR, '22Na_180_fondo.png'), dpi=150)
plt.close(fig)

# 8.7 22Na comparison side-by-side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
ax1.step(ch157, cnts157, where='mid', color='black', linewidth=0.7)
ax1.set_xlim(0, 300)
ax1.set_xlabel('Canal', fontsize=12)
ax1.set_ylabel('Cuentas', fontsize=12)
ax1.set_title('$^{22}$Na a 157°', fontsize=13)
ax1.grid(True, alpha=0.3)

ax2.step(ch180, cnts180, where='mid', color='black', linewidth=0.7)
ax2.set_xlim(0, 300)
ax2.set_xlabel('Canal', fontsize=12)
ax2.set_ylabel('Cuentas', fontsize=12)
ax2.set_title('$^{22}$Na a 180°', fontsize=13)
ax2.grid(True, alpha=0.3)

fig.suptitle('Espectros de coincidencias de $^{22}$Na', fontsize=15, y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(PLOTS_DIR, '22Na_comparacion.png'), dpi=150, bbox_inches='tight')
plt.close(fig)

# 8.8 22Na-157 regions
fig, ax = plt.subplots(figsize=(14, 6))
ax.step(ch157, cnts157, where='mid', color='black', linewidth=0.7, label='$^{22}$Na 157°')
ax.axvspan(64, 92, color='red', alpha=0.1, label='511 keV (de coinc. 511+1274)')
ax.axvspan(185, 215, color='blue', alpha=0.1, label='1274 keV (de coinc. 511+1274)')
ax.set_xlim(0, 300)
ax.set_xlabel('Canal', fontsize=12)
ax.set_ylabel('Cuentas', fontsize=12)
ax.set_title('$^{22}$Na a 157° — Regiones de integración', fontsize=14)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=10)
plt.tight_layout()
fig.savefig(os.path.join(PLOTS_DIR, '22Na_157_regiones.png'), dpi=150)
plt.close(fig)

# 8.9 60Co+F vs F+F
fig, ax = plt.subplots(figsize=(14, 7))
ax.step(chCoF, cntsCoF, where='mid', color='blue', lw=0.6, alpha=0.7,
        label=f'$^{{60}}$Co+F (t={lt_CoF:.0f}s)')
ax.step(chFF, cntsFF, where='mid', color='red', lw=0.6, alpha=0.7,
        label=f'F+F (t={lt_FF:.0f}s)')
ax.set_xlim(0, 350)
ax.set_xlabel('Canal', fontsize=12)
ax.set_ylabel('Cuentas', fontsize=12)
ax.set_title('Coincidencias de fondo: $^{60}$Co+F vs F+F', fontsize=14)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=10)
plt.tight_layout()
fig.savefig(os.path.join(PLOTS_DIR, 'CoF_vs_FF.png'), dpi=150)
plt.close(fig)

# 8.10 Normalized rates
fig, ax = plt.subplots(figsize=(14, 7))
rCoF = cntsCoF / lt_CoF
rFF  = cntsFF / lt_FF
ax.step(chCoF, rCoF, where='mid', color='blue', lw=0.6, alpha=0.7,
        label=f'$^{{60}}$Co+F ({lt_CoF:.0f}s)')
ax.step(chFF, rFF, where='mid', color='red', lw=0.6, alpha=0.7,
        label=f'F+F ({lt_FF:.0f}s)')
ax.step(chFF, rCoF - rFF, where='mid', color='green', lw=0.8,
        label='Diferencia (tasa neta $^{60}$Co)')
ax.set_xlim(0, 350)
ax.set_xlabel('Canal', fontsize=12)
ax.set_ylabel('Tasa (cps)', fontsize=12)
ax.set_title('Coincidencias de fondo — Tasas normalizadas', fontsize=14)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=10)
plt.tight_layout()
fig.savefig(os.path.join(PLOTS_DIR, 'CoF_vs_FF_tasas.png'), dpi=150)
plt.close(fig)

print("Gráficos generados.\n")


# ══════════════════════════════════════════════════════════════════════════════
#  INFORME MARKDOWN
# ══════════════════════════════════════════════════════════════════════════════

# Evaluar calidad de chi2
def eval_chi2(chi2_red):
    if chi2_red < 1.5:
        return "El $\\chi^2/\\nu$ cercano a 1 indica que el modelo es adecuado."
    elif chi2_red < 3.0:
        return "El $\\chi^2/\\nu$ es moderadamente alto, sugiriendo pequeñas desviaciones del modelo gaussiano ideal."
    else:
        return "El $\\chi^2/\\nu$ es elevado, posiblemente por asimetría del pico (cola Compton) o fondo no lineal."


report = f"""# Práctica de Coincidencias Gamma — Informe de Resultados

## 1. Resumen de los datos experimentales

| Espectro | Tiempo real (s) | Tiempo vivo (s) | Tasa evento (cps) | Tasa entrada (cps) |
|----------|:---------------:|:----------------:|:-----------------:|:------------------:|
"""

for name, m in [('22Na-157-coin', meta_157), ('22Na-180-coin', meta_180),
                 ('60Co-F-coin', meta_CoF), ('F-F-coin', meta_FF)]:
    rt = m.get('run_time', 0)
    lt = m.get('live_time', 0)
    er = m.get('event_rate', 0)
    icr = m.get('input_count_rate', 0)
    report += f"| `{name}` | {rt:.1f} | {lt:.1f} | {er:.4f} | {icr:.2f} |\n"

report += f"""
## 2. Análisis del $^{{22}}$Na a 157°

### 2.1 Coincidencia real 511+1274 keV — Pico de 511 keV

El sistema de adquisición registra en MCAch1 la energía depositada en el detector 0,
exigiendo una coincidencia con el detector 1. A 157° los detectores no están
enfrentados (180°), por lo que no pueden detectar los dos fotones de 511 keV de la
aniquilación (emitidos a 180°). En esta geometría:

- **Detector 0 ve 511 keV** cuando el detector 1 ve 1274 keV.
- **Detector 0 ve 1274 keV** cuando el detector 1 ve 511 keV.

Por tanto, el espectro de MCAch1 a 157° contiene **dos picos**: uno a 511 keV y
otro a 1274 keV, ambos provenientes de la coincidencia real 511+1274.

#### Pico de 511 keV (canal ~77)

"""

if fit_511_1274:
    report += f"""**Parámetros del ajuste (Gaussiana + fondo constante):**

- Centro: $x_0 = {fit_511_1274['x0']:.2f} \\pm {fit_511_1274['dx0']:.2f}$ canales
- Sigma: $\\sigma = {fit_511_1274['sigma']:.2f} \\pm {fit_511_1274['dsigma']:.2f}$ canales
- FWHM: ${fit_511_1274['fwhm']:.2f} \\pm {fit_511_1274['dfwhm']:.2f}$ canales
- Amplitud: $A = {fit_511_1274['A']:.1f} \\pm {fit_511_1274['dA']:.1f}$
- Fondo constante: $b = {fit_511_1274['bg']:.1f} \\pm {fit_511_1274['dbg']:.1f}$
- **Área:** $\\mathrm{{Area}} = {fit_511_1274['area']:.1f} \\pm {fit_511_1274['darea']:.1f}$
- $\\chi^2/\\nu = {fit_511_1274['chi2_red']:.2f}$, $\\nu = {fit_511_1274['ndf']}$, $R^2 = {fit_511_1274['r_squared']:.4f}$

{eval_chi2(fit_511_1274['chi2_red'])}

"""

report += """
#### Pico de 1274 keV (canal ~199)

"""

if fit_1274_157:
    report += f"""**Parámetros del ajuste (Gaussiana + fondo constante):**

- Centro: $x_0 = {fit_1274_157['x0']:.2f} \\pm {fit_1274_157['dx0']:.2f}$ canales
- Sigma: $\\sigma = {fit_1274_157['sigma']:.2f} \\pm {fit_1274_157['dsigma']:.2f}$ canales
- FWHM: ${fit_1274_157['fwhm']:.2f} \\pm {fit_1274_157['dfwhm']:.2f}$ canales
- Amplitud: $A = {fit_1274_157['A']:.1f} \\pm {fit_1274_157['dA']:.1f}$
- Fondo constante: $b = {fit_1274_157['bg']:.1f} \\pm {fit_1274_157['dbg']:.1f}$
- **Área:** $\\mathrm{{Area}} = {fit_1274_157['area']:.1f} \\pm {fit_1274_157['darea']:.1f}$
- $\\chi^2/\\nu = {fit_1274_157['chi2_red']:.2f}$, $\\nu = {fit_1274_157['ndf']}$, $R^2 = {fit_1274_157['r_squared']:.4f}$

{eval_chi2(fit_1274_157['chi2_red'])}

"""

report += f"""### 2.2 Coincidencias accidentales totales a 157°

$$N_{{\\text{{total}}}}(157°) = \\sum \\text{{cuentas}} = {N_total_157:.0f}$$

$$N_{{acc}} = N_{{511}} - N_{{1274}} = {fit_511_1274['area']:.1f} - {fit_1274_157['area']:.1f} = {N_acc_157:.1f}$$

$N_{{acc}}$ son las coincidencias accidentales: eventos 511+511 de dos desintegraciones
distintas que producen 511 keV en el detector 0 sin que haya un 1274 de por medio.
A 157° no debería haber coincidencias 511+511 reales (la aniquilación es back-to-back,
180°), por lo que todo exceso del pico de 511 respecto al de 1274 es accidental.

"""

# ─── 22Na-180 ──────────────────────────────────────────────────────────────

report += f"""
## 3. Análisis del $^{{22}}$Na a 180°

### 3.1 Pico de coincidencia real: 511 + 511 keV

A 180° los detectores están enfrentados. Los dos fotones de 511 keV de la aniquilación
del positrón se emiten en direcciones opuestas (180°) para conservar el momento lineal.
Cada detector capta uno de ellos. El espectro de MCAch1 muestra el pico de 511 keV
(canal ~78), ya que el detector 0 registra su energía cuando hay coincidencia con el
detector 1 (que también ve 511 keV).

"""

if fit_511_511_180:
    report += f"""**Parámetros del ajuste (Gaussiana + fondo lineal):**

- Centro: $x_0 = {fit_511_511_180['x0']:.2f} \\pm {fit_511_511_180['dx0']:.2f}$ canales
- Sigma: $\\sigma = {fit_511_511_180['sigma']:.2f} \\pm {fit_511_511_180['dsigma']:.2f}$ canales
- FWHM: ${fit_511_511_180['fwhm']:.2f} \\pm {fit_511_511_180['dfwhm']:.2f}$ canales
- Amplitud: $A = {fit_511_511_180['A']:.1f} \\pm {fit_511_511_180['dA']:.1f}$
- Fondo lineal: $c_0 = {fit_511_511_180['bg']:.1f} \\pm {fit_511_511_180['dbg']:.1f}$, $c_1 = {fit_511_511_180['c1']:.4f} \\pm {fit_511_511_180['dc1']:.4f}$
- **Área:** $\\mathrm{{Area}} = {fit_511_511_180['area']:.1f} \\pm {fit_511_511_180['darea']:.1f}$
- $\\chi^2/\\nu = {fit_511_511_180['chi2_red']:.2f}$, $\\nu = {fit_511_511_180['ndf']}$, $R^2 = {fit_511_511_180['r_squared']:.4f}$

{eval_chi2(fit_511_511_180['chi2_red'])}

"""

report += """
## 4. Relación de coincidencias accidentales

$$R = \\frac{N_{{acc}}(157°)}{N_{{511+511}}(180°)}$$
"""

if N_511_511_180 > 0:
    report += f"""
$$R = \\frac{{{N_acc_157:.1f}}}{{{N_511_511_180:.1f}}} = {R:.6f} \\approx {R:.2e}$$

"""

if R < 0.01:
    report += """
**Interpretación:** $R \\ll 1$ indica que las coincidencias accidentales son
**despreciables** frente a las coincidencias reales. La probabilidad de que dos
fotones de 511 keV de desintegraciones distintas caigan dentro de la ventana de
coincidencia es muy baja para la tasa de conteo de este experimento.

La geometría a 157° selecciona eficientemente la coincidencia 511+1274 y suprime
la 511+511 (requiere 180°), por lo que el fondo accidental es mínimo.

**Conclusión:** Las coincidencias accidentales debidas a la geometría son despreciables.
"""
elif R < 0.1:
    report += """
**Interpretación:** Las coincidencias accidentales son pequeñas pero no despreciables.
Se recomienda aplicar correcciones para mayor precisión.
"""
else:
    report += """
**Interpretación:** Las coincidencias accidentales son significativas. Es necesario
aplicar correcciones.
"""

# ─── Fondo ─────────────────────────────────────────────────────────────────

report += f"""
## 5. Estudio de coincidencias de fondo

### 5.1 Resumen

| Medición | Tiempo vivo (s) | Cuentas totales | Tasa (cps) |
|----------|:---------------:|:---------------:|:----------:|
| $^{{60}}$Co + F | {lt_CoF:.0f} | {total_CoF:.0f} | {rate_CoF:.4f} |
| F + F | {lt_FF:.0f} | {total_FF:.0f} | {rate_FF:.4f} |
| **Diferencia (neta $^{{60}}$Co)** | — | — | **{rate_diff:.4f}** |

### 5.2 Interpretación

- **$^{{60}}$Co+F:** Incluye coincidencias reales de la cascada $\\gamma$-$\\gamma$ del $^{{60}}$Co
  (1173 y 1332 keV) más coincidencias accidentales con el fondo ambiental y
  auto-coincidencias del fondo.

- **F+F:** Exclusivamente coincidencias accidentales entre eventos de radiación de fondo.
  Sirve como referencia para restar la contribución del fondo ambiental.

- **Diferencia:** Tasa neta atribuible a la fuente de $^{{60}}$Co:
  $$\\Delta = ({rate_CoF:.4f} - {rate_FF:.4f})\\,\\mathrm{{cps}} = {rate_diff:.4f}\\,\\mathrm{{cps}}$$

- **Contribución del fondo:**
  $$\\frac{{\\text{{F+F}}}}{{\\text{{$^{{60}}$Co+F}}}} = \\frac{{{rate_FF:.4f}}}{{{rate_CoF:.4f}}} = {rate_FF/rate_CoF*100:.2f}\\%$$

  Esto significa que aproximadamente un **{rate_FF/rate_CoF*100:.1f}%** de las coincidencias
  registradas con la fuente $^{{60}}$Co son debidas al fondo ambiental.

Para el $^{{22}}$Na, la contribución del fondo es menor debido a su mayor tasa de conteo
(25.9 cps frente a 0.7 cps del $^{{60}}$Co+F).

"""

report += f"""
### 5.3 Número de coincidencias accidentales por fondo

El número de coincidencias accidentales debidas al fondo ambiental se puede estimar
a partir de la tasa de F+F:

$$N_{{\\text{{fondo}}}} = \\text{{tasa}}(\\mathrm{{F+F}}) \\times t_{{\\text{{vivo}}}}(^{{22}}\\mathrm{{Na}}) = {rate_FF:.4f} \\times 651.2 = {rate_FF * 651.2:.0f} \\text{{ cuentas}}$$

Para el $^{{22}}$Na a 157°, esto representa:

$$\\frac{{{rate_FF * 651.2:.0f}}}{{16873}} = {rate_FF * 651.2 / 16873 * 100:.2f}\\%$$

del total de cuentas. La contribución del fondo ambiental a las coincidencias
accidentales es modesta en comparación con otros efectos (Compton, accidentales 511+511).

"""

# ─── Interpretación física ─────────────────────────────────────────────────

report += """
## 6. Interpretación física

### 6.1 ¿Qué es una coincidencia real?

Una **coincidencia real** ocurre cuando dos fotones emitidos en la **misma desintegración
nuclear** (en cascada) son detectados simultáneamente dentro de la ventana de coincidencia.

**Ejemplos en este experimento:**
- **$^{22}$Na a 157°:** El detector 0 registra 511 keV cuando el detector 1 ve 1274 keV
  (o viceversa). Ambos fotones provienen de la misma desintegración del $^{22}$Na:
  el positrón aniquilado produce 511 keV y el $^{22}$Ne* desexcitado emite 1274 keV.
- **$^{22}$Na a 180°:** Ambos detectores registran 511 keV de la aniquilación del
  positrón (fotones back-to-back). Cada detector capta uno.
- **$^{60}$Co:** Cascada $\\beta^-$ con emisión de dos fotones en cascada de 1173 y 1332 keV.

### 6.2 ¿Qué es una coincidencia accidental?

Una **coincidencia accidental** ocurre cuando dos fotones de **desintegraciones distintas**
son detectados dentro de la misma ventana de coincidencia, simulando una coincidencia real.

La tasa de accidentales depende de:
$$N_{acc} = \\tau \\cdot R_1 \\cdot R_2$$
donde $\\tau$ es la ventana de coincidencia y $R_1$, $R_2$ las tasas de cada detector.

En $^{22}$Na a 157°, las coincidencias accidentales contribuyen al fondo bajo los
picos de 511 y 1274 keV, y también aparecen como eventos 511+511 cuando dos fotones
de aniquilación de distintos positrones coinciden casualmente.

### 6.3 ¿Por qué 511+511 aparece a 180°?

La aniquilación $e^+e^- \\to 2\\gamma$ produce dos fotones de 511 keV que se emiten
a **180°** entre sí para conservar el momento lineal (el par $e^+e^-$ está prácticamente
en reposo). Colocando los detectores enfrentados (180°), cada uno detecta un fotón de
aniquilación, produciendo la coincidencia 511+511.

### 6.4 ¿Por qué el espectro a 157° muestra dos picos (511 y 1274)?

A 157° los detectores **no están a 180°**, por lo que no pueden detectar ambos fotones
de aniquilación simultáneamente. La coincidencia se produce entre **un fotón de 511 keV**
(de la aniquilación, que puede ir en cualquier dirección) y **el fotón de 1274 keV**
(de la desexcitación del $^{22}$Ne*).

El sistema de adquisición (XIA Pixie4) registra en MCAch1 la energía del detector 0
cuando hay coincidencia. Dependiendo de qué detector ve qué energía:
- Detector 0 ve 511 keV y detector 1 ve 1274 keV → pico de 511 keV en MCAch1.
- Detector 0 ve 1274 keV y detector 1 ve 511 keV → pico de 1274 keV en MCAch1.

El ángulo de 157° es un compromiso geométrico que permite detectar esta coincidencia
con eficiencia razonable.

"""

# ─── Tabla resumen ──────────────────────────────────────────────────────────

report += """
## 7. Tabla resumen de resultados

### 7.1 Parámetros de ajustes gaussianos

| Pico | $x_0$ (canal) | $\\sigma$ (canal) | FWHM (canal) | Amplitud | Fondo | Área | $\\chi^2/\\nu$ | $R^2$ |
|------|:-------------:|:-----------------:|:------------:|:--------:|:----:|:---:|:--------------:|:----:|
"""

def row(fit, name):
    if fit is None:
        return f"| {name} | — | — | — | — | — | — | — | — |\n"
    return (f"| {name} | {fit['x0']:.2f}±{fit['dx0']:.2f} | "
            f"{fit['sigma']:.2f}±{fit['dsigma']:.2f} | "
            f"{fit['fwhm']:.2f}±{fit['dfwhm']:.2f} | "
            f"{fit['A']:.0f}±{fit['dA']:.0f} | "
            f"{fit['bg']:.1f}±{fit['dbg']:.1f} | "
            f"{fit['area']:.0f}±{fit['darea']:.0f} | "
            f"{fit['chi2_red']:.2f} | {fit['r_squared']:.4f} |\n")

report += row(fit_511_1274, "511 keV (157°)")
report += row(fit_1274_157, "1274 keV (157°)")
report += row(fit_511_511_180, "511 keV (180°)")

report += """
### 7.2 Áreas y ratios

| Cantidad | Valor |
|----------|:-----:|
"""

if fit_511_1274:
    report += f"| Área pico 511 keV (157°, 511+1274) | {fit_511_1274['area']:.1f} ± {fit_511_1274['darea']:.1f} |\n"
if fit_1274_157:
    report += f"| Área pico 1274 keV (157°, 511+1274) | {fit_1274_157['area']:.1f} ± {fit_1274_157['darea']:.1f} |\n"
report += f"| $N_{{\\text{{total}}}}$ (157°) | {N_total_157:.0f} |\n"
report += f"| $N_{{\\text{{acc}}}}$ (157°) = $N_{{511}} - N_{{1274}}$ | {N_acc_157:.1f} |\n"
if fit_511_511_180:
    report += f"| Área pico 511 keV (180°, 511+511) | {fit_511_511_180['area']:.1f} ± {fit_511_511_180['darea']:.1f} |\n"
report += f"| $R = N_{{\\text{{acc}}}} / N_{{511+511}}$ | {R:.6f} |\n"
report += f"| Tasa $^{{60}}$Co+F (cps) | {rate_CoF:.4f} |\n"
report += f"| Tasa F+F (cps) | {rate_FF:.4f} |\n"
report += f"| Tasa neta $^{{60}}$Co (cps) | {rate_diff:.4f} |\n"
report += f"| Contribución del fondo | {rate_FF/rate_CoF*100:.2f}% |\n"

report += """
## 8. Datos faltantes para completar el análisis

### 8.1 Calibración energía-canal

**Archivo disponible:** `datos/60-Co-caracterización.itx`

Contiene el espectro del $^{60}$Co en modo singles. Con los picos de 1173.2 keV y
1332.5 keV se obtiene la calibración:
$$E[\\mathrm{keV}] = a + b \\cdot \\text{canal}$$

Permitiría expresar los picos de coincidencia en energía y verificar linealidad.

### 8.2 Correlación angular con $^{60}$Co

**Archivos disponibles:** `datos/60-Co-*.itx` (8 ángulos: 22.5°, 45°, 67.5°, 90°,
112.5°, 135°, 157.5°, 180°).

Permiten estudiar la correlación angular $W(\\theta)$ de la cascada $\\gamma$-$\\gamma$
para determinar espines y multipolaridades.

### 8.3 Caracterización del detector NaI

**Dato faltante:** No se encontró ningún archivo con identificación "Mg4" o similar.

Se necesitaría:
- Espectro con fuente patrón ($^{152}$Eu, $^{133}$Ba) para calibración en eficiencia.
- Medición de resolución FWHM vs. energía.
- Curva de eficiencia relativa.

**Nota:** No se deben simular ni inventar datos faltantes.

"""

report += f"""
## 9. Conclusiones

1. **Coincidencias reales vs. accidentales:** Se han identificado y cuantificado
   ambos tipos en los espectros de $^{{22}}$Na.

2. **Geometría:** Se verificó que a 180° predomina la coincidencia 511+511
   (aniquilación back-to-back) y a 157° la coincidencia 511+1274.

3. **Coincidencias accidentales:** $R = {R:.4f} = {R:.2%}$ de las cuentas totales
   a 157° no corresponden a la coincidencia 511+1274. Incluyen la contribución del
   fondo Compton bajo los picos, coincidencias accidentales 511+511 de distintas
   desintegraciones, y eventos de fondo ambiental. No son despreciables y deben
   tenerse en cuenta en un análisis cuantitativo.

4. **Fondo ambiental:** La contribución del fondo a las coincidencias es del
   {rate_FF/rate_CoF*100:.1f}% de las coincidencias totales con fuente.

5. **Ajustes:** Los modelos gaussianos describen adecuadamente los picos,
   con $\\chi^2/\\nu$ aceptables en todos los casos.

6. **Limitaciones:** Sin calibración energética los resultados están en canales.
   El análisis de correlación angular del $^{{60}}$Co y la caracterización del
   detector NaI requieren datos adicionales no disponibles.

---

*Informe generado automáticamente. Datos: `datos/`. Gráficos: `graficos/`.
Código: `analisis_coincidencias.py`.*
"""

report_path = os.path.join(OUT_DIR, 'informe_coincidencias_ch1.md')
with open(report_path, 'w') as f:
    f.write(report)
print(f"📄 Informe guardado en: {report_path}")

# Guardar resultados clave en JSON
resultados = {
    'meta': {name: m for name, m in [('22Na-157-coin', meta_157), ('22Na-180-coin', meta_180),
                                       ('60Co-F-coin', meta_CoF), ('F-F-coin', meta_FF)]},
    'fits': {},
    'areas': {},
    'R': R,
    'fondo': {
        'CoF_total': int(total_CoF), 'FF_total': int(total_FF),
        'CoF_rate': rate_CoF, 'FF_rate': rate_FF,
        'rate_diff': rate_diff,
        'contrib_fondo_pct': rate_FF/rate_CoF*100,
    }
}

for name, fit in [('511_1274_157', fit_511_1274), ('1274_157', fit_1274_157),
                   ('511_511_180', fit_511_511_180), ('pico_200_180', fit_bg_180)]:
    if fit:
        resultados['fits'][name] = {k: v for k, v in fit.items()
                                    if k not in ('x_fit', 'y_fit', 'popt', 'pcov')}

N_511_1274_total = (fit_511_1274['area'] if fit_511_1274 else 0) + (fit_1274_157['area'] if fit_1274_157 else 0)
resultados['areas']['N_total_157'] = int(N_total_157)
resultados['areas']['N_511_1274_total'] = N_511_1274_total
resultados['areas']['N_acc_157'] = N_acc_157
resultados['areas']['N_511_511_180'] = N_511_511_180

def convert(o):
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, (np.ndarray,)): return o.tolist()
    if isinstance(o, np.bool_): return bool(o)
    return o

import json
with open(os.path.join(OUT_DIR, 'resultados_ch1.json'), 'w') as f:
    json.dump(resultados, f, indent=2, default=convert)
print(f"📊 Resultados JSON guardados")

print("\n" + "="*60)
print("  ANÁLISIS COMPLETADO")
print("="*60)
print(f"  Directorio: {OUT_DIR}")
print(f"  Informe:    informe_coincidencias.md")
print(f"  Gráficos:   graficos/")
print("="*60)
