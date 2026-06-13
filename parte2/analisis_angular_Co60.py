#!/usr/bin/env python3
"""
Análisis de correlación angular del ⁶⁰Co en coincidencia.
Procesa archivos 60-Co-XXX-coin.itx para MCAch1 y genera:
  - Espectros superpuestos normalizados por tiempo
  - Curva de correlación W(θ) con ajuste teórico
  - Tabla resumen para el informe
"""

import os, sys, glob, math, json, re
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_DIR = os.path.abspath(os.path.join(OUT_DIR, '..'))
DATA_DIR = os.path.join(PROJ_DIR, 'datos')
PLOTS_DIR = os.path.join(OUT_DIR, 'graficos_ch1')
os.makedirs(PLOTS_DIR, exist_ok=True)

plt.rcParams.update({'font.size': 11, 'figure.dpi': 150})

# ─── Parse .itx for MCAch1 ──────────────────────────────────────────────

def parse_itx_ch(filepath, ch='MCAch1'):
    meta = {}
    counts = []
    with open(filepath, 'r', errors='replace') as f:
        lines = f.readlines()
    in_ch = False
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith('X //') or s.startswith('X ////'):
            if 'Run Time [s]' in s:
                meta['run_time'] = float(s.split('=')[1].strip().split()[0])
            elif 'Event Rate [cps]' in s:
                meta['event_rate'] = float(s.split('=')[1].strip().split()[0])
            # Track which channel we're in for metadata
            if 'Channel 0' in s:
                meta['_last_ch'] = 'Channel 0'
            elif 'Channel 1' in s:
                meta['_last_ch'] = 'Channel 1'
            if 'Live Time [s]' in s:
                ch_tag = meta.get('_last_ch', '')
                if 'Channel 0' in ch_tag:
                    meta['live_time_ch0'] = float(s.split('=')[1].strip().split()[0])
                elif 'Channel 1' in ch_tag:
                    meta['live_time_ch1'] = float(s.split('=')[1].strip().split()[0])
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
    meta.pop('_last_ch', None)
    return meta, np.array(counts, dtype=float)

def cargar_espectro(nombre):
    path = os.path.join(DATA_DIR, nombre)
    if not os.path.exists(path):
        path = os.path.join(DATA_DIR, nombre + '.itx')
    meta, counts = parse_itx_ch(path, 'MCAch1')
    channels = np.arange(len(counts), dtype=float)
    return meta, channels, counts

# ─── Gaussian fit ───────────────────────────────────────────────────────

def gauss_linear(x, A, x0, s, c0, c1):
    return A * np.exp(-(x - x0)**2 / (2 * s**2)) + c0 + c1 * x

def fit_peak_1332(ch, counts, window=15):
    est = 205
    lo, hi = est - window, est + window
    mask = (ch >= lo) & (ch <= hi)
    x = ch[mask]
    y = counts[mask]
    i_max = np.argmax(y)
    A0 = y[i_max] - np.median(y[:min(5, len(y)//3)] + y[-min(5, len(y)//3):])
    x00 = x[i_max]
    s0 = 4.0
    c00 = np.median(y[:3]) if len(y) > 6 else np.median(y)
    c01 = 0.0
    try:
        popt, _ = curve_fit(
            gauss_linear, x, y,
            p0=[A0, x00, s0, c00, c01],
            bounds=([0, lo-3, 0.5, -np.inf, -np.inf],
                    [np.inf, hi+3, 15, np.inf, np.inf]),
            sigma=np.where(y > 1, np.sqrt(y), 1.0),
            maxfev=10000
        )
        A, x0, sigma, bg, c1 = popt
        y_fit = gauss_linear(x, *popt)
        resid = y - y_fit
        chi2 = np.sum((resid / np.maximum(np.sqrt(y), 0.5))**2)
        ndf = len(x) - 5
        area = A * sigma * math.sqrt(2 * math.pi)
        fwhm = 2.35482 * sigma
        return {
            'x0': x0, 'sigma': sigma, 'A': A, 'bg': bg, 'c1': c1,
            'area': area, 'chi2': chi2, 'ndf': ndf,
            'chi2_red': chi2/ndf if ndf > 0 else 0,
            'x_fit': x, 'y_fit': y_fit,
        }
    except Exception as e:
        return None

# ─── Angle name helpers ────────────────────────────────────────────────

def parse_angle(basename):
    m = re.search(r'([\d,]+)-coin', basename)
    if not m:
        m = re.search(r'(\d+)-coin', basename)
    if m:
        s = m.group(1).replace(',', '.')
        return float(s)
    return None

angle_label = {
    22.5: '22.5°', 45: '45°', 67.5: '67.5°', 90: '90°',
    112.5: '112.5°', 135: '135°', 157.5: '157.5°', 180: '180°',
}

# ─── Main ──────────────────────────────────────────────────────────────

print('='*60)
print('  CORRELACIÓN ANGULAR ⁶⁰Co (MCAch1)')
print('='*60)

itx_files = sorted(glob.glob(os.path.join(DATA_DIR, '60*-coin.itx')))
itx_files = [f for f in itx_files
             if 'caracterizacion' not in f and 'Co-F' not in f]

results = []
for fp in itx_files:
    name = os.path.basename(fp)
    base = os.path.splitext(name)[0]
    ang = parse_angle(base)
    if ang is None:
        print(f'  Skipping {name} (no angle detected)')
        continue

    meta, ch, counts = cargar_espectro(base)
    lt = meta.get('live_time_ch1', meta.get('run_time', 1))
    total_counts = np.sum(counts)
    rate = total_counts / lt

    fit = fit_peak_1332(ch, counts)
    if fit:
        area_norm = fit['area'] / lt
        results.append({
            'angle': ang,
            'label': angle_label.get(ang, f'{ang}°'),
            'basename': base,
            'run_time': meta.get('run_time', 0),
            'live_time': lt,
            'total_counts': int(total_counts),
            'rate': rate,
            'fit': fit,
            'area': fit['area'],
            'area_norm': area_norm,
            'x0': fit['x0'],
            'sigma': fit['sigma'],
        })
        print(f'  {name:28s}  θ={angle_label.get(ang,"?"):>8s}  '
              f'area={fit["area"]:.0f}  '
              f'area/t={area_norm:.4f} cps  '
              f'χ²/ν={fit["chi2_red"]:.1f}')
    else:
        print(f'  {name:28s}  θ={angle_label.get(ang,"?"):>8s}  NO FIT')

print()

# ─── Plots ─────────────────────────────────────────────────────────────

colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(results)))

# 1) Overlay of all normalized spectra
fig, ax = plt.subplots(figsize=(14, 7))
for r, c in zip(results, colors):
    _, ch, counts = cargar_espectro(r['basename'])
    rate_spec = counts / r['live_time']
    ax.step(ch, rate_spec, where='mid', color=c, linewidth=0.5,
            label=f"θ={r['label']} ({r['rate']:.2f} cps)")
ax.set_xlim(100, 250)
ax.set_xlabel('Canal', fontsize=12)
ax.set_ylabel('Tasa (cps)', fontsize=12)
ax.set_title('$^{60}$Co — Espectros de coincidencia vs. ángulo (MCAch1)', fontsize=13)
ax.legend(fontsize=8, ncol=2)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, 'Co60_angular_overlay.png'), dpi=150)
plt.close()
print('  -> Co60_angular_overlay.png')

# 2) Angular correlation: W(θ) = normalized area vs angle
angles = np.array([r['angle'] for r in results])
areas_norm = np.array([r['area_norm'] for r in results])
areas_err = np.array([r['area'] / math.sqrt(max(r['area'], 1)) / r['live_time']
                      for r in results])

# Theoretical: W(θ) = 1 + A2*P2(cosθ) + A4*P4(cosθ)
def P2(x): return (3*x**2 - 1)/2
def P4(x): return (35*x**4 - 30*x**2 + 3)/8

def W(theta_deg, A2, A4, norm):
    th = np.radians(theta_deg)
    ct = np.cos(th)
    return norm * (1 + A2*P2(ct) + A4*P4(ct))

cos_angles = np.cos(np.radians(angles))

# Simple normalization: divide by mean
areas_norm_mean = areas_norm / np.mean(areas_norm)
areas_err_norm = areas_err / np.mean(areas_norm)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Left: W(θ) vs angle
theta_fine = np.linspace(0, 180, 200)
ax1.errorbar(angles, areas_norm_mean, yerr=areas_err_norm,
             fmt='o', color='#d62728', capsize=4, ms=8, label='Datos MCAch1')
try:
    popt, _ = curve_fit(lambda th, A2, A4: W(th, A2, A4, 1.0),
                        angles, areas_norm_mean,
                        p0=[0.1, 0.05], maxfev=10000,
                        sigma=areas_err_norm, absolute_sigma=False)
    A2_fit, A4_fit = popt
    ax1.plot(theta_fine, W(theta_fine, A2_fit, A4_fit, 1.0),
             'r-', lw=2, label=f'Ajuste: A₂={A2_fit:.3f}, A₄={A4_fit:.3f}')
    print(f'  Ajuste W(θ): A₂={A2_fit:.4f}, A₄={A4_fit:.4f}')
except Exception as e:
    print(f'  Ajuste no convergió: {e}')

ax1.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
ax1.set_xlabel('Ángulo θ (grados)', fontsize=12)
ax1.set_ylabel('W(θ) normalizado', fontsize=12)
ax1.set_title('Correlación angular $^{60}$Co — Pico 1332 keV', fontsize=13)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Right: W(θ) vs cosθ
cos_fine = np.cos(np.radians(theta_fine))
ax2.errorbar(cos_angles, areas_norm_mean, yerr=areas_err_norm,
             fmt='o', color='#d62728', capsize=4, ms=8, label='Datos MCAch1')
if 'A2_fit' in dir():
    ax2.plot(cos_fine, W(theta_fine, A2_fit, A4_fit, 1.0),
             'r-', lw=2, label=f'A₂={A2_fit:.3f}, A₄={A4_fit:.3f}')
ax2.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
ax2.set_xlabel('cos θ', fontsize=12)
ax2.set_ylabel('W(θ) normalizado', fontsize=12)
ax2.set_title('Correlación angular — vs. cos θ', fontsize=13)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, 'Co60_angular_correlacion.png'), dpi=150)
plt.close()
print('  -> Co60_angular_correlacion.png')

# 3) Multi-panel: one spectrum per angle
n = len(results)
ncols = 4
nrows = (n + ncols - 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(16, 3*nrows))
axes = axes.flatten()
for i, (r, c) in enumerate(zip(results, colors)):
    ax = axes[i]
    _, ch, counts = cargar_espectro(r['basename'])
    rate_spec = counts / r['live_time']
    ax.step(ch, rate_spec, where='mid', color=c, linewidth=0.6)
    # Overlay fit
    f = r['fit']
    if f is not None:
        x_fit = f['x_fit']
        ax.plot(x_fit, f['y_fit'] / r['live_time'], 'r-', lw=1.5)
    ax.set_xlim(120, 240)
    ax.set_title(f"θ={r['label']}", fontsize=10)
    ax.set_xlabel('Canal', fontsize=8)
    ax.set_ylabel('cps', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(True, alpha=0.3)
for j in range(i+1, len(axes)):
    axes[j].set_visible(False)
fig.suptitle('$^{60}$Co — Espectros de coincidencia a diferentes ángulos (MCAch1)',
             fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, 'Co60_angular_multipanel.png'), dpi=150)
plt.close()
print('  -> Co60_angular_multipanel.png')

# ─── Save results ──────────────────────────────────────────────────────

out_data = []
for r in results:
    f = r['fit']
    out_data.append({
        'angle': r['angle'],
        'label': r['label'],
        'total_counts': r['total_counts'],
        'live_time': r['live_time'],
        'rate_cps': r['rate'],
        'area_1332': f['area'] if f else None,
        'area_norm_cps': r['area_norm'],
        'x0': f['x0'] if f else None,
        'sigma': f['sigma'] if f else None,
        'chi2_red': f['chi2_red'] if f else None,
    })

summary = {
    'A2': float(A2_fit) if 'A2_fit' in dir() else None,
    'A4': float(A4_fit) if 'A4_fit' in dir() else None,
    'data': out_data,
}

with open(os.path.join(OUT_DIR, 'resultados_angular_ch1.json'), 'w') as f:
    json.dump(summary, f, indent=2)
print('\n  Resultados guardados en resultados_angular_ch1.json')

print('\n' + '='*60)
print('  ANÁLISIS ANGULAR COMPLETADO')
print('='*60)
