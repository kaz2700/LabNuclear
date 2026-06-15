#!/usr/bin/env python3
"""
Genera los plots para el informe de correlaci\u00f3n angular del \u2076\u2070Co.
Analiza los picos de 1173 y 1332 keV.
"""

import os, sys, json, re, glob
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(BASE_DIR, 'datos')
SCRIPTS_DIR = os.path.join(BASE_DIR, 'scripts')

sys.path.insert(0, SCRIPTS_DIR)
import _fit_methods as fm

FORTRAN_BIN = os.path.join(BASE_DIR, 'src', 'gaussian_background')

plt.rcParams.update({
    'font.size': 11,
    'figure.dpi': 200,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
})


def parse_itx_ch(filepath, ch='MCAch1'):
    meta = {}
    counts = []
    with open(filepath, 'r', errors='replace') as f:
        lines = f.readlines()
    in_ch = False
    for line in lines:
        s = line.strip()
        if s.startswith('X //') or s.startswith('X ////'):
            if 'Run Time [s]' in s:
                meta['run_time'] = float(s.split('=')[1].strip().split()[0])
            if 'Channel 1' in s:
                meta['_last_ch'] = 'Channel 1'
            if 'Live Time [s]' in s and meta.get('_last_ch') == 'Channel 1':
                meta['live_time_ch1'] = float(s.split('=')[1].strip().split()[0])
        if s.startswith('WAVES') and ch in s:
            in_ch = True
            continue
        if in_ch:
            if s == 'BEGIN':
                continue
            if s == 'END':
                break
            try:
                counts.append(int(s))
            except ValueError:
                pass
    meta.pop('_last_ch', None)
    return meta, np.array(counts, dtype=float)


# ─── Peak configurations ──────────────────────────────────────────────

PEAKS = {
    '1173': {'est': 180, 'window': 15, 'label': '1173 keV'},
    '1332': {'est': 205, 'window': 15, 'label': '1332 keV'},
}

# ─── Load and fit all angles ──────────────────────────────────────────

itx_files = sorted(glob.glob(os.path.join(DATA_DIR, '60*-coin.itx')))
itx_files = [f for f in itx_files if 'caracterizacion' not in f and 'Co-F' not in f]

angles_list = []
peak_data = {'1173': [], '1332': []}

for fp in itx_files:
    name = os.path.basename(fp)
    base = os.path.splitext(name)[0]
    m = re.search(r'([\d,]+)-coin', base)
    if not m:
        continue
    ang = float(m.group(1).replace(',', '.'))
    meta, counts = parse_itx_ch(fp, 'MCAch1')
    ch = np.arange(len(counts), dtype=float)
    lt = meta.get('live_time_ch1', meta.get('run_time', 1))

    for key, cfg in PEAKS.items():
        lo = cfg['est'] - cfg['window']
        hi = cfg['est'] + cfg['window']
        mask = (ch >= lo) & (ch <= hi)
        fit = fm.run_fortran_fit(ch[mask], counts[mask], cfg['est'], FORTRAN_BIN)
        if fit is None:
            peak_data[key].append(None)
            continue
        area = fit['area']
        darea = fit['darea']
        n = area / lt
        dn = darea / lt
        peak_data[key].append({
            'angle': ang,
            'label': f'{ang}°',
            'lt': lt,
            'area': area, 'darea': darea,
            'n': n, 'dn': dn,
            'A': fit['A'], 'dA': fit['dA'],
            'sigma': fit['sigma'], 'dsigma': fit['dsigma'],
            'x0': fit['x0'],
        })

    if ang not in angles_list:
        angles_list.append(ang)

angles_list.sort()

# Remove angles missing any peak
for key in peak_data:
    peak_data[key] = [d for d in peak_data[key] if d is not None]

# Filter to common angles
common_angles = sorted(set(d['angle'] for d in peak_data['1173']) &
                       set(d['angle'] for d in peak_data['1332']))
for key in peak_data:
    peak_data[key] = [d for d in peak_data[key] if d['angle'] in common_angles]
    peak_data[key].sort(key=lambda r: r['angle'])


def compute_epsilon(data_list):
    angles_a = np.array([r['angle'] for r in data_list])
    n_vals = np.array([r['n'] for r in data_list])
    dn_vals = np.array([r['dn'] for r in data_list])
    ref = next(r for r in data_list if abs(r['angle'] - 90) < 1)
    n90, dn90 = ref['n'], ref['dn']
    eps = (n_vals - n90) / n90
    deps = np.sqrt((1 / n90) ** 2 * dn_vals ** 2 +
                   (-n_vals / n90 ** 2) ** 2 * dn90 ** 2)
    return angles_a, n_vals, dn_vals, eps, deps, n90, dn90

eps_data = {}
for key in PEAKS:
    eps_data[key] = compute_epsilon(peak_data[key])

# ─── Solid-angle correction and Legendre fit ─────────────────────────

R = 0.5        # source-detector distance (m)
diam = 81.8e-3 # detector diameter (m)
r_det = diam / 2
beta0 = np.arctan(r_det / R)  # half-angle

Q2 = np.cos(beta0) * (1 + np.cos(beta0)) / 2
Q4 = np.cos(beta0) * (1 + np.cos(beta0)) * (7 * np.cos(beta0)**2 - 3) / 8

def P2(x):
    return (3 * x**2 - 1) / 2

def P4(x):
    return (35 * x**4 - 30 * x**2 + 3) / 8

def eps_model(theta, a2, a4):
    ct = np.cos(np.radians(theta))
    num = 1 + a2 * P2(ct) + a4 * P4(ct)
    den = 1 + a2 * P2(0) + a4 * P4(0)
    return num / den - 1

fit_results = {}
theta_grid = np.linspace(0, 180, 361)
for key in PEAKS:
    angles_a, _, _, eps, deps, _, _ = eps_data[key]
    popt, pcov = curve_fit(eps_model, angles_a, eps, sigma=deps,
                           p0=[0.1, 0.05], absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    a2_obs, a4_obs = popt
    da2_obs, da4_obs = perr
    a2_corr = a2_obs / Q2**2
    a4_corr = a4_obs / Q4**2
    da2_corr = da2_obs / Q2**2
    da4_corr = da4_obs / Q4**2
    # eps_fit is ε(θ) = W(θ)/W(90) - 1
    eps_fit_vals = eps_model(theta_grid, a2_obs, a4_obs)
    fit_results[key] = {
        'a2_obs': a2_obs, 'da2_obs': da2_obs,
        'a4_obs': a4_obs, 'da4_obs': da4_obs,
        'a2_corr': a2_corr, 'da2_corr': da2_corr,
        'a4_corr': a4_corr, 'da4_corr': da4_corr,
        'eps_fit': eps_fit_vals,
    }

# ─── Solid-angle angular uncertainty Δε_ang ──────────────────────────

R_mm = 500.0
r_det_mm = 40.9
dtheta = r_det_mm / R_mm  # angular uncertainty from finite detector size

for key in PEAKS:
    ct = np.cos(np.radians(angles_a))
    st = np.sin(np.radians(angles_a))
    deps_dtheta = 0.25 * ct * st + (1/6) * ct**3 * st
    deps_ang = np.abs(deps_dtheta) * dtheta
    for i, r in enumerate(peak_data[key]):
        r['deps_ang'] = deps_ang[i]

# Save correction results for LaTeX
corr_lines = []
corr_lines.append(f"R = {R} m  d = {diam*1e3:.1f} mm")
corr_lines.append(f"beta0 = {np.degrees(beta0):.3f} deg")
corr_lines.append(f"Q2 = {Q2:.6f}  Q4 = {Q4:.6f}")
corr_lines.append(f"Omega = {np.pi*r_det**2/R**2:.6f} sr")
corr_lines.append(f"dtheta = {dtheta:.6f} rad = {np.degrees(dtheta):.4f} deg")
corr_lines.append("")
corr_lines.append("Pico  a2_obs  da2_obs  a4_obs  da4_obs  a2_corr  da2_corr  a4_corr  da4_corr")
for key in PEAKS:
    r = fit_results[key]
    corr_lines.append(
        f"{key}  {r['a2_obs']:.4f} {r['da2_obs']:.4f} {r['a4_obs']:.4f} {r['da4_obs']:.4f} "
        f"{r['a2_corr']:.4f} {r['da2_corr']:.4f} {r['a4_corr']:.4f} {r['da4_corr']:.4f}")
corr_lines.append("")
corr_lines.append("Delta epsilon_ang (angular uncertainty) for each angle:")
for key in PEAKS:
    corr_lines.append(f"\n{key} keV:")
    corr_lines.append("theta  deps_ang  deps_areas  deps_total")
    _, _, _, _, deps_arr, _, _ = eps_data[key]
    for i, r in enumerate(peak_data[key]):
        deps_areas = deps_arr[i]
        deps_ang = r.get('deps_ang', 0)
        deps_tot = np.sqrt(deps_areas**2 + deps_ang**2)
        corr_lines.append(f"{r['angle']:.1f}  {deps_ang:.6f}  {deps_areas:.6f}  {deps_tot:.6f}")
with open(os.path.join(SCRIPT_DIR, 'correccion_angulo_solido.txt'), 'w') as f:
    f.write('\n'.join(corr_lines))

# ─── Figure 0: Individual spectra per angle (multi-panel) ─────────────

colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(common_angles)))

fig1, ax1 = plt.subplots(figsize=(8, 5))
for i, ang in enumerate(common_angles):
    d = next(r for r in peak_data['1332'] if r['angle'] == ang)
    # Use full counts from any angle's first file
    fp = next(f for f in itx_files if f'60-Co-{ang}'.replace('.0', '').replace('.', ',') in os.path.basename(f) or f'60Co-{ang}'.replace('.0', '') in os.path.basename(f))
    meta, counts = parse_itx_ch(fp, 'MCAch1')
    ch = np.arange(len(counts), dtype=float)
    rate = counts / d['lt']
    ax1.step(ch, rate, where='mid', color=colors[i], linewidth=0.5,
             label=rf'$\theta = {ang}\degree$')
ax1.set_xlim(100, 250)
ax1.set_xlabel('Canal', fontsize=12)
ax1.set_ylabel('Tasa de conteo (cps)', fontsize=12)
ax1.set_title(r'Espectros de coincidencia $^{60}$Co a diferentes \'angulos',
              fontsize=13)
ax1.legend(fontsize=7, ncol=2, framealpha=0.8)
ax1.grid(True, alpha=0.25, linewidth=0.3)
# Mark peak positions
ax1.axvline(180, color='gray', linestyle=':', alpha=0.3, linewidth=0.6)
ax1.axvline(205, color='gray', linestyle=':', alpha=0.3, linewidth=0.6)
ax1.text(180, ax1.get_ylim()[1] * 0.95, '1173', fontsize=7, color='gray',
         ha='center', rotation=90)
ax1.text(205, ax1.get_ylim()[1] * 0.95, '1332', fontsize=7, color='gray',
         ha='center', rotation=90)
plt.tight_layout()
fig1.savefig(os.path.join(SCRIPT_DIR, 'espectros_superpuestos.pdf'))
plt.close()


# ─── Figure 2: Anisotropy ε(θ) — both peaks ───────────────────────────

colors_pk = {'1173': '#2c7bb6', '1332': '#d62728'}
markers_pk = {'1173': 's', '1332': 'o'}

fig2, ax2 = plt.subplots(figsize=(7.5, 5.5))

for key in PEAKS:
    angles_a, n_vals, dn_vals, eps, deps, n90, dn90 = eps_data[key]
    label = rf'$\varepsilon_{{{key}}}(\theta)$'
    ax2.errorbar(angles_a, eps, yerr=deps,
                 fmt=markers_pk[key], color=colors_pk[key],
                 capsize=4, capthick=1.2,
                 ms=7, markerfacecolor='white',
                 markeredgecolor=colors_pk[key], markeredgewidth=1.5,
                 label=label)

for key in PEAKS:
    ax2.plot(theta_grid, fit_results[key]['eps_fit'],
             color=colors_pk[key], linestyle='-', linewidth=0.7, alpha=0.6)

ax2.axhline(0.0, color='gray', linestyle='--', alpha=0.4, linewidth=0.8,
            label='Isotrop\u00eda')

ax2.set_xlabel(r'\'Angulo $\theta$ (grados)', fontsize=12)
ax2.set_ylabel(r'$\varepsilon(\theta) = [n(\theta)-n(90^\circ)]/n(90^\circ)$',
               fontsize=11)
ax2.set_title(r'$^{60}$Co — Correlaci\u00f3n angular $\gamma$-$\gamma$',
              fontsize=13)
ax2.legend(fontsize=10, framealpha=0.8)
ax2.grid(True, alpha=0.25, linewidth=0.3)
ax2.set_xlim(0, 190)
ax2.xaxis.set_major_locator(mticker.MultipleLocator(30))
ax2.xaxis.set_minor_locator(mticker.MultipleLocator(15))

plt.tight_layout()
fig2.savefig(os.path.join(SCRIPT_DIR, 'anisotropia.pdf'))
plt.close()

# ─── Figure 0: Individual spectra per angle (multi-panel) ─────────────

n_angles = len(common_angles)
ncols = 2
nrows = (n_angles + ncols - 1) // ncols

fig0, axes0 = plt.subplots(nrows, ncols, figsize=(7, 2.2 * nrows), sharex=True, sharey=False, constrained_layout=True)
axes0 = axes0.flatten()

for i, ang in enumerate(common_angles):
    d1332 = next(r for r in peak_data['1332'] if r['angle'] == ang)
    fp = next(f for f in itx_files if
              f'60-Co-{ang}'.replace('.0', '').replace('.', ',') in os.path.basename(f)
              or f'60Co-{ang}'.replace('.0', '') in os.path.basename(f))
    meta, counts = parse_itx_ch(fp, 'MCAch1')
    ch = np.arange(len(counts), dtype=float)
    rate = counts / d1332['lt']
    ax = axes0[i]
    ax.step(ch, rate, where='mid', color='#2c3e50', linewidth=0.4)
    ax.axvline(180, color='gray', linestyle=':', alpha=0.4, linewidth=0.5)
    ax.axvline(205, color='gray', linestyle=':', alpha=0.4, linewidth=0.5)
    ax.text(180, ax.get_ylim()[1] * 0.92, '1173', fontsize=6, color='gray',
            ha='center', rotation=90)
    ax.text(205, ax.get_ylim()[1] * 0.92, '1332', fontsize=6, color='gray',
            ha='center', rotation=90)
    ax.set_title(f'$\\theta = {ang}\\degree$', fontsize=10, fontweight='bold')
    ax.set_xlim(0, 512)
    ax.grid(True, alpha=0.15, linewidth=0.3)
    ax.tick_params(labelsize=7)
    if i >= n_angles - ncols:
        ax.set_xlabel('Canal', fontsize=8)
    if i % ncols == 0:
        ax.set_ylabel('Tasa (cps)', fontsize=8)

for j in range(i + 1, len(axes0)):
    axes0[j].set_visible(False)

fig0.savefig(os.path.join(SCRIPT_DIR, 'espectros_individuales.pdf'))
plt.close()

# ─── Figure 3: MCAch1 vs MCAch0 comparison (multi-panel) ──────────

n_angles = len(common_angles)
ncols_c = 2
nrows_c = (n_angles + ncols_c - 1) // ncols_c

fig3, axes3 = plt.subplots(nrows_c, ncols_c, figsize=(8, 2.8 * nrows_c), sharex=True, sharey=False, constrained_layout=True)
axes3 = axes3.flatten()

for i, ang in enumerate(common_angles):
    d = next(r for r in peak_data['1332'] if r['angle'] == ang)
    fp = next(f for f in itx_files if
              f'60-Co-{ang}'.replace('.0', '').replace('.', ',') in os.path.basename(f)
              or f'60Co-{ang}'.replace('.0', '') in os.path.basename(f))
    meta0, counts0 = parse_itx_ch(fp, 'MCAch0')
    meta1, counts1 = parse_itx_ch(fp, 'MCAch1')
    lt0 = meta0.get('live_time_ch1', meta0.get('run_time', 1))
    lt1 = meta1.get('live_time_ch1', meta1.get('run_time', 1))
    ch_arr = np.arange(len(counts0), dtype=float)
    rate0 = counts0 / lt0
    rate1 = counts1 / lt1
    ax = axes3[i]
    ax.step(ch_arr, rate0, where='mid', color='#1f77b4', linewidth=0.35, alpha=0.8, label='MCAch0 (fijo)')
    ax.step(ch_arr, rate1, where='mid', color='#d62728', linewidth=0.35, alpha=0.8, label='MCAch1 (movil)')
    ax.axvline(180, color='gray', linestyle=':', alpha=0.3, linewidth=0.4)
    ax.axvline(205, color='gray', linestyle=':', alpha=0.3, linewidth=0.4)
    ax.set_title(f'$\\theta = {ang}\\degree$', fontsize=10, fontweight='bold')
    ax.set_xlim(0, 512)
    ax.grid(True, alpha=0.15, linewidth=0.3)
    ax.tick_params(labelsize=7)
    if i >= n_angles - ncols_c:
        ax.set_xlabel('Canal', fontsize=8)
    if i % ncols_c == 0:
        ax.set_ylabel('Tasa (cps)', fontsize=8)
    if i == 0:
        ax.legend(fontsize=6, framealpha=0.8, loc='upper right')

for j in range(i + 1, len(axes3)):
    axes3[j].set_visible(False)

fig3.savefig(os.path.join(SCRIPT_DIR, 'comparacion_fijo_movil.pdf'))
plt.close()

print('Plots generados: espectros_individuales.pdf, espectros_superpuestos.pdf, anisotropia.pdf, comparacion_fijo_movil.pdf')

# ─── Save processed data for LaTeX ─────────────────────────────────────

lines = []
header = ('angulo  '
          'A1173 dA1173 sig1173 dsig1173 n1173 dn1173 eps1173 deps1173 deps1173_ang deps1173_tot '
          'A1332 dA1332 sig1332 dsig1332 n1332 dn1332 eps1332 deps1332 deps1332_ang deps1332_tot '
          't_vivo')
lines.append(header)
for i, ang in enumerate(common_angles):
    d3 = peak_data['1332'][i]
    d1 = peak_data['1173'][i]
    _, _, _, e3, de3, _, _ = eps_data['1332']
    _, _, _, e1, de1, _, _ = eps_data['1173']
    de3_ang = d3.get('deps_ang', 0)
    de1_ang = d1.get('deps_ang', 0)
    de3_tot = np.sqrt(de3[i]**2 + de3_ang**2)
    de1_tot = np.sqrt(de1[i]**2 + de1_ang**2)
    lines.append(
        f'{ang:.1f}  '
        f'{d1["A"]:.1f} {d1["dA"]:.1f} {d1["sigma"]:.4f} {d1["dsigma"]:.4f} '
        f'{d1["n"]:.4f} {d1["dn"]:.4f} {e1[i]:.4f} {de1[i]:.4f} {de1_ang:.6f} {de1_tot:.6f}  '
        f'{d3["A"]:.1f} {d3["dA"]:.1f} {d3["sigma"]:.4f} {d3["dsigma"]:.4f} '
        f'{d3["n"]:.4f} {d3["dn"]:.4f} {e3[i]:.4f} {de3[i]:.4f} {de3_ang:.6f} {de3_tot:.6f}  '
        f'{d3["lt"]:.0f}')

with open(os.path.join(SCRIPT_DIR, 'tabla_datos.txt'), 'w') as f:
    f.write('\n'.join(lines))
print('Datos guardados: tabla_datos.txt')
