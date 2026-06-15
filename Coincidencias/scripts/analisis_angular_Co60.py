#!/usr/bin/env python3
"""
Análisis de correlación angular del ⁶⁰Co en coincidencia.
Sigue instrucciones.md paso a paso para ambos picos (1173 y 1332 keV).

Pasos:
  1. Ajuste gaussiano (Fortran) de ambos picos → área, A, σ
  2. n(θ) = área / tiempo_vivo,  dn = d(área) / tiempo_vivo
  3. ε(θ) = n(θ) / n(90°)  con propagación de errores
"""

import os, sys, glob, math, json, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(BASE_DIR, 'datos')
PLOTS_DIR = os.path.join(BASE_DIR, 'graficos', 'ch1')
os.makedirs(PLOTS_DIR, exist_ok=True)

sys.path.insert(0, SCRIPT_DIR)
import _fit_methods as fm

FORTRAN_BIN = os.path.join(BASE_DIR, 'src', 'gaussian_background')

plt.rcParams.update({'font.size': 11, 'figure.dpi': 150})


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


def cargar_espectro(nombre):
    path = os.path.join(DATA_DIR, nombre)
    if not os.path.exists(path):
        path = os.path.join(DATA_DIR, nombre + '.itx')
    meta, counts = parse_itx_ch(path, 'MCAch1')
    channels = np.arange(len(counts), dtype=float)
    return meta, channels, counts


def fit_peak(ch, counts, est, fortran_bin=FORTRAN_BIN):
    window = 15
    lo, hi = est - window, est + window
    mask = (ch >= lo) & (ch <= hi)
    return fm.run_fortran_fit(ch[mask], counts[mask], est, fortran_bin)


def parse_angle(basename):
    m = re.search(r'([\d,]+)-coin', basename)
    if not m:
        m = re.search(r'(\d+)-coin', basename)
    if m:
        return float(m.group(1).replace(',', '.'))
    return None


# Peak definitions
PEAKS = {
    '1173': {'est': 180, 'label': '1173 keV'},
    '1332': {'est': 205, 'label': '1332 keV'},
}

print('=' * 60)
print('  CORRELACIÓN ANGULAR ⁶⁰Co — ambos picos')
print('=' * 60)

itx_files = sorted(glob.glob(os.path.join(DATA_DIR, '60*-coin.itx')))
itx_files = [f for f in itx_files if 'caracterizacion' not in f and 'Co-F' not in f]

all_results = {}
for pk in PEAKS:
    all_results[pk] = []

for fp in itx_files:
    name = os.path.basename(fp)
    base = os.path.splitext(name)[0]
    ang = parse_angle(base)
    if ang is None:
        print(f'  Saltando {name} (no se detectó ángulo)')
        continue

    meta, ch, counts = cargar_espectro(base)
    lt = meta.get('live_time_ch1', meta.get('run_time', 1))

    for pk, cfg in PEAKS.items():
        fit = fit_peak(ch, counts, cfg['est'])
        if fit is None:
            print(f'  {name:28s}  θ={ang:>5.1f}°  {pk:>4s}  AJUSTE FALLIDO')
            continue

        area = fit['area']
        darea = fit['darea']
        n = area / lt
        dn = darea / lt
        A = fit['A']
        dA = fit['dA']
        sigma = fit['sigma']
        dsigma = fit['dsigma']
        x0 = fit['x0']

        all_results[pk].append({
            'angle': ang,
            'label': f'{ang}°',
            'basename': base,
            'live_time': lt,
            'area': area, 'darea': darea,
            'n': n, 'dn': dn,
            'A': A, 'dA': dA,
            'sigma': sigma, 'dsigma': dsigma,
            'x0': x0,
        })

        print(f'  {name:28s}  θ={ang:>5.1f}°  {pk:>4s}  '
              f'área={area:.0f}±{darea:.0f}  '
              f'A={A:.0f}±{dA:.0f}  σ={sigma:.2f}±{dsigma:.2f}  '
              f'n={n:.4f}±{dn:.4f} cps')

print()

# Compute epsilon for each peak
for pk in PEAKS:
    all_results[pk].sort(key=lambda r: r['angle'])
    ref = [r for r in all_results[pk] if abs(r['angle'] - 90) < 1]
    if not ref:
        print(f'ERROR: No se encontró 90° para {pk}')
        continue
    n90 = ref[0]['n']
    dn90 = ref[0]['dn']
    for r in all_results[pk]:
        r['eps'] = (r['n'] - n90) / n90
        r['deps'] = math.sqrt((1 / n90) ** 2 * r['dn'] ** 2 +
                               (-r['n'] / n90 ** 2) ** 2 * dn90 ** 2)

    print(f'  --- {pk} keV ---')
    print(f'  Referencia 90°: n(90°) = {n90:.4f} ± {dn90:.4f} cps')
    angs = [r['angle'] for r in all_results[pk]]
    eps_v = [r['eps'] for r in all_results[pk]]
    deps_v = [r['deps'] for r in all_results[pk]]
    print(f'  {"θ":>7s}  {"n(cps)":>12s}  {"ε":>8s}  {"Δε":>8s}')
    for a, e, de in zip(angs, eps_v, deps_v):
        print(f'  {a:>5.1f}°  {all_results[pk][angs.index(a)]["n"]:.4f}  {e:.4f}  {de:.4f}')
    print()

# Plot combined
fig, ax = plt.subplots(figsize=(10, 6))
colors_pk = {'1173': '#2c7bb6', '1332': '#d62728'}
markers_pk = {'1173': 's', '1332': 'o'}

for pk in PEAKS:
    angs = np.array([r['angle'] for r in all_results[pk]])
    eps_v = np.array([r['eps'] for r in all_results[pk]])
    deps_v = np.array([r['deps'] for r in all_results[pk]])
    ax.errorbar(angs, eps_v, yerr=deps_v,
                fmt=markers_pk[pk], color=colors_pk[pk],
                capsize=4, ms=7, markerfacecolor='white',
                markeredgecolor=colors_pk[pk], markeredgewidth=1.5,
                label=rf'$\varepsilon_{{{pk}}}(\theta)$')

ax.axhline(0.0, color='gray', linestyle='--', alpha=0.4, label='Isotropía')
ax.set_xlabel(r'Ángulo $\theta$ (grados)', fontsize=12)
ax.set_ylabel(r'$\varepsilon(\theta) = [n(\theta)-n(90^\circ)]/n(90^\circ)$', fontsize=12)
ax.set_title(r'$^{60}$Co — Correlación angular $\gamma$-$\gamma$ (picos 1173 y 1332 keV)', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 190)
plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, 'Co60_anisotropia_ambos_picos.png'), dpi=150)
plt.close()
print('  -> Co60_anisotropia_ambos_picos.png')

# Save CSV and JSON
for pk in PEAKS:
    csv_path = os.path.join(BASE_DIR, 'resultados', f'anisotropia_Co60_{pk}.csv')
    with open(csv_path, 'w') as f:
        f.write('angulo,area,darea,n_cps,dn_cps,anisotropia,danisotropia,A,dA,sigma,dsigma,tiempo_vivo\n')
        for r in all_results[pk]:
            f.write(f'{r["angle"]},{r["area"]},{r["darea"]},{r["n"]},{r["dn"]},'
                    f'{r["eps"]},{r["deps"]},{r["A"]},{r["dA"]},'
                    f'{r["sigma"]},{r["dsigma"]},{r["live_time"]}\n')
    print(f'  -> {csv_path}')

out_data_pk = {}
for pk in PEAKS:
    out_data_pk[pk] = []
    for r in all_results[pk]:
        out_data_pk[pk].append({
            'angle': r['angle'],
            'area': r['area'], 'darea': r['darea'],
            'n_cps': r['n'], 'dn_cps': r['dn'],
            'epsilon': r['eps'], 'depsilon': r['deps'],
            'A': r['A'], 'dA': r['dA'],
            'sigma': r['sigma'], 'dsigma': r['dsigma'],
            'live_time': r['live_time'],
        })

summary = {
    'metodo': 'Fortran con propagación de errores (instrucciones.md)',
    'picos_analizados': ['1173', '1332'],
    'data': out_data_pk,
}

json_path = os.path.join(BASE_DIR, 'resultados', 'anisotropia_Co60.json')
with open(json_path, 'w') as f:
    json.dump(summary, f, indent=2)
print(f'  -> {json_path}')

print('\n' + '=' * 60)
print('  ANÁLISIS COMPLETADO')
print('=' * 60)
