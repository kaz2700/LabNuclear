#!/usr/bin/env python3
"""
Shared fitting methods for parte2 analysis.
Fortran-only: uses gaussian_background binary via subprocess.
"""

import math, os, subprocess, re, json
import numpy as np

# ─── Common models ────────────────────────────────────────────────────────

def gauss_const(x, A, x0, s, bg):
    return A * np.exp(-(x - x0)**2 / (2 * s**2)) + bg

def gauss_linear(x, A, x0, s, c0, c1):
    return A * np.exp(-(x - x0)**2 / (2 * s**2)) + c0 + c1 * x

def gauss(x, A, x0, s):
    return A * np.exp(-(x - x0)**2 / (2 * s**2))


# ─── Fortran fit ──────────────────────────────────────────────────────────

def run_fortran_fit(x, y, x0_approx, fortran_bin='bin/gaussian_background'):
    n = len(x)
    idx = int(np.argmin(np.abs(x - x0_approx)))
    target = 15
    shift = target - 1 - idx
    x_reindex = np.arange(1, n + 1) + shift

    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
    os.makedirs(data_dir, exist_ok=True)
    pico_path = os.path.join(data_dir, 'pico.dat')

    with open(pico_path, 'w') as f:
        for xi, yi in zip(x_reindex, y):
            f.write(f'{int(xi)} {int(yi)}\n')

    PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result = subprocess.run([fortran_bin], capture_output=True, text=True, cwd=PROJECT_ROOT)
    out = result.stdout

    try:
        bg_line = [l for l in out.split('\n') if 'Background' in l][0]
        fit_line = [l for l in out.split('\n') if 'Fit:' in l][0]
    except IndexError:
        return None

    bg_m = re.search(r'Background: y =\s+([-\d.]+)\s+\+\s+([-\d.]+)', bg_line)
    if not bg_m:
        return None
    a_bg_f = float(bg_m.group(1))
    b_bg_f = float(bg_m.group(2))

    nums = re.findall(r'[-]?[\d.]+', fit_line.replace('Fit:', ''))
    if len(nums) < 3:
        return None
    a_fit, b_fit, c_fit = map(float, nums[:3])

    params = {}
    for l in out.split('\n'):
        m = re.search(r'Peak position \(x0\)\s*=\s*([-\d.]+)', l)
        if m: params['x0_re'] = float(m.group(1))
        m = re.search(r'Width \(sigma\)\s+=\s*([-\d.]+)', l)
        if m: params['sigma_re'] = float(m.group(1))
        m = re.search(r'Amplitude \(A\)\s+=\s*([-\d.]+)', l)
        if m: params['amp'] = float(m.group(1))

    if 'x0_re' not in params:
        return None

    x0_re = params['x0_re']
    sigma_re = params['sigma_re']
    amp = params['amp']
    x0 = x0_re - shift + x[0] - 1
    sigma = sigma_re
    area = amp * sigma * math.sqrt(2 * math.pi)
    fwhm = 2.3548 * sigma

    x_fit_arr = x
    bg_line_arr = a_bg_f + b_bg_f * (x0_re - shift + x - x[0])
    y_fit_arr = amp * np.exp(-(x - x0)**2 / (2 * sigma**2)) + bg_line_arr
    resid_arr = y - y_fit_arr

    return {
        'method': 'fortran',
        'xmin': x[0], 'xmax': x[-1],
        'A': amp, 'dA': 0,
        'x0': x0, 'dx0': 0,
        'sigma': sigma, 'dsigma': 0,
        'fwhm': fwhm, 'dfwhm': 0,
        'bg': bg_line_arr[int(len(bg_line_arr)/2)], 'dbg': 0,
        'c1': b_bg_f, 'dc1': 0,
        'area': area, 'darea': 0,
        'chi2': 0, 'ndf': 0, 'chi2_red': 0,
        'r_squared': 0,
        'bg_order': 'linear',
        'n_params': 5,
        'x_fit': x_fit_arr, 'y_fit': y_fit_arr,
        'resid': resid_arr,
        'a_bg': a_bg_f, 'b_bg': b_bg_f,
        'x0_re': x0_re, 'sigma_re': sigma_re,
        'shift': shift,
        '_raw_out': out,
    }


# ─── Plot helpers ─────────────────────────────────────────────────────────

METHOD_STYLE = {
    'fortran': {'color': '#d62728', 'linestyle': '-', 'lw': 1.5, 'label': 'Fortran (log+LFIT)'},
}

def plot_fortran_fit(ax, x_fine, fit):
    if fit is None:
        return
    sf = fit.get('shift', 0)
    x0_re = fit.get('x0_re', fit.get('x0', 0))
    bg = fit.get('a_bg', fit.get('bg', 0)) + fit.get('b_bg', 0) * (x0_re - sf + x_fine - x_fine[0])
    y_gauss = gauss(x_fine, fit['A'], fit['x0'], fit['sigma'])
    ax.plot(x_fine, y_gauss + bg, color='#d62728', linestyle='-', linewidth=1.5,
            label='Fortran (log+LFIT)')

def add_params_box(ax, fit, x=0.03, y=0.55, fs=9):
    if fit is None:
        return
    parts = []
    parts.append(f"$x_0 = {fit['x0']:.2f}$")
    parts.append(f"$\\sigma = {fit['sigma']:.2f}$")
    parts.append(f"$\\mathrm{{FWHM}} = {fit['fwhm']:.2f}$")
    parts.append(f"$A = {fit['A']:.0f}$")
    parts.append(f"$\\mathrm{{Area}} = {fit['area']:.0f}$")
    txt = '\n'.join(parts)
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=fs,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))

def make_resid_panel(ax, x, y, fit, xmin, xmax):
    if fit is None:
        return
    f = fit
    if 'x_fit' in f and f['x_fit'] is not None:
        mask = (x >= f['x_fit'][0]) & (x <= f['x_fit'][-1])
    else:
        mask = (x >= xmin) & (x <= xmax)
    x_local = x[mask]
    y_local = y[mask]
    sf = f.get('shift', 0)
    x0_re = f.get('x0_re', f.get('x0', 0))
    bg = f.get('a_bg', f.get('bg', 0)) + f.get('b_bg', 0) * (x0_re - sf + x_local - x_local[0])
    y_model = gauss(x_local, f['A'], f['x0'], f['sigma']) + bg
    resid = y_local - y_model
    ax.step(x_local, resid, where='mid', color='#d62728', linewidth=0.6, alpha=0.7)
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    ax.set_xlim(xmin, xmax)
    ax.set_ylabel('Residuo', fontsize=10)
