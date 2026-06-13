#!/usr/bin/env python3
"""
Shared fitting methods for parte2 analysis.
Contains all 3 methods: scipy curve_fit, scripts-style log-parabola, Fortran.
"""

import math, os, subprocess, re, json
import numpy as np
from scipy.optimize import curve_fit

# ─── Common models ────────────────────────────────────────────────────────

def gauss_const(x, A, x0, s, bg):
    return A * np.exp(-(x - x0)**2 / (2 * s**2)) + bg

def gauss_linear(x, A, x0, s, c0, c1):
    return A * np.exp(-(x - x0)**2 / (2 * s**2)) + c0 + c1 * x

def gauss(x, A, x0, s):
    return A * np.exp(-(x - x0)**2 / (2 * s**2))


# ─── Method B: scipy curve_fit ────────────────────────────────────────────

def fit_scipy(channels, counts, xmin, xmax, bg_order='const',
              sigma_lim=(0.5, 12.0), debug=False):
    mask = (channels >= xmin) & (channels <= xmax)
    x = channels[mask]
    y = counts[mask]
    n_pts = len(x)
    if n_pts < 6:
        return None

    i_max = np.argmax(y)
    x0_guess = x[i_max]
    A_guess = y[i_max] - np.median(y)

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
            print(f"  fit_scipy failed [{xmin},{xmax}]: {e}")
        return None

    y_fit = model(x, *popt)
    resid = y - y_fit
    chi2 = np.sum((resid / np.maximum(np.sqrt(y), 0.5))**2)
    ndf = n_pts - n_params
    chi2_red = chi2 / ndf if ndf > 0 else 999
    ss_res = np.sum(resid**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

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

    if sigma < sigma_lim[0] or sigma > sigma_lim[1]:
        if debug:
            print(f"  fit_scipy: sigma={sigma:.2f} out of bounds, rejecting")
        return None
    if abs(x0_f - x0_guess) > (xmax - xmin) * 0.6:
        return None

    return {
        'method': 'scipy',
        'xmin': xmin, 'xmax': xmax,
        'A': A, 'dA': dA,
        'x0': x0_f, 'dx0': dx0,
        'sigma': sigma, 'dsigma': dsigma,
        'fwhm': fwhm, 'dfwhm': dfwhm,
        'bg': bg, 'dbg': dbg,
        'c1': c1, 'dc1': dc1,
        'area': area, 'darea': d_area,
        'chi2': chi2, 'ndf': ndf, 'chi2_red': chi2_red,
        'r_squared': r2,
        'bg_order': bg_order, 'n_params': n_params,
        'x_fit': x, 'y_fit': y_fit,
        'resid': resid,
        'popt': popt, 'pcov': pcov,
    }


# ─── Method A: scripts-style log-parabola ─────────────────────────────────

def fit_scripts_style(channels, counts, pk_guess, debug=False):
    x = np.asarray(channels, dtype=float)
    y = np.asarray(counts, dtype=float)
    pk_idx = int(np.argmin(np.abs(x - pk_guess)))
    pk = x[pk_idx]

    best_fit = None
    best_score = None

    for hw in [6, 7, 8, 9, 10, 12]:
        for bg_offset in [10, 12, 14, 16, 18, 20]:
            bg_left = list(range(max(0, pk_idx-bg_offset), max(0, pk_idx-bg_offset+6)))
            bg_right = list(range(pk_idx+8, min(len(x), pk_idx+8+bg_offset-2)))
            bg_indices = bg_left + bg_right
            if len(bg_indices) < 6:
                continue

            bg_x = x[bg_indices]
            bg_y = y[bg_indices]
            n = len(bg_x)
            sx = np.sum(bg_x); sy = np.sum(bg_y)
            sxx = np.sum(bg_x*bg_x); sxy = np.sum(bg_x*bg_y)
            denom = n*sxx - sx*sx
            if abs(denom) < 1e-30:
                continue
            b_bg = (n*sxy - sx*sy) / denom
            a_bg = (sy - b_bg*sx) / n

            flo = pk_idx - hw
            fhi = pk_idx + hw
            x_net_list, y_net_list = [], []
            for i in range(flo, fhi+1):
                if i < 0 or i >= len(x):
                    continue
                bg = a_bg + b_bg * x[i]
                net = y[i] - bg
                if net > 0:
                    x_net_list.append(float(x[i]))
                    y_net_list.append(net)
            if len(x_net_list) < 5:
                continue

            x_net = np.array(x_net_list)
            y_net = np.array(y_net_list)
            y_log = np.log(y_net)
            sig_log = 1.0 / np.sqrt(np.maximum(y_net, 0.01))

            w = 1.0 / (sig_log**2)
            s0 = np.sum(w); s1 = np.sum(w*x_net); s2 = np.sum(w*x_net**2)
            s3 = np.sum(w*x_net**3); s4 = np.sum(w*x_net**4)
            t0 = np.sum(w*y_log); t1 = np.sum(w*x_net*y_log)
            t2 = np.sum(w*x_net**2*y_log)

            det = s0*(s2*s4 - s3*s3) - s1*(s1*s4 - s2*s3) + s2*(s1*s3 - s2*s2)
            if abs(det) < 1e-40:
                continue

            a00 = (s2*s4 - s3*s3)/det; a01 = -(s1*s4 - s2*s3)/det; a02 = (s1*s3 - s2*s2)/det
            a11 = (s0*s4 - s2*s2)/det; a12 = -(s0*s3 - s1*s2)/det; a22 = (s0*s2 - s1*s1)/det
            a0 = a00*t0 + a01*t1 + a02*t2
            a1 = a01*t0 + a11*t1 + a12*t2
            a2 = a02*t0 + a12*t1 + a22*t2
            if a2 >= 0:
                continue

            sigma = math.sqrt(-1.0/(2.0*a2))
            x0 = a1 * sigma * sigma
            ampl = math.exp(a0 + x0*x0/(2.0*sigma*sigma))
            if sigma < 1 or sigma > 12 or abs(x0 - pk) > 6:
                continue

            # Compute RMS of residuals
            sum_sq = 0; n_pts = 0
            for i in range(flo, fhi+1):
                if i < 0 or i >= len(x):
                    continue
                bg = a_bg + b_bg*x[i]
                g = ampl * math.exp(-(x[i]-x0)**2/(2.0*sigma*sigma))
                res = y[i] - (bg + g)
                sum_sq += res*res; n_pts += 1
            rms_v = math.sqrt(sum_sq/n_pts) if n_pts > 0 else 0

            chi2_parab = np.sum(((y_log - (a0 + a1*x_net + a2*x_net**2))/sig_log)**2)
            ndf = len(x_net) - 3
            score = abs(x0 - pk) * 3 + chi2_parab/ndf + (rms_v/ampl)*10

            if best_fit is None or score < best_fit[0]:
                best_fit = (score, hw, bg_offset, x0, sigma, ampl, rms_v,
                            chi2_parab, ndf, a_bg, b_bg, pk, flo, fhi,
                            x_net, y_net)

    if best_fit is None:
        return None

    (sc, hw, bg_off, x0, sigma, ampl, rms_v,
     chi2_parab, ndf, a_bg, b_bg, pk, flo, fhi,
     x_net, y_net) = best_fit

    fwhm = 2.3548 * sigma
    area = ampl * sigma * math.sqrt(2 * math.pi)

    # Build full x_fit range, y_fit
    x_fit_arr = x[max(0, flo):min(len(x), fhi+1)]
    y_fit_arr = a_bg + b_bg*x_fit_arr + ampl * np.exp(-(x_fit_arr-x0)**2/(2*sigma**2))
    resid_arr = y[max(0, flo):min(len(x), fhi+1)] - y_fit_arr

    return {
        'method': 'scripts',
        'xmin': x[max(0, flo)], 'xmax': x[min(len(x)-1, fhi)],
        'A': ampl, 'dA': 0,
        'x0': x0, 'dx0': 0,
        'sigma': sigma, 'dsigma': 0,
        'fwhm': fwhm, 'dfwhm': 0,
        'bg': a_bg + b_bg*x0, 'dbg': 0,
        'c1': b_bg, 'dc1': 0,
        'area': area, 'darea': 0,
        'chi2': chi2_parab, 'ndf': ndf, 'chi2_red': chi2_parab/ndf if ndf > 0 else 0,
        'r_squared': 0,
        'bg_order': 'linear',
        'n_params': 2 + 3,
        'x_fit': x_fit_arr, 'y_fit': y_fit_arr,
        'resid': resid_arr,
        'a_bg': a_bg, 'b_bg': b_bg,
    }


# ─── Method C: Fortran ────────────────────────────────────────────────────

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

    # Build fit curve
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


# ─── Plot helper ──────────────────────────────────────────────────────────

METHOD_STYLE = {
    'scipy':   {'color': '#1f77b4', 'linestyle': '--', 'lw': 1.5, 'label': 'Scipy (curve_fit)'},
    'scripts': {'color': '#2ca02c', 'linestyle': ':',  'lw': 1.5, 'label': 'Scripts (log-parabola)'},
    'fortran': {'color': '#d62728', 'linestyle': '-',  'lw': 1.5, 'label': 'Fortran (log+LFIT)'},
}

def plot_all_methods(ax, x_fine, fits, show_params=False):
    for method_key in ['fortran', 'scipy', 'scripts']:
        if method_key not in fits or fits[method_key] is None:
            continue
        f = fits[method_key]
        style = METHOD_STYLE[method_key]
        if method_key == 'scipy':
            bg = f['bg'] + f.get('c1', 0) * x_fine
        elif method_key == 'scripts':
            bg = f.get('a_bg', f.get('bg', 0)) + f.get('b_bg', 0) * x_fine
        elif method_key == 'fortran':
            sf = f.get('shift', 0)
            x0_re = f.get('x0_re', f.get('x0', 0))
            bg = f.get('a_bg', f.get('bg', 0)) + f.get('b_bg', 0) * (x0_re - sf + x_fine - x_fine[0])
        y_gauss = gauss(x_fine, f['A'], f['x0'], f['sigma'])
        ax.plot(x_fine, y_gauss + bg, color=style['color'],
                linestyle=style['linestyle'], linewidth=style['lw'],
                label=style['label'])

def add_params_box(ax, fit, x=0.03, y=0.55, fs=9):
    if fit is None:
        return
    parts = []
    if fit.get('dx0', 0) > 0:
        parts.append(f"$x_0 = {fit['x0']:.2f}\\pm{fit['dx0']:.2f}$")
    else:
        parts.append(f"$x_0 = {fit['x0']:.2f}$")
    parts.append(f"$\\sigma = {fit['sigma']:.2f}$")
    parts.append(f"$\\mathrm{{FWHM}} = {fit['fwhm']:.2f}$")
    parts.append(f"$A = {fit['A']:.0f}$")
    parts.append(f"$\\mathrm{{Area}} = {fit['area']:.0f}$")
    if fit['chi2_red'] > 0:
        parts.append(f"$\\chi^2/\\nu = {fit['chi2_red']:.2f}$")
    txt = '\n'.join(parts)
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=fs,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))

def make_resid_panel(ax, x, y, fits, xmin, xmax):
    for method_key in ['fortran', 'scipy', 'scripts']:
        if method_key not in fits or fits[method_key] is None:
            continue
        f = fits[method_key]
        style = METHOD_STYLE[method_key]
        if 'x_fit' in f and f['x_fit'] is not None:
            mask = (x >= f['x_fit'][0]) & (x <= f['x_fit'][-1])
        else:
            mask = (x >= xmin) & (x <= xmax)
        x_local = x[mask]
        y_local = y[mask]
        # Compute residual at exact data points
        if method_key == 'scipy':
            bg = f['bg'] + f.get('c1', 0) * x_local
        elif method_key == 'scripts':
            bg = f.get('a_bg', f.get('bg', 0)) + f.get('b_bg', 0) * x_local
        elif method_key == 'fortran':
            sf = f.get('shift', 0)
            x0_re = f.get('x0_re', f.get('x0', 0))
            bg = f.get('a_bg', f.get('bg', 0)) + f.get('b_bg', 0) * (x0_re - sf + x_local - x_local[0])
        y_model = gauss(x_local, f['A'], f['x0'], f['sigma']) + bg
        resid = y_local - y_model
        ax.step(x_local, resid, where='mid', color=style['color'],
                linewidth=0.6, alpha=0.7)
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    ax.set_xlim(xmin, xmax)
    ax.set_ylabel('Residuo', fontsize=10)
    ax.legend([METHOD_STYLE[k]['label'] for k in ['fortran', 'scipy', 'scripts']
               if k in fits and fits[k] is not None],
              fontsize=7, ncol=3)
