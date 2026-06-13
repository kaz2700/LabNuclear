import os, sys, math, subprocess, json, re
import numpy as np

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(OUT_DIR), 'datos')
PLOTS_DIR = os.path.join(OUT_DIR, 'graficos_ch1')
os.makedirs(PLOTS_DIR, exist_ok=True)

sys.path.insert(0, OUT_DIR)
import _fit_methods as fm

# ─── Parse .itx ──────────────────────────────────────────────────────────────

def parse_itx(filepath):
    meta = {}
    counts = []
    with open(filepath, 'r', errors='replace') as f:
        lines = f.readlines()
    in_ch0 = False
    for line in lines:
        s = line.strip()
        if s.startswith('X //') or s.startswith('X ////'):
            parts = s.replace('X ////', '').replace('X //', '').strip().split()
            for p in parts:
                if '=' in p:
                    k, v = p.split('=', 1)
                    meta[k.strip()] = v.strip()
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
    if not os.path.exists(path):
        raise FileNotFoundError(f'{path} not found')
    meta, counts = parse_itx(path)
    channels = np.arange(len(counts), dtype=float)
    return meta, channels, counts

# ─── Run Fortran fit for a peak ─────────────────────────────────────────────

def run_fortran_fit(x, y, x0_approx, fortran_bin='bin/gaussian_background'):
    n = len(x)
    idx = int(np.argmin(np.abs(x - x0_approx)))
    target = 15
    shift = target - 1 - idx  # so reindexed[idx] = target
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

    # Parse Fortran output
    bg_line = [l for l in out.split('\n') if 'Background' in l][0]
    fit_line = [l for l in out.split('\n') if 'Fit:' in l][0]

    bg_m = re.search(r'Background: y =\s+([-\d.]+)\s+\+\s+([-\d.]+)', bg_line)
    a_bg = float(bg_m.group(1))
    b_bg = float(bg_m.group(2))

    nums = re.findall(r'[-]?[\d.]+', fit_line.replace('Fit:', ''))
    a_fit, b_fit, c_fit = map(float, nums[:3])

    # chi^2 from log-space fit
    chi2_line = [l for l in out.split('\n') if 'chi^2' in l and 'chi^2/ν' not in l]
    chi2_log = float(re.search(r'=\s*([\d.]+)', chi2_line[0]).group(1)) if chi2_line else 0

    # Parse Gaussian parameters from output
    params = {}
    for l in out.split('\n'):
        m = re.search(r'Peak position \(x0\)\s*=\s*([-\d.]+)', l)
        if m: params['x0_re'] = float(m.group(1))
        m = re.search(r'Width \(sigma\)\s+=\s*([-\d.]+)', l)
        if m: params['sigma_re'] = float(m.group(1))
        m = re.search(r'Amplitude \(A\)\s+=\s*([-\d.]+)', l)
        if m: params['amp'] = float(m.group(1))

    x0_re = params['x0_re']
    sigma_re = params['sigma_re']
    amp = params['amp']

    # Convert back to original channel coordinates
    x0 = x0_re - shift + x[0] - 1
    sigma = sigma_re  # sigma unchanged by translation

    # Background in original coords
    bg_x0 = a_bg + b_bg * (x0 - x[0] + 1 + shift)

    area = amp * sigma * math.sqrt(2 * math.pi)

    bg_line_arr = a_bg + b_bg * (x0_re - shift + x - x[0])
    y_fit_arr = amp * np.exp(-(x - x0)**2 / (2 * sigma**2)) + bg_line_arr
    resid_arr = y - y_fit_arr

    return {
        'x0': x0, 'sigma': sigma, 'A': amp,
        'area': area,
        'a_bg': a_bg, 'b_bg': b_bg, 'bg_x0': bg_x0,
        'chi2_log': chi2_log,
        'x_reindex': x_reindex,
        'shift': shift,
        'x0_re': x0_re, 'sigma_re': sigma_re,
        '_raw_out': out,
        'x_fit': x, 'y_fit': y_fit_arr, 'resid': resid_arr,
    }

# ─── Plot functions ──────────────────────────────────────────────────────────

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 11, 'figure.dpi': 150})

COLORS = {'data':'#1f77b4', 'fortran':'#d62728', 'python':'#2ca02c'}

def gauss(x, A, x0, s):
    return A * np.exp(-(x - x0)**2 / (2 * s**2))

def plot_peak_3methods(ch, counts, xmin, xmax, title, fname, fits):
    mask = (ch >= xmin) & (ch <= xmax)
    x = ch[mask]
    y = counts[mask]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8),
                                     gridspec_kw={'height_ratios': [3, 1]})
    fig.suptitle(title, fontsize=13, fontweight='bold')

    ax1.step(x, y, where='mid', color=COLORS['data'], label='Datos', linewidth=0.8)
    x_fine = np.linspace(xmin, xmax, 500)
    fm.plot_all_methods(ax1, x_fine, fits)
    ax1.legend(fontsize=9)
    ax1.set_ylabel('Cuentas')
    ax1.set_xlim(xmin, xmax)

    fm.make_resid_panel(ax2, ch, counts, fits, xmin, xmax)
    ax2.set_xlabel('Canal')

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, fname), dpi=150)
    plt.close()
    print(f'  -> {fname}')

# ─── Main ────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print('='*60)
    print('  ANÁLISIS CON CÓDIGO FORTRAN (src/)')
    print('='*60)

    fortran_bin = os.path.join(os.path.dirname(OUT_DIR), 'src', 'gaussian_background')
    if not os.path.exists(fortran_bin):
        print('ERROR: Compile Fortran first: cd src && gfortran -o gaussian_background gaussian_background.f subrutina.f')
        sys.exit(1)

    # ── Load all spectra ───────────────────────────────────────────────────
    meta_157, ch_157, counts_157 = cargar_espectro('22Na-157-coin')
    meta_180, ch_180, counts_180 = cargar_espectro('22Na-180-coin')
    meta_cof, ch_cof, counts_cof = cargar_espectro('60Co-F-coin')
    meta_ff, ch_ff, counts_ff = cargar_espectro('F-F-coin')

    # ── Fit all peaks with all 3 methods ──────────────────────────────────
    peaks = [
        ('511_157',  ch_157, counts_157, 66,  90,  78.05),
        ('1274_157', ch_157, counts_157, 188, 210, 196.45),
        ('511_180',  ch_180, counts_180, 66,  92,  78.39),
        ('pico_200', ch_180, counts_180, 190, 215, 196.81),
    ]

    results = {}
    scripts_fits = {}
    scipy_fits = {}
    for label, ch, counts, xmin, xmax, x0_guess in peaks:
        mask = (ch >= xmin) & (ch <= xmax)
        x = ch[mask]
        y = counts[mask]
        fort = run_fortran_fit(x, y, x0_guess, fortran_bin)
        results[label] = fort
        n_pts = len(x)
        print(f'\n{label}:')
        print(f'  FORTRAN: x₀ = {fort["x0"]:.2f}  σ = {fort["sigma"]:.2f}  A = {fort["A"]:.1f}')
        # Scripts-style
        scripts_fits[label] = fm.fit_scripts_style(ch, counts, x0_guess)
        if scripts_fits[label]:
            print(f'  SCRIPTS: x₀ = {scripts_fits[label]["x0"]:.2f}  σ = {scripts_fits[label]["sigma"]:.2f}')

    # Hardcoded scipy results from analisis_coincidencias_ch1.py
    scipy_hardcoded = {
        '511_157':  {'x0': 78.22, 'sigma': 2.41, 'A': 702.4, 'area': 4245.3, 'bg': 35.6, 'c1': 0},
        '1274_157': {'x0': 196.11, 'sigma': 4.30, 'A': 98.1, 'area': 1056.9, 'bg': 0.0, 'c1': 0},
        '511_180':  {'x0': 78.39, 'sigma': 2.43, 'A': 8889.7, 'area': 54247.7, 'bg': 304.5, 'c1': -3.26},
        'pico_200': {'x0': 195.95, 'sigma': 3.88, 'A': 38.1, 'area': 370.9, 'bg': 1.3, 'c1': 0},
    }
    for k, v in scipy_hardcoded.items():
        v['method'] = 'scipy'
        v['chi2_red'] = 0
        v['fwhm'] = 2.3548 * v['sigma']
        v['dstgma'] = 0
        v['bg_order'] = 'const' if v.get('c1', 0) == 0 else 'linear'
    scipy_fits = scipy_hardcoded

    # ── Statistics ──────────────────────────────────────────────────────────
    N_total_157 = int(counts_157.sum())
    N_511 = results['511_157']['area']
    N_1274 = results['1274_157']['area']
    N_511_1274 = N_511 + N_1274
    N_acc = N_total_157 - N_511_1274
    N_511_511 = results['511_180']['area']
    R = N_acc / N_511_511

    print(f'\n{"="*60}')
    print(f'  N_total(157°) = {N_total_157}')
    print(f'  N_511(157°)  = {N_511:.0f}')
    print(f'  N_1274(157°) = {N_1274:.0f}')
    print(f'  N_511+1274   = {N_511_1274:.0f}')
    print(f'  N_acc(157°)  = {N_acc:.0f}')
    print(f'  N_511+511(180°) = {N_511_511:.0f}')
    print(f'  R = {R:.4f}')
    print(f'{"="*60}')

    # Background analysis
    CoF_total = int(counts_cof.sum())
    FF_total = int(counts_ff.sum())
    t_cof = float(meta_cof.get('run_time', 22607))
    t_ff = float(meta_ff.get('run_time', 148117))
    r_cof = CoF_total / t_cof
    r_ff = FF_total / t_ff
    contrib_fondo = r_ff / r_cof * 100

    print(f'\nFondo:')
    print(f'  60Co+F: {CoF_total} ctos, {t_cof:.0f}s, {r_cof:.4f} cps')
    print(f'  F+F:   {FF_total} ctos, {t_ff:.0f}s, {r_ff:.4f} cps')
    print(f'  Contrib. fondo: {contrib_fondo:.1f}%')

    # ── Generate plots ────────────────────────────────────────────────────
    print(f'\nGenerando gráficos...')

    # 1) 22Na-157 completo
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.step(ch_157, counts_157, where='mid', color=COLORS['data'], linewidth=0.6)
    ax.set_xlabel('Canal'); ax.set_ylabel('Cuentas')
    ax.set_title('$^{22}$Na a 157° — Espectro completo')
    ax.set_xlim(0, 250)
    ax.text(78, 600, '511 keV', ha='center', fontsize=10,
            bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
    ax.text(196, 80, '1274 keV', ha='center', fontsize=10,
            bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))

    ax_inset = ax.inset_axes([0.55, 0.55, 0.35, 0.35])
    ax_inset.step(ch_157, counts_157, where='mid', color=COLORS['data'], linewidth=0.6)
    ax_inset.set_xlim(170, 230); ax_inset.set_ylim(0, 90)
    ax_inset.set_title('1274 keV', fontsize=9)
    ax.indicate_inset_zoom(ax_inset, edgecolor='gray')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, '22Na_157_completo.png'), dpi=150); plt.close()
    print('  -> 22Na_157_completo.png')

    # 2-4) Peak fits with all 3 methods
    def build_fits(label):
        d = {}
        if label in results:
            d['fortran'] = results[label]
        if label in scipy_fits:
            d['scipy'] = scipy_fits[label]
        if label in scripts_fits:
            d['scripts'] = scripts_fits[label]
        return d

    plot_peak_3methods(ch_157, counts_157, 66, 90,
              '$^{22}$Na a 157° — Pico 511 keV (comparación de métodos)',
              '22Na_157_pico_511_1274.png', build_fits('511_157'))
    plot_peak_3methods(ch_157, counts_157, 188, 210,
              '$^{22}$Na a 157° — Pico 1274 keV (comparación de métodos)',
              '22Na_157_pico_1274.png', build_fits('1274_157'))
    plot_peak_3methods(ch_180, counts_180, 66, 92,
              '$^{22}$Na a 180° — Pico 511+511 keV (comparación de métodos)',
              '22Na_180_pico_511_511.png', build_fits('511_180'))
    plot_peak_3methods(ch_180, counts_180, 190, 215,
              '$^{22}$Na a 180° — Pico residual ~510 keV (comparación de métodos)',
              '22Na_180_pico_fondo.png', build_fits('pico_200'))

    # 5) 22Na-180 completo
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.step(ch_180, counts_180, where='mid', color=COLORS['data'], linewidth=0.6)
    ax.set_xlabel('Canal'); ax.set_ylabel('Cuentas')
    ax.set_title('$^{22}$Na a 180° — Espectro completo')
    ax.set_xlim(0, 250)

    f511 = results['511_180']
    x_f = np.linspace(66, 92, 200)
    bg_f = f511['a_bg'] + f511['b_bg'] * (f511['x0_re'] - f511['shift'] + x_f - 66)
    ax.plot(x_f, gauss(x_f, f511['A'], f511['x0'], f511['sigma']) + bg_f, '-',
            color=COLORS['fortran'], linewidth=1.5, label='Ajuste Fortran')
    ax.text(78, 7200, '511+511 keV', ha='center', fontsize=10,
            bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, '22Na_180_completo.png'), dpi=150); plt.close()
    print('  -> 22Na_180_completo.png')

    # 6) Comparación 157 vs 180
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    ax1.step(ch_157, counts_157, where='mid', color=COLORS['data'], linewidth=0.6)
    ax1.set_title('$^{22}$Na a 157°', fontsize=11); ax1.set_ylabel('Cuentas')
    ax1.set_xlim(0, 250)
    ax2.step(ch_180, counts_180, where='mid', color=COLORS['data'], linewidth=0.6)
    ax2.set_title('$^{22}$Na a 180°', fontsize=11)
    ax2.set_xlabel('Canal'); ax2.set_ylabel('Cuentas')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, '22Na_comparacion.png'), dpi=150); plt.close()
    print('  -> 22Na_comparacion.png')

    # 7) Fondo
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7))
    ax1.step(ch_cof, counts_cof, where='mid', color='gray', linewidth=0.5, label='$^{60}$Co+F')
    ax1.step(ch_ff, counts_ff, where='mid', color='orange', linewidth=0.5, alpha=0.7, label='F+F')
    ax1.set_xlabel('Canal'); ax1.set_ylabel('Cuentas')
    ax1.set_title('Espectros de fondo')
    ax1.set_xlim(0, 300); ax1.legend(fontsize=9)

    ax2.bar(0, r_cof, width=0.3, color='gray', label=f'$^{60}$Co+F ({r_cof:.3f} cps)')
    ax2.bar(1, r_ff, width=0.3, color='orange', alpha=0.7, label=f'F+F ({r_ff:.3f} cps)')
    ax2.set_xticks([0, 1]); ax2.set_xticklabels(['$^{60}$Co+F', 'F+F'])
    ax2.set_ylabel('Tasa (cps)'); ax2.set_title('Tasas de coincidencia')
    ax2.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'CoF_vs_FF.png'), dpi=150); plt.close()
    print('  -> CoF_vs_FF.png')

    print('\nGráficos generados.')

    # ── Save results ──────────────────────────────────────────────────────
    out = {
        'metodo': 'Fortran (src/gaussian_background.f + subrutina.f)',
        'descripcion': 'Ajuste gaussiano con fondo lineal mediante transformación logarítmica y ajuste cuadrático',
        'fortran_results': {
            k: {'x0': v['x0'], 'sigma': v['sigma'], 'A': v['A'], 'area': v['area']}
            for k, v in results.items()
        },
        'estadisticas': {
            'N_total_157': N_total_157,
            'N_511_157': N_511,
            'N_1274_157': N_1274,
            'N_511_1274': N_511_1274,
            'N_acc_157': N_acc,
            'N_511_511_180': N_511_511,
            'R': R,
        },
        'fondo': {
            'CoF_total': CoF_total, 'FF_total': FF_total,
            'CoF_rate': r_cof, 'FF_rate': r_ff,
            'contrib_fondo_pct': contrib_fondo,
        }
    }
    with open(os.path.join(OUT_DIR, 'resultados_fortran_ch1.json'), 'w') as f:
        json.dump(out, f, indent=2)
    print('\nResultados guardados en resultados_fortran_ch1.json')
