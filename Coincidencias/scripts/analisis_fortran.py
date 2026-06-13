import os, sys, math, subprocess, json, re
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(BASE_DIR, 'datos')
PLOTS_DIR = os.path.join(BASE_DIR, 'graficos', 'ch0')
os.makedirs(PLOTS_DIR, exist_ok=True)

sys.path.insert(0, SCRIPT_DIR)
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
        if s.startswith('WAVES') and 'MCAch0' in s:
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

# ─── Plot functions ──────────────────────────────────────────────────────────

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 11, 'figure.dpi': 150})

COLORS = {'data':'#1f77b4', 'fortran':'#d62728'}

def gauss(x, A, x0, s):
    return A * np.exp(-(x - x0)**2 / (2 * s**2))

def plot_peak_fortran(ch, counts, xmin, xmax, title, fname, fort_fit):
    mask = (ch >= xmin) & (ch <= xmax)
    x = ch[mask]
    y = counts[mask]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8),
                                     gridspec_kw={'height_ratios': [3, 1]})
    fig.suptitle(title, fontsize=13, fontweight='bold')

    ax1.step(x, y, where='mid', color=COLORS['data'], label='Datos', linewidth=0.8)
    x_fine = np.linspace(xmin, xmax, 500)
    if fort_fit is not None:
        fm.plot_fortran_fit(ax1, x_fine, fort_fit)
    ax1.legend(fontsize=9)
    ax1.set_ylabel('Cuentas')
    ax1.set_xlim(xmin, xmax)

    fm.make_resid_panel(ax2, ch, counts, fort_fit, xmin, xmax)
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

    fortran_bin = os.path.join(BASE_DIR, 'src', 'gaussian_background')
    if not os.path.exists(fortran_bin):
        print('ERROR: Compile Fortran first: cd src && gfortran -o gaussian_background gaussian_background.f subrutina.f')
        sys.exit(1)

    # ── Load all spectra ───────────────────────────────────────────────────
    meta_157, ch_157, counts_157 = cargar_espectro('22Na-157-coin')
    meta_180, ch_180, counts_180 = cargar_espectro('22Na-180-coin')
    meta_cof, ch_cof, counts_cof = cargar_espectro('60Co-F-coin')
    meta_ff, ch_ff, counts_ff = cargar_espectro('F-F-coin')

    # ── Fit all peaks with Fortran ─────────────────────────────────────────
    peaks = [
        ('511_157',  ch_157, counts_157, 66,  90,  78.05),
        ('1274_157', ch_157, counts_157, 188, 210, 196.45),
        ('511_180',  ch_180, counts_180, 66,  92,  78.39),
        ('pico_200', ch_180, counts_180, 190, 215, 196.81),
    ]

    results = {}
    for label, ch, counts, xmin, xmax, x0_guess in peaks:
        mask = (ch >= xmin) & (ch <= xmax)
        x = ch[mask]
        y = counts[mask]
        fort = fm.run_fortran_fit(x, y, x0_guess, fortran_bin)
        results[label] = fort
        print(f'\n{label}:')
        if fort:
            print(f'  FORTRAN: x₀ = {fort["x0"]:.2f}  σ = {fort["sigma"]:.2f}  A = {fort["A"]:.1f}')

    # ── Statistics ──────────────────────────────────────────────────────────
    N_total_157 = int(counts_157.sum())
    N_511 = results['511_157']['area']
    N_1274 = results['1274_157']['area']
    N_acc = N_511 - N_1274
    N_511_511 = results['511_180']['area']
    R = N_acc / N_511_511

    print(f'\n{"="*60}')
    print(f'  N_total(157°) = {N_total_157}')
    print(f'  N_511(157°)  = {N_511:.0f}')
    print(f'  N_1274(157°) = {N_1274:.0f}')
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

    # 2-4) Peak fits with Fortran
    plot_peak_fortran(ch_157, counts_157, 66, 90,
              '$^{22}$Na a 157° — Pico 511 keV (ajuste Fortran)',
              '22Na_157_pico_511_1274.png', results['511_157'])
    plot_peak_fortran(ch_157, counts_157, 188, 210,
              '$^{22}$Na a 157° — Pico 1274 keV (ajuste Fortran)',
              '22Na_157_pico_1274.png', results['1274_157'])
    plot_peak_fortran(ch_180, counts_180, 66, 92,
              '$^{22}$Na a 180° — Pico 511+511 keV (ajuste Fortran)',
              '22Na_180_pico_511_511.png', results['511_180'])
    plot_peak_fortran(ch_180, counts_180, 190, 215,
              '$^{22}$Na a 180° — Pico residual ~510 keV (ajuste Fortran)',
              '22Na_180_pico_fondo.png', results['pico_200'])

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
    with open(os.path.join(BASE_DIR, 'resultados', 'resultados_fortran.json'), 'w') as f:
        json.dump(out, f, indent=2)
    print('\nResultados guardados en resultados_fortran.json')
