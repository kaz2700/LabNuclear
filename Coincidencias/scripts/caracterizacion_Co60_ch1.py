import os, sys, math, json
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
    meta, counts = parse_itx(path)
    channels = np.arange(len(counts), dtype=float)
    return meta, channels, counts

def gauss_const(x, A, x0, s, bg):
    return A * np.exp(-(x - x0)**2 / (2 * s**2)) + bg

print('='*60)
print('  CARACTERIZACIÓN DEL DETECTOR CON ⁶⁰Co (MCAch1)')
print('='*60)

meta, ch, counts = cargar_espectro('60-Co-caracterización')
print(f'  Tiempo: {meta.get("run_time", "?")} s')
print(f'  Tasa: {meta.get("event_rate", "?")} cps')

print(f'  Canales del espectro: {len(counts)} (0-{len(counts)-1})')

# Smooth with numpy moving average and find peaks
window = np.ones(20) / 20.0
smoothed = np.convolve(counts, window, mode='same')

peaks = []
for i in range(30, len(counts)-30):
    if smoothed[i] > 50 and smoothed[i] == smoothed[i-10:i+11].max():
        refined = np.sum(ch[i-2:i+3] * counts[i-2:i+3]) / max(np.sum(counts[i-2:i+3]), 1)
        peaks.append((refined, counts[i], smoothed[i]))

peaks.sort(key=lambda x: -x[2])
print(f'\n  Picos principales:')
for pos, raw, sm in peaks[:10]:
    print(f'    canal {pos:.1f}: {int(raw)} cuentas')

# The two highest peaks are the 1173 and 1332
if len(peaks) >= 2:
    top_peaks = sorted([p for p in peaks if p[0] > 50], key=lambda x: -x[1])
    c1, c2 = int(round(top_peaks[0][0])), int(round(top_peaks[1][0]))
    ch1, ch2 = min(c1, c2), max(c1, c2)
    print(f'\n  Pico 1 (1173 keV): canal ~{ch1}')
    print(f'  Pico 2 (1332 keV): canal ~{ch2}')

    all_fits = {}
    for label, E, ch_est in [('1173', 1173.2, ch1), ('1332', 1332.5, ch2)]:
        dw = 15
        mask = (ch >= ch_est-dw) & (ch <= ch_est+dw)
        fort_fit = fm.run_fortran_fit(ch[mask], counts[mask], ch_est, FORTRAN_BIN)
        all_fits[label] = {'fortran': fort_fit} if fort_fit else {}
        if fort_fit:
            print(f'\n  --- Pico {label} keV ---')
            print(f'    x₀ = {fort_fit["x0"]:.2f}')
            print(f'    σ = {fort_fit["sigma"]:.2f}')
            print(f'    FWHM = {fort_fit["fwhm"]:.2f} canales')
            print(f'    A = {fort_fit["A"]:.0f}')
            print(f'    Fondo = {fort_fit["bg"]:.1f}')

            if label == '1173':
                x0_1, E1 = fort_fit['x0'], E
            else:
                x0_2, E2 = fort_fit['x0'], E

    # Energy calibration using Fortran centroids
    b_cal = (E2 - E1) / (x0_2 - x0_1)
    a_cal = E1 - b_cal * x0_1
    print(f'\n  --- CALIBRACIÓN ENERGÍA-CANAL ---')
    print(f'    E(keV) = {a_cal:.4f} + {b_cal:.4f} × canal')
    print(f'    Factor: {b_cal:.4f} keV/canal')

    # Compton edge and back-scattering
    for E_gamma, label_g in [(1173.2, '1173'), (1332.5, '1332')]:
        E_c = E_gamma / (1 + 2*E_gamma/511)
        ch_c = (E_c - a_cal) / b_cal
        print(f'\n  --- Borde Compton ({label_g} keV) ---')
        print(f'    E_Compton = {E_c:.1f} keV → canal ~{ch_c:.0f}')

    ch_bs_lo = (180 - a_cal) / b_cal
    ch_bs_hi = (260 - a_cal) / b_cal
    print(f'\n  --- Back-scattering (~200-250 keV) ---')
    print(f'    Región: canales {ch_bs_lo:.0f} - {ch_bs_hi:.0f}')

    # Full spectrum overview plot
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.plot(ch, counts, color='#1f77b4', linewidth=0.4)
    ax.set_xlabel('Canal')
    ax.set_ylabel('Cuentas')
    ax.set_title('$^{60}$Co — Espectro completo (modo singles)')
    ax.set_xlim(0, 500)
    peak_info = [
        (x0_1, 1173.2, '1173 keV'),
        (x0_2, 1332.5, '1332 keV'),
    ]
    for x0, E, label in peak_info:
        y_max = counts[int(round(x0))]
        ax.axvline(x0, color='red', linestyle='--', alpha=0.5, linewidth=0.8)
        ax.text(x0, y_max*1.05, label, rotation=90, fontsize=8,
                color='red', ha='right', va='bottom')
    for E_gamma, label_g in [(1173.2, '1173'), (1332.5, '1332')]:
        E_c = E_gamma / (1 + 2*E_gamma/511)
        ch_c = (E_c - a_cal) / b_cal
        ax.axvline(ch_c, color='green', linestyle=':', alpha=0.5, linewidth=0.8)
        ax.text(ch_c, ax.get_ylim()[1]*0.7, f'Compton\n{label_g}', fontsize=7,
                color='green', ha='center', rotation=90)
    ch_bs_lo = (180 - a_cal) / b_cal
    ch_bs_hi = (260 - a_cal) / b_cal
    ax.axvspan(ch_bs_lo, ch_bs_hi, alpha=0.08, color='orange')
    ax.text((ch_bs_lo+ch_bs_hi)/2, ax.get_ylim()[1]*0.5, 'Back-\nscattering',
            fontsize=7, color='orange', ha='center', va='center')
    cal_text = f'$E = {a_cal:.2f} + {b_cal:.2f}\\times$canal\n({b_cal:.2f} keV/canal)'
    ax.text(0.97, 0.95, cal_text, transform=ax.transAxes, fontsize=9,
            va='top', ha='right', bbox=dict(boxstyle='round,pad=0.3',
            facecolor='lightyellow', alpha=0.8))
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'Co60_caracterizacion_completo.png'), dpi=150)
    plt.close()
    print('  -> Co60_caracterizacion_completo.png')

    # Plot fitted peaks
    peak_conf = [
        ('1173 keV', ch1, 15, '1173'),
        ('1332 keV', ch2, 15, '1332'),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    for i, (label, ch_est, dw, key) in enumerate(peak_conf):
        ax1 = axes[i, 0]
        ax2 = axes[i, 1]
        mask = (ch >= ch_est-dw) & (ch <= ch_est+dw)
        x_p = ch[mask]
        y_p = counts[mask]
        ax1.step(x_p, y_p, where='mid', color='#1f77b4', linewidth=0.8, label='Datos')
        xf = np.linspace(ch_est-dw, ch_est+dw, 500)
        fits_dict = all_fits.get(key, {})
        if 'fortran' in fits_dict and fits_dict['fortran'] is not None:
            fm.plot_fortran_fit(ax1, xf, fits_dict['fortran'])
        ax1.set_ylabel('Cuentas')
        ax1.set_title(label)
        ax1.legend(fontsize=8)
        fm.make_resid_panel(ax2, ch, counts, fits_dict.get('fortran'), ch_est-dw, ch_est+dw)
        if i == len(peak_conf) - 1:
            ax2.set_xlabel('Canal')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'Co60_picos_caracterizacion.png'), dpi=150)
    plt.close()
    print('  -> Co60_picos_caracterizacion.png')

    # Save calibration results
    cal = {
        'metodo': 'Fortran (src/gaussian_background.f)',
        'a_keV': a_cal, 'b_keV_per_ch': b_cal,
        'E1_keV': E1, 'ch1': x0_1,
        'E2_keV': E2, 'ch2': x0_2,
    }
    with open(os.path.join(BASE_DIR, 'resultados', 'calibracion_ch1.json'), 'w') as f:
        json.dump(cal, f, indent=2)
    print('\n  Calibración guardada en calibracion_ch1.json')

    # Also plot with energy axis
    E_axis = a_cal + b_cal * ch
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(E_axis, counts, color='#1f77b4', linewidth=0.5)
    ax.set_xlabel('Energía (keV)')
    ax.set_ylabel('Cuentas')
    ax.set_title('$^{60}$Co — Espectro con calibración en energía')
    ax.set_xlim(0, 3000)
    for E_line, label in [(1173.2, '1173'), (1332.5, '1332')]:
        ax.axvline(E_line, color='red', linestyle='--', alpha=0.4)
        ax.text(E_line, ax.get_ylim()[1]*0.9, label, rotation=90, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'Co60_calibrado.png'), dpi=150)
    plt.close()
    print('  -> Co60_calibrado.png')

print('\n' + '='*60)
print('  CARACTERIZACIÓN COMPLETADA')
print('='*60)
