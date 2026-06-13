import os, sys, math, json
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2 as chi2_dist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_DIR = os.path.abspath(os.path.join(OUT_DIR, '..'))
DATA_DIR = os.path.join(PROJ_DIR, 'datos')
PLOTS_DIR = os.path.join(OUT_DIR, 'graficos_ch1')
os.makedirs(PLOTS_DIR, exist_ok=True)

sys.path.insert(0, OUT_DIR)
import _fit_methods as fm
FORTRAN_BIN = os.path.join(PROJ_DIR, 'src', 'gaussian_background')

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

def fit_peak(ch, counts, xmin, xmax):
    mask = (ch >= xmin) & (ch <= xmax)
    x = ch[mask]
    y = counts[mask]
    i_max = np.argmax(y)
    A0 = y[i_max] - np.median(y)
    x00 = x[i_max]
    s0 = 3.0
    bg0 = np.median(np.concatenate([y[:3], y[-3:]]))
    try:
        popt, pcov = curve_fit(gauss_const, x, y,
                               p0=[A0, x00, s0, bg0],
                               bounds=([0, xmin-5, 0.5, 0], [np.inf, xmax+5, 20, np.inf]),
                               sigma=np.where(y > 1, np.sqrt(y), 1.0), maxfev=10000)
        perr = np.sqrt(np.diag(pcov))
        chi2 = np.sum(((y - gauss_const(x, *popt)) / np.maximum(np.sqrt(y), 0.5))**2)
        ndf = len(x) - 4
        return popt, perr, chi2, ndf
    except:
        return None, None, None, None

print('='*60)
print('  CARACTERIZACIÓN DEL DETECTOR CON ⁶⁰Co')
print('='*60)

meta, ch, counts = cargar_espectro('60-Co-caracterización')
print(f'  Tiempo: {meta.get("run_time", "?")} s')
print(f'  Tasa: {meta.get("event_rate", "?")} cps')

from scipy.ndimage import uniform_filter1d
print(f'  Canales del espectro: {len(counts)} (0-{len(counts)-1})')

# Smooth and find all significant peaks
smoothed = uniform_filter1d(counts, size=20)
peaks = []
for i in range(30, len(counts)-30):
    if smoothed[i] > 50 and smoothed[i] == smoothed[i-10:i+11].max():
        refined = np.sum(ch[i-2:i+3] * counts[i-2:i+3]) / max(np.sum(counts[i-2:i+3]), 1)
        peaks.append((refined, counts[i], smoothed[i]))

peaks.sort(key=lambda x: -x[2])
print(f'\n  Picos principales:')
for pos, raw, sm in peaks[:10]:
    print(f'    canal {pos:.1f}: {int(raw)} cuentas')

# The two highest peaks are the 1173 and 1332 (channels ~181 and ~206)
if len(peaks) >= 2:
    # Take the two highest peaks by raw counts
    top_peaks = sorted([p for p in peaks if p[0] > 50], key=lambda x: -x[1])
    c1, c2 = int(round(top_peaks[0][0])), int(round(top_peaks[1][0]))
    ch1, ch2 = min(c1, c2), max(c1, c2)  # lower channel = lower energy = 1173
    print(f'\n  Pico 1 (1173 keV): canal ~{ch1}')
    print(f'  Pico 2 (1332 keV): canal ~{ch2}')

    # Fit windows: narrow to avoid peak overlap (peaks separated by ~52 channels)
    all_fits = {}
    for label, E, ch_est in [('1173', 1173.2, ch1), ('1332', 1332.5, ch2)]:
        dw = 15
        popt, perr, chi2, ndf = fit_peak(ch, counts, ch_est-dw, ch_est+dw)
        if popt is not None:
            A, x0, s, bg = popt
            dA, dx0, ds, dbg = perr
            fwhm = 2.355 * s
            area = A * s * math.sqrt(2*math.pi)
            mask = (ch >= ch_est-dw) & (ch <= ch_est+dw)
            x_local = ch[mask]
            y_local = counts[mask]
            y_fit_local = gauss_const(x_local, *popt)
            scipy_fit = {
                'method': 'scipy', 'x0': x0, 'sigma': s, 'A': A, 'area': area,
                'bg': bg, 'c1': 0, 'fwhm': fwhm, 'chi2_red': chi2/ndf if ndf > 0 else 0,
                'x_fit': x_local, 'y_fit': y_fit_local, 'resid': y_local - y_fit_local,
            }
            all_fits[label] = {'scipy': scipy_fit}
            print(f'\n  --- Pico {label} keV ---')
            print(f'    x₀ = {x0:.1f} ± {dx0:.1f}')
            print(f'    σ = {s:.2f} ± {ds:.2f}')
            print(f'    FWHM = {fwhm:.2f} canales')
            print(f'    A = {A:.0f} ± {dA:.0f}')
            print(f'    Fondo = {bg:.1f} ± {dbg:.1f}')
            print(f'    χ²/ν = {chi2/ndf:.2f}, ν = {ndf}')
            
            # Scripts-style and Fortran
            scripts_fit = fm.fit_scripts_style(ch, counts, ch_est)
            if scripts_fit:
                all_fits[label]['scripts'] = scripts_fit
                print(f'    Scripts: x₀={scripts_fit["x0"]:.2f} σ={scripts_fit["sigma"]:.2f}')
            mask = (ch >= ch_est-dw) & (ch <= ch_est+dw)
            fort_fit = fm.run_fortran_fit(ch[mask], counts[mask], ch_est, FORTRAN_BIN)
            if fort_fit:
                all_fits[label]['fortran'] = fort_fit
                print(f'    Fortran: x₀={fort_fit["x0"]:.2f} σ={fort_fit["sigma"]:.2f}')
            
            # Energy calibration: E[keV] = a + b * channel
            if label == '1173':
                x0_1, E1 = x0, E
            else:
                x0_2, E2 = x0, E
    
    # Energy calibration
    b_cal = (E2 - E1) / (x0_2 - x0_1)
    a_cal = E1 - b_cal * x0_1
    print(f'\n  --- CALIBRACIÓN ENERGÍA-CANAL ---')
    print(f'    E(keV) = {a_cal:.4f} + {b_cal:.4f} × canal')
    print(f'    Factor: {b_cal:.4f} keV/canal')
    
    # Find sum peak (1173+1332 = 2505.7 keV)
    x0_sum = (2505.7 - a_cal) / b_cal
    print(f'\n  --- PICO SUMA (1173+1332 = 2505.7 keV) ---')
    print(f'    Canal esperado: {x0_sum:.0f}')
    
    # Look for peak near x0_sum
    dw = 40
    x0_sum_int = int(round(x0_sum))
    popt_s, perr_s, chi2_s, ndf_s = fit_peak(ch, counts, x0_sum_int-dw, x0_sum_int+dw)
    all_fits['sum'] = {}
    if popt_s is not None:
        A, x0, s, bg = popt_s
        mask_s = (ch >= x0_sum_int-dw) & (ch <= x0_sum_int+dw)
        x_local_s = ch[mask_s]
        y_local_s = counts[mask_s]
        y_fit_s = gauss_const(x_local_s, *popt_s)
        all_fits['sum']['scipy'] = {
            'method': 'scipy', 'x0': x0, 'sigma': s, 'A': A, 'area': A * s * math.sqrt(2*math.pi),
            'bg': bg, 'c1': 0, 'fwhm': 2.355*s, 'chi2_red': chi2_s/ndf_s if ndf_s > 0 else 0,
            'x_fit': x_local_s, 'y_fit': y_fit_s, 'resid': y_local_s - y_fit_s,
        }
        print(f'    Pico encontrado en canal {x0:.1f}')
        print(f'    Energía: {a_cal + b_cal*x0:.1f} keV')
        print(f'    σ = {s:.2f}, FWHM = {2.355*s:.2f}')
        # Scripts-style and Fortran for sum peak
        scripts_s = fm.fit_scripts_style(ch, counts, x0_sum)
        if scripts_s:
            all_fits['sum']['scripts'] = scripts_s
            print(f'    Scripts: x₀={scripts_s["x0"]:.2f} σ={scripts_s["sigma"]:.2f}')
        fort_s = fm.run_fortran_fit(ch[mask_s], counts[mask_s], x0_sum, FORTRAN_BIN)
        if fort_s:
            all_fits['sum']['fortran'] = fort_s
            print(f'    Fortran: x₀={fort_s["x0"]:.2f} σ={fort_s["sigma"]:.2f}')
    else:
        print(f'    No se encontró pico suma visible')
    
    # Compton edge and back-scattering
    # For 1173 keV: Compton edge at E_c = E / (1 + 2E/mc²) where mc²=511 keV
    # E_c = 1173 / (1 + 2*1173/511) = 1173 / (1 + 4.59) = 1173/5.59 = 209.8 keV
    # Back-scattering peak at ~200-250 keV
    # Channel for Compton edge:
    for E_gamma, label_g in [(1173.2, '1173'), (1332.5, '1332')]:
        E_c = E_gamma / (1 + 2*E_gamma/511)
        ch_c = (E_c - a_cal) / b_cal
        print(f'\n  --- Borde Compton ({label_g} keV) ---')
        print(f'    E_Compton = {E_c:.1f} keV → canal ~{ch_c:.0f}')
    
    # Back-scattering: typically around 200-250 keV
    ch_bs_lo = (180 - a_cal) / b_cal
    ch_bs_hi = (260 - a_cal) / b_cal
    print(f'\n  --- Back-scattering (~200-250 keV) ---')
    print(f'    Región: canales {ch_bs_lo:.0f} - {ch_bs_hi:.0f}')
    
    # Plot fitted peaks: 3 columns (peaks) x 2 rows (fit + resid)
    peak_conf = [
        ('1173 keV', ch1, 15, '1173'),
        ('1332 keV', ch2, 15, '1332'),
        ('Pico suma 2505 keV', int(round(x0_sum)), 40, 'sum'),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(16, 14))
    for i, (label, ch_est, dw, key) in enumerate(peak_conf):
        ax1 = axes[i, 0]
        ax2 = axes[i, 1]
        mask = (ch >= ch_est-dw) & (ch <= ch_est+dw)
        x_p = ch[mask]
        y_p = counts[mask]
        ax1.step(x_p, y_p, where='mid', color='#1f77b4', linewidth=0.8, label='Datos')
        xf = np.linspace(ch_est-dw, ch_est+dw, 500)
        fits_dict = all_fits.get(key, {})
        fm.plot_all_methods(ax1, xf, fits_dict)
        ax1.set_ylabel('Cuentas')
        ax1.set_title(label)
        ax1.legend(fontsize=8)
        fm.make_resid_panel(ax2, ch, counts, fits_dict, ch_est-dw, ch_est+dw)
        if i == len(peak_conf) - 1:
            ax2.set_xlabel('Canal')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'Co60_picos_caracterizacion.png'), dpi=150)
    plt.savefig(os.path.join(PLOTS_DIR, 'Co60_caracterizacion_completo.png'), dpi=150)
    plt.close()
    print('  -> Co60_picos_caracterizacion.png')
    print('  -> Co60_caracterizacion_completo.png')
    
    # Save calibration results
    cal = {
        'a_keV': a_cal, 'b_keV_per_ch': b_cal,
        'E1_keV': E1, 'ch1': x0_1,
        'E2_keV': E2, 'ch2': x0_2,
        'pico_suma_canal': x0 if popt_s is not None else None,
        'pico_suma_keV': a_cal + b_cal*x0 if popt_s is not None else None,
    }
    with open(os.path.join(OUT_DIR, 'calibracion_ch1.json'), 'w') as f:
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
    for E_line, label in [(1173.2, '1173'), (1332.5, '1332'), (2505.7, 'Suma')]:
        ax.axvline(E_line, color='red', linestyle='--', alpha=0.4)
        ax.text(E_line, ax.get_ylim()[1]*0.9, label, rotation=90, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'Co60_calibrado.png'), dpi=150)
    plt.close()
    print('  -> Co60_calibrado.png')

print('\n' + '='*60)
print('  CARACTERIZACIÓN COMPLETADA')
print('='*60)
