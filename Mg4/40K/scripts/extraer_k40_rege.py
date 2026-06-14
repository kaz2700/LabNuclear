#!/usr/bin/env python3
"""
Extraer el pico de 40K del espectro de fondo del detector REGe
(fondo_GeHP_REGe.dat) usando ajuste gaussiano + fondo lineal.

Formato del .dat:
  Linea 1: tiempo vivo (en ticks de 20 ms)
  Linea 2: numero de canales
  Lineas 3+: cuentas por canal (una por linea)
"""

import math
import json
import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATOS = ROOT / "datos"
OUTDIR = HERE.parent
DATA_DIR = OUTDIR / "data"
PLOTS_DIR = OUTDIR / "plots"
FORTRAN_BIN = ROOT.parent / "src" / "gaussian_background"

SPECTRUM_PATH = DATOS / "fondo_GeHP_REGe.dat"
SPECTRUM_TXT = DATA_DIR / "fondo_rege.txt"

K40_REF_KEV = 1460.82


def parse_rege_dat(path):
    lines = path.read_text().splitlines()
    vals = [int(l.strip()) for l in lines]
    live_time_ticks = vals[0]
    num_channels = vals[1]
    counts = vals[2:]
    if len(counts) != num_channels:
        raise ValueError(
            f"Expected {num_channels} channels, got {len(counts)}"
        )
    return live_time_ticks, counts


def save_spectrum_txt(counts, path):
    with path.open("w") as f:
        f.write("#canal\tcuentas\n")
        for i, c in enumerate(counts, 1):
            f.write(f"{i}\t{c}\n")


def find_peaks(counts, threshold_sigma=3.0, window_peak=7, window_bg=101):
    nch = len(counts)
    peaks = []
    for i in range(window_peak, nch - window_peak):
        c = counts[i]
        lo = max(0, i - window_bg // 2)
        hi = min(nch, i + window_bg // 2 + 1)
        bg_data = counts[lo:i - window_peak] + counts[i + window_peak + 1:hi]
        if not bg_data:
            continue
        bg_mean = sum(bg_data) / len(bg_data)
        bg_var = max(0, sum((x - bg_mean) ** 2 for x in bg_data) / len(bg_data))
        bg_std = math.sqrt(bg_var)

        is_local_max = True
        for j in range(i - window_peak, i + window_peak + 1):
            if j != i and counts[j] >= c:
                is_local_max = False
                break

        if not is_local_max:
            continue

        excess = c - bg_mean
        significance = excess / max(bg_std, 1.0)
        if significance > threshold_sigma and c >= 5:
            peaks.append((i + 1, c, bg_mean, excess, significance))
    peaks.sort(key=lambda x: x[4], reverse=True)
    return peaks


def linear_calibration(anchors):
    x = [a[0] for a in anchors]
    y = [a[1] for a in anchors]
    n = len(x)
    sx = sum(x)
    sy = sum(y)
    sxx = sum(xi * xi for xi in x)
    sxy = sum(xi * yi for xi, yi in zip(x, y))
    denom = n * sxx - sx * sx
    if denom == 0:
        return None, None
    b = (n * sxy - sx * sy) / denom
    a = (sy - b * sx) / n
    return a, b


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    # --- Parse ---
    live_ticks, counts = parse_rege_dat(SPECTRUM_PATH)
    live_s = live_ticks * 0.020
    nch = len(counts)

    # Save 2-column txt for gnuplot
    save_spectrum_txt(counts, SPECTRUM_TXT)

    print(f"Espectro: {SPECTRUM_PATH.name}")
    print(f"  Canales: {nch}")
    print(f"  Live ticks: {live_ticks}  ->  {live_s:.1f} s ({live_s/3600:.3f} h)")
    print(f"  Cuentas totales: {sum(counts):,}")
    print(f"  Tasa total: {sum(counts)/live_s:.1f} cps")
    print()

    # --- Find peaks for calibration ---
    peaks = find_peaks(counts, threshold_sigma=4.0)
    print(f"Picos detectados (significancia > 4 sigma):")
    for ch, cnt, bg, excess, sig in peaks[:20]:
        print(f"  Ch {ch:5d}  cnt={cnt:>6}  bg={bg:>7.1f}  excess={excess:>+6.0f}  sig={sig:.1f}")
    print()

    # --- Calibration using known background lines ---
    # The two most prominent peaks:
    #   ch 2160 -> 609.31 keV (Bi-214)
    #   ch 6179 -> 1764.49 keV (Bi-214)
    cal_anchors = [(2160, 609.31), (6179, 1764.49)]
    a_kev, b_kev_ch = linear_calibration(cal_anchors)
    if a_kev is None:
        print("ERROR: calibracion fallo")
        return

    print(f"Calibracion lineal (2 puntos):")
    print(f"  E(keV) = {a_kev:.4f} + {b_kev_ch:.6f} * ch")
    print(f"  Anclas: ch 2160 -> 609.31 keV (Bi-214), ch 6179 -> 1764.49 keV (Bi-214)")
    print()

    # Verify calibration with other visible peaks
    print("Verificacion de calibracion:")
    verify_peaks = [
        (788, "Pb-212", 238.63),
        (1010, "Pb-214", 295.22),
        (1187, "Pb-214", 351.93),
        (1491, "varios", 417),
        (3854, "Bi-214", 1120.29),
        (5078, "candidato", 1448),
        (6179, "Bi-214", 1764.49),
    ]
    for ch, name, expected in verify_peaks:
        e_calc = a_kev + b_kev_ch * ch
        print(f"  Ch {ch:5d}  {name:>10s}  E_calc={e_calc:.2f}  E_ref={expected:.2f}  delta={e_calc - expected:+.2f}")

    print()

    # --- 40K region ---
    ch_k40_est = (K40_REF_KEV - a_kev) / b_kev_ch
    print(f"40K: energia de referencia = {K40_REF_KEV} keV")
    print(f"     canal estimado = {ch_k40_est:.1f}")

    # Find local maximum in the region around estimated position
    search_lo = max(0, int(ch_k40_est) - 100)
    search_hi = min(nch, int(ch_k40_est) + 100)
    search_region = counts[search_lo:search_hi]
    max_region_val = max(search_region)
    max_region_idx = search_region.index(max_region_val) + search_lo + 1
    e_at_max = a_kev + b_kev_ch * max_region_idx
    print(f"     maximo local en region: ch {max_region_idx} ({max_region_val} cnt, {e_at_max:.2f} keV)")

    # Decide seed channel for fit
    if abs(max_region_idx - ch_k40_est) < 30:
        seed_ch = max_region_idx
    else:
        seed_ch = round(ch_k40_est)
    print(f"     semilla para ajuste: ch {seed_ch}")

    # ROI for fit: 27 channels for Fortran program
    fit_half = 13
    fit_start = seed_ch - fit_half
    fit_end = seed_ch + fit_half

    # Fit channels list
    fit_channels_raw = list(range(max(1, fit_start), min(nch + 1, fit_end + 1)))
    fit_counts_raw = [counts[i - 1] for i in fit_channels_raw]
    n_fit = len(fit_channels_raw)

    # Ensure we have exactly 27 channels
    if n_fit > 27:
        excess = n_fit - 27
        trim = excess // 2
        fit_channels = fit_channels_raw[trim:trim + 27]
        fit_counts = fit_counts_raw[trim:trim + 27]
    elif n_fit < 27:
        pad_before = (27 - n_fit) // 2
        pad_after = 27 - n_fit - pad_before
        ch_first = fit_channels_raw[0]
        ch_last = fit_channels_raw[-1]
        for _ in range(pad_before):
            ch_first -= 1
            if ch_first >= 1:
                fit_channels_raw.insert(0, ch_first)
                fit_counts_raw.insert(0, counts[ch_first - 1])
        for _ in range(pad_after):
            ch_last += 1
            if ch_last <= nch:
                fit_channels_raw.append(ch_last)
                fit_counts_raw.append(counts[ch_last - 1])
        fit_channels = fit_channels_raw
        fit_counts = fit_counts_raw
    else:
        fit_channels = fit_channels_raw
        fit_counts = fit_counts_raw

    print(f"     ROI ajuste: canales {fit_channels[0]}-{fit_channels[-1]} ({len(fit_channels)} canales)")
    print(f"     cuentas en ROI: {fit_counts}")

    # --- Fit with Fortran ---
    fit_ok = False
    if len(fit_channels) == 27 and FORTRAN_BIN.exists():
        pico_path = DATA_DIR / "pico.dat"
        with pico_path.open("w") as f:
            for i, y in enumerate(fit_counts, start=1):
                f.write(f"{i} {y}\n")

        print(f"\n  Ejecutando ajuste Fortran (gaussian_background)...")
        run = subprocess.run(
            [str(FORTRAN_BIN)],
            cwd=str(OUTDIR),
            capture_output=True, text=True, check=False,
        )
        print(f"  stdout:\n{run.stdout}")
        if run.stderr:
            print(f"  stderr:\n{run.stderr}")

        if run.returncode == 0:
            # Parse Fortran output
            f_x0 = re.search(r"Peak position\s*\(x0\)\s*=\s*([\-0-9.]+)", run.stdout)
            f_sigma = re.search(r"Width\s*\(sigma\)\s*=\s*([\-0-9.]+)", run.stdout)
            f_amp = re.search(r"Amplitude\s*\(A\)\s*=\s*([\-0-9.]+)", run.stdout)
            f_bg = re.search(r"Background:\s*y\s*=\s*([\-0-9.]+)\s*\+\s*([\-0-9.]+)", run.stdout)
            f_chi2 = re.search(r"chi\^2\s*=\s*([\-0-9.]+)", run.stdout)

            if f_x0 and f_sigma and f_amp:
                fit_ok = True
                x0_local = float(f_x0.group(1))
                x0_ch = fit_channels[0] + (x0_local - 1.0)
                sigma_ch = float(f_sigma.group(1))
                amp = float(f_amp.group(1))
                fwhm_ch = 2.3548 * sigma_ch
                area_net = amp * sigma_ch * math.sqrt(2.0 * math.pi)
                rate_cps = area_net / live_s if live_s > 0 else 0
                e_k40 = a_kev + b_kev_ch * x0_ch
                delta_e = e_k40 - K40_REF_KEV
                bg_a = float(f_bg.group(1)) if f_bg else 0
                bg_b = float(f_bg.group(2)) if f_bg else 0
                chi2 = float(f_chi2.group(1)) if f_chi2 else None

                print(f"\n  === RESULTADOS DEL AJUSTE ===")
                print(f"    x0 (local)   = {x0_local:.4f}")
                print(f"    x0 (canal)   = {x0_ch:.4f}")
                print(f"    sigma        = {sigma_ch:.4f} canales")
                print(f"    FWHM         = {fwhm_ch:.4f} canales")
                print(f"    Amplitud     = {amp:.2f} cuentas")
                print(f"    Area neta    = {area_net:.2f} cuentas")
                print(f"    Tasa neta    = {rate_cps:.6f} cps")
                print(f"    E(40K)       = {e_k40:.3f} keV")
                print(f"    Delta 40K    = {delta_e:.3f} keV")
                print(f"    Fondo: y = {bg_a:.4f} + {bg_b:.4f} * x")
                if chi2 is not None:
                    print(f"    chi^2        = {chi2:.6f}")

                # Generate fit curve
                n_points = 200
                x_local_vals = [1 + i * (27 - 1) / (n_points - 1) for i in range(n_points)]
                fit_curve = []
                for xl in x_local_vals:
                    bg = bg_a + bg_b * xl
                    g = amp * math.exp(-0.5 * ((xl - x0_local) / sigma_ch) ** 2)
                    ch_real = fit_channels[0] + (xl - 1.0)
                    fit_curve.append((ch_real, bg + g, bg, g))

                # Save fit curve
                curve_path = DATA_DIR / "k40_fit_curve.dat"
                with curve_path.open("w") as f:
                    for cr, ytot, ybg, yg in fit_curve:
                        f.write(f"{cr:.4f} {ytot:.6f} {ybg:.6f} {yg:.6f}\n")

                # Save results JSON
                results = {
                    "input_spectrum": str(SPECTRUM_PATH),
                    "live_time_s": live_s,
                    "live_time_ticks": live_ticks,
                    "num_channels": nch,
                    "total_counts": sum(counts),
                    "total_rate_cps": sum(counts) / live_s if live_s > 0 else 0,
                    "calibration": {
                        "a_keV": a_kev,
                        "b_keV_per_channel": b_kev_ch,
                        "anchors": [
                            {"channel": a[0], "energy_keV": a[1], "nuclide": "Bi-214"}
                            for a in cal_anchors
                        ],
                    },
                    "k40": {
                        "ref_energy_keV": K40_REF_KEV,
                        "seed_channel": seed_ch,
                        "roi_channels": [fit_channels[0], fit_channels[-1]],
                        "x0_local": x0_local,
                        "x0_channel": x0_ch,
                        "sigma_channel": sigma_ch,
                        "fwhm_channel": fwhm_ch,
                        "fwhm_keV": fwhm_ch * b_kev_ch,
                        "amplitude_counts": amp,
                        "area_net_counts": area_net,
                        "rate_net_cps": rate_cps,
                        "energy_estimated_keV": e_k40,
                        "energy_delta_keV": delta_e,
                        "fit_background_a": bg_a,
                        "fit_background_b": bg_b,
                        "chi2": chi2,
                        "fortran_stdout": run.stdout,
                    },
                }

                out_json = OUTDIR / "resultados_k40_rege.json"
                out_json.write_text(json.dumps(results, indent=2))
                print(f"\n  Resultados guardados: {out_json}")

                # --- Gnuplot scripts ---
                cal_str = f"E(keV) = {a_kev:.6f} + {b_kev_ch:.6f} * channel"
                e_str = f"{e_k40:.2f}"
                x0_str = f"{x0_ch:.4f}"

                # 1. Fit detail
                gp_fit = DATA_DIR / "plot_k40_fit.gp"
                png_fit = PLOTS_DIR / "k40_rege_fit.png"
                gp_fit.write_text(
                    "set terminal pngcairo size 1400,850 enhanced font 'sans,13'\n"
                    f"set output '{png_fit}'\n"
                    "set title 'Fondo REGe: pico 40K (ajuste gaussiano + fondo lineal)'\n"
                    "set xlabel 'Canal'\n"
                    "set ylabel 'Cuentas'\n"
                    "set grid lw 0.5\n"
                    "set key top left\n"
                    "plot 'data/k40_fit_curve.dat' using 1:2 with lines lw 2 lc rgb '#d62728' title 'Ajuste (gauss + fondo)', \\\n"
                    "     'data/k40_fit_curve.dat' using 1:3 with lines lw 1.5 lc rgb '#2ca02c' title 'Fondo lineal', \\\n"
                    "     'data/k40_fit_curve.dat' using 1:4 with lines lw 2 lc rgb '#ff7f0e' title 'Gaussiana (pico)', \\\n"
                    "'data/pico.dat' using ($1 + " + str(fit_channels[0] - 1) + "):2 with points pt 7 ps 1.2 lc rgb '#1f77b4' title 'Datos'\n"
                )

                # 2. Full spectrum (channel space)
                gp_full = DATA_DIR / "plot_k40_full.gp"
                png_full = PLOTS_DIR / "k40_rege_full.png"
                gp_full.write_text(
                    "set terminal pngcairo size 1600,900 enhanced font 'sans,13'\n"
                    f"set output '{png_full}'\n"
                    "set title 'Fondo REGe: espectro completo e identificacion del 40K'\n"
                    "set xlabel 'Canal'\n"
                    "set ylabel 'Cuentas'\n"
                    "set grid lw 0.5\n"
                    "set key top right\n"
                    "set logscale y\n"
                    f"set arrow from {x0_str},graph 0 to {x0_str},graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'\n"
                    f"set label 1 '40K ~ {e_str} keV' at {x0_ch + 120:.1f},graph 0.85 tc rgb '#cc0000'\n"
                    "plot 'data/fondo_rege.txt' using 1:2 with lines lw 1 lc rgb '#003366' title 'Espectro REGe', \\\n"
                    f"     '-' using 1:2 with points pt 7 ps 1.5 lc rgb '#cc0000' title 'Pico 40K'\n"
                    f"{x0_str} {amp + bg_a + bg_b * 14:.0f}\n"
                    "e\n"
                )

                # 3. Energy-calibrated spectrum
                gp_ecal = DATA_DIR / "plot_k40_energy.gp"
                png_ecal = PLOTS_DIR / "k40_rege_energy.png"
                gp_ecal.write_text(
                    "set terminal pngcairo size 1600,900 enhanced font 'sans,13'\n"
                    f"set output '{png_ecal}'\n"
                    "set title 'Fondo REGe: espectro completo calibrado en energia'\n"
                    "set xlabel 'Energia (keV)'\n"
                    "set ylabel 'Cuentas'\n"
                    "set grid lw 0.5\n"
                    "set key top right\n"
                    "set logscale y\n"
                    f"set arrow from {e_str},graph 0 to {e_str},graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'\n"
                    f"set label 1 '40K: {e_str} keV' at {e_k40 + 30:.1f},graph 0.85 tc rgb '#cc0000'\n"
                    f"E(ch) = {a_kev:.6f} + {b_kev_ch:.6f} * ch\n"
                    "plot 'data/fondo_rege.txt' using (E($1)):2 with lines lw 1 lc rgb '#003366' title 'Espectro REGe'\n"
                )

                # Run gnuplot scripts
                for gp_script in [gp_fit, gp_full, gp_ecal]:
                    if gp_script.exists():
                        gp_run = subprocess.run(
                            ["gnuplot", str(gp_script)],
                            cwd=str(OUTDIR),
                            capture_output=True, text=True, check=False,
                        )
                        if gp_run.returncode == 0:
                            print(f"  Grafico: {gp_script.parent.name}/{gp_script.name} -> {gp_script.stem}.png")
                        else:
                            print(f"  Error gnuplot {gp_script.name}: {gp_run.stderr}")

                print(f"\n  Todos los graficos en: {PLOTS_DIR}/")

    if not fit_ok:
        if not FORTRAN_BIN.exists():
            print(f"  Binario Fortran no encontrado: {FORTRAN_BIN}")
            print("  Realizando ajuste Python como fallback...")
            _fit_python(counts, fit_channels, fit_counts, live_s, a_kev, b_kev_ch,
                       seed_ch, OUTDIR, DATA_DIR, PLOTS_DIR)
        else:
            print("  Error: el ajuste Fortran no produjo resultados")

    print("\nHecho.")


def _fit_python(counts, fit_channels, fit_counts, live_s, a_kev, b_kev_ch,
                seed_ch, OUTDIR, DATA_DIR, PLOTS_DIR):
    """Fallback: fit using Python (scipy or manual simplex)."""
    try:
        import numpy as np
        from scipy.optimize import curve_fit

        x_data = np.arange(1, len(fit_counts) + 1, dtype=float)
        y_data = np.array(fit_counts, dtype=float)

        def model(x, a, x0, sigma, bg0, bg1):
            return bg0 + bg1 * x + a * np.exp(-0.5 * ((x - x0) / sigma) ** 2)

        x0_guess = fit_counts.index(max(fit_counts)) + 1
        a_guess = max(fit_counts) - np.median(fit_counts)
        sigma_guess = 3.0
        bg0_guess = np.median(fit_counts[:5])
        bg1_guess = 0.0

        p0 = [a_guess, x0_guess, sigma_guess, bg0_guess, bg1_guess]
        popt, pcov = curve_fit(
            model, x_data, y_data, p0=p0,
            sigma=np.sqrt(np.maximum(y_data, 1)),
            absolute_sigma=True,
            maxfev=10000,
        )
        perr = np.sqrt(np.diag(pcov))

        a, x0_l, sigma, bg0, bg1 = popt
        a_err, x0_l_err, sigma_err, bg0_err, bg1_err = perr

        x0_ch = fit_channels[0] + (x0_l - 1.0)
        fwhm_ch = 2.3548 * sigma
        area_net = a * sigma * np.sqrt(2.0 * np.pi)
        rate_cps = area_net / live_s if live_s > 0 else 0
        e_k40 = a_kev + b_kev_ch * x0_ch
        delta_e = e_k40 - 1460.82

        # Compute chi2
        y_fit = model(x_data, *popt)
        chi2 = np.sum((y_data - y_fit) ** 2 / np.maximum(y_data, 1))
        ndf = len(y_data) - len(popt)

        print(f"\n  === RESULTADOS DEL AJUSTE (Python/scipy) ===")
        print(f"    x0 (local)   = {x0_l:.4f} +/- {x0_l_err:.4f}")
        print(f"    x0 (canal)   = {x0_ch:.4f}")
        print(f"    sigma        = {sigma:.4f} +/- {sigma_err:.4f} canales")
        print(f"    FWHM         = {fwhm_ch:.4f} canales")
        print(f"    Amplitud     = {a:.2f} +/- {a_err:.2f}")
        print(f"    Area neta    = {area_net:.2f}")
        print(f"    Tasa neta    = {rate_cps:.6f} cps")
        print(f"    E(40K)       = {e_k40:.3f} keV")
        print(f"    Delta 40K    = {delta_e:.3f} keV")
        print(f"    Fondo: y = {bg0:.4f} + {bg1:.4f} * x")
        print(f"    chi^2/ndf    = {chi2:.2f}/{ndf}")

        # Save fit curve
        n_pts = 200
        x_curve = np.linspace(1, len(fit_counts), n_pts)
        y_curve = model(x_curve, *popt)
        y_bg = bg0 + bg1 * x_curve
        y_g = a * np.exp(-0.5 * ((x_curve - x0_l) / sigma) ** 2)
        curve_path = DATA_DIR / "k40_fit_curve.dat"
        with curve_path.open("w") as f:
            for i in range(n_pts):
                ch_real = fit_channels[0] + (x_curve[i] - 1.0)
                f.write(f"{ch_real:.4f} {y_curve[i]:.6f} {y_bg[i]:.6f} {y_g[i]:.6f}\n")

        # Save results
        results = {
            "input_spectrum": str(SPECTRUM_PATH),
            "live_time_s": live_s,
            "num_channels": 8192,
            "total_counts": int(sum(counts)),
            "total_rate_cps": sum(counts) / live_s if live_s > 0 else 0,
            "calibration": {
                "a_keV": a_kev,
                "b_keV_per_channel": b_kev_ch,
                "anchors": [{"channel": 2160, "energy_keV": 609.31},
                           {"channel": 6179, "energy_keV": 1764.49}],
            },
            "k40": {
                "ref_energy_keV": 1460.82,
                "seed_channel": seed_ch,
                "roi_channels": [fit_channels[0], fit_channels[-1]],
                "x0_local": float(x0_l),
                "x0_local_err": float(x0_l_err),
                "x0_channel": float(x0_ch),
                "sigma_channel": float(sigma),
                "sigma_channel_err": float(sigma_err),
                "fwhm_channel": float(fwhm_ch),
                "fwhm_keV": float(fwhm_ch * b_kev_ch),
                "amplitude_counts": float(a),
                "amplitude_err": float(a_err),
                "area_net_counts": float(area_net),
                "rate_net_cps": float(rate_cps),
                "energy_estimated_keV": float(e_k40),
                "energy_delta_keV": float(delta_e),
                "fit_background_a": float(bg0),
                "fit_background_b": float(bg1),
                "chi2": float(chi2),
                "ndf": ndf,
                "method": "python_scipy",
            },
        }
        out_json = OUTDIR / "resultados_k40_rege.json"
        out_json.write_text(json.dumps(results, indent=2))
        print(f"\n  Resultados guardados: {out_json}")

        # Gnuplot scripts
        cal_str = f"E(keV) = {a_kev:.6f} + {b_kev_ch:.6f} * channel"
        e_str = f"{e_k40:.2f}"
        x0_str = f"{x0_ch:.4f}"

        gp_fit = DATA_DIR / "plot_k40_fit.gp"
        png_fit = PLOTS_DIR / "k40_rege_fit.png"
        gp_fit.write_text(
            "set terminal pngcairo size 1400,850 enhanced font 'sans,13'\n"
            f"set output '{png_fit}'\n"
            "set title 'Fondo REGe: pico 40K (ajuste gaussiano + fondo lineal)'\n"
            "set xlabel 'Canal'\n"
            "set ylabel 'Cuentas'\n"
            "set grid lw 0.5\n"
            "set key top left\n"
            "plot 'data/k40_fit_curve.dat' using 1:2 with lines lw 2 lc rgb '#d62728' title 'Ajuste (gauss + fondo)', \\\n"
            "     'data/k40_fit_curve.dat' using 1:3 with lines lw 1.5 lc rgb '#2ca02c' title 'Fondo lineal', \\\n"
            "     'data/k40_fit_curve.dat' using 1:4 with lines lw 2 lc rgb '#ff7f0e' title 'Gaussiana (pico)', \\\n"
            f"'data/pico.dat' using ($1 + {fit_channels[0] - 1}):2 with points pt 7 ps 1.2 lc rgb '#1f77b4' title 'Datos'\n"
        )

        gp_full = DATA_DIR / "plot_k40_full.gp"
        png_full = PLOTS_DIR / "k40_rege_full.png"
        gp_full.write_text(
            "set terminal pngcairo size 1600,900 enhanced font 'sans,13'\n"
            f"set output '{png_full}'\n"
            "set title 'Fondo REGe: espectro completo e identificacion del 40K'\n"
            "set xlabel 'Canal'\n"
            "set ylabel 'Cuentas'\n"
            "set grid lw 0.5\n"
            "set key top right\n"
            "set logscale y\n"
            f"set arrow from {x0_str},graph 0 to {x0_str},graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'\n"
            f"set label 1 '40K ~ {e_str} keV' at {x0_ch + 120:.1f},graph 0.85 tc rgb '#cc0000'\n"
            "plot 'data/fondo_rege.txt' using 1:2 with lines lw 1 lc rgb '#003366' title 'Espectro REGe', \\\n"
            f"     '-' using 1:2 with points pt 7 ps 1.5 lc rgb '#cc0000' title 'Pico 40K'\n"
            f"{x0_str} {a + bg0 + bg1 * 14:.0f}\n"
            "e\n"
        )

        gp_ecal = DATA_DIR / "plot_k40_energy.gp"
        png_ecal = PLOTS_DIR / "k40_rege_energy.png"
        gp_ecal.write_text(
            "set terminal pngcairo size 1600,900 enhanced font 'sans,13'\n"
            f"set output '{png_ecal}'\n"
            "set title 'Fondo REGe: espectro completo calibrado en energia'\n"
            "set xlabel 'Energia (keV)'\n"
            "set ylabel 'Cuentas'\n"
            "set grid lw 0.5\n"
            "set key top right\n"
            "set logscale y\n"
            f"set arrow from {e_str},graph 0 to {e_str},graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'\n"
            f"set label 1 '40K: {e_str} keV' at {e_k40 + 30:.1f},graph 0.85 tc rgb '#cc0000'\n"
            f"E(ch) = {a_kev:.6f} + {b_kev_ch:.6f} * ch\n"
            "plot 'data/fondo_rege.txt' using (E($1)):2 with lines lw 1 lc rgb '#003366' title 'Espectro REGe'\n"
        )

        for gp_script in [gp_fit, gp_full, gp_ecal]:
            if gp_script.exists():
                gp_run = subprocess.run(
                    ["gnuplot", str(gp_script)],
                    cwd=str(OUTDIR),
                    capture_output=True, text=True, check=False,
                )
                if gp_run.returncode == 0:
                    print(f"  Grafico: {gp_script.stem}")
                else:
                    print(f"  Error gnuplot {gp_script.name}: {gp_run.stderr}")

    except ImportError:
        print("  scipy no disponible. No se puede hacer ajuste.")

if __name__ == "__main__":
    main()
