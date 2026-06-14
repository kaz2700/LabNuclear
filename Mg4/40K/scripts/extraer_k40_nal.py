#!/usr/bin/env python3
"""
Extraer picos de 40K y aniquilacion del espectro NaI(Tl)
usando el programa Fortran gaussian_background (modificado).

Procedimiento:
  1. Lee fondo_nal.txt (512 canales)
  2. Define ROIs para aniquilacion (~ch 145-195) y 40K (~ch 435-505)
  3. Ejecuta Fortran para cada ROI -> obtiene x0, sigma, amplitud
  4. Calcula calibracion E(ch) = a + b*ch
  5. Guarda resultados en resultados_k40_nal.json
  6. Genera curvas de ajuste y scripts gnuplot
  7. Ejecuta gnuplot
"""

import math
import json
import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTDIR = HERE.parent
DATA_DIR = OUTDIR / "data"
PLOTS_DIR = OUTDIR / "plots"
SPECTRUM_PATH = DATA_DIR / "fondo_nal.txt"
FORTRAN_BIN = ROOT.parent / "src" / "gaussian_background"
PICO_PATH = DATA_DIR / "pico.dat"

K40_REF_KEV = 1460.82
ANN_REF_KEV = 511.0


def read_spectrum(path):
    lines = path.read_text().splitlines()
    data = []
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        ch = int(parts[0])
        cnt = int(parts[1])
        data.append((ch, cnt))
    return data


def run_fortran_on_roi(channels):
    """
    Run Fortran with the given list of (global_channel, counts).
    Returns parsed results dict.
    """
    nf = len(channels)
    with PICO_PATH.open("w") as f:
        for i, (_, cnt) in enumerate(channels, 1):
            f.write(f"{i} {cnt}\n")

    print(f"  ROI: ch {channels[0][0]}-{channels[-1][0]} ({nf} channels)")
    print(f"  Running Fortran...")
    run = subprocess.run(
        [str(FORTRAN_BIN), "0.60"],
        cwd=str(OUTDIR),
        capture_output=True, text=True, check=False,
    )
    stdout = run.stdout
    if run.stderr:
        print(f"  stderr: {run.stderr}")

    print(stdout)

    f_x0 = re.search(r"Peak position\s*\(x0\)\s*=\s*([\d.-]+)", stdout)
    f_sigma = re.search(r"Width\s*\(sigma\)\s*=\s*([\d.-]+)", stdout)
    f_amp = re.search(r"Amplitude\s*\(A\)\s*=\s*([\d.-]+)", stdout)
    f_bg = re.search(r"Background:\s*y\s*=\s*([\d.-]+)\s*\+\s*([\d.-]+)", stdout)
    f_chi2 = re.search(r"chi\^2\s*=\s*([\d.-]+)", stdout)

    if not (f_x0 and f_sigma and f_amp):
        print("  ERROR: Could not parse Fortran output")
        return None

    x0_local = float(f_x0.group(1))
    sigma = float(f_sigma.group(1))
    amplitude = float(f_amp.group(1))
    chi2 = float(f_chi2.group(1)) if f_chi2 else None
    bg_a = float(f_bg.group(1)) if f_bg else 0
    bg_b = float(f_bg.group(2)) if f_bg else 0

    first_ch = channels[0][0]
    x0_global = first_ch + (x0_local - 1.0)
    fwhm_ch = 2.3548 * sigma

    return {
        "x0_local": x0_local,
        "x0_global": x0_global,
        "sigma": sigma,
        "fwhm_channel": fwhm_ch,
        "amplitude": amplitude,
        "chi2": chi2,
        "bg_a": bg_a,
        "bg_b": bg_b,
        "nf": nf,
    }


def generate_fit_curve(first_ch, nf, x0_local, sigma, amplitude, bg_a, bg_b):
    n_pts = 200
    curve = []
    for i in range(n_pts):
        x_local = 1 + i * (nf - 1) / (n_pts - 1)
        bg = bg_a + bg_b * x_local
        g = amplitude * math.exp(-0.5 * ((x_local - x0_local) / sigma) ** 2)
        ch_global = first_ch + (x_local - 1.0)
        curve.append((ch_global, bg + g, bg, g))
    return curve


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    # --- Read spectrum ---
    data = read_spectrum(SPECTRUM_PATH)
    total_counts = sum(c for _, c in data)
    print(f"Espectro NaI(Tl): {len(data)} canales, {total_counts:,} cuentas")
    print()

    # --- ROI definitions (from visual inspection of the spectrum) ---
    # Annihilation peak: shoulder around ch 170
    # K-40 peak: main peak around ch 468 (max at 468)
    ann_roi = [(ch, cnt) for ch, cnt in data if 140 <= ch <= 195]
    k40_roi = [(ch, cnt) for ch, cnt in data if 430 <= ch <= 502]

    print("=== Pico de aniquilacion (511 keV) ===")
    ann_result = run_fortran_on_roi(ann_roi)
    if ann_result is None:
        print("FAILED")
        return
    x0_ann = ann_result["x0_global"]
    print(f"  x0 (global)  = {x0_ann:.4f}")
    print(f"  sigma        = {ann_result['sigma']:.4f} canales")
    print(f"  FWHM         = {ann_result['fwhm_channel']:.4f} canales")
    print(f"  Amplitud     = {ann_result['amplitude']:.2f}")
    print(f"  chi^2        = {ann_result['chi2']:.4f}")
    print()

    print("=== Pico de 40K (1460.82 keV) ===")
    k40_result = run_fortran_on_roi(k40_roi)
    if k40_result is None:
        print("FAILED")
        return
    x0_k40 = k40_result["x0_global"]
    print(f"  x0 (global)  = {x0_k40:.4f}")
    print(f"  sigma        = {k40_result['sigma']:.4f} canales")
    print(f"  FWHM         = {k40_result['fwhm_channel']:.4f} canales")
    print(f"  Amplitud     = {k40_result['amplitude']:.2f}")
    print(f"  chi^2        = {k40_result['chi2']:.4f}")
    print()

    # --- Calibration ---
    b_kev = (K40_REF_KEV - ANN_REF_KEV) / (x0_k40 - x0_ann)
    a_kev = ANN_REF_KEV - b_kev * x0_ann
    print(f"=== Calibracion ===")
    print(f"  Anchors: ch {x0_ann:.2f} -> 511.0 keV, ch {x0_k40:.2f} -> 1460.82 keV")
    print(f"  E(keV) = {a_kev:.4f} + {b_kev:.6f} * ch")
    print()

    # --- Net K-40 parameters ---
    live_s = 6375.0
    sigma_k40 = k40_result["sigma"]
    amp_k40 = k40_result["amplitude"]
    area_net = amp_k40 * sigma_k40 * math.sqrt(2.0 * math.pi)
    rate_net = area_net / live_s
    fwhm_kev = k40_result["fwhm_channel"] * b_kev
    energy_est = a_kev + b_kev * x0_k40
    delta_e = energy_est - K40_REF_KEV
    nf_k40 = k40_result["nf"]
    chi2_k40 = k40_result["chi2"]
    ndf_k40 = k40_result["nf"] - 3

    print(f"=== Parametros netos 40K ===")
    print(f"  x0            = {x0_k40:.4f} canales")
    print(f"  sigma         = {sigma_k40:.4f} canales ({sigma_k40 * b_kev:.2f} keV)")
    print(f"  FWHM          = {k40_result['fwhm_channel']:.4f} canales ({fwhm_kev:.2f} keV)")
    print(f"  Amplitud      = {amp_k40:.2f} cuentas")
    print(f"  Area neta     = {area_net:.2f} cuentas")
    print(f"  Tasa neta     = {rate_net:.4f} cps")
    print(f"  E estimada    = {energy_est:.3f} keV (delta = {delta_e:.3f} keV)")
    print(f"  chi^2 / ndf   = {chi2_k40:.2f} / {ndf_k40}")
    print()

    # --- Generate fit curves ---
    print("Generando curvas de ajuste...")

    ann_curve = generate_fit_curve(
        ann_roi[0][0], len(ann_roi),
        ann_result["x0_local"], ann_result["sigma"],
        ann_result["amplitude"], ann_result["bg_a"], ann_result["bg_b"],
    )
    curve_path = DATA_DIR / "nal_ann_fit_curve.dat"
    with curve_path.open("w") as f:
        f.write("#ch_global y_total y_bg y_gaussian\n")
        for ch, ytot, ybg, yg in ann_curve:
            f.write(f"{ch:.4f} {ytot:.6f} {ybg:.6f} {yg:.6f}\n")

    k40_curve = generate_fit_curve(
        k40_roi[0][0], len(k40_roi),
        k40_result["x0_local"], k40_result["sigma"],
        k40_result["amplitude"], k40_result["bg_a"], k40_result["bg_b"],
    )
    curve_path = DATA_DIR / "nal_fit_curve.dat"
    with curve_path.open("w") as f:
        f.write("#ch_global y_total y_bg y_gaussian\n")
        for ch, ytot, ybg, yg in k40_curve:
            f.write(f"{ch:.4f} {ytot:.6f} {ybg:.6f} {yg:.6f}\n")

    # --- Save results JSON ---
    results = {
        "detector": "NaI(Tl)",
        "live_time_s": live_s,
        "num_channels": len(data),
        "total_counts": total_counts,
        "total_rate_cps": total_counts / live_s if live_s > 0 else 0,
        "calibration": {
            "a_keV": a_kev,
            "b_keV_per_channel": b_kev,
            "anchors": [
                {
                    "channel": x0_ann,
                    "fitted_x0": x0_ann,
                    "energy_keV": ANN_REF_KEV,
                    "nuclide": "e+e-",
                },
                {
                    "channel": x0_k40,
                    "fitted_x0": x0_k40,
                    "energy_keV": K40_REF_KEV,
                    "nuclide": "K-40",
                },
            ],
            "formula": f"E = {a_kev:.4f} + {b_kev:.6f}*ch",
        },
        "annihilation": {
            "x0_channel": x0_ann,
            "sigma_channel": ann_result["sigma"],
            "fwhm_channel": ann_result["fwhm_channel"],
            "amplitude_counts": ann_result["amplitude"],
            "chi2": ann_result["chi2"],
            "ndf": len(ann_roi) - 3,
            "background_a": ann_result["bg_a"],
            "background_b": ann_result["bg_b"],
        },
        "k40": {
            "x0_channel": x0_k40,
            "sigma_channel": sigma_k40,
            "fwhm_channel": k40_result["fwhm_channel"],
            "fwhm_keV": fwhm_kev,
            "amplitude_counts": amp_k40,
            "area_net_counts": area_net,
            "rate_net_cps": rate_net,
            "energy_estimated_keV": energy_est,
            "energy_delta_keV": delta_e,
            "chi2": chi2_k40,
            "chi2_red": chi2_k40 / ndf_k40 if ndf_k40 > 0 else 0,
            "ndf": ndf_k40,
            "roi_channels": [k40_roi[0][0], k40_roi[-1][0]],
            "background_a": k40_result["bg_a"],
            "background_b": k40_result["bg_b"],
        },
    }

    out_json = OUTDIR / "resultados_k40_nal.json"
    out_json.write_text(json.dumps(results, indent=2))
    print(f"Resultados guardados: {out_json}")
    print()

    # --- Generate gnuplot scripts ---
    print("Generando scripts gnuplot...")

    # 1. Full spectrum
    gp = DATA_DIR / "plot_nal_full.gp"
    gp.write_text(
        "set terminal pngcairo size 1800,1000 enhanced font 'sans,14'\n"
        f"set output '{PLOTS_DIR / 'nal_full.png'}'\n"
        "set title 'Espectro de fondo detector NaI(Tl) (512 canales, t_l = 6375 s)'\n"
        "set xlabel 'Canal'\n"
        "set ylabel 'Cuentas'\n"
        "set grid lw 0.5\n"
        "set logscale y\n"
        f"set arrow from {x0_k40:.1f},graph 0 to {x0_k40:.1f},graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'\n"
        f"set label 1 '^{{40}}K: ch {x0_k40:.1f}' at {x0_k40 + 10:.1f},graph 0.85 tc rgb '#cc0000'\n"
        "plot 'data/fondo_nal.txt' using 1:2 with lines lw 1 lc rgb '#003366' title 'Fondo NaI(Tl)'\n"
    )

    # 2. Energy-calibrated
    gp = DATA_DIR / "plot_nal_energy.gp"
    ann_e = a_kev + b_kev * x0_ann
    k40_e = a_kev + b_kev * x0_k40
    gp.write_text(
        "set terminal pngcairo size 1800,1000 enhanced font 'sans,14'\n"
        f"set output '{PLOTS_DIR / 'nal_energy.png'}'\n"
        f"set title 'Espectro NaI(Tl) calibrado en energia'\n"
        "set xlabel 'Energia (keV)'\n"
        "set ylabel 'Cuentas'\n"
        "set grid lw 0.5\n"
        "set logscale y\n"
        f"E(ch) = {a_kev:.6f} + {b_kev:.6f} * ch\n"
        f"set arrow from {k40_e:.1f},graph 0 to {k40_e:.1f},graph 1 nohead dt 2 lw 1.5 lc rgb '#cc0000'\n"
        f"set label 1 '^{{40}}K: {k40_e:.1f} keV' at {k40_e + 30:.1f},graph 0.85 tc rgb '#cc0000'\n"
        "plot 'data/fondo_nal.txt' using (E($1)):2 with lines lw 1 lc rgb '#003366' title 'Fondo NaI(Tl)'\n"
    )

    # 3. Calibration verification
    gp = DATA_DIR / "plot_nal_calibracion.gp"
    lines_ref = [
        ("511.0", "e+e-", x0_ann),
        ("583.2", "Tl-208", (583.2 - a_kev) / b_kev),
        ("609.3", "Bi-214", (609.3 - a_kev) / b_kev),
        ("911.2", "Ac-228", (911.2 - a_kev) / b_kev),
        ("969.0", "Ac-228", (969.0 - a_kev) / b_kev),
        ("1120.3", "Bi-214", (1120.3 - a_kev) / b_kev),
        ("1460.8", "K-40", x0_k40),
        ("1764.5", "Bi-214", (1764.5 - a_kev) / b_kev),
    ]
    arrow_lines = ""
    for ename, enuc, ech in lines_ref:
        if 1 <= ech <= 512:
            arrow_lines += (
                f"set arrow from {ech:.1f},graph 0 to {ech:.1f},graph 1 "
                f"nohead dt 3 lw 1 lc rgb '#888888'\n"
            )
    anchors_str = f"set arrow from {x0_ann:.1f},graph 0 to {x0_ann:.1f},graph 1 nohead dt 2 lw 2 lc rgb '#cc0000'\n"
    anchors_str += f"set arrow from {x0_k40:.1f},graph 0 to {x0_k40:.1f},graph 1 nohead dt 2 lw 2 lc rgb '#cc0000'\n"
    gp.write_text(
        "set terminal pngcairo size 2000,1200 enhanced font 'sans,13'\n"
        f"set output '{PLOTS_DIR / 'nal_calibracion.png'}'\n"
        "set title 'Verificacion calibracion NaI(Tl): lineas de fondo esperadas'\n"
        "set xlabel 'Canal'\n"
        "set ylabel 'Cuentas'\n"
        "set grid lw 0.5\n"
        "set logscale y\n"
        f"{arrow_lines}"
        f"{anchors_str}"
        "plot 'data/fondo_nal.txt' using 1:2 with lines lw 1 lc rgb '#003366' title 'Fondo NaI(Tl)'\n"
    )

    # 4. Annihilation fit
    gp = DATA_DIR / "plot_nal_ann_fit.gp"
    gp.write_text(
        "set terminal pngcairo size 1200,700 enhanced font 'sans,13'\n"
        f"set output '{PLOTS_DIR / 'nal_ann_fit.png'}'\n"
        f"set title 'Pico de aniquilacion (511 keV) en NaI(Tl)'\n"
        "set xlabel 'Canal'\n"
        "set ylabel 'Cuentas'\n"
        "set grid lw 0.5\n"
        "set key top left\n"
        f"set arrow from {x0_ann:.1f},graph 0 to {x0_ann:.1f},graph 1 "
        f"nohead dt 2 lw 1.5 lc rgb '#cc0000'\n"
        "plot 'data/nal_ann_fit_curve.dat' using 1:2 with lines lw 2 lc rgb '#d62728' "
        "title 'Ajuste (gauss+fondo)', \\\n"
        "     'data/nal_ann_fit_curve.dat' using 1:3 with lines lw 1.5 lc rgb '#2ca02c' "
        "title 'Fondo lineal', \\\n"
        f"     'data/fondo_nal.txt' using 1:2 every ::{ann_roi[0][0]-1}::{ann_roi[-1][0]-1} "
        "with points pt 7 ps 1.2 lc rgb '#1f77b4' title 'Datos'\n"
    )

    # 5. K-40 fit
    gp = DATA_DIR / "plot_nal_k40_fit.gp"
    gp.write_text(
        "set terminal pngcairo size 1400,850 enhanced font 'sans,13'\n"
        f"set output '{PLOTS_DIR / 'nal_k40_fit.png'}'\n"
        "set title 'Fondo NaI(Tl): pico K-40 (ajuste gaussiano + fondo lineal)'\n"
        "set xlabel 'Canal'\n"
        "set ylabel 'Cuentas'\n"
        "set grid lw 0.5\n"
        "set key top left\n"
        f"set arrow from {x0_k40:.1f},graph 0 to {x0_k40:.1f},graph 1 "
        f"nohead dt 2 lw 1.5 lc rgb '#cc0000'\n"
        "plot 'data/nal_fit_curve.dat' using 1:2 with lines lw 2 lc rgb '#d62728' "
        "title 'Ajuste (gauss+fondo)', \\\n"
        "     'data/nal_fit_curve.dat' using 1:3 with lines lw 1.5 lc rgb '#2ca02c' "
        "title 'Fondo lineal', \\\n"
        f"     'data/fondo_nal.txt' using 1:2 every ::{k40_roi[0][0]-1}::{k40_roi[-1][0]-1} "
        "with points pt 7 ps 1.2 lc rgb '#1f77b4' title 'Datos'\n"
    )

    # 6. K-40 region zoom
    gp = DATA_DIR / "plot_nal_k40_region.gp"
    gp.write_text(
        "set terminal pngcairo size 1400,700 enhanced font 'sans,14'\n"
        f"set output '{PLOTS_DIR / 'nal_k40_region.png'}'\n"
        "set title 'Region del K-40 (canales "
        f"{k40_roi[0][0]}-{k40_roi[-1][0]})'\n"
        "set xlabel 'Canal'\n"
        "set ylabel 'Cuentas'\n"
        "set grid lw 0.5\n"
        "set key top left\n"
        f"set arrow from {x0_k40:.1f},graph 0 to {x0_k40:.1f},graph 1 "
        f"nohead dt 2 lw 1.5 lc rgb '#cc0000'\n"
        f"plot [{k40_roi[0][0]}:{k40_roi[-1][0]}] 'data/fondo_nal.txt' using 1:2 "
        "with lines lw 1.5 lc rgb '#003366' title 'Datos', \\\n"
        "     'data/nal_fit_curve.dat' using 1:2 with lines lw 2 lc rgb '#d62728' "
        "title 'Ajuste'\n"
    )

    # 7. Calibration verification (energy-calibrated)
    gp = DATA_DIR / "plot_nal_calibracion_energy.gp"
    gp.write_text(
        "set terminal pngcairo size 2000,1200 enhanced font 'sans,13'\n"
        f"set output '{PLOTS_DIR / 'nal_calibracion_energy.png'}'\n"
        "set title 'Verificacion calibracion NaI(Tl): lineas de fondo en energia'\n"
        "set xlabel 'Energia (keV)'\n"
        "set ylabel 'Cuentas'\n"
        "set grid lw 0.5\n"
        "set logscale y\n"
        f"E(ch) = {a_kev:.6f} + {b_kev:.6f} * ch\n"
        "plot 'data/fondo_nal.txt' using (E($1)):2 with lines lw 1 lc rgb '#003366' "
        "title 'Fondo NaI(Tl)'\n"
    )

    # --- Run gnuplot ---
    print("Ejecutando gnuplot...")
    for gp_name in [
        "plot_nal_full.gp",
        "plot_nal_energy.gp",
        "plot_nal_ann_fit.gp",
        "plot_nal_k40_fit.gp",
        "plot_nal_k40_region.gp",
        "plot_nal_calibracion.gp",
        "plot_nal_calibracion_energy.gp",
    ]:
        gp_path = DATA_DIR / gp_name
        if gp_path.exists():
            result = subprocess.run(
                ["gnuplot", str(gp_path)],
                cwd=str(OUTDIR),
                capture_output=True, text=True, check=False,
            )
            if result.returncode == 0:
                print(f"  OK: {gp_name} -> {gp_path.stem}.png")
            else:
                print(f"  ERROR: {gp_name}: {result.stderr}")

    print("\nHecho.")


if __name__ == "__main__":
    main()
