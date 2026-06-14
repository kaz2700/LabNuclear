#!/usr/bin/env python3
"""Calibracion multipunto usando picos de fondo identificados manualmente.

Puntos de calibracion (identificados por inspeccion visual del espectro):
  114.673 keV -> ch 468  (pico intenso, origen: muestra Mg4 + fondo)
  186.211 keV -> ch 691  (226Ra, pico claro)
  238.632 keV -> ch 905  (212Pb, pico pequeno pero claro)
  1460.820 keV -> ch 5189 (40K, pico claro)
  2614.511 keV -> ch 9205 (208Tl, pico claro)

Se usa regresion lineal con los 5 puntos.
Luego se verifica cada linea de fondo contra el maximo local mas cercano.
"""
import json
import math
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK_DIR = ROOT / "caracterizacion_sin_fondo"
DATA_DIR = WORK_DIR / "data"
PLOTS_DIR = WORK_DIR / "plots"
SPECTRUM_PATH = ROOT / "datos" / "Mg4_GeHP.txt"
FONDO_PATH = ROOT / "datos" / "fondo_GeHP_BEGe.dat"
FORTRAN_BIN = ROOT.parent / "src" / "gaussian_background"

# Puntos de calibracion identificados manualmente (canal, energia)
ANCHORS = [
    (468, 114.673),
    (691, 186.211),
    (905, 238.632),
    (5189, 1460.820),
    (9205, 2614.511),
]


def parse_hpge_txt(path):
    live_s = None
    counts = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("#"):
            m = re.search(r"Live:\s*([0-9]+(?:\.[0-9]+)?)\s*s", line)
            if m:
                live_s = float(m.group(1))
            continue
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        counts.append(int(parts[1]))
    return live_s, counts


def parse_fondo_dat(path):
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    lines = [l.rstrip() for l in lines if l.strip() and not l.startswith("Energ")]
    entries = []
    for line in lines:
        parts = line.split("\t")
        parts = [p for p in parts if p]
        if len(parts) < 2:
            continue
        e_str = parts[0].strip().replace(",", ".")
        rest = parts[1].strip()
        rate_match = re.search(r"(\d+)\s*\(\s*(\d+)\s*\)\s*$", rest)
        if not rate_match:
            continue
        rate_val = float(rate_match.group(1))
        rate_err = float(rate_match.group(2))
        reaction = rest[:rate_match.start()].strip()
        try:
            energy = float(e_str)
        except ValueError:
            continue
        entries.append({
            "energy_keV": energy,
            "reaction": reaction,
            "rate_cps_e6": rate_val,
            "rate_err_cps_e6": rate_err,
        })
    return entries


def find_local_max(counts, center_ch, half_window):
    lo = max(0, center_ch - half_window - 1)
    hi = min(len(counts), center_ch + half_window)
    region = counts[lo:hi]
    if not region:
        return None, None
    max_val = max(region)
    max_idx = region.index(max_val) + lo + 1
    return max_idx, max_val


def is_peak(counts, ch_obs, half=5):
    lo = max(0, ch_obs - half - 1)
    hi = min(len(counts), ch_obs + half)
    near = [counts[i] for i in range(lo, hi)]
    if not near:
        return False
    mx_idx = near.index(max(near))
    if mx_idx < half or mx_idx >= len(near) - half:
        return False
    for k in range(1, half + 1):
        if not (near[mx_idx] > near[mx_idx - k] and near[mx_idx] > near[mx_idx + k]):
            return False
    return True


def is_significant_peak(counts, ch, half=5, min_prominence=1.5):
    if not is_peak(counts, ch, half):
        return False
    lo = max(0, ch - 2 * half - 1)
    hi = min(len(counts), ch + 2 * half)
    region = [counts[i] for i in range(lo, hi)]
    bg = sorted(region)[len(region) // 4]
    cnt = counts[ch - 1]
    if cnt < bg * min_prominence:
        return False
    return cnt >= 15


def linear_regression(points):
    n = len(points)
    if n < 2:
        return None, None, None
    sx = sum(p[0] for p in points)
    sy = sum(p[1] for p in points)
    sxx = sum(p[0] ** 2 for p in points)
    sxy = sum(p[0] * p[1] for p in points)
    denom = n * sxx - sx * sx
    if denom == 0:
        return None, None, None
    b = (n * sxy - sx * sy) / denom
    a = (sy - b * sx) / n
    if n > 2:
        ss_res = sum((p[1] - (a + b * p[0])) ** 2 for p in points)
        ss_tot = sum((p[1] - sy / n) ** 2 for p in points)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    else:
        r2 = 1.0
    return a, b, r2


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    live_s, counts = parse_hpge_txt(SPECTRUM_PATH)
    n = len(counts)

    fondo_entries = parse_fondo_dat(FONDO_PATH)

    # --- Calibracion lineal con los 5 anclas ---
    a_cal, b_cal, r2 = linear_regression(ANCHORS)

    # --- Verificar cada linea de fondo ---
    verification_lines = []
    problems = []

    for entry in fondo_entries:
        e_ref = entry["energy_keV"]
        ch_pred = (e_ref - a_cal) / b_cal
        if ch_pred < 1 or ch_pred > n:
            continue

        ch_obs, cnt_obs = find_local_max(counts, round(ch_pred), half_window=30)
        if ch_obs is None:
            problems.append(f"  {e_ref:.1f} keV: no se encontro maximo local cerca de ch {ch_pred:.1f}")
            continue

        e_obs = a_cal + b_cal * ch_obs
        delta_e = e_obs - e_ref
        es_peak = is_peak(counts, ch_obs)
        significativo = cnt_obs >= 20 and is_significant_peak(counts, ch_obs, half=5, min_prominence=1.3)

        # Ver si este canal es uno de los anclas
        in_cal = False
        for ach, ae in ANCHORS:
            if abs(ae - e_ref) < 0.1:
                in_cal = True
                break

        if abs(delta_e) > 6:
            problems.append(
                f"  {e_ref:.1f} keV: discrepancia grande dE={delta_e:+.1f} keV "
                f"(ch pred {ch_pred:.1f}, ch obs {ch_obs}, cnts {cnt_obs})"
            )

        verification_lines.append({
            "energy_ref_keV": e_ref,
            "reaction": entry["reaction"],
            "channel_predicted": ch_pred,
            "channel_observed": ch_obs,
            "counts_observed": cnt_obs,
            "energy_observed_keV": e_obs,
            "delta_keV": delta_e,
            "is_peak": es_peak,
            "significant": significativo,
            "in_calibration": in_cal,
        })

    # --- Print results ---
    print("=" * 78)
    print("CALIBRACION CON LINEAS DE FONDO (5 anclas identificadas manualmente)")
    print("=" * 78)
    print(f"\nEspectro: Mg4_GeHP.txt ({live_s} s, {n} canales)")
    print(f"\nAnclas de calibracion:")
    for ch, e in ANCHORS:
        print(f"  ch {ch:>5d} -> {e:.3f} keV")
    print()
    print(f"Calibracion lineal ({len(ANCHORS)} puntos, R² = {r2:.6f}):")
    print(f"  E(keV) = {a_cal:.4f} + {b_cal:.6f} * canal")
    print()

    print(f"Verificacion de {len(verification_lines)} lineas de fondo:")
    print(f"{'E_ref':>8s} {'Reaccion':<22s} {'Ch_pred':>7s} {'Ch_obs':>6s} {'Cnts':>5s} {'dE(keV)':>8s} {'Pico':>5s} {'Cal':>3s}")
    print("-" * 70)
    for v in sorted(verification_lines, key=lambda x: x["energy_ref_keV"]):
        e_ref = v["energy_ref_keV"]
        rxn = v["reaction"][:20]
        ch_pred = v["channel_predicted"]
        ch_obs = v["channel_observed"]
        cnt = v["counts_observed"]
        de = v["delta_keV"]
        pk = "SI" if v["significant"] else ("no" if v["is_peak"] else "--")
        cal = "SI" if v["in_calibration"] else ""
        print(f"{e_ref:>8.1f} {rxn:<22s} {ch_pred:>7.1f} {ch_obs:>6d} {cnt:>5d} {de:>+8.3f} {pk:>5s} {cal:>3s}")

    if problems:
        print()
        print("PROBLEMAS DETECTADOS:")
        for p in problems:
            print(p)

    # --- 40K Fortran fit (detailed) ---
    half = 13
    ch_start = 5189 - half
    ch_end = 5189 + half
    roi_counts = counts[ch_start - 1:ch_end]
    pico_path = DATA_DIR / "pico.dat"
    with pico_path.open("w", encoding="utf-8") as f:
        for i, y in enumerate(roi_counts, start=1):
            f.write(f"{i} {y}\n")

    run = subprocess.run(
        [str(FORTRAN_BIN)],
        cwd=str(WORK_DIR),
        capture_output=True, text=True, check=False,
    )
    if run.returncode == 0:
        m_x0 = re.search(r"Peak position \(x0\) =\s*([\-0-9.]+)", run.stdout)
        m_s = re.search(r"Width \(sigma\)\s*=\s*([\-0-9.]+)", run.stdout)
        m_a = re.search(r"Amplitude \(A\)\s*=\s*([\-0-9.]+)", run.stdout)
        if m_x0:
            x0_local = float(m_x0.group(1))
            x0_ch = ch_start + (x0_local - 1.0)
            print(f"\nAjuste Fortran 40K:")
            print(f"  Centroide: ch {x0_ch:.4f} ({a_cal + b_cal * x0_ch:.3f} keV)")
            if m_s:
                sigma = float(m_s.group(1))
                fwhm = 2.3548 * sigma
                print(f"  Sigma: {sigma:.3f} ch, FWHM: {fwhm:.3f} ch = {fwhm * b_cal:.3f} keV")
            if m_a:
                amp = float(m_a.group(1))
                area = amp * sigma * math.sqrt(2.0 * math.pi)
                net_rate = area / live_s
                print(f"  Area neta: {area:.1f} cnt, Tasa: {net_rate:.5f} cps")

    # --- Generate plots ---
    colors = ["#888888", "#006600", "#cc0000", "#0066cc", "#cc6600", "#990099", "#669900"]

    # Full spectrum plot with predicted channel markers
    full_script = DATA_DIR / "plot_cal_full.gp"
    full_png = PLOTS_DIR / "cal_hpge_full.png"
    arrow_cmds = []
    label_cmds = []
    label_n = 1
    for v in verification_lines:
        ch_pred = v["channel_predicted"]
        e_ref = v["energy_ref_keV"]
        rxn = v["reaction"]
        c = colors[(label_n - 1) % len(colors)]
        arrow_cmds.append(
            f"set arrow {label_n} from {ch_pred},graph 0 to {ch_pred},graph 1 "
            f"nohead dt 2 lw 1.5 lc rgb '{c}'"
        )
        lbl_txt = rxn.split(",")[0] if "," in rxn else rxn
        y_pos = 0.95 - 0.040 * ((label_n - 1) % 5)
        x_off = 60 + 300 * ((label_n - 1) // 5)
        label_cmds.append(
            f"set label {label_n} '{e_ref:.1f} keV ({lbl_txt})' at "
            f"{ch_pred + x_off},graph {y_pos} "
            f"tc rgb '{c}' font 'sans,9'"
        )
        label_n += 1
    full_script.write_text(
        "set terminal pngcairo size 1600,900 enhanced font 'sans,13'\n"
        f"set output '{full_png.as_posix()}'\n"
        "set title 'Espectro HPGe Mg4 con lineas de fondo (calibracion 5 puntos)'\n"
        "set xlabel 'Canal'\n"
        "set ylabel 'Cuentas'\n"
        "set grid lw 0.5\n"
        "set key top right\n"
        "set logscale y\n"
        f"{chr(10).join(arrow_cmds)}\n"
        f"{chr(10).join(label_cmds)}\n"
        f"plot '{SPECTRUM_PATH.as_posix()}' every ::4 using 1:2 with lines lw 1 "
        f"lc rgb '#003366' title 'Espectro HPGe Mg4'\n",
        encoding="utf-8",
    )
    subprocess.run(["gnuplot", str(full_script)], cwd=str(WORK_DIR), capture_output=True, check=False)

    # Energy spectrum verification plot
    ver_script = DATA_DIR / "plot_cal_verificacion.gp"
    ver_png = PLOTS_DIR / "cal_hpge_verificacion.png"
    ver_arrows = ""
    ver_labels = ""
    vn = 1
    for v in verification_lines:
        e_ref = v["energy_ref_keV"]
        rxn = v["reaction"]
        c = colors[(vn - 1) % len(colors)]
        ver_arrows += (
            f"set arrow {vn + 10} from {e_ref},graph 0 to {e_ref},graph 1 "
            f"nohead dt 2 lw 1.5 lc rgb '{c}'\n"
        )
        lbl_v = rxn.split(",")[0]
        y_pos = 0.90 - 0.05 * vn
        x_pos = e_ref + (5 if e_ref < 1800 else -130)
        x_align = "left" if e_ref < 1800 else "right"
        ver_labels += (
            f"set label {vn + 10} '{e_ref:.1f} keV ({lbl_v})' at "
            f"{x_pos},graph {y_pos} "
            f"tc rgb '{c}' font 'sans,9' {x_align}\n"
        )
        vn += 1
    ver_script.write_text(
        "set terminal pngcairo size 1600,900 enhanced font 'sans,13'\n"
        f"set output '{ver_png.as_posix()}'\n"
        "set title 'Verificacion en energia (calibracion 5 puntos con lineas de fondo)'\n"
        "set xlabel 'Energia (keV)'\n"
        "set ylabel 'Cuentas'\n"
        "set grid lw 0.5\n"
        "set key top right\n"
        "set logscale y\n"
        "set xrange [0:2000]\n"
        f"E(ch) = {a_cal:.6f} + {b_cal:.6f} * ch\n"
        f"{ver_arrows}\n"
        f"{ver_labels}\n"
        f"plot '{SPECTRUM_PATH.as_posix()}' every ::4 using "
        f"(E($1)):2 with lines lw 1 lc rgb '#003366' "
        f"title 'Espectro HPGe Mg4'\n",
        encoding="utf-8",
    )
    subprocess.run(["gnuplot", str(ver_script)], cwd=str(WORK_DIR), capture_output=True, check=False)

    # --- Save verification table ---
    table_lines = [
        f"{'E_ref(keV)':<12s} {'Reaccion':<28s} {'Ch_pred':>7s} {'Ch_obs':>7s} "
        f"{'Cuentas':>7s} {'dE(keV)':>8s} {'Pico':>5s}"
    ]
    table_lines.append("-" * 75)
    for v in sorted(verification_lines, key=lambda x: x["energy_ref_keV"]):
        e_ref = v["energy_ref_keV"]
        rxn = v["reaction"][:26]
        ch_pred = v["channel_predicted"]
        ch_obs = v["channel_observed"]
        cnt = v["counts_observed"]
        de = v["delta_keV"]
        pk = "SI" if v["significant"] else ("no" if v["is_peak"] else "--")
        table_lines.append(
            f"{e_ref:<12.1f} {rxn:<28s} {ch_pred:>7.1f} "
            f"{ch_obs:>7d} {cnt:>7d} {de:>+8.3f} {pk:>5s}"
        )
    table_path = DATA_DIR / "tabla_verificacion.txt"
    table_path.write_text("\n".join(table_lines) + "\n", encoding="utf-8")

    # --- Save JSON with calibration results ---
    json_path = DATA_DIR / "resultados_k40_hpge.json"
    json_data = {
        "calibracion": {
            "pendiente_keV_por_ch": b_cal,
            "offset_keV": a_cal,
            "R2": r2,
            "n_puntos": len(ANCHORS),
            "anclas": [{"canal": ch, "energia_keV": e} for ch, e in ANCHORS],
        },
        "espectro": {
            "archivo": "Mg4_GeHP.txt",
            "tiempo_vivo_s": live_s,
            "canales": n,
        },
        "pico_K40": {
            "canal_aproximado": 5189,
            "energia_keV": 1460.820,
        },
        "verificacion": [
            {
                "energia_ref_keV": v["energy_ref_keV"],
                "reaccion": v["reaction"],
                "canal_predicho": v["channel_predicted"],
                "canal_observado": v["channel_observed"],
                "cuentas": v["counts_observed"],
                "delta_keV": v["delta_keV"],
                "es_pico": v["significant"],
            }
            for v in verification_lines
        ],
        "problemas": problems,
    }
    with json_path.open("w", encoding="utf-8") as f:
        json.dump(json_data, f, indent=2, ensure_ascii=False)

    print()
    print("Archivos generados:")
    print(f"  Espectro completo: {full_png}")
    print(f"  Verificacion: {ver_png}")
    print(f"  Tabla: {table_path}")
    print(f"  JSON: {json_path}")


if __name__ == "__main__":
    main()
