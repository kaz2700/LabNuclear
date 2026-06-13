import zipfile
import sys
from pathlib import Path

from parse_chn import parse_chn
from parse_rpt import parse_rpt
from to_itx import rpt_to_itx


def find_zip():
    candidates = [
        Path("cosas_chnquenosabemoscmoabrir.zip"),
        Path(__file__).resolve().parent.parent / "cosas_chnquenosabemoscmoabrir.zip",
        Path.cwd() / "cosas_chnquenosabemoscmoabrir.zip",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


ZIP_PATH = find_zip()


def fmt_dt(dt, val, unit=""):
    if val is None:
        return "?"
    if isinstance(val, float) and val == int(val):
        val = int(val)
    return f"{val}{unit}"


def show_chn_analysis(name):
    if not ZIP_PATH:
        return {}, []
    meta, counts = parse_chn(name, from_zip=str(ZIP_PATH))

    print(f"  File:        {name}")
    print(f"  Size:        {meta['file_size']} bytes")
    if meta.get("empty"):
        print(f"  Status:      EMPTY (all zeros)")
        return meta, counts

    print(f"  Device:      {meta['device'] or 'unknown'}")
    print(f"  Detector:    {meta['detector'] or 'unknown'}")
    print(f"  Date:        {meta['date_str']}  →  {meta['date'] or '?'}")
    print(f"  Channels:    {meta['num_channels']}")
    print(f"  Live time:   {fmt_dt(meta['date'], meta['live_time_raw'])} (raw units)")
    print(f"  Real time:   {fmt_dt(meta['date'], meta['real_time_raw'])} (raw units)")
    print(f"  Data range:  {min(counts)} – {max(counts)}")

    return meta, counts


def show_rpt_analysis(name):
    if not ZIP_PATH:
        return {}, []
    meta, counts = parse_rpt(name, from_zip=str(ZIP_PATH))

    lv = meta["live_time"]
    rv = meta["real_time"]
    print(f"  File:        {name}")
    print(f"  Channels:    {meta['num_channels']} ({meta['first_channel']}-based)")
    print(f"  Live time:   {lv} s")
    print(f"  Real time:   {rv} s")
    if lv and rv:
        print(f"  Dead time:   {100 * (1 - lv / rv):.1f}%")

    total = sum(counts)
    rate = total / lv if lv else 0
    print(f"  Total:       {total} counts  ({rate:.1f} cps)")

    peaks = find_peaks(counts, min_height=200)
    print(f"  Peaks:       {len(peaks)} found (≥200 cnt)")
    for ch, cnt in peaks[:8]:
        print(f"    Ch {ch:>5d}  {cnt:>6d}")
    if len(peaks) > 8:
        print(f"    ... +{len(peaks) - 8} more")

    return meta, counts


def find_peaks(counts, min_height=0, window=5):
    n = len(counts)
    peaks = []
    for i in range(window, n - window):
        c = counts[i]
        if c < min_height:
            continue
        section = counts[i - window:i + window + 1]
        if c == max(section) and c > counts[i - 1] and c > counts[i + 1]:
            peaks.append((i + 1, c))
    peaks.sort(key=lambda x: x[1], reverse=True)
    return peaks


def main():
    if not ZIP_PATH:
        print("ERROR: Zip file not found.")
        sys.exit(1)

    print(f"Zip: {ZIP_PATH.resolve()}\n")

    with zipfile.ZipFile(str(ZIP_PATH)) as zf:
        for info in zf.infolist():
            print(f"  {info.filename:<30s} {info.file_size:>8,} B")

    print("\n" + "─" * 50)
    print("  .Chn ─ binary files (Easy-MCA-8k)")
    print("─" * 50)
    for name in ["fondodia3.Chn", "dia2.Chn"]:
        show_chn_analysis(name)
        print()

    print("─" * 50)
    print("  .RPT ─ text report")
    print("─" * 50)
    show_rpt_analysis("raulceciliapablo.RPT")
    print()

    meta_ch2, cnt_ch2 = parse_chn("dia2.Chn", from_zip=str(ZIP_PATH))
    _, cnt_rpt = parse_rpt("raulceciliapablo.RPT", from_zip=str(ZIP_PATH))

    if cnt_ch2 and cnt_rpt:
        ch2_range = max(cnt_ch2) - min(cnt_ch2)
        rpt_range = max(cnt_rpt) - min(cnt_rpt)

        match_ratio = sum(
            1 for i, c in enumerate(cnt_ch2)
            if i < len(cnt_rpt) and abs(c - cnt_rpt[i]) < 0.1 * max(cnt_rpt + [1])
        ) / len(cnt_ch2) if cnt_ch2 else 0

        print(f"  NOTE: The binary dia2.Chn data values ({meta_ch2['num_channels']} × 32-bit)")
        print(f"        range {min(cnt_ch2)}–{max(cnt_ch2)}, which is incompatible")
        print(f"        with the RPT counts (range {min(cnt_rpt)}–{max(cnt_rpt)}).")
        print(f"        These are likely DIFFERENT measurements or the binary stores")
        print(f"        calibration/config data rather than raw channel counts.")
        print()

    print("─" * 50)
    print("  Conversion to .itx")
    print("─" * 50)
    itx_path = Path(__file__).parent / "raulceciliapablo_converted.itx"
    rpt_to_itx("raulceciliapablo.RPT", outpath=str(itx_path),
               source_name="raulceciliapablo", from_zip=str(ZIP_PATH))
    print(f"  RPT  →  {itx_path.name}  ({itx_path.stat().st_size:,} bytes)")
    print(f"  Compatible with existing scripts in ../parte2/")

    print("\nDone.")


if __name__ == "__main__":
    main()
