from pathlib import Path
from parse_chn import parse_chn
from parse_rpt import parse_rpt

HERE = Path(__file__).parent
ZIP = HERE / "cosas_chnquenosabemoscmoabrir.zip"


def chn_to_txt(name, outname):
    meta, counts = parse_chn(name, from_zip=str(ZIP))
    if meta["empty"]:
        out = HERE / outname
        out.write_text("# vacio (todo ceros)\n")
        print(f"  {outname}  — vacio")
        return

    cal = meta["calibration"]
    slope = cal["slope"] if cal else 1.0
    zero = cal["zero"] if cal else 0.0
    livetime_s = meta["livetime_s"]
    device = meta["device"] or "desconocido"

    lines = [
        f"# {name}  |  {device}  |  {meta['date_str']}",
        f"# Live: {livetime_s:.0f} s  |  Real: {meta['realtime_s']:.0f} s",
        f"# Cal: {zero:.2f} + {slope:.4f} * canal  + {cal['quad']} * canal^2" if cal else "# Cal: (sin calibracion)",
        f"# Total cuentas: {sum(counts):,}  ({sum(counts)/livetime_s:.1f} cps)" if livetime_s else "",
        "#canal\tcuentas\tenergia_keV",
    ]

    for i, c in enumerate(counts, 1):
        e = zero + slope * i + (cal["quad"] if cal else 0) * i * i
        lines.append(f"{i}\t{c}\t{e:.2f}")

    out = HERE / outname
    out.write_text("\n".join(lines) + "\n")
    print(f"  {outname}  — {len(counts)} canales, {sum(counts):,} cuentas")


def rpt_to_txt(name, outname):
    meta, counts = parse_rpt(name, from_zip=str(ZIP))
    lv = meta["live_time"] or 0
    lines = [
        f"# {name}  |  Easy-MCA-8k (RPT)",
        f"# Live: {lv:.0f} s  |  Real: {meta['real_time']:.0f} s" if meta["real_time"] else "# Live/Real: ?",
        f"# Total cuentas: {sum(counts):,}  ({sum(counts)/lv:.1f} cps)" if lv else "",
        "#canal\tcuentas",
    ]
    for i, c in enumerate(counts, meta["first_channel"] or 1):
        lines.append(f"{i}\t{c}")

    out = HERE / outname
    out.write_text("\n".join(lines) + "\n")
    print(f"  {outname}  — {len(counts)} canales, {sum(counts):,} cuentas")


def main():
    print("Convirtiendo .Chn → .txt ...")
    chn_to_txt("dia2.Chn", "dia2_Chn_256ch_espectro.txt")
    chn_to_txt("fondodia3.Chn", "fondodia3_Chn_vacio.txt")

    print("Convirtiendo .RPT → .txt ...")
    rpt_to_txt("raulceciliapablo.RPT", "raulceciliapablo_RPT_16384ch_espectro.txt")

    print("\nListo.")


if __name__ == "__main__":
    main()
