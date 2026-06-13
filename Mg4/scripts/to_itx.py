from datetime import datetime
from pathlib import Path
from parse_rpt import parse_rpt
from parse_chn import parse_chn

def counts_to_itx(counts, live_time=None, real_time=None, label="MCAch0",
                  source_name="unknown", device="Easy-MCA-8k", outpath=None):
    now = datetime.now().strftime("%a, %d de %b de %Y, %H:%M:%S")
    nch = len(counts)

    lines = ["IGOR"]
    lines.append(f"X // {device} data converted {now}")
    if live_time is not None:
        lines.append(f"X // Live Time [s]= {live_time}")
    if real_time is not None:
        lines.append(f"X // Real Time [s]= {real_time}")
    lines.append(f"X // Source: {source_name}")
    lines.append("")
    lines.append(f"WAVES\t{label}")
    lines.append("BEGIN")
    for c in counts:
        lines.append(f"\t{int(c)}")
    lines.append("END")

    text = "\n".join(lines) + "\n"

    if outpath:
        Path(outpath).write_text(text)

    return text


def rpt_to_itx(rpt_path, outpath=None, label="MCAch0", source_name="unknown",
               from_zip=None):
    meta, counts = parse_rpt(rpt_path, from_zip=from_zip)
    return counts_to_itx(
        counts,
        live_time=meta["live_time"],
        real_time=meta["real_time"],
        label=label,
        source_name=source_name,
        device="Easy-MCA-8k",
        outpath=outpath,
    )


def chn_to_itx(chn_path, outpath=None, label="MCAch0", source_name="unknown",
               from_zip=None):
    meta, counts = parse_chn(chn_path, from_zip=from_zip)
    return counts_to_itx(
        counts,
        live_time=meta["live_time_raw"],
        real_time=meta["real_time_raw"],
        label=label,
        source_name=source_name,
        device=meta.get("device", "Easy-MCA-8k"),
        outpath=outpath,
    )
