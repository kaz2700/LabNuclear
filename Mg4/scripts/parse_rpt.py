import zipfile
from pathlib import Path

def parse_rpt(filepath, from_zip=None):
    if from_zip:
        with zipfile.ZipFile(from_zip) as zf:
            text = zf.read(filepath).decode('latin-1')
    else:
        text = Path(filepath).read_text(encoding='latin-1')

    lines = text.splitlines()

    meta = {
        "live_time": None,
        "real_time": None,
        "comment": None,
    }

    channels = []
    data_start = False

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("Live time"):
            val = stripped.split(":", 1)[1].strip()
            try:
                meta["live_time"] = float(val)
            except ValueError:
                pass
            continue

        if stripped.startswith("Real Time"):
            val = stripped.split(":", 1)[1].strip()
            try:
                meta["real_time"] = float(val)
            except ValueError:
                pass
            continue

        if "Canal" in stripped or "Cuentas" in stripped:
            data_start = True
            meta["comment"] = stripped.strip()
            continue

        if not data_start:
            continue

        parts = stripped.split()
        if len(parts) < 2:
            continue

        try:
            ch = int(parts[0])
            cnt = int(parts[1])
            channels.append((ch, cnt))
        except ValueError:
            continue

    counts = [c for _, c in channels]
    first_ch = channels[0][0] if channels else None

    meta["first_channel"] = first_ch
    meta["num_channels"] = len(counts)

    return meta, counts
