import struct
import zipfile
from pathlib import Path
from datetime import datetime

HEADER_SIZE = 32
MAGIC = 0xFFFF
CAL_MARKER = 0xFF9A


def _parse_header(raw):
    magic, mca_num, segment = struct.unpack_from('<HHH', raw, 0)
    if magic != MAGIC:
        return None

    secs_bytes = raw[6:8]
    secs_str = secs_bytes.decode('ascii', errors='replace')
    realtime_ticks, livetime_ticks = struct.unpack_from('<II', raw, 8)
    date_str = raw[16:24].decode('ascii', errors='replace').strip()
    time_str = raw[24:28].decode('ascii', errors='replace').strip()
    chan_offset, num_chans = struct.unpack_from('<hh', raw, 28)

    full_date_str = f"{date_str} {time_str}:{secs_str}" if secs_str.isdigit() else f"{date_str} {time_str}"

    dt = None
    try:
        dt = datetime.strptime(full_date_str, "%d%b%Y %H:%M:%S")
    except ValueError:
        try:
            dt = datetime.strptime(full_date_str, "%d%b%y %H:%M:%S")
        except ValueError:
            pass

    hdr = {
        "mca_num": mca_num,
        "segment": segment,
        "date_str": full_date_str,
        "date": dt,
        "realtime_s": realtime_ticks * 0.020,
        "livetime_s": livetime_ticks * 0.020,
        "realtime_ticks": realtime_ticks,
        "livetime_ticks": livetime_ticks,
        "chan_offset": chan_offset,
        "num_chans": num_chans,
    }
    return hdr


def _parse_calibration(raw, cal_offset):
    marker = struct.unpack_from('<H', raw, cal_offset)[0]
    if marker != CAL_MARKER:
        return None
    zero, slope, quad = struct.unpack_from('<fff', raw, cal_offset + 4)
    return {"zero": zero, "slope": slope, "quad": quad}


def _find_device_string(raw, start):
    idx = raw.find(b"'55438C76E7C84EC", start)
    if idx < 0:
        idx = raw.find(b"Easy-MCA", start)
    if idx < 0:
        idx = raw.find(b"C060", start)
    if idx < 0:
        return ""
    return raw[idx:].split(b'\x00')[0].decode('ascii', errors='replace').strip()


def parse_chn(filepath, from_zip=None):
    if from_zip:
        with zipfile.ZipFile(from_zip) as zf:
            raw = zf.read(filepath)
    else:
        raw = Path(filepath).read_bytes()

    if len(raw) < HEADER_SIZE:
        raise ValueError(f"File too small: {len(raw)} bytes")

    hdr = _parse_header(raw)
    if hdr is None:
        is_empty = all(b == 0 for b in raw)
        meta = {
            "mca_num": 0, "segment": 0,
            "date": None, "date_str": "", "realtime_s": 0, "livetime_s": 0,
            "device": "", "calibration": None,
            "num_channels": 0, "file_size": len(raw), "empty": is_empty,
        }
        return meta, []

    num_chans = hdr["num_chans"]
    data_start = HEADER_SIZE
    data_end = data_start + num_chans * 4

    data = list(struct.iter_unpack('<I', raw[data_start:data_end]))
    counts = [v[0] for v in data]

    cal = _parse_calibration(raw, data_end)

    device = _find_device_string(raw, data_end)

    meta = {
        "mca_num": hdr["mca_num"],
        "segment": hdr["segment"],
        "date": hdr["date"],
        "date_str": hdr["date_str"],
        "realtime_s": hdr["realtime_s"],
        "livetime_s": hdr["livetime_s"],
        "realtime_ticks": hdr["realtime_ticks"],
        "livetime_ticks": hdr["livetime_ticks"],
        "chan_offset": hdr["chan_offset"],
        "device": device,
        "calibration": cal,
        "num_channels": num_chans,
        "file_size": len(raw),
        "empty": False,
    }

    return meta, counts
