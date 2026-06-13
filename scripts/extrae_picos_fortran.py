import os
import numpy as np

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'datos')

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
        if s.startswith('WAVES') and 'MCAch0' in s:
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
    if not os.path.exists(path):
        raise FileNotFoundError(f'{path} not found')
    meta, counts = parse_itx(path)
    channels = np.arange(len(counts), dtype=float)
    return meta, channels, counts

# Peak centers from Python analysis
peak_info = [
    ('22Na-157-coin', '511_157',  66,   90,  78.05),
    ('22Na-157-coin', '1274_157', 188, 210, 196.45),
    ('22Na-180-coin', '511_180',  66,   92,  78.39),
    ('22Na-180-coin', 'pico_200_180', 190, 215, 196.81),
]

os.makedirs('data', exist_ok=True)
for fname, label, xmin, xmax, x0 in peak_info:
    meta, ch, counts = cargar_espectro(fname)
    mask = (ch >= xmin) & (ch <= xmax)
    x = ch[mask]
    y = counts[mask]
    # Find index of peak center, reindex so peak is at position 15
    idx = np.argmin(np.abs(x - x0))
    target = 15
    # Shift so peak lands at position 15
    x_reindex = np.arange(1, len(x) + 1) + (target - 1 - idx)
    outpath = os.path.join('data', f'pico_{label}.dat')
    with open(outpath, 'w') as f:
        for xi, yi in zip(x_reindex, y):
            f.write(f'{xi} {int(yi)}\n')
    print(f'{outpath}: {len(x)} pts, ymax={y.max()}, reindexed_peak_pos={x_reindex[idx]}')
