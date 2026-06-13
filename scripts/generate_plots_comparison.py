#!/usr/bin/env python3
"""
Generate comparison plots: MCAch0 (fixed detector, red) vs MCAch1 (mobile, green).
Output: plots/comparison/
"""
import os, sys, math, glob
import numpy as np

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots', 'comparison')
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'datos')
os.makedirs(OUT_DIR, exist_ok=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10, 'figure.dpi': 150})

COLORS = {'ch0': '#d62728', 'ch1': '#2ca02c'}

def parse_itx_ch(filepath, ch='MCAch0'):
    counts = []
    with open(filepath, 'r', errors='replace') as f:
        lines = f.readlines()
    in_ch = False
    for line in lines:
        s = line.strip()
        if s.startswith('WAVES') and ch in s:
            in_ch = True
            continue
        if in_ch:
            if s == 'BEGIN': continue
            if s == 'END': break
            try:
                counts.append(int(s))
            except ValueError:
                pass
    return np.array(counts, dtype=float)

def make_comparison_plot(counts0, counts1, base, label):
    fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(14, 8),
                                             gridspec_kw={'height_ratios': [3, 1]})

    ax_top.step(np.arange(len(counts0)), counts0, where='mid',
                color=COLORS['ch0'], linewidth=0.6, alpha=0.8,
                label='MCAch0 (fijo, rojo)')
    ax_top.step(np.arange(len(counts1)), counts1, where='mid',
                color=COLORS['ch1'], linewidth=0.6, alpha=0.8,
                label='MCAch1 (móvil, verde)')

    nonzero = np.where(counts0 + counts1 > 0)[0]
    if len(nonzero) > 0:
        ax_top.set_xlim(max(0, nonzero[0] - 10), min(len(counts0), nonzero[-1] + 10))
    ax_top.set_ylabel('Cuentas', fontsize=11)
    ax_top.set_title(f'{label} — Comparación MCAch0 vs MCAch1', fontsize=12)
    ax_top.grid(True, alpha=0.3)
    ax_top.legend(fontsize=9)

    ratio = np.where(counts0 > 0, counts1 / np.maximum(counts0, 1), 0)
    ax_bottom.step(np.arange(len(ratio)), ratio, where='mid',
                   color='purple', linewidth=0.6, alpha=0.8)
    ax_bottom.axhline(1.0, color='gray', linestyle='--', linewidth=0.5)
    ax_bottom.set_xlim(ax_top.get_xlim())
    ax_bottom.set_xlabel('Canal (raw ADC)', fontsize=11)
    ax_bottom.set_ylabel('Ratio ch1/ch0', fontsize=11)
    ax_bottom.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, f'{base}_comparison.png'), dpi=150)
    plt.close()

def make_overlay_comparison_plot(counts0, counts1, base, label):
    """Overlay plot with zoom on peak region for Co-60 spectra."""
    fig, ax = plt.subplots(figsize=(14, 7))

    ax.step(np.arange(len(counts0)), counts0, where='mid',
            color=COLORS['ch0'], linewidth=0.6, alpha=0.8,
            label='MCAch0 (fijo, rojo — peor resolución)')
    ax.step(np.arange(len(counts1)), counts1, where='mid',
            color=COLORS['ch1'], linewidth=0.6, alpha=0.8,
            label='MCAch1 (móvil, verde — mejor resolución)')

    # Zoom on the peak region (~ch 180-230 for 1332 keV)
    ax.set_xlim(170, 240)
    ax.set_xlabel('Canal (raw ADC)', fontsize=11)
    ax.set_ylabel('Cuentas', fontsize=11)
    ax.set_title(f'{label} — Comparación directa, región del pico de 1332 keV', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, f'{base}_overlay.png'), dpi=150)
    plt.close()

# ─── Main ──────────────────────────────────────────────────────────────────

itx_files = sorted(glob.glob(os.path.join(DATA_DIR, '*.itx')))

for itx_path in itx_files:
    base = os.path.splitext(os.path.basename(itx_path))[0]
    print(f'Processing {base}...')

    counts0 = parse_itx_ch(itx_path, 'MCAch0')
    counts1 = parse_itx_ch(itx_path, 'MCAch1')

    if len(counts0) != len(counts1):
        print(f'  Warning: length mismatch ch0={len(counts0)}, ch1={len(counts1)}')
        continue

    # Full spectrum comparison
    make_comparison_plot(counts0, counts1, base, base)

    # Overlay zoom for Co-60 spectra (have clear 1332 peak)
    if 'Co' in base and 'caracterización' in base.lower() or (
        'Co' in base and any(a in base for a in ['22,5', '45', '67,5', '90', '112,5', '135', '157,5', '180'])
    ):
        make_overlay_comparison_plot(counts0, counts1, base, base)

print(f'\nAll comparison plots saved to: {OUT_DIR}')
