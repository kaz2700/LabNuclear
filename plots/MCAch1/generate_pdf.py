#!/usr/bin/env python3
"""Generate PDF with all MCAch1 plots."""
import os, glob
os.environ['MPLBACKEND'] = 'Agg'
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

DIR = os.path.dirname(os.path.abspath(__file__))
pngs = sorted(glob.glob(os.path.join(DIR, '*.png')))

with PdfPages(os.path.join(DIR, 'plots_MCAch1.pdf')) as pdf:
    for path in pngs:
        img = plt.imread(path)
        fig, ax = plt.subplots(figsize=(14, 8))
        ax.imshow(img)
        ax.axis('off')
        fig.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
        print(f'  Added: {os.path.basename(path)}')

print(f'\nPDF saved: {os.path.join(DIR, "plots_MCAch1.pdf")}')
