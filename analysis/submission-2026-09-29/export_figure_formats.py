#!/usr/bin/env python3
"""Recreate the numerical plots and export lossless 600-dpi TIFFs for submission.

Use the analysis environment. This invokes the existing figure-producing
drivers, with their fixed analysis specifications, and does not edit pixels.
"""
from pathlib import Path
import os
import runpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
os.chdir(ROOT)
import sys
sys.path.insert(0, str(ROOT))   # drivers import local modules such as hgvs
for number, driver in (("Schematic", "run_schematic.py"), (1, "run_scn5a.py"), (2, "run_figure2.py"), (3, "run_bounds.py"), ("Titration", "run_titration.py"), ("TitrationMain", "run_titration_main.py")):
    figures = []
    original_close = plt.close
    plt.close = lambda fig=None: figures.append(fig) if hasattr(fig, "savefig") else None
    try:
        namespace = runpy.run_path(str(ROOT / driver), run_name="__main__")
        fig = namespace.get("fig") or figures[-1]
        stem = number if isinstance(number, str) else f"Figure{number}"
        target = ROOT / f"{stem}.tiff"
        fig.savefig(target, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
        # PLOS requires flattened RGB TIFFs without an alpha channel.
        from PIL import Image
        with Image.open(target) as image:
            rgba = image.convert("RGBA")
        flat = Image.new("RGB", rgba.size, "white")
        flat.paste(rgba, mask=rgba.split()[3])
        flat.save(target, compression="tiff_lzw", dpi=(600, 600))
    finally:
        plt.close = original_close
        plt.close("all")
    print(f"Exported {stem}.tiff at 600 dpi (RGB, LZW) from plotting code")
