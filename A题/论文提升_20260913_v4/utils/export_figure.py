from __future__ import annotations

import argparse
import os
import sys
from typing import Iterable

import matplotlib.pyplot as plt


VECTOR_FORMATS = {"pdf", "svg", "eps"}
RASTER_FORMATS = {"png", "tiff", "tif", "jpg", "jpeg"}
SUPPORTED_FORMATS = VECTOR_FORMATS | RASTER_FORMATS


def _ensure_parent(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent and not os.path.exists(parent):
        os.makedirs(parent, exist_ok=True)


def export_figure(
    fig,
    basename: str,
    formats: Iterable[str] | None = None,
    dpi: int = 300,
    size_inches: tuple[float, float] | None = None,
    grayscale_preview: bool = False,
    tight: bool = True,
    pad_inches: float = 0.05,
    transparent: bool = False,
) -> list[str]:
    if formats is None:
        formats = ("pdf", "svg", "png")
    formats = [f.lower().lstrip(".") for f in formats]
    unknown = [f for f in formats if f not in SUPPORTED_FORMATS]
    if unknown:
        raise ValueError(f"Unsupported formats: {unknown}. "
                         f"Supported: {sorted(SUPPORTED_FORMATS)}")

    if size_inches is not None:
        if len(size_inches) != 2:
            raise ValueError("size_inches must be (width, height)")
        fig.set_size_inches(*size_inches)

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"  

    saved: list[str] = []
    for fmt in formats:
        if fmt in {"jpg", "jpeg"}:
            print(f"[scipilot-figure-skill] WARNING: skipping {fmt} — "
                  "JPEG is lossy and unsuitable for line/text figures.",
                  file=sys.stderr)
            continue
        path = f"{basename}.{fmt}"
        _ensure_parent(path)
        kwargs: dict = {
            "bbox_inches": "tight" if tight else None,
            "pad_inches": pad_inches,
            "transparent": transparent,
        }
        if fmt in RASTER_FORMATS:
            kwargs["dpi"] = dpi
        fig.savefig(path, **kwargs)
        saved.append(path)
        print(f"[scipilot-figure-skill] wrote {path}")

    if grayscale_preview:
        gray_path = _grayscale_from(fig, basename, dpi=dpi)
        if gray_path:
            saved.append(gray_path)
    return saved


def _grayscale_from(fig, basename: str, dpi: int) -> str | None:
    try:
        from PIL import Image
    except ImportError:
        print("[scipilot-figure-skill] Pillow not available; "
              "grayscale preview skipped.", file=sys.stderr)
        return None

    png_path = f"{basename}.png"
    _ensure_parent(png_path)
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")

    gray_path = f"{basename}_grayscale.png"
    Image.open(png_path).convert("L").save(gray_path)
    print(f"[scipilot-figure-skill] wrote {gray_path} (grayscale preview)")
    return gray_path


def _demo(out_basename: str) -> None:
    import numpy as np
    rng = np.random.default_rng(7)
    x = np.linspace(0, 10, 50)
    y1 = np.sin(x) + rng.normal(0, 0.1, x.size)
    y2 = np.cos(x) + rng.normal(0, 0.1, x.size)

    fig, ax = plt.subplots(figsize=(3.5, 2.625))
    ax.plot(x, y1, label="sin", marker="o", markersize=3)
    ax.plot(x, y2, label="cos", marker="s", markersize=3, linestyle="--")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend(frameon=False)

    paths = export_figure(
        fig, out_basename,
        formats=["pdf", "svg", "png"],
        size_inches=(3.5, 2.625),
        dpi=300,
        grayscale_preview=True,
    )
    print("\nDemo done. Files:")
    for p in paths:
        print(f"  {p}")


def _cli() -> int:
    p = argparse.ArgumentParser(description="scipilot-figure-skill figure exporter")
    p.add_argument("cmd", choices=["demo"], help="`demo`: 跑一张演示图导出 4 种格式")
    p.add_argument("--out", default="./scipilot_demo",
                   help="输出 basename (默认 ./scipilot_demo)")
    args = p.parse_args()
    if args.cmd == "demo":
        _demo(args.out)
    return 0


if __name__ == "__main__":
    sys.exit(_cli())
