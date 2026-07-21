#!/usr/bin/env python3
"""Embed the three generated figures into the standalone report HTML."""

import base64
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FIGURES = {
    "{{FIG1}}": ROOT / "figs" / "fig1_mechanism.png",
    "{{FIG2}}": ROOT / "figs" / "fig2_retention.png",
    "{{FIG_R1}}": ROOT / "figs" / "figR1_a549_entanglement.png",
}


def main() -> None:
    html = (ROOT / "qc_page.html").read_text()
    for token, figure_path in FIGURES.items():
        encoded = base64.b64encode(figure_path.read_bytes()).decode()
        html = html.replace(token, f"data:image/png;base64,{encoded}")

    unresolved = [token for token in FIGURES if token in html]
    if unresolved:
        raise RuntimeError(f"Unresolved figure placeholders: {unresolved}")

    output = ROOT / "pairwise-qc-recommendation.html"
    output.write_text(html)
    print(f"Wrote {output} ({output.stat().st_size / 1024:.0f} KiB)")


if __name__ == "__main__":
    main()
