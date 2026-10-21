"""Report the printed size of each figure at \\linewidth, against the text block."""
from pathlib import Path

from PIL import Image

TEXT_W_CM = 21.0 - 2 * 2.4      # a4 width minus the 2.4 cm margins
TEXT_H_CM = 29.7 - 2 * 2.4      # a4 height minus the 2.4 cm margins
LINE_W_PT = TEXT_W_CM / 2.54 * 72.27
PAGE_H_PT = (TEXT_H_CM - 3.0) / 2.54 * 72.27   # room for the caption

for path in sorted(Path("bench/figures").glob("*.png")):
    with Image.open(path) as im:
        w, h = im.size
        dpi = im.info.get("dpi", (None, None))
    height_pt = LINE_W_PT * h / w
    # the tallest width at which the figure still fits on one page with its caption
    fit = min(1.0, 0.9 * PAGE_H_PT * w / (LINE_W_PT * h))
    print(f"{path.name:28s} {w:5d}x{h:5d} dpi={dpi}  height at linewidth = "
          f"{height_pt:6.1f} pt (page allows {PAGE_H_PT:.0f}); max width = {fit:.2f}")
