"""Draw Figure 5: a screenshot of the web interface, produced by the browser.

    python bench/fig_ui_screenshot.py --pdb 5I6X

The guidelines ask for a figure of the interface, and a photograph of a browser window is
not reproducible: it depends on the fonts that happen to be installed and on the width of the
window.  This script drives a headless Chromium through Playwright against the live
application, waits until the topology it asked for has actually arrived, and saves the page at
a fixed width.  The screenshot is taken from the live site because the figure is a picture of
the deployed software, and the caption says so.

Two waits are used, and both are needed.  The page reports "Computing ..." while a request is
in flight, but the text alone is not enough: the status line and the spinner keep their
wording for a moment after the response lands, and the viewer below them is drawn afterwards.
The figure waits, first, for the "Computing ..." label to disappear and, second, until the
page has stopped growing, which is when the 3-D view and the summary have been laid out.

The interface is dark; the figure is therefore not re-coloured.  A white margin is added and
a thin rule is drawn round it, so that the panel reads as a screenshot on a page rather than
as a printing fault.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image, ImageDraw

WIDTH = 1440
SCALE = 2                       # the device scale factor of the browser
MARGIN = 28                     # white margin added round the screenshot, at 2x


def viewer_drawn(page) -> bool:
    """Is the three-dimensional view actually painted?

    The ribbon is drawn by 3dmol.js on a WebGL canvas.  A headless browser can report the
    canvas as present and still leave it blank, so the check is made on the pixels: the
    viewer panel is white, and a ribbon of a few thousand atoms always changes that.
    """
    js = ("() => { const c = [...document.querySelectorAll('canvas')]"
          ".find(el => el.getBoundingClientRect().height > 150);"
          " if (!c) return 0;"
          " try { const t = document.createElement('canvas');"
          "   t.width = 120; t.height = 60;"
          "   const g = t.getContext('2d');"
          "   g.drawImage(c, 0, 0, 120, 60);"
          "   const d = g.getImageData(0, 0, 120, 60).data;"
          "   let n = 0; for (let i = 0; i < d.length; i += 4)"
          "     if (d[i] < 240 || d[i + 1] < 240 || d[i + 2] < 240) n++;"
          "   return n / (120 * 60); } catch (e) { return -1; } }")
    return (page.evaluate(js) or 0) > 0.02


def settle(page) -> None:
    """Wait until the requested topology has been drawn."""
    page.wait_for_function(
        "() => /Computing/i.test(document.body.innerText)", timeout=60_000)
    page.wait_for_function(
        "() => !/Computing/i.test(document.body.innerText)", timeout=300_000)
    # the viewer and the summary are laid out after the response: the page stops growing
    page.wait_for_timeout(2_000)
    height = page.evaluate("document.body.scrollHeight")
    for _ in range(20):
        page.wait_for_timeout(1_000)
        grown = page.evaluate("document.body.scrollHeight")
        if grown == height:
            break
        height = grown


def shots(browser, url: str, query: str, out: Path, tries: int = 4) -> bool:
    """Capture until the three-dimensional view is painted; return whether it was."""
    for attempt in range(1, tries + 1):
        page = browser.new_page(viewport={"width": WIDTH, "height": 900},
                                device_scale_factor=SCALE)
        page.goto(url, wait_until="networkidle", timeout=120_000)
        page.locator("input").first.fill(query)
        page.get_by_role("button", name="Fetch").click()
        settle(page)
        page.wait_for_timeout(6_000)     # the ribbon is drawn asynchronously
        page.screenshot(path=str(out), full_page=True)
        painted = viewer_drawn(page)
        page.close()
        if painted:
            return True
        print(f"  attempt {attempt}: the 3-D view came back blank; retrying")
    return False


def capture(url: str, query: str, out: Path) -> None:
    from playwright.sync_api import sync_playwright

    with sync_playwright() as p:
        browser = p.chromium.launch()
        drawn = shots(browser, url, query, out)
        browser.close()
    if not drawn:
        raise SystemExit(
            "the three-dimensional view never painted: run the script again, or take the\n"
            "screenshot of the 2-D map alone and say so in the legend")


def frame(path: Path) -> None:
    """White margin and a hairline round the screenshot."""
    shot = Image.open(path).convert("RGB")
    canvas = Image.new("RGB", (shot.width + 2 * MARGIN, shot.height + 2 * MARGIN),
                       (255, 255, 255))
    canvas.paste(shot, (MARGIN, MARGIN))
    draw = ImageDraw.Draw(canvas)
    draw.rectangle([MARGIN - 1, MARGIN - 1, MARGIN + shot.width, MARGIN + shot.height],
                   outline=(190, 196, 204), width=1)
    canvas.save(path, dpi=(300, 300))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default="https://protein-visualize.vercel.app")
    ap.add_argument("--pdb", default="5I6X")
    ap.add_argument("--out", default="bench/figures/figure5_interface.png")
    args = ap.parse_args()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    capture(args.url, args.pdb, out)
    frame(out)
    print(f"  {out}: {Image.open(out).size[0]}x{Image.open(out).size[1]} px, "
          f"screenshot of {args.url} with PDB {args.pdb}")


if __name__ == "__main__":
    main()
