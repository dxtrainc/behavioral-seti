#!/bin/bash
# Build the submission PDF from the paper's HTML rendering.
#
#   ./build.sh <source-index.html> [output.pdf]
#
# There is no LaTeX or pandoc anywhere on the fleet, so the PDF is produced by
# printing the HTML rendering through headless Chromium or Chrome.  Needs one of them
# and the four Charis TTFs in fonts/ (see fetch-fonts.sh).
#
# Two snap-confinement rules bite here, and both fail as a one-page "file
# couldn't be accessed" PDF rather than as an error:
#   * Chromium cannot read outside $HOME      -> a /tmp source gives ERR_FILE_NOT_FOUND
#   * nor can it read dot-directories in $HOME -> a work dir like ~/.build gives
#     ERR_ACCESS_DENIED
# So the work directory is under $HOME and is NOT hidden. Keep it that way.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="${1:?usage: build.sh <source-index.html> [output.pdf]}"
OUT="${2:-$HOME/beacon-paper-submission.pdf}"
WORK="${WORK:-$HOME/beacon-pdfbuild-work}"   # must be under $HOME and not hidden

# Chromium was renamed/removed more than once on this fleet; take whichever
# Blink-based browser is actually installed. All of these accept the same
# --headless=old --print-to-pdf flags and render identically.
BROWSER=""
for B in chromium-browser chromium google-chrome-stable google-chrome; do
    if command -v "$B" >/dev/null 2>&1; then BROWSER="$B"; break; fi
done
[ -n "$BROWSER" ] || { echo "no chromium/chrome found (tried chromium-browser, chromium, google-chrome-stable, google-chrome)"; exit 1; }
echo "  browser: $BROWSER"

mkdir -p "$WORK"
cp "$SRC" "$WORK/index.html"
for f in charts.js plotly.min.js; do
  [ -f "$(dirname "$SRC")/$f" ] && cp "$(dirname "$SRC")/$f" "$WORK/"
done
[ -d "$HERE/fonts" ] && mkdir -p "$WORK/fonts" && cp "$HERE/fonts"/Charis-*.ttf "$WORK/fonts/" 2>/dev/null || true

if [ ! -f "$WORK/fonts/Charis-Regular.ttf" ]; then
  echo "  Charis TTFs missing -- run $HERE/fetch-fonts.sh first"; exit 1
fi

echo "tighten figures"
python3 "$HERE/tighten.py" "$WORK/charts.js" "$WORK/charts-print.js"

echo "prepare html"
python3 "$HERE/mkpdf.py" "$WORK/index.html" "$WORK/build.html" "$HERE/print.css"
python3 "$HERE/inline_fonts.py" "$WORK/build.html" "$WORK/fonts"

echo "render"
rm -f "$OUT"
# Chrome treats every file:// URL as its OWN opaque origin, so the @font-face
# rules are same-origin-blocked and the faces silently do not load -- the PDF
# renders in a fallback and pdffonts shows Type 3 only. Chromium did not need
# this flag; Chrome does. The page-count and font checks below catch it either way.
timeout 600 "$BROWSER" --headless=old --disable-gpu --no-sandbox \
  --allow-file-access-from-files \
  --virtual-time-budget=60000 --no-pdf-header-footer \
  --print-to-pdf="$OUT" "file://$WORK/build.html" 2>&1 \
  | grep -iE "written|ERR_" || true

[ -f "$OUT" ] || { echo "  render produced nothing"; exit 1; }

PAGES=$(pdfinfo "$OUT" | awk '/^Pages/{print $2}')
if [ "${PAGES:-0}" -lt 5 ]; then
  echo "  FAILED: $PAGES page(s) -- $BROWSER could not read the work dir."
  pdftotext "$OUT" - 2>/dev/null | head -3 | sed 's/^/    /'
  exit 1
fi

echo
echo "  pages : $(pdfinfo "$OUT" | awk '/^Pages/{print $2}')"
echo "  size  : $(stat -c%s "$OUT") bytes"
echo "  serif : $(pdffonts "$OUT" | awk 'NR>2{print $1}' | sed 's/^[A-Z]*+//' \
            | grep -ci charis) Charis faces embedded (expect 4)"
echo "  out   : $OUT"
