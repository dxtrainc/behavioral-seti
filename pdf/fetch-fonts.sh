#!/bin/bash
# Fetch the four Charis TTFs the PDF build needs.
#
# They are not vendored into the repo: the release zip is ~10 MB and the licence
# (OFL 1.1) is better satisfied by pointing at upstream than by copying the
# binaries in.  The four faces land in fonts/, which build.sh copies into its
# work directory so Chromium can load them by @font-face.
#
# Why these and not the system font: the Bitstream Charter that ships with X11
# is Type 1 (.pfb), which Skia cannot load -- it falls back to Liberation Serif
# without reporting anything.  And Google's webfont CSS splits a face into
# unicode-range subsets, so Greek and maths characters fall through the stack.
# The full TTFs avoid both failure modes.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="$HERE/fonts"
mkdir -p "$DEST"

if [ -f "$DEST/Charis-Regular.ttf" ]; then
  echo "  already present in $DEST"; exit 0
fi

URL=$(curl -fsSL --max-time 30 \
  "https://api.github.com/repos/silnrsi/font-charis/releases/latest" \
  | python3 -c "import json,sys;print([a['browser_download_url'] for a in json.load(sys.stdin)['assets'] if a['name'].endswith('.zip')][0])")
echo "  $URL"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
curl -fsSL --max-time 300 "$URL" -o "$TMP/charis.zip"
unzip -q -o "$TMP/charis.zip" -d "$TMP"
for f in Regular Bold Italic BoldItalic; do
  src=$(find "$TMP" -name "Charis-$f.ttf" | head -1)
  [ -n "$src" ] || { echo "  missing Charis-$f.ttf in the release"; exit 1; }
  cp "$src" "$DEST/"
done
find "$TMP" -iname "OFL*" -o -iname "*LICENSE*" | head -1 | xargs -r -I{} cp {} "$DEST/OFL.txt"
ls -la "$DEST"
