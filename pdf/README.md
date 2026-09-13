# PDF build

Builds the submission PDF from the paper's HTML rendering.

```sh
./fetch-fonts.sh                              # once
./build.sh /path/to/site/index.html out.pdf
```

Produces a 74-page Letter PDF with the four Charis faces embedded. Requires
`chromium-browser`, `python3`, and poppler (`pdfinfo`, `pdffonts`).

## Why print HTML instead of typesetting it

There is no LaTeX or pandoc anywhere on the fleet. The paper's HTML rendering is
the same content as the markdown draft, so the PDF is produced by printing that
through headless Chromium. The HTML route also carries the nine figures and the
full table content, which the earlier LaTeX build did not.

## Files

| file | what it does |
|---|---|
| `build.sh` | driver: copies the source to a work dir, runs the two Python steps, renders, sanity-checks the result |
| `tighten.py` | rewrites `charts.js` → `charts-print.js` with figure boxes scaled ×0.62 |
| `mkpdf.py` | injects `print.css` and the `@font-face` block, rewrites Unicode exponents as `<sup>`, repoints the page at the tightened charts |
| `print.css` | the `@media print` rules: page size, margins, type, break behaviour |
| `fetch-fonts.sh` | downloads the four Charis TTFs from SIL upstream into `fonts/` |

`fonts/` is not committed: the release zip is ~10 MB, and the OFL is better
satisfied by pointing at upstream than by vendoring the binaries.

## Four traps, each of which fails silently

Every one of these produces a plausible-looking PDF rather than an error, so
they are worth knowing before changing anything here.

**1. The system Charter cannot be used.** The Bitstream Charter that ships with
X11 is Type 1 (`.pfb`), which Skia declines to load. It falls back to Liberation
Serif — a Times clone, thinner than what we were trying to get away from —
without reporting anything. The full TTFs are loaded from `fonts/` instead.

**2. Google's webfont CSS splits faces into `unicode-range` subsets.** Every
Greek letter and maths sign then drops through to the next font in the stack.
This is why `fetch-fonts.sh` pulls whole TTFs rather than linking the Google CSS.
Source Serif 4 was evaluated and rejected for a related reason: Google's delivery
would not yield a bold instance, so bold text reverted to Liberation Serif.

**3. The document writes exponents as literal Unicode superscripts.** 588 of
them — `1.4×10` followed by U+207B U+2076. Charis has no U+207x block, so each
one fell back; that alone put a foreign face on 69 of 83 pages. `mkpdf.py`
rewrites them as `<sup>` markup, which Charis sets from its own digits and which
is better typography than the Unicode glyphs in any case.

**4. Snap confinement.** Chromium cannot read outside `$HOME` (a `/tmp` source
gives `ERR_FILE_NOT_FOUND`) and cannot read dot-directories inside it (a work dir
like `~/.build` gives `ERR_ACCESS_DENIED`). Both surface as a one-page "file
couldn't be accessed" PDF. The work directory is therefore under `$HOME` and not
hidden; `build.sh` fails loudly if the render comes back under five pages.

## Typeface

Charis (SIL v7), Matthew Carter's Charter redrawn: sturdy stems, low stroke
contrast, drawn to hold colour in output that thins type. It replaces TeX Gyre
Pagella, a Palatino clone whose high contrast and light stems read as washed out
at text sizes, and Newsreader, the web page's face, which is thinner still.

## Page count

83 → 74, measured at each step rather than estimated:

| change | pages |
|---|---|
| figure boxes ×0.62 (3,580 px → 2,218 px) | 83 → 80 |
| leading 1.52→1.44, margins 0.85→0.78 in, tables 8.7→8.4 pt | → 75 |
| long note boxes allowed to break across pages | → 74 |

The figures were only three pages of it. The larger cost was orphaning: a
half-page `.note` carrying `break-inside:avoid` that will not fit pushes
wholesale to the next page and leaves the rest of the current one blank — one
page held a two-inch figure and 73% white. Mean trailing whitespace is now 12.4%,
down from 18%, with only three pages over 30% blank (title, last, and one other).

Body type is deliberately left at 10.4 pt. Shrinking it would recover more pages
and undo the reason for re-rendering.
