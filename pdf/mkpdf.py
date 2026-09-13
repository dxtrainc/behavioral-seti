#!/usr/bin/env python3
"""Prepare the paper's HTML for print: fonts, exponents, tightened figures.

Run by build.sh; see the README in this directory for the whole pipeline.

Three things have to be done to the source HTML or the PDF comes out wrong in a
way nothing warns you about:

1. FONT.  The original PDF was set in TeX Gyre Pagella, a Palatino clone whose
   high stroke contrast and light stems read as washed out at text sizes; the
   web page uses Newsreader, which is thinner still.  Both are replaced by
   Charis (SIL v7) -- Charter redrawn, sturdy stems, low contrast, drawn to hold
   colour in output that thins type.  The Bitstream Charter that ships with X11
   cannot be used: it is Type 1 (.pfb), which Skia silently declines, falling
   back to Liberation Serif.  Google's webfont CSS is no good either -- it
   splits a face into unicode-range subsets, so every Greek letter and maths
   sign drops through to the next font in the stack.  The full TTFs are loaded
   from fonts/ instead, which avoids both.

2. EXPONENTS.  The document writes them as literal Unicode superscripts (588 of
   them: "1.4x10" then U+207B U+2076).  Charis has no U+207x block, so every one
   fell back -- that alone put a foreign face on 69 of 83 pages.  They are
   rewritten as <sup> markup, which Charis sets from its own digits, and which
   is better typography than the Unicode glyphs regardless.

3. FIGURES.  The script points the page at charts-print.js, the tightened chart
   definitions produced by tighten.py.
"""
import io, os, re, sys

SUP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
       "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9",
       "⁻": "−", "⁺": "+", "ⁿ": "n", "ⁱ": "i"}
SUB = {"₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4",
       "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9",
       "ₖ": "k", "ₐ": "a", "ₙ": "n"}
EMOJI = {"❌": "✗", "✅": "✓"}          # colour emoji print badly

FACES = """<style id="CHARIS-FACES">
@font-face{font-family:"Charis";src:url("fonts/Charis-Regular.ttf") format("truetype");font-weight:400;font-style:normal;font-display:block}
@font-face{font-family:"Charis";src:url("fonts/Charis-Bold.ttf") format("truetype");font-weight:700;font-style:normal;font-display:block}
@font-face{font-family:"Charis";src:url("fonts/Charis-Italic.ttf") format("truetype");font-weight:400;font-style:italic;font-display:block}
@font-face{font-family:"Charis";src:url("fonts/Charis-BoldItalic.ttf") format("truetype");font-weight:700;font-style:italic;font-display:block}
@media print{
  sup,sub{font-size:.68em;line-height:0;position:relative;vertical-align:baseline}
  sup{top:-.46em}
  sub{bottom:-.22em}
}
</style>
"""

STACK = '"Charis",Georgia,serif'
# the placeholder stacks print.css ships with, rewritten to the real one
PLACEHOLDERS = [
    '"Charter","Bitstream Charter","Source Serif 4","Charis SIL",Georgia,serif',
    '"Charter","Bitstream Charter",Georgia,serif',
]


def map_text_nodes(html, fn):
    """Apply fn to text outside tags, <script> and <style>."""
    skip = [(m.start(), m.end()) for m in
            re.finditer(r"<(script|style)[^>]*>[\s\S]*?</\1>", html, re.I)]

    def guarded(i):
        return any(a <= i < b for a, b in skip)

    out, pos = [], 0
    for m in re.finditer(r"<[^>]+>", html):
        if m.start() > pos:
            seg = html[pos:m.start()]
            out.append(seg if guarded(pos) else fn(seg))
        out.append(m.group(0))
        pos = m.end()
    if pos < len(html):
        seg = html[pos:]
        out.append(seg if guarded(pos) else fn(seg))
    return "".join(out)


def rewrite(seg):
    for bad, good in EMOJI.items():
        seg = seg.replace(bad, good)
    for tbl, tag in ((SUP, "sup"), (SUB, "sub")):
        rx = re.compile("[" + "".join(re.escape(c) for c in tbl) + "]+")
        seg = rx.sub(lambda m: "<%s>%s</%s>"
                     % (tag, "".join(tbl[c] for c in m.group(0)), tag), seg)
    return seg


def build(src, dst, css_path):
    s = io.open(src, encoding="utf-8").read()

    before = sum(s.count(c) for c in list(SUP) + list(SUB))
    s = map_text_nodes(s, rewrite)
    after = sum(s.count(c) for c in list(SUP) + list(SUB))

    css = io.open(css_path, encoding="utf-8").read()
    for ph in PLACEHOLDERS:
        css = css.replace(ph, STACK)
    if 'id="PRINT-STYLESHEET"' not in s:
        # after the last author style block, so it wins at equal specificity
        i = s.rfind("</style>") + len("</style>")
        s = s[:i] + '\n<style id="PRINT-STYLESHEET">\n' + css + '\n</style>\n' + s[i:]
    if "CHARIS-FACES" not in s:
        s = s.replace("</head>", FACES + "</head>", 1)

    if os.path.exists(os.path.join(os.path.dirname(dst) or ".", "charts-print.js")):
        s = s.replace('src="charts.js"', 'src="charts-print.js"')

    io.open(dst, "w", encoding="utf-8").write(s)
    print("  superscript/subscript codepoints: %d -> %d" % (before, after))
    print("  wrote %s" % dst)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("usage: mkpdf.py <src.html> <dst.html> [print.css]")
    css = sys.argv[3] if len(sys.argv) > 3 else \
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "print.css")
    build(sys.argv[1], sys.argv[2], css)
