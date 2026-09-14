#!/usr/bin/env python3
"""Inline the Charis faces as data: URIs, and bind them OUTSIDE @media print.

Two separate Chrome behaviours had to be worked around here, and only the
second one is obvious in hindsight.

1. Chrome renders every file:// document in its own opaque origin, so an
   @font-face pointing at fonts/Charis-Regular.ttf is same-origin-blocked and
   silently does not load. --allow-file-access-from-files does NOT lift it
   (tested: byte-identical output). A data: URI is not a subresource fetch, so
   there is nothing left to block.

2. A font-family declared ONLY inside @media print is not applied by Chrome's
   headless --print-to-pdf. Bisected against the real stylesheet: four faces,
   font-display, font-synthesis, font-optical-sizing, text-rendering and
   print-color-adjust all embed correctly; moving the same body{} rule inside
   @media print is the single change that breaks it. Chromium honoured it,
   which is why the build worked until chromium was removed from the box, and
   why print.css binds the family there. The binding is therefore repeated
   unconditionally here, after the page's own stylesheet so it wins the
   cascade. print.css keeps its @media print rules for sizing and colour.

Without this the PDF renders in Georgia, every serif glyph is emitted as an
outlined Type 3 font, and the result is the washed-out rendering the Charis
work existed to fix.

    inline_fonts.py <build.html> <fonts-dir>
"""
import base64, io, os, re, sys

BIND = (
    '<style id="CHARIS-BIND">'
    'body,p,li,td,th,dt,dd,figcaption,blockquote,h1,h2,h3,h4,h5,h6,'
    'caption,em,strong,i,b,sup,sub'
    '{font-family:"Charis","Charis SIL",Georgia,serif !important}'
    'code,pre,kbd,samp,tt{font-family:ui-monospace,Menlo,Consolas,monospace !important}'
    '</style>'
)


def main(html, fontdir):
    s = io.open(html, encoding="utf-8").read()
    n, missing = 0, []

    def repl(m):
        nonlocal n
        name = m.group(1)
        p = os.path.join(fontdir, name)
        if not os.path.exists(p):
            missing.append(name)
            return m.group(0)
        n += 1
        return 'url("data:font/ttf;base64,%s")' % base64.b64encode(
            open(p, "rb").read()).decode("ascii")

    out = re.sub(r'url\("fonts/([A-Za-z0-9._-]+\.ttf)"\)', repl, s)
    if missing:
        print("  MISSING: %s" % ", ".join(sorted(set(missing))))
        return 1
    if not n:
        print("  no fonts/*.ttf url() references found -- nothing inlined")
        return 1

    if "CHARIS-BIND" not in out:
        i = out.rfind("</head>")
        if i < 0:
            print("  no </head> -- cannot bind the family")
            return 1
        out = out[:i] + BIND + out[i:]

    io.open(html, "w", encoding="utf-8").write(out)
    print("    inlined %d font faces as data: URIs  (%.1f MB -> %.1f MB)"
          % (n, len(s) / 1e6, len(out) / 1e6))
    print("    bound the family outside @media print (Chrome ignores it inside)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
