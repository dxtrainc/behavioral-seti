#!/usr/bin/env python3
"""Produce a print variant of charts.js with the figures tightened.

The nine Plotly figures are authored at screen heights (520, 520, 480, 520, 520,
420, 420, 340, 360 px = 4,100 px total, about 4.6 pages of pure figure).  Two
costs follow in print: the height itself, and orphaning -- each figure carries
break-inside:avoid, so a 520 px block that does not fit the remaining space
pushes to the next page and leaves up to a full page blank behind it.  Measured
trailing whitespace across the document averaged 18% of page height.

Scaling the plot area down cuts both.  Type inside the plots is scaled far less
than the box, so labels stay readable: at SCALE 0.62 a 520 px figure becomes
3.4 in tall, which is a normal single-column figure, while tick labels go from
13 px to 11 px and axis titles track them.
"""
import io, re, sys

SCALE = 0.62          # plot box
FONT = 11             # was 13
AXIS_TITLE = 11       # was 13
LEGEND = 10           # was 12
TICK = 10

def tighten(js):
    n = [0]

    def h(m):
        n[0] += 1
        return "height:%d" % max(210, round(int(m.group(1)) * SCALE))
    js = re.sub(r"height:(\d+)", h, js)

    # base font, axis titles, legend
    js = js.replace('font:{family:"Inter, system-ui, sans-serif", size:13, color:C.ink}',
                    'font:{family:"Inter, system-ui, sans-serif", size:%d, color:C.ink}' % FONT)
    js = js.replace("title:{font:{size:13,color:C.slate}}",
                    "title:{font:{size:%d,color:C.slate}}" % AXIS_TITLE)
    js = js.replace("legend:{orientation:\"h\", y:-0.22, x:0, font:{size:12}}",
                    "legend:{orientation:\"h\", y:-0.20, x:0, font:{size:%d}}" % LEGEND)
    js = js.replace("font:{size:12}", "font:{size:%d}" % LEGEND)

    # trim the generous plot margins; the long category labels on the ranking and
    # receiver charts need their left gutter, so those are reduced proportionally
    js = js.replace("margin:{l:70,r:30,t:30,b:60}", "margin:{l:58,r:20,t:16,b:46}")
    js = js.replace("margin:{l:290,r:30,t:20,b:70}", "margin:{l:232,r:20,t:12,b:52}")
    js = js.replace("margin:{l:320,r:70,t:20,b:60}", "margin:{l:256,r:48,t:12,b:46}")
    return js, n[0]


if __name__ == "__main__":
    src, dst = sys.argv[1], sys.argv[2]
    js = io.open(src, encoding="utf-8").read()
    out, count = tighten(js)
    io.open(dst, "w", encoding="utf-8").write(out)
    before = [int(x) for x in re.findall(r"height:(\d+)", js)]
    after = [int(x) for x in re.findall(r"height:(\d+)", out)]
    print("  %d figure heights scaled x%.2f" % (count, SCALE))
    print("    before: %s  = %d px" % (before, sum(before)))
    print("    after : %s  = %d px" % (after, sum(after)))
    print("    saved : %d px = %.1f in = %.1f pages of figure box"
          % (sum(before) - sum(after), (sum(before) - sum(after)) / 96.0,
             (sum(before) - sum(after)) / 96.0 / 9.28))
