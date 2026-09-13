# web

Standalone pages that accompany the paper.

| file | what it is |
|---|---|
| `summary.html` | a five-page summary of the paper, readable on its own |

## summary.html

Covers the premise, where the work sits against the 2018 NASA Technosignatures
Workshop taxonomy, the twelve-step gate chain, gate 0 and what recognition costs
in compute, the search stated as a cryptanalyst's procedure, the 29 searches and
their outcomes, future research directions, and what the hypothesis does and does
not do to the Fermi debate.

Self-contained apart from two things: webfonts from Google (Newsreader and Inter,
for screen) and the Charis TTFs (for print). Serve it from any static server:

```sh
cd web && python3 -m http.server 8111
```

### Fonts

The page reads on screen in Newsreader and **prints in Charis**, which needs the
TTFs present alongside it:

```sh
../pdf/fetch-fonts.sh ./fonts
```

Without them the print view falls back to Georgia — legible, but thinner on paper
than intended; the reasoning is the same as in `pdf/README.md`. `fonts/` is
gitignored for the same reason it is there: the release is ~10 MB and the OFL is
better served by pointing at upstream than by vendoring binaries.

### Print

The page carries an `@media print` block that takes it from ten pages to five
without cutting content — screen type is set at 18 px for a monitor, which prints
at roughly 13.5 pt. Two things in that block are load-bearing and were arrived at
by measuring rather than by taste:

- **Tables and asides must be allowed to break.** Holding a half-page block
  together orphans whatever sits above it. With `break-inside: avoid` on them the
  page content totalled 4.81 pages but spanned six; letting them flow, with
  `thead` repeating on each fragment, put pages 1–5 at 90–94% full.
- **The last 5% comes from page furniture, not type.** Margins, cell padding and
  paragraph spacing were tightened to land on five pages so that the body could
  stay at 9.4 pt.

To check a change has not pushed it over:

```sh
chromium-browser --headless=old --disable-gpu --no-sandbox \
  --virtual-time-budget=18000 --no-pdf-header-footer \
  --print-to-pdf=/tmp/summary.pdf "file://$PWD/summary.html"
pdfinfo /tmp/summary.pdf | grep ^Pages      # expect 5
```

Chromium cannot read the file from `/tmp` or from a dot-directory under `$HOME`;
both fail as a one-page "file couldn't be accessed" PDF rather than as an error.
