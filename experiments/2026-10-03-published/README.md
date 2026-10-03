# LOSS: no explicit long Loeschian arithmetic progression is publicly available anywhere as of 2026-10-03 03:50 CDT; the IBM September 2026 solution is still an unpublished commented-out placeholder in the blog post, and no (a, d) pair for any Loeschian AP of 20+ terms exists in any source checked.

Re-check of the 2026-10-03 earlier negative. Budget ~40 min, no computation needed (nothing to verify).

## Decisive finding (READ, not inferred)

Ponder This solutions are published **inline in the blog post** under a `## Solution` heading,
not at `/haifa/ponderthis/solutions/<Month>2026.html` (that path is 404 for every recent month,
including July 2026 whose solution *is* live). The September 2026 blog post's raw HTML contains
the solution section still wrapped in an HTML comment, i.e. the editorial template, unfilled:

```
<!--## Solution

:::details

::summary[Click here to view the solution]

Insert the solution text here.

:::
--->

## Solvers
- **\*Paul Lupascu** (1/9/2026 9\:45 AM IDT)
...
```

Source (raw HTML, `curl -sL https://research.ibm.com/blog/ponder-this-september-2026`).
For contrast, the August 2026 post at https://research.ibm.com/blog/ponder-this-august-2026
has a live, uncommented `## Solution` block ("The numerical solutions are: * n = 8: 62...",
and contains the 23-digit integer 27636190239652591799943). So August is published and
September is not; the September solution is pending, not hidden.

## What the September page DOES publish (verbatim)

The `**` holders, in reverse order of appearance with their original `**`-worthy record:

1. Jackson La Vallee (57)
2. Daniel Bitin (55)
3. Henk Waßmann (49)
4. Bertram Felgenhauer (45)
5. Paul Lupascu (44)

Final standings by length:

- n = 57: Jackson La Vallee
- n = 55: Daniel Bitin, Hubert Puszklewicz, Henk Waßmann, Achyuth Jayadevan, Emin-Ali CALYAKA
- n = 51: Daniel Chong Jyh Tar, Alper Halbutogullari
- n = 50: Bertram Felgenhauer
- n = 48: Peyo, Paul Lupascu, Jean-François Hermant
- n = 47: Lawrence Hon, George Jiri Spitalsky, Prashant Wankhede, Tamir Ganor & Shouky Dan
- n = 46: King Pig, Daniel Cahill, Nadir S.
- n = 45: Kang Jin Cho, Lazar Ilic, Vladimir Volevich, Stéphane Higueret
- n = 44: Aayaam Panigrahi

No starting value `a` and no step size `d` is given for any of these, anywhere on the page.
The only multi-digit integers in the September page's HTML are site-wide chrome numbers
(130431, 146887, 393803, 451973, 5067096, 57495365, 59855420, 6485594, 713848, 787157217,
8180510, 950799339, 9588225) — the identical set appears in the August post, so they are
template/analytics values, not puzzle data. INFERRED from the fact that the two pages share
the exact same number set.

Note the `**` ladder (44 -> 45 -> 49 -> 55 -> 57) is new information relative to the earlier
check in that it shows the records fell in that order; it still carries no (a, d).

## Every URL checked and what it returned

IBM:

| URL | Status / result |
|---|---|
| https://research.ibm.com/blog/ponder-this-september-2026 | 200, 60565 bytes. Challenge text + solver list + `**` ladder. Solution section present but **HTML-commented placeholder** ("Insert the solution text here."). No (a, d). |
| https://research.ibm.com/blog/ponder-this-august-2026 | 200, 116226 bytes. Solution section **live** (uncommented). Confirms the publish mechanism and that Sept is simply not out yet. |
| https://research.ibm.com/blog/ponder-this-september-2026-solution | 404 (no such URL pattern) |
| https://research.ibm.com/haifa/ponderthis/solutions/September2026.html | 404 |
| https://research.ibm.com/haifa/ponderthis/solutions/August2026.html | 404 |
| https://research.ibm.com/haifa/ponderthis/solutions/July2026.html | 404 (yet July's solution IS live in the blog — so this whole directory is dead, not evidence of anything) |
| https://research.ibm.com/haifa/ponderthis/challenges/September2026.html | 404 |
| https://research.ibm.com/haifa/ponderthis/index.shtml | 200. Lists September 2026 as current challenge, solution pending; August 2026 "The Wheel of Buttons"; July 2026 and earlier with solutions. |
| https://research.ibm.com/labs/israel/ponder-this | 200. September 2026 listed as current challenge. No September solution link. |

Wayback Machine:

| URL | Result |
|---|---|
| http://archive.org/wayback/available?url=research.ibm.com/haifa/ponderthis/solutions/September2026.html | `{"url": "...", "archived_snapshots": {}}` — never captured, consistent with the page never having existed. |

GitHub (via `gh search repos ponderthis` and the contents API):

| Target | Result |
|---|---|
| https://github.com/marty777/ponderthis (repo pushed 2026-09-08) | `2026/` contains only `01, 02, 04, 05, 06, 07, 08, README.md`. **No `09` directory.** Nothing for September 2026. |
| https://github.com/AlexFleischerParis/ponderthis | Flat file list; newest files are `challenge2025*`. Nothing from 2026 at all. |
| https://github.com/Lazar-Ilic/Lazar | Root contains only `Notes/` and `Resume/`. No Ponder This September 2026 content surfaced. |
| `gh search repos ponderthis --limit 30` | 30 repos; the only 2026-dated one besides marty777 is `CoinGH/IBM-PonderThis-June-2026` (June, irrelevant). No September 2026 repo. |
| GitHub code search API for "Loeschian" | 401 Requires authentication — NOT checked. Gap. |
| https://gist.github.com/danielchong/080207c606de85dd03298434c0b8fe03 (Daniel Chong's solved-challenges gist, surfaced by search) | Listed by search as a solved-challenge index only; no progression data surfaced. |

OEIS (JSON API):

| Query | Result |
|---|---|
| `https://oeis.org/search?q=Loeschian+arithmetic+progression` | Only A301430 (Landau–Ramanujan analog constant, Waldschmidt 2018). **No Loeschian-AP sequence exists.** |
| `https://oeis.org/search?q=A003136+arithmetic+progression` | Only A301430. |
| `https://oeis.org/search?q=keyword:new+Loeschian` | `null` (no results). Nothing submitted since 2026-10-03. |
| A003136 itself | Last modified 2026-08-22, revision 344. No Ponder This link, no AP comment. |

Hugo Pfoertner's sums-of-two-squares analogues (the one hypothesis that could have moved):

| Sequence | State as of this check |
|---|---|
| A398603 "least number which is the end of an AP of n numbers that are sums of two squares (A001481)" | 33 terms, last = 33344305. Keyword `more`. Author Hugo Pfoertner, created 2026-08-04. |
| A398604 (corresponding steps) | 33 terms, last = 825132. Created 2026-08-04, revision 5, last touched 2026-08-04T15:09:14-04:00. |
| A093365 "least number which is the end of an AP of n numbers that are sums of two *nonzero* squares" | **37 terms**, data ends `..., 852688273, 1138267241, 3106000349, 3131438993`. Comments: "a(34) > 225000000" and "a(38) > 10^10. - Hugo Pfoertner, Aug 04 2026". |

**Not extended past n = 37.** The known record a = 2215647809, d = 25438644, n = 37 has last term
2215647809 + 36*25438644 = 3131438993, which is exactly A093365(37) — consistent, and a(38) is
still only bounded below (> 10^10), not found. So nothing new there either.

Web/forum searches run (all returned only the IBM challenge page itself plus generic
math references; no (a, d) anywhere):

- "Ponder This September 2026 solution Loeschian arithmetic progression"
- "github ponder this 2026 09 Loeschian progression 57 terms"
- "\"Loeschian\" arithmetic progression 57 terms starting value step size solution"
- "\"Jackson La Vallee\" Loeschian OR \"ponder this\" 2026" — surfaced only IBM solver lists; no personal page, repo, blog, or writeup for La Vallee exists in the index.
- "mersenneforum OR mersenne.ca Loeschian arithmetic progression ponder this September 2026" — only old (2023) Ponder This threads. `mersenneforum.org/search.php` returns 301 to a login/redirect, so the forum's internal search was NOT exercised directly. Gap.
- "OEIS \"Loeschian\" arithmetic progression sequence 2026 Pfoertner longest"
- "Loeschian numbers long arithmetic progression construction CRT \"x^2+xy+y^2\" record 2026 blog"

## Verification

Nothing to verify. `src/verify.py` and `src/crosscheck.py` were not run because no candidate
(a, d, n) was found. Zero external examples exist.

## When to re-check

August 2026's solution is live and September's is a filled-in-later placeholder in the same
blog post, so the September solution will appear by editing
**https://research.ibm.com/blog/ponder-this-september-2026** in place — the URL will not change
and no new page will appear. Re-check that single URL's raw HTML for the string
`Insert the solution text here` (present => still unpublished). The IBM archive page states
solutions are posted "after the beginning of the next month"; the challenge closed 2026-09-30
and this check is only ~3 days later. **Worth re-checking 2026-10-06 and again 2026-10-10.**
Historically the solution text also need not contain the record (a, d) at all — the August post's
solution discusses method and gives only the requested numerical answers — so even once published
it may only restate n = 57 without La Vallee's progression. Do not plan around it.

## Residual gaps (honest)

1. GitHub *code* search for "Loeschian" requires auth and was not run; only repo-name search and
   directory listings were.
2. mersenneforum.org internal search redirects (301) and was not driven; only web-index search.
3. Private/unindexed solver correspondence (the IBM mailing list) is unreachable by construction.
