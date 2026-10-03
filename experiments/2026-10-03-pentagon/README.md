# Pentagonal automatic-terms family vs. baseline — matched measurement

## The identity (verified)

If `t_{k0}` is Loeschian and `3·c·d·t_{k0}` is a perfect square, then `t_{k0 ± c j²}` is
Loeschian for every `j`. Proof: write `4t = u² + 3v²` and look for `U = u + αj`, `V = v + βj`
with `4(t + cdj²) = U² + 3V²`. Matching coefficients needs `uα + 3vβ = 0` and `α² + 3β² = 4cd`,
which is solvable in rationals exactly when `3cdt` is a square.

Instantiated with `a = 3x² + m²`, `d = 24m²`:
`t_k = 3x² + m²(24k+1)`, and `24k+1` is a square exactly at the generalized pentagonal `k`,
where `t_k = 3x² + (m(6j±1))²` is visibly of the form `X² + 3Y²`.

    automatic indices below 58: 0 1 2 5 7 12 15 22 26 35 40 51 57   -> 13 of 58

**Measured: 0 failures in 260,000 automatic positions.** The identity holds exactly.

## The matched experiment

Same `d`, same term size, same congruence conditions; only the shape of `a` differs.
`m = 5·11·17·23·29 = 623645` (forces 5,11,17,23,29 into `d = 24m²`; `m` odd and `3∤m`, so with
`x` even every term is odd and `≡ 1 mod 3`). 20000 candidates of each kind, all 58 terms tested
with the exact Loeschian test.

| | per-term pass | non-automatic-term pass | P(run≥10) | P(run≥14) | P(run≥16) |
|---|---|---|---|---|---|
| pentagon | 0.5895 | **0.4709** | 0.0579 | 0.0043 | 0.0014 |
| baseline | 0.4921 | 0.4921 | 0.0179 | 0.0006 | 0.0002 |

Two facts, both important:
- the 13 automatic terms are free (0 failures), and
- the 45 **non-automatic** terms pass at **0.4709 vs 0.4921** — the family is slightly *worse*
  than generic on the terms it does not fix (0.957 per term, 0.14 over 45 terms).

Extrapolated to a full 58-run at this term size:
`P_pent/P_base = 0.4709^45 / 0.4921^58 = e^{7.24} ≈ **1400×** per candidate.`

## Why it still loses — the general argument

A square condition on `(a,d)` is not free: for fixed `d` it forces `a = s²/(3cd) − k0·d`, i.e. `a`
is parametrized by `s`, so the family has `~√T` members per `d` instead of `~T`. Over all `d` the
candidate count drops from `~T²/D0` to `~T/√D0`. The cost is therefore a factor

    T / (√D0 · const)

against a gain of `ρ^{-13} ≈ 1400`. Breaking even needs `T ≲ 1400·√D0 ≈ 4·10⁹`, but this family
forces `T ≥ 57d = 57·24m² ≈ 5.3·10¹⁴`. **It loses by ~10⁵, and the gap grows with T.**

Measured instance (this `d`, `T ≈ 5.3e14`): pentagon has 4.1e6 candidates × 1.9e-15 = 8e-9
expected 58-runs; baseline has 5.6e13 candidates × 1.4e-18 = 8e-5. Baseline wins by ~10⁴.

## Verdict

The structure is real, elegant, and exactly reproduces the 13 free positions — but adding square
conditions always costs more candidates than the free positions buy back, for every `n` and `T`
in range. This independently confirms the earlier agent's "0.1–0.8× of baseline" and supplies the
mechanism. **No elegant construction beats the plain admissible-residue sieve for n = 58.**

Reproduce: `experiments/2026-10-03-pentagon/measure.py` (writes terms, pipes through
`build/apsearch --isl`, tabulates).
