# Caputo Fractional Differential Equation with Bounded Delay — Monotone Iterative Method

This repository contains the Python/SymPy code accompanying the paper:

> J. D. Ramírez, *"Existence of Extremal Solutions for Caputo Fractional Differential
> Equations with Bounded Delay,"* Progress in Fractional Differentiation and
> Applications, **9**(3), 487–498 (2023).
> http://dx.doi.org/10.18576/pfda/090311

The notebook `Caputo_FDE_with_Delay.ipynb` reproduces **Example 1 / Section 4** of the
paper: it builds the coupled monotone (lower/upper) sequences from **Theorem 3** and
shows numerically that they squeeze down to the unique solution of a Caputo fractional
delay IVP.

## What the code does

The paper studies the delay IVP

```
 cD^q u(t) = f(t, u(t), u_t(1)) + g(t, u(t), u_t(1)),   t ∈ J = [0, 1]
 u(s) = φ(s),                                          s ∈ [-1, 0]
```

with Caputo derivative order `q = 1/2`, and in the worked example:

```
f(t, u(t), u_t(1)) =  (1/10) u(t) + (1/6) u^2(t-1)        (increasing in u, u_t)
g(t, u(t), u_t(1)) = -(1/11) u^2(t) - (1/4) u(t-1)        (decreasing in u, u_t)
φ(s) = 1,  s ∈ [-1, 0]
```

Because `f` is increasing and `g` is decreasing in both arguments, `v⁰ ≡ 0` and
`w⁰ ≡ 2` form a pair of **Type II coupled lower/upper solutions** (Definition 5,
eq. (12)). Theorem 3 then defines two intertwined sequences via the equivalent
Volterra fractional integral equation (eq. (6)):

```
 v^{n+1}(t) = φ(0) + (1/Γ(q)) ∫₀ᵗ (t-x)^{q-1} [ f(x, w^n(x), w^n_x(1)) + g(x, v^n(x), v^n_x(1)) ] dx
 w^{n+1}(t) = φ(0) + (1/Γ(q)) ∫₀ᵗ (t-x)^{q-1} [ f(x, v^n(x), v^n_x(1)) + g(x, w^n(x), w^n_x(1)) ] dx
```

which converge monotonically and uniformly to the coupled minimal solution `ρ`
and maximal solution `r`. Under the one-sided Lipschitz condition (Remark after
Theorem 3, eq. (18)) `ρ = u = r`, i.e. the sequences squeeze onto the unique
solution `u`.

Since `q = 1/2`, `Γ(1/2) = √π`, so every kernel `(t-x)^{q-1}` becomes
`(t-x)^{-1/2}` and the prefactor `1/Γ(q)` becomes `1/√π` — this is exactly what
appears in the notebook's `integrate(...)` calls.

## Notebook structure

| Cell(s) | Purpose |
|---|---|
| 1 | Imports: `sympy` (symbolic integration/`Piecewise`/`plot`), `numpy`, `matplotlib`, `torch` (unused device stub), `fractions.Fraction` |
| 3 | Symbols `t, x, s`; an earlier scratch definition of `v0`, `w0` (superseded later) |
| 4–8 | Scratch/earlier-attempt versions of `v1`–`w6` for a *different* toy equation — kept for reference but **not** the ones used in the final Example 1 computation/plot |
| 9–10 | Markdown restating Example 1 and why `v⁰=0, w⁰=2` are Type II lower/upper solutions |
| 11 | Sets the actual `v0 = 0`, `w0 = 2` used in the example |
| 12–18 | Recursively build `v1..v7` and `w1..w7` as `sympy.Piecewise` functions: `1` for `t ≤ 0` (the history `φ(s)=1`) and the fractional Volterra integral for `t > 0`, using the previous iterate's opposite-partner in the `f`/`g` split (Type II coupling) |
| 19 | Plots all 16 functions (`v0..v7`, `w0..w7`) on `t ∈ [-1, 1]` — reproduces Fig. 1 of the paper ("Monotone method showing fourteen iterates") |

Each iterate is computed **symbolically** with SymPy's `integrate`, so the code
is exact/closed-form (not a numerical quadrature scheme) — this mirrors the
paper's remark that the method avoids needing an explicit solution formula and
instead works directly from the Volterra integral equation (6).

## Requirements

```
python >= 3.9
sympy
numpy
matplotlib
torch          # imported but not required for the computation (device line is commented out)
```

Install with:

```bash
pip install sympy numpy matplotlib torch
```

(`torch` can be omitted if you delete the `import torch` line — it is not used
in any calculation.)

## Running

Open and run all cells in order:

```bash
jupyter notebook Caputo_FDE_with_Delay.ipynb
```

Cells 12–18 must run sequentially since each iterate depends on the previous
`v_n`, `w_n` pair. **Note:** symbolic integration of `(t-x)^{-1/2}` times
increasingly complicated polynomial/`Piecewise` expressions gets slower with
each iterate — computing `v7`/`w7` (cell 18) can take a noticeably longer time
than the earlier cells.

## Output

- A plot (`p`, cell 19) equivalent to **Figure 1** in the paper: `v⁰ ≤ v¹ ≤ … ≤ v⁷
  ≤ u ≤ w⁷ ≤ … ≤ w¹ ≤ w⁰` visibly closing in on the true solution `u` on
  `t ∈ [-1, 1]`.
- Evaluating `v7`, `w7` at `t = 0, 0.1, …, 1.0` reproduces **Table 1** of the
  paper (values ≈ 0.9235–1.0 across the two sequences, agreeing to 3–4 decimal
  places by `t = 1`).

## Relating code objects to the paper's notation

| Code | Paper |
|---|---|
| `v0, w0` | `v⁰(t), w⁰(t)` — Type II lower/upper solutions, eq. (12) |
| `vN, wN` (`N=1..7`) | `v^N(t), w^N(t)` — intertwined sequences from Theorem 3, eqs. (16)-(17) |
| `q = 1/2` (hardcoded as `-1/2` exponent, `1/√π` prefactor) | fractional order `q` |
| `Piecewise((1, t<=0), (..., t>0))` | the history segment `φ(t-t0)` for `t0-τ ≤ t < t0` vs. the Volterra integral for `t0 ≤ t ≤ t0+a` in eq. (6) |
| `1/10`, `1/6`, `1/11`, `1/4` coefficients | the coefficients in `f` and `g` of eq. (19) |

## Caveats / things to know before extending this code

1. Cells 4–8 are leftover scratch work for a different (non-delay-specific,
   earlier draft) equation and are **not** wired into the Example 1 pipeline —
   don't mix them into the final computation.
2. The recursion is unrolled by hand (`v1`…`v7` are separate named cells rather
   than a loop). To go beyond 7 iterates or change the equation/coefficients,
   generalize this into a `for n in range(N)` loop that stores each
   `(v_n, w_n)` pair, substituting the new `f`, `g`, `q`, and history function
   `φ`.
3. Symbolic (`sympy`) evaluation is exact but does not scale well; for higher
   iterate counts or more complex `f`, `g`, consider numerically discretizing
   the Volterra integral (e.g., a fractional Adams/product-trapezoidal scheme)
   instead of symbolic integration.
