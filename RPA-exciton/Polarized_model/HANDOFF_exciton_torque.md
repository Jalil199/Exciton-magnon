# Handoff — Exciton-induced torque & damping on classical spins (CrSBr)

**Date:** 2026-07-09 · **Dir:** `RPA-exciton/Polarized_model/`

## Goal
How a photoexcited **exciton population** exerts **torque and non-Gilbert damping** on the classical
layer magnetizations in CrSBr, reproducing the user's own paper *"Landau-Lifshitz equation with
nonlocal damping due to exciton bath"* (Garcia-Gaitan, Varela-Manjarres, Nikolić —
`LLG_for_exciton_bath.pdf`). Target = generalized LLG:
`∂_t M_L = M_L × [B_L^cl + j_θ^L M_L̄ + Σ η_θ^{LL'} ∂_t M_L']`, with torque `j_θ ∝ N` (exciton number)
and nonlocal, non-Gilbert damping `η^{LL'}` (N-independent).

⚠️ **`LLG_for_exciton_bath.pdf` is only a DRAFT and may still contain incorrect results/formulas**
(user, 2026-07-24). The point of this whole notebook line is to independently re-derive and
numerically check each piece from the microscopic BSE/Green's-function machinery — not to assume
the draft is already right. When a notebook result disagrees with the draft, flag it; don't force
a reconciliation.

**Experimental cross-check for the damping (new, 2026-07-24):** Brennan, Tang, Varela-Manjarres,
Chang, Chica, Zhu, Roy, Nikolić, Ren, Bae, *"Excitonic spin torque in a magnetic semiconductor,"*
Nature Materials, 2026-06-15 (same collaboration; the user is a co-author). CrSBr, canted AFM,
2.4 eV optical pump. Reports a **nonlinear, angle-dependent damping**: "the same exciton density
exerts a stronger spin torque when spins are more canted," from misalignment of the
non-equilibrium spin polarization **s**(t) with **M**(t) (**s**×**M** torque). Phenomenological
form `X'' + ω₀²X + γX' + (μ²/9)X³ − μXX' = 0`; LLG terms
`τ_DL = A_d[M_i·M_j/(M_i·M_i)](M_i·dM_i/dt)` (damping-like, μ≈0.56 roughly material-intrinsic),
`τ_FL = g_c μ₀ e^(−t/τ_X) â·M_i` (field-like). Energy balance `dE/dt = -(γ-μX)X'²`: for `X>γ/μ`
the effective damping goes **negative** (exciton bath pumps energy into the spins). The `−μXX'`
term is the experimental target the Step 10 microscopic damping derivation (below) should be
checked against — better than the draft's own DMRG numbers, since this is a published result.

## Model (4-band minimal, θ = angle between layer moments; FM θ=0, canted π/2, AFM π)
- Two layers × (valence + conduction). Interlayer tunneling ∝ `cos(θ/2)=√((1+M1·M2)/2)`.
- `H_valence(k,θ)`, `H_conduction(k,θ)` (see any notebook cell). Params: `L=200`, `γ_c=γ_v=0.40`,
  `Δ_0=2.0`, `Jpd_v=0.069`, `Jpd_c=0.031`, `V_intra=1.708`, `V_inter=0.691`. Bright exciton
  calibrated to CrSBr A exciton: `E_X=1.3521 eV`, binding `E_bind≈0.68 eV` (QP gap 2.0 vs optical 1.35).
- Excited state: FD occupations with separate quasi-Fermi `μ_c=Δ_0/2+s`, `μ_v=Δ_0/2−s`, hot `Tx=0.10`;
  working point `s=0.7` → carrier density `n_c≈8.6e-3`/cell.
- Same model as the LLG paper (params differ slightly; we keep our calibration).

## Notebooks
- **`ExcitonBathKernel.ipynb`** (MAIN, 61 cells) — step-by-step reproduction of Perfetto-Stefanucci
  PRB 94 245303 (2016) incoherent-exciton quantities, toward the torque. See status below.
- `ExcitonTorque.ipynb` — first-variation single-D torque; per-exciton field `∂E_X/∂θ ∝ sinθ`
  (**7.5 meV/rad at θ=π/2**, HF-validated). **This is the correct torque** (see finding below).
- `ExcitonSpinSusceptibility.ipynb`, `DrivenExcitonTorque.ipynb`, `BrightExcitonMode.ipynb` — earlier.
- Reference PDFs: `Perfetto et al. - 2016 ...pdf` (the 2016 paper), `LLG_for_exciton_bath.pdf` (user paper).

## ExcitonBathKernel status
Steps 1–6 (**verified**, standalone Julia before appending; PyPlot cells run in user's kernel):
1. Π^R bubble (Eq 15). 2. χ_Φ=Π+ΠDΠ, pole `E_X=1.3521`. 3. Excited FD occupations, `E_X` blueshifts
   with n_c. 4. Lesser χ^< and exciton weight `Z_X=F^bright` (Eq 21/37); **mass action `Z_X∝n_c²`**
   (neutral photoexcitation n_c=n_h). Coherent contraction `coh=|Σ_p Y_p|²` (Eq 65) — bright dominates.
5. Finite-q BSE (`transitions_q`, `bse_q_vertex` → `Ω,Fλ,Fbλ,coh,vtx`), `Σ^R_cc, Σ^<_cc`, sideband at
   `ε_v(P-q)+Ω_λ(q)`. 6. `G^R_cc, A_cc, G^<_cc`. **Vertex EXACT (fixed 2026-07-24):** weight `vtx/L²`
   with `vtx=|Σ_I K_I|²`, `K_I=(ΔE_I−Ω_λ)Y_I` (2026 Eq.6) — no free factor, replaces the old
   `Ccal=V̄²` calibration. At q=0 bright mode, exact `vtx` is **~2.0×** the old `Ccal·coh` (measured).
   Side effect: spectral sum rule moves `∫A/2π: 0.999→0.87` — the old value was only ≈1 because
   `Ccal` had been reverse-tuned to hit it; the ~13% deficit with the untuned exact vertex reflects
   other truncations already present (finite `Nq=48`, 2 bound modes/q, finite `η`, single sector).

Added 2026-07-09:
- θ-swept spectra `A_{P=0}(ω,θ)`, `-iG^<_{P=0}(ω,θ)` (FM→AFM); gap edge opens ~95 meV, sideband ~12 meV.
- **N_bright figures**: `N_bright=Σ_q F^bright_q` vs n_c (slope-2 mass action) and vs θ. At s=0.7:
  `N_bright≈9e-5`/cell (per-CELL density; L=200 is the k-grid NOT sample size; macroscopic sample has
  10²–10⁵ excitons; bound fraction ~1%). n² is because photon makes e+h (n_c=n_h); F^λ∝f_c(1-f_v) bilinear.
- **Step 7 valence `Σ_vv` (cells 50-51) — ⚠️ WRONG, DO NOT TRUST (see finding).**
- **Step 8 conduction-only torque (cells 53-55, added by user)** — uses `∫(-iG^<_cc)dω/2π = n_c(P)`
  weighted by `X_c`. The **occupation route is correct**, but it is conduction-only ⇒ **overestimates
  the torque ~3.6×** (missing valence, see finding). Superseded by Step 9.
- **Step 9 FULL cc+vv torque (cells 56-57)** — per-exciton torque `∂E_X/∂θ` (benchmark-verified),
  decomposed `cc(+) − vv(+)`, ×`N_bright(θ)`, LLG `j_θ=B_e/(4cos(θ/2))`, `τ_1=j_θ sinθ`. θ-sweep shows
  cc/vv cancellation 60–93% (stronger toward AFM); net `∂E_X/∂θ ∝ sinθ` peaks 7.50 meV/rad at π/2.
  ⚠️ Uses the `_top valence branch is θ-independent` lesson: the hole lives in BOTH branches, weighted
  by `|Y_I|²` in the transition basis (do NOT use `εv_top` for ∂ε_v/∂θ — it gives 0). Right panel plots
  `τ_1=N∂E_X/∂θ` (physical) AND raw `B_e=N∂E_X/∂c`. **Units: meV** (per-exciton coeff meV/rad; total
  torque tiny ~1e-3..1e-4 meV = per-cell, dilute).
- **Step 9b density family (cells 59-60)** — same torque panel for several `n_c` (s=0.4→0.9). `∂E_X/∂θ`
  is density-independent, so `B_e,τ_1(θ;n_c) = N_bright(θ;n_c)·∂E_X/∂θ` form a **∝n_c² family with an
  ~universal θ-shape** (log scale, curves span ~4 decades). Reuses `bright_split_t9` (run Step 9 first).

## ★ KEY FINDING TODAY — single-exciton benchmark (settles the valence factors)
Exciton state `|X⟩=Σ_k Y_k c†_k v_k|GS⟩` gives, **exactly**:
- `n_c(k) = +|Y_k|²`  (conduction occupation added by the exciton)
- `n_v(k) = 1 − |Y_k|²`  (valence depleted by `|Y_k|²`)

**The torque is an EQUAL-TIME quantity:** `B_e = i Tr[G^< X̂] = Σ_k [n_c(k) ∂ε_c/∂θ + n_v(k) ∂ε_v/∂θ]`
— it needs the **occupations**, NOT the spectral satellite. Therefore:
```
B_e^exc / N = Σ_I |Y_I|² (∂ε_c − ∂ε_v)/∂θ = ∂E_X/∂θ   [Hellmann-Feynman, V θ-independent]
            = 7.4952 meV/rad at θ=π/2   (ratio to finite-diff ∂E_X/∂θ = 1.0000, EXACT)
```
So the physical torque `τ_1 = N ∂E_X/∂θ ∝ N sinθ` — matches `ExcitonTorque.ipynb` (7.5 meV) and LLG `j_θ∝N`.

**⚠️ c vs θ chain rule (fixed 2026-07-09, Step 9/9b):** H depends on θ ONLY through `c=cos(θ/2)`.
The operator in `B_e=iTr[G^<X̂]` is `X̂=∂H/∂c` ⇒ `B_e = N ∂E_X/∂c` (the raw ⟨G^<X⟩), NOT `N∂E_X/∂θ`.
They differ by `dc/dθ=-½sin(θ/2)`. The LLG torque `j_θ=B_e/(4c)`, `τ_1=j_θ sinθ`, and the identity
`sinθ/(4c)=½sin(θ/2)=-dc/dθ` gives cleanly **`τ_1 = N ∂E_X/∂θ`** (physical, finite). Step 8 (c-operator
+ /(4c)·sinθ) was self-consistent; my first Step 9 wrongly fed `∂E_X/∂θ` into the `/(4c)·sinθ` pipeline
(double-counted the chain rule, spurious `sin(θ/2)/2`). Now: `τ_1=N∂E_X/∂θ` (0.68 µeV at π/2, s=0.7),
raw `B_e=N∂E_X/∂c` (diverges at FM as `1/sin(θ/2)`).

**cc/vv split** (per exciton, of the 7.50 meV/rad):
- conduction `+Σ|Y|²∂ε_c/∂θ = +26.72 meV/rad`
- valence `−Σ|Y|²∂ε_v/∂θ = −19.22 meV/rad`  (`|vv/cc|=0.72`, since `Jpd_v>Jpd_c`)
- **VALENCE IS NOT NEGLIGIBLE** — partial cancellation gives net 7.50. Conduction-only (Step 8)
  gives 26.7 ⇒ ~3.6× too large.

**Why Step 7 (spectral Σ_vv) is wrong:** I built it as a naive c↔v mirror with an in-gap satellite at
`ε_c(P+q)−Ω ≈ +0.65 eV` and guessed weights (`f_c F^λ`). The exact benchmark shows the exciton-induced
valence weight is a **two-hole shake-up BELOW the valence qp** (`≈ε_v−E_bind`), and — more importantly —
**the satellite is irrelevant for the torque** (equal-time ⇒ occupations only). Discard the Step 7
spectral picture; use the occupation route.

## Scope / approximation (asked 2026-07-09)
NOT assuming N=1. The BSE describes the **structure** of one e-h pair (the exciton wavefunction Y, E_X);
the **population** N enters linearly (torque ∝ N) via the occupations `n_c=N|Y|²`, `n_v=−N|Y|²`. The real
approximation is **dilute / independent excitons** (inter-exciton distance ≫ Bohr radius; our n_c~1e-4..1e-2/cell,
~1% bound ⇒ well satisfied). Population feedback is included ONLY through phase-space filling (occupations in
the excited BSE, Step 3 → blueshift, N_bright(n_c)). NEGLECTED: exciton-exciton scattering, biexcitons,
dynamical screening from the population, Mott transition (high density). `∂E_X/∂θ` uses the equilibrium
(N→0) exciton wavefunction — mild simplification (could recompute Y at each n_c for higher-order density effect).
Why `n_c=|Y|²`: from `|X⟩=Σ_k Y_k c†_k v_k|GS⟩`, `⟨c†_k c_k⟩=|Y_k|²` exactly (definition of the exciton amplitude).
Consistent with the 2016/2026 incoherent-exciton framework (`N_λQ=δ_Q0|ρ_λ|²+N^inc`).

## NEXT STEPS
1. ✅ DONE (Step 9): full cc+vv torque `B_e=N∂E_X/∂θ`, LLG `j_θ`,`τ_1`. Step 7 spectral discarded,
   Step 8 superseded. (Could still delete/annotate the wrong Step 7 cells 50-51 to avoid confusion.)
2. **Damping `η^{LL'}`** — the real remaining physics, NOT YET STARTED. Genuinely needs the
   dynamical/retarded part: `η ~ ∂_ω Σ^R` (non-adiabatic correction ∝ ∂_t M), N-independent, nonlocal
   non-Gilbert. NOT equal-time (unlike the torque) → the Σ^R/G^R spectral machinery (Steps 5-6, now
   using the exact vertex, item 3) IS needed here. This is where the LLG paper's `η_θ^{LL'}` comes
   from. Cross-check target: Brennan et al. 2026 Nature Materials (see top of file) — expect an
   angle-dependent, population-nonlinear damping (`−μXX'`-like), not a plain N-independent Gilbert
   term; reconcile with `η^{LL'}` being N-independent but θ-dependent.
5. **Sawtooth/limit-cycle search (2026-07-24, extensive, all 4 mechanisms closed):** full
   writeup in `scratchpad/SAWTOOTH_SEARCH_NOTES.md`. Population-lag (adiabatic + non-adiabatic),
   q-resolved dispersion, bright-CT Landau-Zener, AND bright-DARK Landau-Zener all ruled out
   with concrete evidence. The dark-exciton lead looked promising (splitting closes to exact
   zero at θ=π) but closed MORE definitively than bright-CT: `V_diag` (interaction in the
   {++,--,+-,-+} channel basis) is EXACTLY block-diagonal between bright{++,--} and dark{+-,-+}
   (checked numerically, cross-block terms = 0.0) — a symmetry-protected true crossing with zero
   coupling, not just a small one. Conclusion: no mechanism within this 4-band minimal model
   produces the sign-changing structure needed for a sawtooth limit cycle; would need physics
   beyond the model (more bands, spin-orbit, lower symmetry) to have any chance.
3. ✅ DONE (2026-07-24): exact vertex from Stefanucci-Perfetto 2026 Eq.6 (`K_I=(ΔE_I−Ω_λ)Y_I`,
   `vtx=|Σ_I K_I|²`) implemented in `bse_q_vertex`, replacing the calibrated `Ccal=V̄²` throughout
   Steps 5b, 6, 8 (cells 27,29,32,34,37,40,43,53). Does NOT change the Step 9 torque number itself
   (that route never used `Ccal`/`coh`, only `Fλ`) — it matters for the self-energy/spectral-function
   family and will matter for the Step 10 damping. Verified standalone
   (`scratchpad/vertex_check.jl`-style check): exact `vtx` ≈2.0× old `Ccal·coh` at q=0 bright mode;
   sum rule `∫A/2π` moves from the old (tuned) 0.999 to an untuned 0.87 — flagged as a real, expected
   residual from other truncations, not a new bug. Notebook cells re-run pending (outputs cleared).
4. Connect the absolute `j_θ(θ)`, `τ_1(θ)` numbers to LLG Eq 1 / the paper's DMRG results.

## Scratchpad (verification scripts, `.../scratchpad/`)
`calib_check.jl` (Σ_cc/G_cc + `exc`,`qs`,`μc`,`μv`,`Ccal`,`η`,`Lη`), `single_exc_valence.jl`
(THE benchmark: n_c=|Y|², n_v=1−|Y|², HF torque ratio 1.0000, cc/vv split +26.7/−19.2), `verify2.jl`
(full torque θ-sweep, cc/vv split + N_bright), `verify3.jl` (multi-density torque, ∝n_c² check),
`valence_check2.jl` (the wrong spectral Σ_vv — kept for reference), `bse2016_check.jl`, `step2/3_check.jl`.
NOTE: matplotlib mathtext does NOT support `\tfrac`/`\dfrac` (only `\frac`); OK in Jupyter markdown (MathJax).
