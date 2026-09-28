# Exciton-induced torque and damping on layer magnetizations in CrSBr

Theory project on how a photoexcited **exciton population** acts back on the
classical magnetic moments of a layered magnetic semiconductor (CrSBr, A-type
antiferromagnet). The target is a generalized Landau-Lifshitz equation for the
macrospin of each layer, derived from the microscopic electron-hole problem:

```
∂_t M_L = M_L × [ B_L^cl + j_θ^L M_L̄ + Σ_L' η_θ^{LL'} ∂_t M_L' ]
```

with a torque `j_θ` set by the number of excitons and a nonlocal, non-Gilbert
damping `η^{LL'}`. It accompanies the draft *Landau-Lifshitz equation with
nonlocal damping due to exciton bath* (Garcia-Gaitan, Varela-Manjarres,
Nikolić) and the experiment *Excitonic spin torque in a magnetic semiconductor*
(Brennan et al., Nature Materials 2026). The draft is still being validated:
the notebooks here re-derive and check each piece independently.

## Physical idea

- The two layers carry classical moments `M_1`, `M_2` at relative angle θ.
  Interlayer electron hopping is spin-filtered, so it scales as `cos(θ/2)`
  (Anderson-Hasegawa): allowed in the FM configuration, blocked in the AFM one.
- Excitons are bound electron-hole pairs obtained from the Bethe-Salpeter
  equation, `D = (V⁻¹ − Π)⁻¹`. Their energy `E_X(θ)` depends on the angle
  through the interlayer hopping. This is the 15 meV redshift of the CrSBr
  A exciton between AFM and FM.
- Integrating out the electrons gives an effective action for the moments.
  Its first variation is the torque; the second variation is the exciton
  contribution to the spin susceptibility (exchange, inertia, damping).
- With a real exciton population `N`, the torque is an equal-time quantity:
  `τ = N ∂E_X/∂θ ∝ N sin θ`, zero at FM and AFM, maximal at θ = π/2.
- Damping needs the dynamical (retarded) part of the exciton kernel and is the
  main open problem.

## Model

Minimal four-band model of Heissenbüttel et al. (PRB 111, 075107): two layers
× (valence, conduction), one-dimensional bands, intralayer and interlayer
electron-hole attraction, interlayer hopping `Jpd cos(θ/2)`. Parameters are
calibrated to the CrSBr A exciton (1.36 eV in AFM, 15 meV redshift to FM,
`V_intra/V_inter` matching ab initio). The symmetric layer sector is the bright
exciton, the antisymmetric one the dark exciton. Excited-state populations
follow the nonequilibrium BSE of Perfetto-Stefanucci (quasi-Fermi occupations,
exciton weights, self-energies, dressed Green's functions).

## Repository layout

- `RPA-exciton/Polarized_model/` — main line of work (Julia notebooks, calibrated model).
  - `Static_D.ipynb` — bare bubble, dressed exciton propagator, calibration.
  - `BrightExcitonMode.ipynb` — single-mode reduction of the bright exciton, `E_X(θ)`, `E_X(Q)`.
  - `ExcitonSpinSusceptibility.ipynb` — equilibrium exciton contribution to the spin susceptibility (one-D and two-D terms, slow-spin moments J, α, inertia).
  - `DrivenExcitonTorque.ipynb` — coherently driven exciton: detuning-controlled damping / anti-damping.
  - `ExcitonTorque.ipynb` — first-variation torque, per-exciton coefficient `∂E_X/∂θ`.
  - `ExcitonBathKernel.ipynb` — step-by-step incoherent-population framework (Perfetto-Stefanucci 2016 / Stefanucci-Perfetto 2026) up to the full torque `N ∂E_X/∂θ`.
  - `HANDOFF_exciton_torque.md` — detailed status, verified numbers, known wrong cells.
  - `scratchpad/` — standalone verification scripts, two-macrospin LLG integrators, and `SAWTOOTH_SEARCH_NOTES.md` (search for a self-sustained limit cycle; all mechanisms tried within this model are ruled out).
- `RPA-exciton/*.ipynb` — first generation (spin-channel BSE with s-d coupling, static susceptibility, light-induced torque). Superseded by the polarized model.
- `DMRG/tJ_DMRG.jl`, `DMRG/Dyn_tJ_DMRG.jl` — original bilayer t-J model with DMRG (ITensors) and classical spins. Starting point of the project, no longer active.

## Status (September 2026)

Verified: torque `τ = N ∂E_X/∂θ ∝ N sin θ` (7.5 meV/rad per exciton at θ = π/2,
valence and conduction contributions partially cancel); exciton number
`N ∝ n_c²` (mass action); no intrinsic Gilbert damping from a gapped virtual
exciton in equilibrium; exact sideband vertex implemented.

Open: the microscopic damping `η^{LL'}` (the draft's derivation is unfinished;
the population-lag mechanism only gives ordinary damping, never the
anti-damping seen in experiment); the four-band model has no channel that
produces the sign-changing damping needed for the experimental limit cycle.

## Main references

- Heißenbüttel et al., PRB 111, 075107 (2025) — four-band CrSBr model.
- Perfetto, Sangalli, Marini, Stefanucci, PRB 94, 245303 (2016); Stefanucci, Perfetto, arXiv:2601.16786 (2026) — nonequilibrium BSE, excited Green's functions.
- Reyes-Osorio, Nikolić, PRB 109, 024413 (2024); PRL 135, 246701 (2025) — LLG from Schwinger-Keldysh field theory.
- Bae et al., Nature 609, 282 (2022); Diederich et al., Nat. Nanotechnol. (2022, 2025); Datta et al., Nat. Mater. 24, 1027 (2025); Brennan et al., Nat. Mater. (2026) — CrSBr exciton-magnon experiments.
- Scheie et al., Adv. Sci. 9, 2202467 (2022) — magnetic Hamiltonian of CrSBr.
