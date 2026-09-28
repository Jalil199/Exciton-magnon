# Exciton-magnon coupling in CrSBr

Theory project on the interplay between excitons and magnetic order in the
layered antiferromagnetic semiconductor CrSBr. The question is how a
photoexcited exciton population acts back on the magnetization of the layers,
and how that feedback can be written as effective torques and damping in the
equation of motion of the classical moments.

## Physical picture

- Each layer carries a classical moment; the relative angle between layers is
  the relevant magnetic variable.
- Interlayer electron hopping is spin-filtered, so the electronic structure,
  and therefore the exciton energy, depends on that angle.
- Excitons are obtained from the Bethe-Salpeter equation on top of a minimal
  multiband model calibrated to the CrSBr A exciton.
- Integrating out the electronic degrees of freedom gives an effective action
  for the moments. Its variations yield the exciton contribution to the torque
  and to the spin susceptibility, including the dynamical (damping-like) parts.
- Nonequilibrium exciton populations are treated with the excited-state
  Green's-function framework of Perfetto and Stefanucci.

## Repository layout

- `RPA-exciton/Polarized_model/` — current line of work. Julia notebooks
  (ITensors-free, plain linear algebra + PyPlot) for the calibrated multiband
  model: exciton propagator, bright/dark modes, spin susceptibility, torque,
  and the incoherent-population kernel. `HANDOFF_exciton_torque.md` keeps the
  running notes; `scratchpad/` holds standalone verification scripts and
  classical-spin dynamics tests.
- `RPA-exciton/*.ipynb` — first generation of the calculation (spin-channel
  BSE with s-d coupling, static susceptibility, light-induced torque).
  Superseded by the polarized model.
- `DMRG/` — original bilayer t-J model with DMRG (ITensors) and classical
  spins. Starting point of the project, no longer active.

## Main references

- Heißenbüttel et al., PRB 111, 075107 (2025) — multiband CrSBr model.
- Perfetto, Sangalli, Marini, Stefanucci, PRB 94, 245303 (2016); Stefanucci,
  Perfetto, arXiv:2601.16786 (2026) — nonequilibrium BSE and excited-state
  Green's functions.
- Reyes-Osorio, Nikolić, PRB 109, 024413 (2024); PRL 135, 246701 (2025) — LLG
  from Schwinger-Keldysh field theory.
- Bae et al., Nature 609, 282 (2022); Diederich et al., Nat. Nanotechnol.
  (2022, 2025); Datta et al., Nat. Mater. 24, 1027 (2025); Brennan et al.,
  Nat. Mater. (2026) — CrSBr exciton-magnon experiments.
- Scheie et al., Adv. Sci. 9, 2202467 (2022) — magnetic Hamiltonian of CrSBr.
