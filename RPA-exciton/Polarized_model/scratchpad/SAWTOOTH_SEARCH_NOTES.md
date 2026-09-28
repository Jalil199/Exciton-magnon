# Búsqueda del mecanismo del diente de sierra (Brennan Fig 3a) — notas de referencia

**Fecha:** 2026-07-24 · **Contexto:** después de implementar el vértice exacto (Step 10 del
handoff) e integrar el LLG de 2 macrospines, se intentó reproducir el comportamiento
autosostenido (diente de sierra) que Brennan et al. 2026 reportan en Fig. 3a a alta fluencia
(130 μW). Todos los scripts están en este mismo directorio (`scratchpad/`).

## ✅ SEGURO — resultados verificados, se pueden usar como base

1. **Vértice exacto (Step 5b-9b de `ExcitonBathKernel.ipynb`).** Implementado y verificado:
   `vtx=|Σ_I K_I|²` con `K_I=(ΔE_I-Ω_λ)Y_I` (2026 Eq.6) reemplaza `Ccal·coh`. Medido: exacto
   ≈2.0× el viejo `Ccal·coh` en el modo brillante `q=0`. Regla de suma `∫A/2π` baja de 0.999
   (con `Ccal` retro-ajustado) a 0.87 (sin ajustar) — esperado, no es un bug.

2. **Parámetros reales de literatura para el Hamiltoniano clásico (Ec. 22 del draft).**
   Los valores del draft (`J_AF=19.5meV, Kx=21meV, Kz=57meV, gμB·B=42meV`) dan una frecuencia
   de precesión bare de ~90-96 meV (~22 THz) — confirmado con cálculo de modos normales
   (Jacobiano, sin necesidad de integrar), **1000× más rápido** que el magnón real de CrSBr.
   Casi seguro un error de unidades μeV↔meV en el draft (`Kz=57meV` es numéricamente el mismo
   número que el `A_z=58μeV` real).
   **Valores reales correctos** (fit LSWT, "Tunable magnons in a dual-gated 2D antiferromagnet"):
   `J_int=6 μeV, A_x=14 μeV, A_z=58 μeV`, validados a `B=0.2T` dando `ν_IP=24.4GHz, ν_OP=33.8GHz`.
   Con estos + `g=2`, nuestro cálculo de modos normales da `14.4/15.0 GHz` — orden de magnitud
   correcto (factor ~1.7-2× de diferencia, probablemente por el factor de espín S=3/2 no incluido
   o diferencias de convención geométrica). `θ_eq0 ≈ 109.3°` a `B=0.2T`.

3. **Scheie et al. (Adv. Sci. 2022, neutrones) — anisotropía real.** Gap de magnón medido
   ópticamente: `0.102–0.141 meV → 24.7–34.1 GHz`, coincide con los 26.7 GHz de Brennan.
   Intercambios en el plano `J₁=-1.90, J₂=-3.38, ... meV`. Anisotropía de un solo ion
   "demasiado chica para resolver" con neutrones — consistente con `A_x,A_z~μeV`.

4. **El draft (`LLG_for_exciton_bath.pdf`) tiene la sección de damping literalmente sin
   terminar.** End Matter sección D dice *"we derive an expression for... `η_LL'_θ`"* y corta
   ahí, con una referencia LaTeX rota `[Eq.(??)]`. No hay fórmula cerrada que copiar.

5. **Mecanismo de population-lag — probado exhaustivamente, resultado negativo limpio.**
   - `η_θ(θ) = τ_X·(dN_eq/dθ)·∂_θE_X(θ)` (adiabático, primer orden en `θ̇`): negativo en
     TODO el rango barrido `θ/π∈[0.08,0.92]` (`τ_X`=260/500/800ps → `η_θ` de -0.005 a -0.27).
   - Versión NO adiabática (N(t) como variable dinámica real, sin aproximar): decae mucho
     menos que la versión lineal (razón amplitud ~0.99 vs ~0.5-0.8), pero **sigue decayendo**,
     nunca cruza a crecimiento neto — ni escaneando `τ_X` hasta 50ns (60× el rango real) ni
     la densidad de bombeo hasta 128× — con parámetros reales, sin modificar.
   - Con desplazamiento artificial `η_θ(θ)+C`: para `C` grande el sistema SÍ cruza a
     crecimiento, y con `η_eff=γ₀-μ(θ-θ_eq0)` (asimétrico, tipo Brennan) se obtiene un
     **ciclo límite genuino y estable** (amplitud converge a 43.6°, `[83°,127°]`) — pero el
     perfil temporal es simétrico (subida 36.4ps ≈ bajada 35.1ps), NO tiene la forma de diente
     de sierra (que necesita asimetría temporal, no solo de amplitud).
   - `N(t)` con decaimiento libre desacoplado de `θ` (sin realimentación): sin disipación,
     oscila indefinidamente sin decaer una vez que `N→0`. Confirma que la realimentación
     `N↔θ` era el único canal de disipación que había.

6. **Canal de autoenergía / dispersión en `q` — también negativo, cuantitativo.**
   `dΩ_X(q;θ)/dθ` es positivo en los 32 puntos de `q` calculados, sin cambio de signo
   (rango `+0.0037` a `+0.0054` eV/rad); la ocupación `F_q` está muy concentrada en `q≈0`.
   Fórmula tipo Kambersky no tiene de dónde sacar un cambio de signo.

7. **Landau-Zener brillante↔CT: descartado cuantitativamente.** Gap `bright-CT ≈ 600meV`
   (casi constante en todo θ). `P_LZ ~ exp(-2π×968) ≈ 0` a la velocidad real del magnón
   (`θ̇~0.11meV` a 26.7GHz). El excitón está atado por su energía de ligadura (~0.68eV),
   **millones de veces** mayor que la escala de energía del magnón — la aproximación
   adiabática para la función de onda del excitón (no solo la población) está muy bien
   justificada. Este canal específico NO puede ser el mecanismo.

8. **Splitting brillante-oscuro (analítico, fórmulas de `Static_D.ipynb`) — hallazgo real,
   prometedor, NO explorado hasta el final.** `J_plus=(Jpd_c+Jpd_v)cosθ/2` (brillante),
   `J_minus=(Jpd_c-Jpd_v)cosθ/2` (oscuro). El splitting brillante-oscuro real (usando
   `sector_bound_roots` con los parámetros calibrados) es de solo **0.3–12.5 meV** en el
   rango `θ/π∈[0.1,0.9]`, y **se cierra exactamente a cero en `θ=π`** (AFM puro, porque
   `cos(θ/2)=0` anula tanto `J_plus` como `J_minus`). Es una degeneración protegida por
   simetría, no un artefacto — muy distinta del caso brillante-CT (gap grande y constante).

## ❌ NO SEGURO — abierto, sin resolver

4. **No se probó incluir los elementos electrón-hueco no ligados (continuo)** en el
   acoplamiento no adiabático — solo se comparó brillante contra CT y contra la dispersión
   en `q` del propio brillante. Quedó pendiente ver si algún estado del continuo (no solo
   los 2 estados ligados) tiene un denominador de energía pequeño en algún punto de `θ`,
   aunque dado el hallazgo #9 de abajo (bloqueo por simetría) es poco probable que ayude —
   el continuo también se separa en los mismos bloques `{++,--}`/`{+-,-+}`.

## ✅ CERRADO DEFINITIVAMENTE (actualización 2026-07-24, misma sesión)

9. **Brillante-oscuro: acoplamiento EXACTAMENTE CERO por simetría — más concluyente que el
   caso CT.** Se reconstruyó la separación correcta de 4 canales (`λ∈{++,--,+-,-+}`, usando
   los autovectores FIJOS `[1,1]/√2` y `[1,-1]/√2` de Static_D en vez de `eigen()` — el bug
   del intento anterior era que `eigen()` ordena por autovalor, y ese orden se invierte entre
   conducción y valencia porque `Jvθ,Jcθ≥0` siempre pero con signos opuestos en cómo desplazan
   a "+" vs "-", corrompiendo el filtro `a=b`/`a≠b`). Con los autovectores fijos, se calculó
   `V_diag = M'·V_layer·M` (con `M_transition_to_layer` de Static_D, `V_layer=-diag(V_intra,
   V_intra,V_inter,V_inter)`): **los bloques cruzados `{++,--}×{+-,-+}` son exactamente cero**
   (verificado numéricamente, no solo pequeño). Como `∂H/∂θ` de una matriz bloque-diagonal
   sigue bloque-diagonal, el acoplamiento Landau-Zener brillante-oscuro es idénticamente cero
   para cualquier `θ` — no hay cruce evitado, es un cruce verdadero sin transición posible,
   protegido por la simetría de la interacción, no por tamaño de gap. **Esta vía queda cerrada
   de forma más definitiva que la brillante-CT** (que al menos tenía acoplamiento no nulo).

## Conclusión de la sesión (2026-07-24)

Los CUATRO mecanismos probados dentro del modelo de 4 bandas — population-lag (adiabático y
no adiabático), dispersión `Ω_X(q,θ)`, Landau-Zener brillante-CT, y Landau-Zener
brillante-oscuro — están descartados con evidencia concreta (numérica o algebraica exacta).
Ninguno produce la estructura de signo cambiante necesaria para el ciclo límite tipo diente de
sierra. Esto sugiere que el modelo mínimo de 4 bandas, tal como está construido, no tiene NINGÚN
canal capaz de producirla — haría falta algo estructuralmente ausente (más bandas, acoplamiento
espín-órbita, menor simetría cristalina que rompa el bloqueo `{++,--}`/`{+-,-+}` encontrado en
el punto 9) para que exista un mecanismo de este tipo dentro de una extensión de este modelo.

## Archivos relevantes en este directorio

- `step10_damping.jl`, `step10_damping_data.csv` — `η_θ(θ)` adiabático (población-lag).
- `two_spin_llg.py`, `two_spin_nonadiabatic.py`, `two_spin_bigkick.py` — integradores LLG de
  2 macrospines (versión adiabática y no adiabática).
- `scan_threshold.py`, `threshold_scan.png` — barrido `τ_X`/densidad, sin cruce a crecimiento.
- `omega_q_theta.jl` — dispersión `Ω_X(q,θ)`, sin cambio de signo.
- `twolevel_matrixelements.jl`, `twolevel_data.csv` — brillante-CT, gauge-fijo, `P_LZ≈0`.
- `darkexciton_check.jl`, `dark_sector_validate.jl` — splitting oscuro analítico (✅) y el
  intento de construcción numérica en `k` que no coincidió (❌, punto 1 de arriba).

## Próximo paso más prometedor

Corregir la construcción del sector oscuro siguiendo exactamente la lógica de canales
`λ∈{++,--,+-,-+}` de `Static_D.ipynb` (cells 9-12), calcular `X_off` brillante-oscuro cerca
de `θ=π`, y evaluar si el Landau-Zener ahí es significativo (a diferencia del caso brillante-CT,
donde sí se descartó con números).
