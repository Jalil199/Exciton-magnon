# Step 10 (standalone) — population-lag damping η_θ(θ) from a finite exciton lifetime τ_X.
#
# Physical picture: N does not track N_eq(θ(t)) instantaneously; it relaxes with rate 1/τ_X
# (a real, dissipative lifetime — phonons/recombination/exciton-exciton scattering — NOT part
# of the closed Ĥ). Adiabatic elimination (valid since τ_X ~ 260-800 ps ≪ magnon period ~37 ns,
# Brennan 2026) gives, to first order in θ̇:
#
#   N(t) ≈ N_eq(θ(t)) − τ_X · (dN_eq/dθ) · θ̇(t)
#   τ_1(t) = N(t)·∂_θE_X(θ) = [N_eq·∂_θE_X](θ)  −  θ̇(t)·[τ_X·(dN_eq/dθ)·∂_θE_X](θ)
#                              \_______________/          \_______________________/
#                              reactive (Step 9)                 η_θ(θ)  ≡ damping
#
# In natural units (ħ=1): τ_X[1/eV] = τ_X[s]/ħ[eV·s]; η_θ = τ_X[1/eV]·(dN_eq/dθ)·∂_θE_X[eV]
# comes out DIMENSIONLESS — directly analogous to the standard Gilbert α.
#
# Reuses exactly the verified Step 9 (bright_split_t9, ∂_θE_X) and Step 9b (N_bright(θ)) machinery.

using LinearAlgebra, Printf
BLAS.set_num_threads(1)   # avoid oversubscription: parallelize at the Julia-thread (θ-grid) level instead
println("Julia threads available: ", Threads.nthreads())

# ---------------- model setup (same params as ExcitonBathKernel.ipynb) ----------------
L = 200
γ_c = 0.40; γ_v = 0.40
Δ_0 = 2.0
Jpd_v = 0.069; Jpd_c = 0.031
V_intra = 1.708; V_inter = 0.691
ks = collect(range(-π, stop=π - 2π/L, length=L))

σ_0 = ComplexF64[1 0; 0 1]; σ_x = ComplexF64[0 1; 1 0]
ϵk(k, γ) = -2γ*cos(k)
H_valence(k, θ)    = (Jv=Jpd_v*cos(θ/2); -(ϵk(k,γ_v)+2γ_v+Jv)*σ_0 + Jv*σ_x)
H_conduction(k, θ) = (c=cos(θ/2); Jc=Jpd_c*c;
                      (ϵk(k,γ_c)+2γ_c+Δ_0+(Jpd_c+Jpd_v)*(1-c)+Jc)*σ_0 - Jc*σ_x)
FD(E, μ, Tx) = 1/(exp(clamp((E-μ)/Tx,-60,60))+1)

Γ_sector = [σ_0/√2, σ_x/√2]
Vα = [V_intra, V_inter]

function transitions_q(θ, q; μc, μv, Tx, ks=ks)
    nT=4*length(ks); ΔE=zeros(nT); fc=zeros(nT); fv=zeros(nT); M=zeros(ComplexF64,2,nT); t=0
    for k in ks
        Fv=eigen(Hermitian(H_valence(k,θ))); Fc=eigen(Hermitian(H_conduction(k+q,θ)))
        mα=[Fc.vectors'*Γ*Fv.vectors for Γ in Γ_sector]
        for a in 1:2, b in 1:2
            t+=1; ΔE[t]=Fc.values[a]-Fv.values[b]
            fc[t]=FD(Fc.values[a],μc,Tx); fv[t]=FD(Fv.values[b],μv,Tx)
            for α in 1:2; M[α,t]=mα[α][a,b]; end
        end
    end
    return (; ΔE, fc, fv, M)
end

function bse_q_vertex(θ, q; μc, μv, Tx)
    d=transitions_q(θ,q;μc=μc,μv=μv,Tx=Tx); fI=d.fv.-d.fc; sf=sqrt.(abs.(fI)); σz=sign.(fI)
    Kt=zeros(ComplexF64,length(fI),length(fI)); for α in 1:2; u=sf.*d.M[α,:]; Kt.+=(Vα[α]/L).*(u*u'); end
    H=Diagonal(ComplexF64.(d.ΔE)).-Diagonal(ComplexF64.(σz))*Kt; F=eigen(H); Ω=real.(F.values)
    msk=abs.(fI).>1e-4
    wF =[msk[i] ? d.fc[i]*(1-d.fv[i])/abs(fI[i]) : 0.0 for i in eachindex(fI)]
    Fλ=real.((abs2.(F.vectors))'*wF)
    return (; Ω, Fλ)
end

# ---------------- Step 9: q=0 bright exciton E_X(θ), Hellmann-Feynman ∂_θE_X ----------------
function bright_split_t9(θ)
    nT=4L; EC=zeros(nT); EV=zeros(nT); Mm=zeros(ComplexF64,2,nT); t=0
    for k in ks
        Fv=eigen(Hermitian(H_valence(k,θ))); Fc=eigen(Hermitian(H_conduction(k,θ)))
        mα=[Fc.vectors'*Γ*Fv.vectors for Γ in Γ_sector]
        for a in 1:2,b in 1:2; t+=1; EC[t]=Fc.values[a]; EV[t]=Fv.values[b]
            for α in 1:2; Mm[α,t]=mα[α][a,b]; end; end
    end
    Kt=zeros(ComplexF64,nT,nT); for α in 1:2; u=Mm[α,:]; Kt.+=(Vα[α]/L).*(u*u'); end
    F=eigen(Diagonal(ComplexF64.(EC.-EV)).-Kt); i0=argmin(real.(F.values))
    (; EX=real(F.values[i0]))
end

# ---------------- N_bright(θ) at the working excited state (s=0.7, matches Step 9b) ----------
Tx = 0.10; s = 0.7; μc = Δ_0/2 + s; μv = Δ_0/2 - s
NqD = 18; qsD = collect(range(-π, stop=π-2π/NqD, length=NqD))
Nbr(θ) = sum(begin e=bse_q_vertex(θ,q; μc=μc, μv=μv, Tx=Tx); e.Fλ[argmin(e.Ω)] end for q in qsD)/NqD

# ---------------- θ sweep: E_X(θ), N_eq(θ), and both θ-derivatives ----------------
θg = collect(range(0.08π, 0.92π, length=18)); dθ = 1e-4

nθ = length(θg)
EXg  = zeros(nθ); dEXg = zeros(nθ); Ng = zeros(nθ); dNg = zeros(nθ)

println("θ sweep: E_X(θ) (parallel over $(nθ) points, $(Threads.nthreads()) threads) ..."); flush(stdout)
Threads.@threads for i in 1:nθ
    θ = θg[i]
    EXg[i]  = bright_split_t9(θ).EX
    dEXg[i] = (bright_split_t9(θ+dθ).EX - bright_split_t9(θ-dθ).EX)/(2dθ)
end
println("  E_X(θ) done")

println("θ sweep: N_eq(θ) (parallel over $(nθ) points) ..."); flush(stdout)
Threads.@threads for i in 1:nθ
    θ = θg[i]
    Ng[i]  = Nbr(θ)
    dNg[i] = (Nbr(θ+dθ) - Nbr(θ-dθ))/(2dθ)
end
println("  N_eq(θ) done")

τ1g = Ng .* dEXg     # reactive torque (Step 9 result, sanity cross-check)

# ---------------- τ_X → natural units, η_θ(θ) ----------------
ħ_eVs = 6.582119569e-16      # eV·s
τX_ps_list = [260.0, 500.0, 800.0]   # Brennan's measured range

println("\nθ/π    E_X(eV)   N_eq       ∂θE_X(meV)  dN/dθ        τ1(meV)")
for i in eachindex(θg)
    @printf("%.3f  %.4f   %.3e   %8.4f   %+10.3e   %+8.4f\n",
            θg[i]/π, EXg[i], Ng[i], 1e3*dEXg[i], dNg[i], 1e3*τ1g[i])
end

println("\nη_θ(θ) = τ_X[1/eV]·(dN/dθ)·∂θE_X   (dimensionless, natural units ħ=1)")
println("θ/π   " * join([@sprintf("η(τX=%.0fps)", t) for t in τX_ps_list], "   "))
ηg = Dict(τX => Float64[] for τX in τX_ps_list)
for i in eachindex(θg)
    row = @sprintf("%.3f", θg[i]/π)
    for τX in τX_ps_list
        τX_natural = (τX*1e-12) / ħ_eVs         # 1/eV
        η = τX_natural * dNg[i] * dEXg[i]
        push!(ηg[τX], η)
        row *= @sprintf("   %+10.4e", η)
    end
    println(row)
end

# ---------------- full-precision CSV (avoid the earlier hand-transcription rounding bug) ----------------
csvpath = joinpath(@__DIR__, "step10_damping_data.csv")
open(csvpath, "w") do io
    println(io, "theta_over_pi,EX_eV,Neq,dEXdtheta_eV,dNdtheta,tau1_eV,eta_260ps,eta_500ps,eta_800ps")
    for i in eachindex(θg)
        etas = [ (τX*1e-12/ħ_eVs) * dNg[i] * dEXg[i] for τX in τX_ps_list ]
        println(io, join([θg[i]/π, EXg[i], Ng[i], dEXg[i], dNg[i], τ1g[i], etas...], ","))
    end
end
println("\nwrote full-precision data to $csvpath")

# sign check
for τX in τX_ps_list
    v = ηg[τX]
    if any(v .> 0) && any(v .< 0)
        i = findfirst(k -> sign(v[k]) != sign(v[1]), eachindex(v))
        @printf("\nτ_X=%.0fps: η_θ CHANGES SIGN near θ/π≈%.3f (η=%.3e → %.3e)\n",
                τX, θg[i-1]/π, v[i-1], v[i])
    else
        @printf("\nτ_X=%.0fps: η_θ keeps constant sign (%s) over the whole sweep, max|η|=%.3e\n",
                τX, v[1]>0 ? "positive" : "negative", maximum(abs.(v)))
    end
end
