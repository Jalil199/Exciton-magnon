# Precompute, on a theta grid, the reduced 2-level (bright, CT exciton) Hamiltonian data needed
# for time-dependent Schrodinger + classical-spin coupled dynamics (Landau-Zener-like):
#   E_bright(theta), E_CT(theta)                      -- diagonal energies (Hellmann-Feynman, exact)
#   Xdiag_bright(theta) = dE_bright/dtheta,  Xdiag_CT  -- diagonal torque matrix elements
#   Xoff(theta) = <bright|dH/dtheta|CT>                -- OFF-diagonal torque/coupling matrix element
# (real, since H(theta) is real-symmetric-representable here with a consistent eigenvector phase
# choice -- checked numerically below).
using LinearAlgebra, Printf
BLAS.set_num_threads(1)
println("Julia threads: ", Threads.nthreads())

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
Γ_sector = [σ_0/√2, σ_x/√2]
Vα = [V_intra, V_inter]

# Build the FULL 800x800 (q=0, equilibrium bright-sector) BSE Hamiltonian matrix at given theta
function H_BSE(θ)
    nT=4L; EC=zeros(nT); EV=zeros(nT); Mm=zeros(ComplexF64,2,nT); t=0
    for k in ks
        Fv=eigen(Hermitian(H_valence(k,θ))); Fc=eigen(Hermitian(H_conduction(k,θ)))
        mα=[Fc.vectors'*Γ*Fv.vectors for Γ in Γ_sector]
        for a in 1:2,b in 1:2; t+=1; EC[t]=Fc.values[a]; EV[t]=Fv.values[b]
            for α in 1:2; Mm[α,t]=mα[α][a,b]; end; end
    end
    Kt=zeros(ComplexF64,nT,nT); for α in 1:2; u=Mm[α,:]; Kt.+=(Vα[α]/L).*(u*u'); end
    return Diagonal(ComplexF64.(EC.-EV)).-Kt
end

θg = collect(range(0.10π, 0.90π, length=17))
dθ = 1e-4
nθ = length(θg)

Ebright = zeros(nθ); ECT = zeros(nθ)
dEb = zeros(nθ); dECT = zeros(nθ)
Xoff = zeros(ComplexF64, nθ)
Ybright_all = Vector{Vector{ComplexF64}}(undef, nθ)
YCT_all     = Vector{Vector{ComplexF64}}(undef, nθ)

println("Phase 1: diagonalizing H_BSE(θ) over $(nθ) theta points (parallel) ..."); flush(stdout)
Threads.@threads for i in 1:nθ
    θ = θg[i]
    F = eigen(Hermitian(H_BSE(θ)))
    Ybright_all[i] = F.vectors[:,1]; YCT_all[i] = F.vectors[:,2]
    Ebright[i] = F.values[1]; ECT[i] = F.values[2]
end
println("  done")

println("Phase 2: fixing eigenvector gauge/phase along theta (serial, sequential overlap) ..."); flush(stdout)
for i in 2:nθ
    if real(dot(Ybright_all[i-1], Ybright_all[i])) < 0
        Ybright_all[i] .*= -1
    end
    if real(dot(YCT_all[i-1], YCT_all[i])) < 0
        YCT_all[i] .*= -1
    end
end
println("  done")

println("Phase 3: computing dH/dtheta and matrix elements with gauge-fixed vectors (parallel) ..."); flush(stdout)
Threads.@threads for i in 1:nθ
    θ = θg[i]
    Hp = H_BSE(θ+dθ); Hm = H_BSE(θ-dθ)
    dH = (Hp .- Hm) ./ (2dθ)
    Yb = Ybright_all[i]; Yc = YCT_all[i]
    dEb[i]  = real(Yb' * dH * Yb)
    dECT[i] = real(Yc' * dH * Yc)
    Xoff[i] = Yb' * dH * Yc
end
println("  done")

@printf("\n%6s  %10s %10s  %10s %10s  %12s %12s  %8s\n",
        "θ/π","E_bright","E_CT","dEb/dθ","dECT/dθ","Re(Xoff)","Im(Xoff)","gap(meV)")
for i in 1:nθ
    gap = 1e3*(ECT[i]-Ebright[i])
    @printf("%6.3f  %10.4f %10.4f  %+10.5f %+10.5f  %+12.5e %+12.5e  %8.2f\n",
            θg[i]/π, Ebright[i], ECT[i], dEb[i], dECT[i], real(Xoff[i]), imag(Xoff[i]), gap)
end

println("\nmax |Im(Xoff)| = ", maximum(abs.(imag.(Xoff))), "  (should be ~0 if real-representable)")

# write CSV for the Python time-evolution
open(joinpath(@__DIR__, "twolevel_data.csv"), "w") do io
    println(io, "theta,Ebright_eV,ECT_eV,dEb_dtheta_eV,dECT_dtheta_eV,Xoff_real_eV,Xoff_imag_eV")
    for i in 1:nθ
        println(io, join([θg[i], Ebright[i], ECT[i], dEb[i], dECT[i], real(Xoff[i]), imag(Xoff[i])], ","))
    end
end
println("wrote twolevel_data.csv")
