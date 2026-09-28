using LinearAlgebra, Printf
L = 200
γ_c=0.40; γ_v=0.40; Δ_0=2.0; Jpd_v=0.069; Jpd_c=0.031; V_intra=1.708; V_inter=0.691
ks = collect(range(-π, stop=π - 2π/L, length=L))
σ_0 = ComplexF64[1 0; 0 1]; σ_x = ComplexF64[0 1; 1 0]
ϵk(k, γ) = -2γ*cos(k)
H_valence(k, θ)    = (Jv=Jpd_v*cos(θ/2); -(ϵk(k,γ_v)+2γ_v+Jv)*σ_0 + Jv*σ_x)
H_conduction(k, θ) = (c=cos(θ/2); Jc=Jpd_c*c;
                      (ϵk(k,γ_c)+2γ_c+Δ_0+(Jpd_c+Jpd_v)*(1-c)+Jc)*σ_0 - Jc*σ_x)
Γ_sector = [σ_0/√2, σ_x/√2]
Vα = [V_intra, V_inter]

# BRIGHT: same-index (a=b) pairs only.  DARK: crossed-index (a!=b) pairs only.
function H_BSE_sector(θ; crossed::Bool)
    idxs = Tuple{Int,Int}[]
    for k in ks, a in 1:2, b in 1:2
        if crossed
            a != b && push!(idxs, (a,b))
        else
            a == b && push!(idxs, (a,b))
        end
    end
    nT = length(ks) * (crossed ? 2 : 2)   # 2 (a,b) combos per k in each case
    ΔE=zeros(nT); Mm=zeros(ComplexF64,2,nT); t=0
    for k in ks
        Fv=eigen(Hermitian(H_valence(k,θ))); Fc=eigen(Hermitian(H_conduction(k,θ)))
        mα=[Fc.vectors'*Γ*Fv.vectors for Γ in Γ_sector]
        for a in 1:2, b in 1:2
            keep = crossed ? (a != b) : (a == b)
            keep || continue
            t+=1; ΔE[t]=Fc.values[a]-Fv.values[b]
            for α in 1:2; Mm[α,t]=mα[α][a,b]; end
        end
    end
    Kt=zeros(ComplexF64,nT,nT); for α in 1:2; u=Mm[α,:]; Kt.+=(Vα[α]/L).*(u*u'); end
    return Diagonal(ComplexF64.(ΔE)).-Kt
end

println("θ/π   bright_numeric   dark_numeric(crossed)   [analytic dark target ~1.357-1.360]")
for θ in range(0.1π, 0.9π, length=9)
    Hb = H_BSE_sector(θ; crossed=false)
    Hd = H_BSE_sector(θ; crossed=true)
    Eb = minimum(real(eigvals(Hermitian(Hb))))
    Ed = minimum(real(eigvals(Hermitian(Hd))))
    @printf("%.3f   %.4f          %.4f\n", θ/π, Eb, Ed)
end
