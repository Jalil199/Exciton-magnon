# Check whether the finite-q exciton dispersion Omega_X(q,theta) has a SIGN-CHANGING
# d(Omega_X)/d(theta) across q -- the missing ingredient for a Kambersky-like damping formula
# from the dressed exciton self-energy route (as opposed to the population-lag ansatz).
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

function bse_q_Omega(θ, q; μc, μv, Tx)
    d=transitions_q(θ,q;μc=μc,μv=μv,Tx=Tx); fI=d.fv.-d.fc; sf=sqrt.(abs.(fI)); σz=sign.(fI)
    Kt=zeros(ComplexF64,length(fI),length(fI)); for α in 1:2; u=sf.*d.M[α,:]; Kt.+=(Vα[α]/L).*(u*u'); end
    H=Diagonal(ComplexF64.(d.ΔE)).-Diagonal(ComplexF64.(σz))*Kt; F=eigen(H); Ω=real.(F.values)
    msk=abs.(fI).>1e-4
    wF =[msk[i] ? d.fc[i]*(1-d.fv[i])/abs(fI[i]) : 0.0 for i in eachindex(fI)]
    Fλ=real.((abs2.(F.vectors))'*wF)
    i0 = argmin(Ω)                # lowest = bright exciton branch at this q
    return (Ω=Ω[i0], Fλ=Fλ[i0])
end

Tx=0.10; s=0.7; μc=Δ_0/2+s; μv=Δ_0/2-s
θ0 = 0.6070*π    # theta_eq0 from the real-parameter macrospin calc
dθ = 0.02
Nq = 32; qs = collect(range(-π, stop=π-2π/Nq, length=Nq))

Om_p = zeros(Nq); Om_m = zeros(Nq); F_0 = zeros(Nq)
print("computing Omega_X(q, theta0 +- dtheta) over $(Nq) q-points ... "); flush(stdout)
Threads.@threads for i in 1:Nq
    q = qs[i]
    rp = bse_q_Omega(θ0+dθ, q; μc=μc, μv=μv, Tx=Tx)
    rm = bse_q_Omega(θ0-dθ, q; μc=μc, μv=μv, Tx=Tx)
    r0 = bse_q_Omega(θ0, q; μc=μc, μv=μv, Tx=Tx)
    Om_p[i] = rp.Ω; Om_m[i] = rm.Ω; F_0[i] = r0.Fλ
end
println("done")

dOmdtheta = (Om_p .- Om_m) ./ (2dθ)
println("\n q/π      Omega_X(q,θ0)   dΩ/dθ (eV/rad)   F_q (occupation)")
for i in 1:Nq
    Om0 = (Om_p[i]+Om_m[i])/2
    @printf("%+.3f    %.4f          %+.5f          %.4e\n", qs[i]/π, Om0, dOmdtheta[i], F_0[i])
end

pos = count(>(0), dOmdtheta); neg = count(<(0), dOmdtheta)
@printf("\ndΩ/dθ: %d positive, %d negative out of %d q-points\n", pos, neg, Nq)
if pos>0 && neg>0
    println("SIGN CHANGE CONFIRMED across q -- Kambersky-like formula is viable")
else
    println("NO sign change across q -- same conclusion as q=0 alone")
end
