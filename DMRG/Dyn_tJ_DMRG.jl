using ITensors
using ITensorMPS
using DelimitedFiles
using LinearAlgebra

ITensors.enable_debug_checks()

function ITensors.space(::SiteType"2tJ";
                        conserve_qns=false, conserve_sz=false)
    if conserve_qns
        if conserve_sz
            return [QN(("Na", 0, -1), ("Nb", 0, -1), ("Sz",  0)) => 1,
                    QN(("Na", 1, -1), ("Nb", 0, -1), ("Sz", +1)) => 1, 
                    QN(("Na", 1, -1), ("Nb", 0, -1), ("Sz", -1)) => 1,
                    QN(("Na", 0, -1), ("Nb", 1, -1), ("Sz", +1)) => 1,
                    QN(("Na", 0, -1), ("Nb", 1, -1), ("Sz", -1)) => 1, 
                    QN(("Na", 1, -1), ("Nb", 1, -1), ("Sz", +2)) => 1,
                    QN(("Na", 1, -1), ("Nb", 1, -1), ("Sz",  0)) => 2, 
                    QN(("Na", 1, -1), ("Nb", 1, -1), ("Sz", -2)) => 1,
                   ]
        else
            return [QN(("Na", 0, -1), ("Nb", 0, -1)) => 1,
                    QN(("Na", 1, -1), ("Nb", 0, -1)) => 2, 
                    QN(("Na", 0, -1), ("Nb", 1, -1)) => 2,
                    QN(("Na", 1, -1), ("Nb", 1, -1)) => 4
                   ]
        end
    else
        return 9
    end
end

ITensors.val(::ValName"AEmpBEmp", ::SiteType"2tJ") = 1
ITensors.val(::ValName"AupBEmp", ::SiteType"2tJ")  = 2
ITensors.val(::ValName"AdnBEmp", ::SiteType"2tJ")  = 3
ITensors.val(::ValName"AEmpBup", ::SiteType"2tJ")  = 4
ITensors.val(::ValName"AEmpBdn", ::SiteType"2tJ")  = 5
ITensors.val(::ValName"AupBup",  ::SiteType"2tJ")  = 6
ITensors.val(::ValName"AupBdn",  ::SiteType"2tJ")  = 7
ITensors.val(::ValName"AdnBup",  ::SiteType"2tJ")  = 8
ITensors.val(::ValName"AdnBdn",  ::SiteType"2tJ")  = 9

ITensors.state(::StateName"AEmpBEmp", ::SiteType"2tJ") = [1.0, 0, 0, 0, 0, 0, 0, 0, 0]
ITensors.state(::StateName"AupBEmp",  ::SiteType"2tJ") = [0, 1.0, 0, 0, 0, 0, 0, 0, 0]
ITensors.state(::StateName"AdnBEmp",  ::SiteType"2tJ") = [0, 0, 1.0, 0, 0, 0, 0, 0, 0]
ITensors.state(::StateName"AEmpBup",  ::SiteType"2tJ") = [0, 0, 0, 1.0, 0, 0, 0, 0, 0]
ITensors.state(::StateName"AEmpBdn",  ::SiteType"2tJ") = [0, 0, 0, 0, 1.0, 0, 0, 0, 0]
ITensors.state(::StateName"AupBup",   ::SiteType"2tJ") = [0, 0, 0, 0, 0, 1.0, 0, 0, 0]
ITensors.state(::StateName"AupBdn",   ::SiteType"2tJ") = [0, 0, 0, 0, 0, 0, 1.0, 0, 0]
ITensors.state(::StateName"AdnBup",   ::SiteType"2tJ") = [0, 0, 0, 0, 0, 0, 0, 1.0, 0]
ITensors.state(::StateName"AdnBdn",   ::SiteType"2tJ") = [0, 0, 0, 0, 0, 0, 0, 0, 1.0]

function ITensors.op!(Op::ITensor, ::OpName"Nt", ::SiteType"2tJ", s::Index)
    Op[s' => 2, s => 2] = 1
    Op[s' => 3, s => 3] = 1
    Op[s' => 4, s => 4] = 1
    Op[s' => 5, s => 5] = 1
    Op[s' => 6, s => 6] = 2
    Op[s' => 7, s => 7] = 2
    Op[s' => 8, s => 8] = 2
    return Op[s' => 9, s => 9] = 2
end

function ITensors.op!(Op::ITensor, ::OpName"Caup", ::SiteType"2tJ", s::Index)
    Op[s' => 1, s => 2] = 1
    Op[s' => 4, s => 6] = -1
    return Op[s' => 5, s => 7] = -1
end

function ITensors.op!(Op::ITensor, ::OpName"Cadagup", ::SiteType"2tJ", s::Index)
    Op[s' => 2, s => 1] = 1
    Op[s' => 6, s => 4] = -1
    return Op[s' => 7, s => 5] = -1
end

function ITensors.op!(Op::ITensor, ::OpName"Cadn", ::SiteType"2tJ", s::Index)
    Op[s' => 1, s => 3] = 1
    Op[s' => 4, s => 8] = -1
    return Op[s' => 5, s => 9] = -1
end

function ITensors.op!(Op::ITensor, ::OpName"Cadagdn", ::SiteType"2tJ", s::Index)
    Op[s' => 3, s => 1] = 1
    Op[s' => 8, s => 4] = -1
    return Op[s' => 9, s => 5] = -1
end

function ITensors.op!(Op::ITensor, ::OpName"Cbup", ::SiteType"2tJ", s::Index)
    Op[s' => 1, s => 4] = 1
    Op[s' => 2, s => 6] = 1
    return Op[s' => 3, s => 8] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"Cbdagup", ::SiteType"2tJ", s::Index)
    Op[s' => 4, s => 1] = 1
    Op[s' => 6, s => 2] = 1
    return Op[s' => 8, s => 3] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"Cbdn", ::SiteType"2tJ", s::Index)
    Op[s' => 1, s => 5] = 1
    Op[s' => 2, s => 7] = 1
    return Op[s' => 3, s => 9] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"Cbdagdn", ::SiteType"2tJ", s::Index)
    Op[s' => 5, s => 1] = 1
    Op[s' => 7, s => 2] = 1
    return Op[s' => 9, s => 3] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"F", ::SiteType"2tJ", s::Index)
    Op[s' => 1, s => 1] = 1
    Op[s' => 2, s => 2] = -1
    Op[s' => 3, s => 3] = -1
    Op[s' => 4, s => 4] = -1
    Op[s' => 5, s => 5] = -1
    Op[s' => 6, s => 6] = 1
    Op[s' => 7, s => 7] = 1
    Op[s' => 8, s => 8] = 1
    return Op[s' => 9, s => 9] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"Naup", ::SiteType"2tJ", s::Index)
    Op[s' => 2, s => 2] = 1
    Op[s' => 6, s => 6] = 1
    return Op[s' => 7, s => 7] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"Nadn", ::SiteType"2tJ", s::Index)
    Op[s' => 3, s => 3] = 1
    Op[s' => 8, s => 8] = 1
    return Op[s' => 9, s => 9] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"Na", ::SiteType"2tJ", s::Index)
    Op[s' => 2, s => 2] = 1
    Op[s' => 3, s => 3] = 1
    Op[s' => 6, s => 6] = 1
    Op[s' => 7, s => 7] = 1
    Op[s' => 8, s => 8] = 1
    return Op[s' => 9, s => 9] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"Nbup", ::SiteType"2tJ", s::Index)
    Op[s' => 4, s => 4] = 1
    Op[s' => 6, s => 6] = 1
    return Op[s' => 8, s => 8] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"Nbdn", ::SiteType"2tJ", s::Index)
    Op[s' => 5, s => 5] = 1
    Op[s' => 7, s => 7] = 1
    return Op[s' => 9, s => 9] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"Nb", ::SiteType"2tJ", s::Index)
    Op[s' => 4, s => 4] = 1
    Op[s' => 5, s => 5] = 1
    Op[s' => 6, s => 6] = 1
    Op[s' => 7, s => 7] = 1
    Op[s' => 8, s => 8] = 1
    return Op[s' => 9, s => 9] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"Sza", ::SiteType"2tJ", s::Index)
    Op[s' => 2, s => 2] = 0.5
    Op[s' => 3, s => 3] = -0.5
    Op[s' => 6, s => 6] = 0.5
    Op[s' => 7, s => 7] = 0.5
    Op[s' => 8, s => 8] = -0.5
    return Op[s' => 9, s => 9] = -0.5
end

function ITensors.op!(Op::ITensor, ::OpName"Szb", ::SiteType"2tJ", s::Index)
    Op[s' => 4, s => 4] = 0.5
    Op[s' => 5, s => 5] = -0.5
    Op[s' => 6, s => 6] = 0.5
    Op[s' => 7, s => 7] = -0.5
    Op[s' => 8, s => 8] = 0.5
    return Op[s' => 9, s => 9] = -0.5
end

function ITensors.op!(Op::ITensor, ::OpName"Sz", ::SiteType"2tJ", s::Index)
    Op[s' => 2, s => 2] = 1
    Op[s' => 3, s => 3] = -1
    Op[s' => 4, s => 4] = 1
    Op[s' => 5, s => 5] = -1
    Op[s' => 6, s => 6] = 2
    return Op[s' => 9, s => 9] = 2
end

function ITensors.op!(Op::ITensor, ::OpName"Spa", ::SiteType"2tJ", s::Index)
    Op[s' => 2, s => 3] = 1
    Op[s' => 6, s => 8] = 1
    return Op[s' => 7, s => 9] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"Sma", ::SiteType"2tJ", s::Index)
    Op[s' => 3, s => 2] = 1
    Op[s' => 8, s => 6] = 1
    return Op[s' => 9, s => 7] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"Sxa", ::SiteType"2tJ", s::Index)
    Op[s' => 2, s => 3] = 0.5
    Op[s' => 3, s => 2] = 0.5
    Op[s' => 6, s => 8] = 0.5
    Op[s' => 8, s => 6] = 0.5
    Op[s' => 7, s => 9] = 0.5
    return Op[s' => 9, s => 7] = 0.5
end

function ITensors.op!(Op::ITensor, ::OpName"Sya", ::SiteType"2tJ", s::Index)
    complex!(Op)
    Op[s' => 2, s => 3] =  0.5*1im
    Op[s' => 3, s => 2] = -0.5*1im
    Op[s' => 6, s => 8] =  0.5*1im
    Op[s' => 8, s => 6] = -0.5*1im
    Op[s' => 7, s => 9] =  0.5*1im
    return Op[s' => 9, s => 7] = -0.5*1im
end

function ITensors.op!(Op::ITensor, ::OpName"Spb", ::SiteType"2tJ", s::Index)
    Op[s' => 4, s => 5] = 1
    Op[s' => 6, s => 7] = 1
    return Op[s' => 8, s => 9] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"Smb", ::SiteType"2tJ", s::Index)
    Op[s' => 5, s => 4] = 1
    Op[s' => 7, s => 6] = 1
    return Op[s' => 9, s => 8] = 1
end

function ITensors.op!(Op::ITensor, ::OpName"Sxb", ::siteType"2tJ", s::Index)
    Op[s' => 4, s => 5] = 0.5
    Op[s' => 5, s => 4] = 0.5
    Op[s' => 6, s => 7] = 0.5
    Op[s' => 7, s => 6] = 0.5
    Op[s' => 8, s => 9] = 0.5
    return Op[s => 9, s => 8] = 0.5
end

function ITensors.op!(Op::ITensor, ::OpName"Syb", ::siteType"2tJ", s::Index)
    complex!(Op)
    Op[s' => 4, s => 5] =  0.5*1im
    Op[s' => 5, s => 4] = -0.5*1im
    Op[s' => 6, s => 7] =  0.5*1im
    Op[s' => 7, s => 6] = -0.5*1im
    Op[s' => 8, s => 9] =  0.5*1im
    return Op[s' => 9, s => 8] = -0.5*1im
end

function ITensors.op!(Op::ITensor, ::OpName"bdagup", ::SiteType"2tJ", s::Index)
    return Op[s' => 4, s => 2] = -1
end

function ITensors.op!(Op::ITensor, ::OpName"bup", ::SiteType"2tJ", s::Index)
    return Op[s' => 2, s => 4] = -1
end

function ITensors.op!(Op::ITensor, ::OpName"bdagdn", ::SiteType"2tJ", s::Index)
    return Op[s' => 5, s => 3] = -1
end

function ITensors.op!(Op::ITensor, ::OpName"bdn", ::SiteType"2tJ", s::Index)
    return Op[s' => 3, s => 5] = -1
end

ITensors.has_fermion_string(::OpName"Cadagup", ::SiteType"2tJ") = true
ITensors.has_fermion_string(::OpName"Caup", ::SiteType"2tJ") = true
ITensors.has_fermion_string(::OpName"Cadagdn", ::SiteType"2tJ") = true
ITensors.has_fermion_string(::OpName"Cadn", ::SiteType"2tJ") = true
ITensors.has_fermion_string(::OpName"Cbdagup", ::SiteType"2tJ") = true
ITensors.has_fermion_string(::OpName"Cbup", ::SiteType"2tJ") = true
ITensors.has_fermion_string(::OpName"Cbdagdn", ::SiteType"2tJ") = true
ITensors.has_fermion_string(::OpName"Cbdn", ::SiteType"2tJ") = true

function Build_H(N, t, J, U)
    os = OpSum()
    for i in 1:N-1
        # Hopping
        os += -t, "Cadagup", i, "Caup", i+1
        os += -t, "Cadagup", i+1, "Caup", i
        os += -t, "Cadagdn", i, "Cadn", i+1
        os += -t, "Cadagdn", i+1, "Cadn", i

        os += -t, "Cbdagup", i, "Cbup", i+1
        os += -t, "Cbdagup", i+1, "Cbup", i
        os += -t, "Cbdagdn", i, "Cbdn", i+1
        os += -t, "Cbdagdn", i+1, "Cbdn", i

        # Exchange
        os += J, "Sza", i, "Sza", i+1
        os += J/2, "Spa", i, "Sma", i+1
        os += J/2, "Spa", i+1, "Sma", i

        os += -0.25*J, "Na", i, "Na", i+1
        os += -0.25*J, "Nb", i, "Nb", i+1

        os += J, "Szb", i, "Szb", i+1
        os += J/2, "Spb", i, "Smb", i+1
        os += J/2, "Spb", i+1, "Smb", i

    end

    for i in 1:N
        os += U, "Na", i, "Nb", i
    end
    return os
end

function initial_state(N, Na, Nb; Sza_target=nothing, Szb_target=nothing)
    @assert 0 ≤ Na ≤ N
    @assert 0 ≤ Nb ≤ N

    Sza_target = isnothing(Sza_target) ? (isodd(Na) ? 1 : 0) : Sza_target
    Szb_target = isnothing(Szb_target) ? (isodd(Nb) ? 1 : 0) : Szb_target

    @assert abs(Sza_target) ≤ Na && (Na - abs(Sza_target)) % 2 == 0
    @assert abs(Szb_target) ≤ Nb && (Nb - abs(Szb_target)) % 2 == 0

    Na_up = (Na + Sza_target) ÷ 2
    Na_dn = Na - Na_up
    Nb_up = (Nb + Szb_target) ÷ 2
    Nb_dn = Nb - Nb_up

    a_occ = distribute_electrons(N, Na_up, Na_dn)
    b_occ = distribute_electrons(N, Nb_up, Nb_dn)

    state = Vector{String}(undef, N)
    for i in 1:N
        a = a_occ[i]
        b = b_occ[i]
        state[i] = "A" * a * "B" * b
    end
    return state
end

function distribute_electrons(N, n_up, n_dn)
    n_tot = n_up + n_dn
    occ = fill("Emp", N)
    n_tot ==  0 && return occ
    sites = round.(Int, LinRange(1, N, Int(n_tot+2)))[2:end-1]
    sites = unique(sites)
    while length(sites) < n_tot
        for j in 1:N
            if !(j in sites)
                push!(sites, j)
                length(sites) == n_tot && break
            end
        end
    end
    sort!(sites)
    spins = Vector{String}(undef, n_tot)
    iu, id = 0, 0
    for k in 1:n_tot
        if isodd(k) && iu<n_up
            spins[k] = "up"; iu += 1
        elseif iseven(k) && id < n_dn
            spins[k] = "dn"; id += 1
        elseif iu < n_up
            spins[k] = "up"; iu += 1
        else 
            spins[k] = "dn"; id += 1
        end
    end

    for (s, σ) in zip(sites, spins)
        occ[s] = σ
    end
    return occ
end

function Find_E(N, Na, Nb, t, J, U)
    sites = siteinds("2tJ", N, conserve_qns=true)
    H = MPO(Build_H(N, t, J, U), sites)
    state = initial_state(N, Na, Nb)

    ψ = randomMPS(sites, state, 10)

    sweeps = Sweeps(20)
    maxdim!(sweeps, 20, 20, 30, 30, 50, 50, 100, 100, 100, 200)
    noise!(sweeps, 1E-6, 1E-7, 1E-8, 1E-9)
    cutoff!(sweeps, 1E-8)
    E, ϕ = dmrg(H, ψ, sweeps)
    return E
end

function Compute_single_exciton_distribution(N, Na, Nb, t, J, U)
    sites = siteinds("2tJ", N, conserve_qns=true)
    H = MPO(Build_H(N, t, J, U), sites)
    state = initial_state(N, Na, Nb)

    ψ = randomMPS(sites, state, 10)

    sweeps = Sweeps(20)
    maxdim!(sweeps, 20, 20, 30, 30, 50, 50, 100, 100, 100, 200)
    noise!(sweeps, 1E-6, 1E-7, 1E-8, 1E-9)
    cutoff!(sweeps, 1E-8)
    E, ϕ = dmrg(H, ψ, sweeps)
    C = correlation_matrix(ϕ, "Bdagup", "Bup")
    k = range(-π, stop=π, N+1)[1:end-1]
    Cex = zeros(N)
    for kval in 1:N
        for x in 1:N
            for y in 1:N
                Cex[kval] += 1/N * exp(1im*k[kval]*(x-y)) * C[x,y]
            end
        end
    end
    return Cex
end

function Find_Dex(N, Na, Nb, t, J, U)
    E1 = Find_E(N, Na, Nb, t, J, U)
    E2 = Find_E(N, Na-1, Nb, t, J, U)
    E3 = Find_E(N, Na, Nb+1, t, J, U)
    E4 = Find_E(N, Na-1, Nb+1, t, J, U)
    println(E1 - E2 - E3 + E4)
    return E1 - E2 - E3 + E4
end

# LLG functions

function Compute_Heff(Spins, Qspins, Jex, Jsd, Jani, eani, Hext)
    # Computes the effective field in the LLG equation
    # For simplicity, I'm only gonna consider two spins, one per layer, 
    # and coupled to the middle of the chain
    # Spins is a 2x3 array with the spin directions
    # Qspins is the vector of expectation valuues of quantum spins
    # Jsd coupling strength between quantum and classical spins
    # Jex is the Heisenberg exchange interaction between classical spins
    # Jani is the easy axis anisotropy
    # Hext is the external magnetic field
    # eani anisotropy unitary vector direction

    Heff = zeros(2, 3)

    # Classical interaction
    Heff[1,:] += Jex * Spins[2,:] + Hext[1,:] + 2*Jani * (eani[1,:] ⋅ Spins[1,:]) * eani[1,:]
    Heff[2,:] += Jex * Spins[1,:] + Hext[1,:] + 2*Jani * (eani[2,:] ⋅ Spins[2,:]) * eani[2,:]

    # Interaction with exciton via mean field coupling
    Heff[1,:] += Jsd * Qspins[1,:]
    Heff[2,:] += Jsd * Qspins[2,:]
    return Heff
end

function Spin_change(Spins, gamma, lambda, dt, Heff)
    # Computes the change for a single LLG discrete time-evolution
    # gamma: gyromagnetic ratio
    # lambda: Gilbert damping
    # dt: time step
    Spinsp = zeros(2,3)

    for i in 1:2
        A = -gamma*(cross(Spins[i,:], Heff[i,:])+lambda*(cross(Spins[i,:], cross(Spins[i,:],Heff[i,:]))))
        Spinsp[i,:] += A
    end
    return Spinsp
end

function Normalize(Spin)
    NSpins = zeros(2,3)
    for i in 1:2
        NSpins[i,:] = Spin[i,:]/norm(Spin[i,:])
    end
    return NSpins
end

function Integrate_LLG(Spins, Qspins, gamma, Lambda, dt, Jex, Jsd, Jani, eani, Hext)
    Heff = Compute_Heff(Spins, Qspins, Jex, Jsd, Jani, eani, Hext)
    Spinsp1 = Spin_change(Spins, gamma, Lambda, dt, Heff)
    New_spins = Normalize(Spins + Spinsp1)
    Heff = Compute_Heff(New_spins, Qspins, Jex, Jsd, Jani, eani, Heff)
    Spinsp2 = Spin_change(New_spins, gamma, Lambda, dt, Heff)
    Spin = Normalize(Spins + 0.5*(Spinsp1+Spinsp2))
    return Spin
end

function get_QSpins(ϕ, Nhalf)
    Qspins = zeros(2,3)
    Qspins[1,1] = real(expect(ϕ, "Sxa"; sites=Nhalf))
    Qspins[1,2] = real(expect(complex(ϕ), "Sya"; sites=Nhalf))
    Qspins[1,3] = real(expect(ϕ, "Sza"; sites=Nhalf))

    Qspins[2,1] = real(expect(ϕ, "Sxb"; sites=Nhalf))
    Qspins[2,2] = real(expect(complex(ϕ), "Syb"; sites=Nhalf))
    Qspins[2,3] = real(expect(ϕ, "Szb"; sites=Nhalf))
    return Qspins
end

N = 24
Na = 18
Nb = 6

t = 1
J = 1.2
U = 3

Dex = zeros(6)

Dex[1] = Find_Dex(N, Na, Nb, t, J, U)

N = 32
Na = 24
Nb = 8

Dex[2] = Find_Dex(N, Na, Nb, t, J, U)

N = 40
Na = 30
Nb = 10

Dex[3] = Find_Dex(N, Na, Nb, t, J, U)

N = 48
Na = 36
Nb = 12

Dex[4] = Find_Dex(N, Na, Nb, t, J, U)

N = 56
Na = 42
Nb = 14

Dex[5] = Find_Dex(N, Na, Nb, t, J, U)

N = 64
Na = 48
Nb = 16

Dex[6] = Find_Dex(N, Na, Nb, t, J, U)

Cex = zeros(64)
Cex = Compute_single_exciton_distribution(N, Na, Nb, t, J, U)

writedlm("Data/Cex.txt", data)
