using Printf
γ_c=0.40; γ_v=0.40; Δ_0=2.0; Jpd_v=0.069; Jpd_c=0.031; V_intra=1.708; V_inter=0.691

E0_pair_analytic() = Δ_0 + Jpd_c + Jpd_v
J_plus_analytic(θ) = (Jpd_c+Jpd_v)*cos(θ/2)
J_minus_analytic(θ) = (Jpd_c-Jpd_v)*cos(θ/2)
pair_dispersion_scale() = 2*(γ_c+γ_v)

function continuum_g_below(Ω, Eedge; γc=γ_c, γv=γ_v)
    A = 2*(γc+γv)
    Ω < Eedge || return NaN
    radicand = (Eedge+A-Ω)^2 - A^2
    radicand > 0 || return NaN
    return -1/sqrt(radicand)
end

function sector_pole_equation(Ω, Jsector; V_intra=V_intra, V_inter=V_inter)
    E0=E0_pair_analytic(); Jabs=abs(Jsector); Vavg=0.5*(V_intra+V_inter)
    gminus=continuum_g_below(Ω, E0-Jabs); gplus=continuum_g_below(Ω, E0+Jabs)
    (isfinite(gminus) && isfinite(gplus)) || return NaN
    return 1 + Vavg*(gminus+gplus) + V_intra*V_inter*gminus*gplus
end

function bisect_root(f,a,b;iterations=80)
    fa=f(a); fb=f(b)
    @assert isfinite(fa) && isfinite(fb) && fa*fb<=0
    lo,hi=a,b; flo=fa
    for _ in 1:iterations
        mid=0.5*(lo+hi); fmid=f(mid)
        if !isfinite(fmid); hi=mid
        elseif flo*fmid<=0; hi=mid
        else; lo=mid; flo=fmid; end
    end
    return 0.5*(lo+hi)
end

function sector_bound_roots(Jsector; V_intra=V_intra, V_inter=V_inter, nscan=2500, margin=1e-7)
    E0=E0_pair_analytic(); Jabs=abs(Jsector); lower_edge=E0-Jabs
    A=pair_dispersion_scale()
    width = 4A + 4max(abs(V_intra),abs(V_inter)) + 2Jabs + 2
    xs = collect(range(lower_edge-width, stop=lower_edge-margin, length=nscan))
    f(Ω)=sector_pole_equation(Ω,Jsector;V_intra=V_intra,V_inter=V_inter)
    roots=Float64[]
    prev=f(xs[1])
    for i in 2:length(xs)
        cur=f(xs[i])
        if isfinite(prev) && isfinite(cur) && prev*cur<0
            push!(roots, bisect_root(f,xs[i-1],xs[i]))
        end
        prev=cur
    end
    return roots
end

function bright_dark_poles(θ)
    s=sector_bound_roots(J_plus_analytic(θ)); a=sector_bound_roots(J_minus_analytic(θ))
    return (isempty(s) ? NaN : minimum(s)), (isempty(a) ? NaN : minimum(a))
end

println("θ/π    bright(eV)   dark(eV)   splitting(meV)")
for θ in range(0.1π, 0.9π, length=9)
    Eb, Ed = bright_dark_poles(θ)
    @printf("%.3f   %.4f      %.4f     %.2f\n", θ/π, Eb, Ed, 1e3*(Ed-Eb))
end

println("\nzoom near AFM (θ→π):")
for θ in range(0.90π, 0.999π, length=15)
    Eb, Ed = bright_dark_poles(θ)
    @printf("%.4f   %.5f      %.5f     %.4f\n", θ/π, Eb, Ed, 1e3*(Ed-Eb))
end
