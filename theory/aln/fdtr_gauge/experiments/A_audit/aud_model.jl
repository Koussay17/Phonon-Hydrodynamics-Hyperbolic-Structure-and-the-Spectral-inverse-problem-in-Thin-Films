# Independent FDTR forward model for the audit of workstream A (no code shared with src/twin_fdtr.py or
# experiments/A/indep.py; physics derived afresh).
#
# Stack: Au (Fourier) / R1 / AlN (anisotropic Cattaneo, tau_R d q/dt + q = -Lambda grad T) / R2 / sapphire (half-space).
# Hankel transform in r (T ~ J0(kr)), exp(+i omega t), p = i omega. AlN at wavenumber k:
#     Q_z' = -(C p + lam_r k^2 / g) T,  T' = -g Q_z / lam_z,  g = 1 + tau_R p
#     => gamma^2 = (C p g + lam_r k^2)/lam_z,  admittance Y = lam_z gamma / g.
# Layers combined as ABCD (thermal quadrupole) matrices [T; Q]_top = M [T; Q]_bottom, each matrix scaled by
# exp(-gamma d) (cancels in the ratio), substrate Q = Y_s T, Z = (M11 + M12 Ys)/(M21 + M22 Ys).
# Gaussian pump/probe (1/e^2 radii wp, wr): H = (1/2pi) int k exp(-eta k^2) Z dk, eta = (wp^2 + wr^2)/8,
# evaluated as H = 1/(4 pi eta) int_0^inf exp(-s) Z(k^2 = s/eta) ds by a vector-valued globally adaptive
# Gauss-Kronrod (7, 15) rule written here, with componentwise relative tolerance on H and on every derivative.
# Derivatives: forward-mode dual numbers (ForwardDiff.Dual), one pass for all 13 thermal/beam parameters.
# Fisher: whitened real Jacobian; marginal sigma of ln tau_R by orthogonal projection (Householder QR of the nuisance
# block, re-orthogonalised), cross-checked by the BigFloat inverse of the normal matrix.

using LinearAlgebra, ForwardDiff, Printf
BLAS.set_num_threads(1)

const HERE = @__DIR__
const OUT = joinpath(HERE, "out")
include(joinpath(OUT, "aud_inputs.jl"))          # defines INP (written by aud0_inputs.py)

# minimal JSON writer (the project environment has no JSON package as a direct dependency)
jsonstr(x::AbstractString) = "\"" * replace(x, "\\" => "\\\\", "\"" => "\\\"") * "\""
jsonstr(x::Bool) = x ? "true" : "false"
jsonstr(::Nothing) = "null"
jsonstr(x::Integer) = string(x)
jsonstr(x::Real) = isfinite(x) ? repr(Float64(x)) : "null"
jsonstr(x::Complex) = jsonstr([real(x), imag(x)])
jsonstr(x::Symbol) = jsonstr(string(x))
jsonstr(x::Union{AbstractVector,Tuple}) = "[" * join((jsonstr(v) for v in x), ",") * "]"
jsonstr(x::AbstractMatrix) = jsonstr([x[i, :] for i in 1:size(x, 1)])
jsonstr(x::AbstractDict) = "{" * join((jsonstr(string(k)) * ":" * jsonstr(v) for (k, v) in x), ",") * "}"
jsonstr(x::NamedTuple) = jsonstr(Dict(string(k) => v for (k, v) in pairs(x)))
writejson(path, x) = open(io -> write(io, jsonstr(x)), path, "w")

# ------------------------------------------------------------------ complex helpers valid for dual numbers
@inline function csqrt(z::Complex{T}) where {T<:Real}
    a, b = reim(z)
    r = sqrt(a * a + b * b)
    if a >= zero(T)
        u = sqrt((r + a) / 2)
        return Complex(u, b / (2u))
    else
        v = sqrt((r - a) / 2)
        v = b >= zero(T) ? v : -v
        return Complex(b / (2v), v)
    end
end
@inline cexp(z::Complex) = exp(real(z)) * Complex(cos(imag(z)), sin(imag(z)))

@inline function layer(γ, Y, d)
    E = cexp(-2 * γ * d)
    c = (1 + E) / 2
    s = (1 - E) / 2
    return (c, s / Y, Y * s, c)
end
@inline mmul(A, B) = (A[1] * B[1] + A[2] * B[3], A[1] * B[2] + A[2] * B[4],
                      A[3] * B[1] + A[4] * B[3], A[3] * B[2] + A[4] * B[4])
@inline resist(R) = (one(R), R, zero(R), one(R))

# ------------------------------------------------------------------ parameters
const PN = ["ln_tauR", "ln_lamz", "ln_lamr", "ln_C", "ln_R1", "ln_R2", "ln_lamAu", "ln_CAu", "ln_dAu",
            "ln_dAlN", "ln_lamS", "ln_CS", "ln_wscale"]
const NP = length(PN)
const ALLN = vcat(PN, ["g", "phi0", "td_ns"])

function base_for(Lkey::String; over = Dict{String,Any}())
    s = INP["stack"]; b = INP["base"][Lkey]
    L = parse(Float64, split(Lkey, "_")[1]) * 1e-9
    d = Dict{String,Float64}("tauR" => b["tauR"], "lamz" => b["lamz"], "lamr" => b["lamr"], "C" => s["C"],
        "R1" => s["R1"], "R2" => s["R2"], "lamAu" => s["lamAu"], "CAu" => s["CAu"], "dAu" => s["dAu"],
        "dAlN" => L, "lamS" => s["lamS"], "CS" => s["CS"], "wp" => s["wp"], "wr" => s["wr"])
    for (k, v) in over
        d[k] = v
    end
    return (; (Symbol(k) => v for (k, v) in d)...)
end

function phys(θ, b)
    e(i) = exp(θ[i])
    return (tauR = b.tauR * e(1), lamz = b.lamz * e(2), lamr = b.lamr * e(3), C = b.C * e(4),
            R1 = b.R1 * e(5), R2 = b.R2 * e(6), lamAu = b.lamAu * e(7), CAu = b.CAu * e(8), dAu = b.dAu * e(9),
            dAlN = b.dAlN * e(10), lamS = b.lamS * e(11), CS = b.CS * e(12), wp = b.wp * e(13), wr = b.wr * e(13))
end

function impedance(k2, p, P)
    γs = csqrt((P.CS * p + P.lamS * k2) / P.lamS)
    Ys = P.lamS * γs
    g = 1 + P.tauR * p
    γn = csqrt((P.C * p * g + P.lamr * k2) / P.lamz)
    Yn = P.lamz * γn / g
    γa = csqrt((P.CAu * p + P.lamAu * k2) / P.lamAu)
    Ya = P.lamAu * γa
    M = mmul(mmul(mmul(layer(γa, Ya, P.dAu), resist(P.R1)), layer(γn, Yn, P.dAlN)), resist(P.R2))
    return (M[1] + M[2] * Ys) / (M[3] + M[4] * Ys)
end

struct Tg end
dualθ(θ0) = [ForwardDiff.Dual{Tg}(θ0[i], ntuple(j -> j == i ? 1.0 : 0.0, NP)) for i in 1:NP]

# vector integrand [Re v, Im v, Re dv/dθ1, Im dv/dθ1, ...] with v = exp(-s) Z(s/eta) / (4 pi eta)
function integrand_vec(s, θd, b, f)
    P = phys(θd, b)
    η = (P.wp^2 + P.wr^2) / 8
    p = Complex(0.0, 2π * f)
    v = exp(-s) * impedance(s / η, p, P) / (4π * η)
    re, im = reim(v)
    out = Vector{Float64}(undef, 2 + 2NP)
    out[1] = ForwardDiff.value(re); out[2] = ForwardDiff.value(im)
    pr = ForwardDiff.partials(re); pim = ForwardDiff.partials(im)
    @inbounds for j in 1:NP
        out[1 + 2j] = pr[j]; out[2 + 2j] = pim[j]
    end
    return out
end

# value only (Float64 parameters), for perturbed evaluations
function integrand_val(s, θ, b, f)
    P = phys(θ, b)
    η = (P.wp^2 + P.wr^2) / 8
    p = Complex(0.0, 2π * f)
    v = exp(-s) * impedance(s / η, p, P) / (4π * η)
    return [real(v), imag(v)]
end

# ------------------------------------------------------------------ Gauss-Kronrod (7,15), QUADPACK nodes
const XGK = [0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
             0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
             0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
             0.207784955007898467600689403773245, 0.0]
const WGK = [0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
             0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
             0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
             0.204432940075298892414161999234649, 0.209482141084727828012999174891714]
const WG = [0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
            0.381830050505118944950369775488975, 0.417959183673469387755102040816327]

function gk15(f, a, b)
    c = (a + b) / 2; h = (b - a) / 2
    fc = f(c)
    rk = WGK[8] .* fc
    rg = WG[4] .* fc
    ra = WGK[8] .* abs.(fc)
    for j in 1:7
        x = h * XGK[j]
        f1 = f(c - x); f2 = f(c + x)
        s = f1 .+ f2
        rk .+= WGK[j] .* s
        ra .+= WGK[j] .* (abs.(f1) .+ abs.(f2))
        if iseven(j)
            rg .+= WG[j ÷ 2] .* s
        end
    end
    return rk .* h, abs.(rk .- rg) .* h, ra .* h
end

const BREAKS = [0.0, 1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 60.0, 120.0]

# Componentwise criterion (each complex component c = (re, im) pair):
#   err_c <= rtol * max(|I_c|, afloor * int |f_c|)
# The floor only matters for components whose integral nearly cancels (substrate parameters at 200 MHz); the
# error estimate is |K15 - G7| (pessimistic for analytic integrands).
function adapt(f, breaks; rtol = 1e-13, afloor = 1e-6, maxint = 3000)
    ints = Vector{Any}()
    It = nothing; Et = nothing; At = nothing
    for i in 1:length(breaks)-1
        I, E, A = gk15(f, breaks[i], breaks[i+1])
        push!(ints, (breaks[i], breaks[i+1], I, E, A))
        It = It === nothing ? copy(I) : It .+ I
        Et = Et === nothing ? copy(E) : Et .+ E
        At = At === nothing ? copy(A) : At .+ A
    end
    tolf(It, At) = begin
        t = similar(It)
        for c in 1:2:length(It)
            m = max(hypot(It[c], It[c+1]), afloor * hypot(At[c], At[c+1]))
            t[c] = rtol * m; t[c+1] = rtol * m
        end
        max.(t, 1e-300)
    end
    tol = tolf(It, At)
    while true
        if all(Et .<= tol) || length(ints) >= maxint
            return It, Et, length(ints), all(Et .<= tol), tol
        end
        sc = [maximum(x[4] ./ tol) for x in ints]
        i = argmax(sc)
        a, b, Io, Eo, Ao = ints[i]
        m = (a + b) / 2
        I1, E1, A1 = gk15(f, a, m); I2, E2, A2 = gk15(f, m, b)
        ints[i] = (a, m, I1, E1, A1); push!(ints, (m, b, I2, E2, A2))
        It .+= I1 .+ I2 .- Io; Et .+= E1 .+ E2 .- Eo; At .+= A1 .+ A2 .- Ao
        tol = tolf(It, At)
    end
end

# ln H and its Jacobian (complex, nf x NP) for the Gaussian average; geometry "1d": k = 0 limit, same normalisation
function lnH_jac(b, freqs; θ0 = zeros(NP), rtol = 1e-13, geometry = "beam")
    nf = length(freqs)
    lnH = Vector{ComplexF64}(undef, nf)
    J = Matrix{ComplexF64}(undef, nf, NP)
    stats = Dict("max_intervals" => 0, "all_converged" => true, "max_rel_err_est" => 0.0)
    θd = dualθ(θ0)
    for (i, f) in enumerate(freqs)
        if geometry == "1d"
            P = phys(θd, b)
            v = impedance(0.0, Complex(0.0, 2π * f), P) * 2 / (π * (P.wp^2 + P.wr^2))
            re, im = reim(v)
            I = vcat([ForwardDiff.value(re), ForwardDiff.value(im)],
                     vec(permutedims(hcat(collect(ForwardDiff.partials(re)), collect(ForwardDiff.partials(im))))))
        else
            I, E, n, ok, _ = adapt(s -> integrand_vec(s, θd, b, f), BREAKS; rtol = rtol)
            stats["max_intervals"] = max(stats["max_intervals"], n)
            stats["all_converged"] &= ok
            for c in (1, 3)    # H and the tau_R derivative: relative error estimates
                m = hypot(I[c], I[c+1])
                key = c == 1 ? "max_rel_err_est" : "max_rel_err_est_tau"
                stats[key] = max(get(stats, key, 0.0), hypot(E[c], E[c+1]) / m)
            end
        end
        H = complex(I[1], I[2])
        lnH[i] = log(H)
        for j in 1:NP
            J[i, j] = complex(I[1+2j], I[2+2j]) / H
        end
    end
    return lnH, J, stats
end

function lnH_val(b, freqs, θ; rtol = 1e-13, geometry = "beam")
    out = Vector{ComplexF64}(undef, length(freqs))
    for (i, f) in enumerate(freqs)
        if geometry == "1d"
            P = phys(θ, b)
            H = impedance(0.0, Complex(0.0, 2π * f), P) * 2 / (π * (P.wp^2 + P.wr^2))
        else
            I, _, _, _, _ = adapt(s -> integrand_val(s, θ, b, f), BREAKS; rtol = rtol)
            H = complex(I[1], I[2])
        end
        out[i] = log(H)
    end
    return out
end

# full complex design matrix: 13 thermal/beam columns + g, phi0 (rad), td (ns)
function full_design(J, freqs)
    nf = length(freqs)
    return hcat(J, ones(ComplexF64, nf), fill(1.0im, nf), -1.0im * 2π .* freqs .* 1e-9)
end

whiten(D, sa, sp_deg) = vcat(real.(D) ./ sa, imag.(D) ./ deg2rad(sp_deg))

const SETS = Dict(
    "N0" => ["ln_tauR"],
    "N1" => ["ln_tauR", "g", "phi0", "td_ns"],
    "N2" => ["ln_tauR", "g", "phi0", "td_ns", "ln_R1", "ln_R2"],
    "N3" => ["ln_tauR", "g", "phi0", "td_ns", "ln_R1", "ln_R2", "ln_lamz", "ln_C"],
    "N4" => ["ln_tauR", "g", "phi0", "td_ns", "ln_R1", "ln_R2", "ln_lamz", "ln_C", "ln_lamAu", "ln_CAu"],
    "N5" => ["ln_tauR", "g", "phi0", "td_ns", "ln_R1", "ln_R2", "ln_lamz", "ln_C", "ln_lamAu", "ln_CAu", "ln_wscale"],
)
const SETNAMES = ["N0", "N1", "N2", "N3", "N4", "N5"]

function sigma_proj(A, t::Int, nuis::Vector{Int})
    at = A[:, t]
    isempty(nuis) && return 1 / norm(at)
    B = A[:, nuis]
    F = qr(B)
    Q = F.Q * Matrix{Float64}(I, size(B, 1), size(B, 2))
    r = at - Q * (Q' * at)
    r = r - Q * (Q' * r)
    return 1 / norm(r)
end

function sigma_big(A, t::Int, cols::Vector{Int})
    setprecision(BigFloat, 256) do
        Ab = BigFloat.(A[:, cols])
        nrm = [sqrt(sum(Ab[:, j] .^ 2)) for j in 1:size(Ab, 2)]
        An = Ab ./ permutedims(nrm)
        Fm = An' * An
        Ci = inv(Fm)
        k = findfirst(==(t), cols)
        return Float64(sqrt(Ci[k, k]) / nrm[k])
    end
end

function cond_norm(A, cols)
    B = A[:, cols]
    B = B ./ permutedims([norm(B[:, j]) for j in 1:size(B, 2)])
    sv = svdvals(B)
    return sv[1] / sv[end]
end

function fisher_table(J, freqs; sa = 0.01, sigs = (0.01, 0.03, 0.1), mask = trues(length(freqs)), big = false)
    D = full_design(J, freqs)[mask, :]
    res = Dict{String,Any}()
    for N in SETNAMES
        cols = [findfirst(==(n), ALLN) for n in SETS[N]]
        res[N] = Dict{String,Any}()
        for sp in sigs
            A = whiten(D, sa, sp)
            s = sigma_proj(A, cols[1], cols[2:end])
            e = Dict{String,Any}("sigma" => s, "cond" => cond_norm(A, cols))
            if big
                e["sigma_bigfloat"] = sigma_big(A, cols[1], cols)
            end
            res[N][string(sp)] = e
        end
    end
    return res
end

# ------------------------------------------------------------------ main computations
function main()
    freqs = Float64.(INP["freq80"])
    t0 = time()
    out = Dict{String,Any}("freq80" => freqs)

    # 0. rule check: GK15 integrates polynomials of degree <= 22 exactly; G7 degree <= 13
    pe = Float64[]
    for n in (0, 5, 13, 22)
        I, E = gk15(x -> [x^n], 0.0, 1.0)
        push!(pe, abs(I[1] - 1 / (n + 1)))
    end
    out["gk15_polynomial_errors_deg_0_5_13_22"] = pe
    @printf("GK15 polynomial errors %s\n", string(pe))

    # 1. Jacobians at all four thicknesses, finite beams, tolerance study
    cases = Dict{String,Any}()
    for Lk in ("30", "100", "300", "500", "500_rta")
        b = base_for(Lk)
        lnH13, J13, st13 = lnH_jac(b, freqs; rtol = 1e-13)
        lnH10, J10, st10 = lnH_jac(b, freqs; rtol = 1e-10)
        lnH1d, J1d, _ = lnH_jac(b, freqs; geometry = "1d")
        tab = fisher_table(J13, freqs; big = true)
        tab10 = fisher_table(J10, freqs)
        tab1d = fisher_table(J1d, freqs)
        # bands as subsets of the 80 points
        m1 = (freqs .>= 1e6 * (1 - 1e-12)) .& (freqs .<= 2e8 * (1 + 1e-12))
        m2 = (freqs .>= 1e4 * (1 - 1e-12)) .& (freqs .<= 5e7 * (1 + 1e-12))
        tabb1 = fisher_table(J13, freqs; mask = m1)
        tabb2 = fisher_table(J13, freqs; mask = m2)
        # quadrature tolerance effect on sigma
        qrel = maximum(abs(tab10[N][sp]["sigma"] / tab[N][sp]["sigma"] - 1) for N in SETNAMES for sp in ("0.01", "0.03", "0.1"))
        bigrel = maximum(abs(tab[N][sp]["sigma_bigfloat"] / tab[N][sp]["sigma"] - 1) for N in SETNAMES for sp in ("0.01", "0.03", "0.1"))
        # phase sensitivity and exact phase change at the admissible w
        dphi = rad2deg.(imag.(J13[:, 1]))
        imax = argmax(abs.(dphi))
        exact = Dict{String,Any}()
        for w in (-0.08, -0.052, -0.037, -0.028, 0.028, 0.037, 0.052, 0.08)
            θ = zeros(NP); θ[1] = log1p(w)
            lw = lnH_val(b, freqs, θ)
            d = rad2deg.(imag.(lw .- lnH13))
            k = argmax(abs.(d))
            exact[@sprintf("%+.3f", w)] = Dict("max_abs_dphi_deg" => abs(d[k]), "f_Hz" => freqs[k],
                "max_abs_dlnA" => maximum(abs.(real.(lw .- lnH13))),
                "linear_pred_deg" => abs(w) * maximum(abs.(dphi)))
        end
        # N0 information per frequency (share of the total at 0.01 deg)
        A0 = whiten(full_design(J13, freqs)[:, 1:1], 0.01, 0.01)
        nf = length(freqs)
        info = A0[1:nf, 1] .^ 2 .+ A0[nf+1:end, 1] .^ 2
        cases[Lk] = Dict("tauR_ps" => b.tauR * 1e12, "lamz" => b.lamz, "lamr" => b.lamr,
            "lnH" => [[real(z), imag(z)] for z in lnH13],
            "J_tau" => [[real(z), imag(z)] for z in J13[:, 1]],
            "J_all" => [[[real(z), imag(z)] for z in J13[:, j]] for j in 1:NP],
            "quad_stats_1e-13" => st13, "quad_stats_1e-10" => st10,
            "lnH_diff_rtol_1e-10_vs_1e-13" => maximum(abs.(lnH10 .- lnH13)),
            "Jtau_rel_diff_rtol_1e-10_vs_1e-13" => maximum(abs.(J10[:, 1] .- J13[:, 1])) / maximum(abs.(J13[:, 1])),
            "sigma_rel_change_rtol_1e-10" => qrel, "sigma_bigfloat_rel_diff" => bigrel,
            "table" => tab, "table_1M_200M" => tabb1, "table_10k_50M" => tabb2, "table_1d" => tab1d,
            "max_abs_dphi_dlntau_deg" => abs(dphi[imax]), "f_at_max_Hz" => freqs[imax],
            "max_abs_dlnA_dlntau" => maximum(abs.(real.(J13[:, 1]))),
            "exact_dphi" => exact,
            "N0_info_fraction_above_10MHz" => sum(info[freqs .> 1e7]) / sum(info),
            "N0_info_fraction_phase" => sum(A0[nf+1:end, 1] .^ 2) / sum(info))
        @printf("L=%s tau=%.4f ps  N0/N1/N2/N3 0.01deg: %.4f %.4f %.4f %.4f | N2 0.1deg %.4f | quad rtol effect %.1e | big %.1e | %.0fs\n",
                Lk, b.tauR * 1e12, tab["N0"]["0.01"]["sigma"], tab["N1"]["0.01"]["sigma"], tab["N2"]["0.01"]["sigma"],
                tab["N3"]["0.01"]["sigma"], tab["N2"]["0.1"]["sigma"], qrel, bigrel, time() - t0)
        flush(stdout)
    end
    out["cases"] = cases
    writejson(joinpath(OUT, "aud_julia_core.json"), out)
    println("done ", time() - t0)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
