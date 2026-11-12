# Audit: 1D Guyer-Krumhansl rows of the branch (k = 0, local temperature convention) with the independent ABCD model.
# At k = 0 a GK film is the Cattaneo film with lam/g replaced by lam_eff = lam (1 + tau_l p)/(1 + tau_R p),
# gamma = sqrt(C p / lam_eff), Y = lam_eff gamma. Parameters as in aud_model.jl plus ln tau_l (14th).
include("aud_model.jl")

function impedance_gk1d(p, P, taul)
    γs = csqrt(P.CS * p / P.lamS); Ys = P.lamS * γs
    leff = P.lamz * (1 + taul * p) / (1 + P.tauR * p)
    γn = csqrt(P.C * p / leff); Yn = leff * γn
    γa = csqrt(P.CAu * p / P.lamAu); Ya = P.lamAu * γa
    M = mmul(mmul(mmul(layer(γa, Ya, P.dAu), resist(P.R1)), layer(γn, Yn, P.dAlN)), resist(P.R2))
    return (M[1] + M[2] * Ys) / (M[3] + M[4] * Ys)
end

struct Tg2 end
function jac_gk1d(b, taul0, freqs)
    np_ = NP + 1
    θd = [ForwardDiff.Dual{Tg2}(0.0, ntuple(j -> j == i ? 1.0 : 0.0, np_)) for i in 1:np_]
    J = Matrix{ComplexF64}(undef, length(freqs), np_)
    for (i, f) in enumerate(freqs)
        P = phys(θd[1:NP], b)
        v = log(impedance_gk1d(Complex(0.0, 2π * f), P, taul0 * exp(θd[np_])) * 2 / (π * (P.wp^2 + P.wr^2)))
        # log of a dual complex: log|v| + i arg v, done through real parts
        re, im = reim(v)
        for j in 1:np_
            J[i, j] = complex(ForwardDiff.partials(re)[j], ForwardDiff.partials(im)[j])
        end
    end
    return J
end
Base.log(z::Complex{<:ForwardDiff.Dual}) = Complex(log(sqrt(real(z)^2 + imag(z)^2)), atan(imag(z), real(z)))

const ALLN_GK = vcat(PN, ["ln_taul", "g", "phi0", "td_ns"])
const SETS_GK = Dict(
    "N0" => ["ln_tauR"],
    "N1" => ["ln_tauR", "g", "phi0", "td_ns"],
    "N2" => ["ln_tauR", "g", "phi0", "td_ns", "ln_R1", "ln_R2"],
    "N3" => ["ln_tauR", "g", "phi0", "td_ns", "ln_R1", "ln_R2", "ln_lamz", "ln_C"],
    "N4" => ["ln_tauR", "g", "phi0", "td_ns", "ln_R1", "ln_R2", "ln_lamz", "ln_C", "ln_taul", "ln_lamAu", "ln_CAu"],
)

function main_gk()
    freqs = Float64.(INP["freq80"])
    out = Dict{String,Any}()
    for Lk in ("30", "100", "300", "500"), tf in (0.5, 2.0)
        b = base_for(Lk; over = Dict{String,Any}("lamr" => INP["base"][Lk]["lamz"]))
        J = jac_gk1d(b, tf * b.tauR, freqs)
        nf = length(freqs)
        D = hcat(J, ones(ComplexF64, nf), fill(1.0im, nf), -1.0im * 2π .* freqs .* 1e-9)
        row = Dict{String,Any}()
        for N in ("N0", "N1", "N2", "N3", "N4"), sp in (0.01, 0.1)
            cols = [findfirst(==(n), ALLN_GK) for n in SETS_GK[N]]
            A = whiten(D[:, cols], 0.01, sp)
            row["$N|$sp"] = sigma_proj(A, 1, collect(2:length(cols)))
        end
        out["$Lk|gk$tf"] = row
        @printf("L=%s tau_l=%.1f tau_R: N0 %.4g N2 %.4g N3 %.4g N4 %.4g (0.01 deg)\n", Lk, tf, row["N0|0.01"], row["N2|0.01"],
                row["N3|0.01"], row["N4|0.01"])
    end
    writejson(joinpath(OUT, "aud_gk1d.json"), out)
end
main_gk()
