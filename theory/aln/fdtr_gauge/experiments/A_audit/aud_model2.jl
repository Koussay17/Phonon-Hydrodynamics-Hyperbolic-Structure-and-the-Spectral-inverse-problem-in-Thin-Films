# Audit, Julia part 2 (uses the independent model of aud_model.jl):
#  1. noiseless profile likelihood (nuisances refitted by Gauss-Newton with dual-number Jacobians) against the
#     Fisher prediction (w/sigma)^2: linear response in tau_R and the reading w* = sigma(ln tau_R);
#  2. linearity of H in tau_R over 0 <= tau_R <= 4 tau_R,0;
#  3. frequency-grid study sigma*sqrt(N), N = 20 ... 1280, endpoint and cell-centre grids;
#  4. unit test of the phase weighting (degrees/radians), phase/amplitude split;
#  5. scenario sensitivities that could lower w* (not part of the audited claim).
include("aud_model.jl")

function design_cols(J, freqs, names)
    D = full_design(J, freqs)
    return D[:, [findfirst(==(n), ALLN) for n in names]]
end

# ---------------------------------------------------------------- 1. profile likelihood
function profile(Lk, setname, sp_deg, ws; sa = 0.01)
    freqs = Float64.(INP["freq80"])
    nf = length(freqs)
    b0 = base_for(Lk)
    d0 = lnH_val(b0, freqs, zeros(NP); rtol = 1e-12)
    names = SETS[setname][2:end]                         # nuisances
    thn = [n for n in names if n in PN]
    ins = [n for n in names if !(n in PN)]
    Dins = full_design(zeros(ComplexF64, nf, NP), freqs)[:, [NP + findfirst(==(n), ["g", "phi0", "td_ns"]) for n in ins]]
    instr(ν) = isempty(ins) ? zeros(ComplexF64, nf) : Dins * ν
    out = Any[]
    for w in ws
        bw = base_for(Lk; over = Dict{String,Any}("tauR" => b0.tauR * (1 + w)))
        θ = zeros(NP)
        ν = zeros(length(ins))
        chi2_prev = Inf
        iters = 0
        if !isempty(names)
            for it in 1:8
                iters = it
                lm, Jm, _ = lnH_jac(bw, freqs; θ0 = θ, rtol = 1e-10)
                R = whiten(log.(exp.(d0 .- lm .- instr(ν))), sa, sp_deg)     # principal branch of ln(H_d/H_m)
                chi2 = sum(abs2, R)
                A = whiten(hcat(Jm[:, [findfirst(==(n), PN) for n in thn]], Dins), sa, sp_deg)
                Δ = A \ R
                for (k, n) in enumerate(thn)
                    θ[findfirst(==(n), PN)] += Δ[k]
                end
                ν .+= Δ[length(thn)+1:end]
                (maximum(abs.(Δ)) < 1e-9 || abs(chi2 - chi2_prev) <= 1e-10 * max(chi2, 1e-12)) && break
                chi2_prev = chi2
            end
        end
        # final misfit with a tighter quadrature
        lm = lnH_val(bw, freqs, θ; rtol = 1e-12)
        chi2 = sum(abs2, whiten(log.(exp.(d0 .- lm .- instr(ν))), sa, sp_deg))
        push!(out, Dict("w" => w, "dchi2" => chi2, "gn_iterations" => iters,
                        "nuisance_thermal" => Dict(n => θ[findfirst(==(n), PN)] for n in thn),
                        "nuisance_instr" => Dict(ins[k] => ν[k] for k in 1:length(ins))))
    end
    _, J0, _ = lnH_jac(b0, freqs)
    A = whiten(design_cols(J0, freqs, SETS[setname]), sa, sp_deg)
    σF = sigma_proj(A, 1, collect(2:size(A, 2)))
    for e in out
        e["fisher_pred_w"] = (e["w"] / σF)^2
        e["fisher_pred_lnw"] = e["w"] > -1 ? (log1p(e["w"]) / σF)^2 : nothing
    end
    return Dict("sigma_F" => σF, "profile" => out)
end

function main2()
    t0 = time()
    freqs = Float64.(INP["freq80"])
    res = Dict{String,Any}()

    # 1. profile likelihood
    prof = Dict{String,Any}()
    ws = [-1.0, -0.5, -0.08, -0.028, 0.028, 0.08, 1.0, 3.0]
    for (Lk, N, sp) in (("500", "N0", 0.01), ("500", "N1", 0.01), ("500", "N2", 0.01), ("500", "N2", 0.1),
                        ("500", "N3", 0.01), ("300", "N0", 0.01), ("300", "N2", 0.01), ("100", "N2", 0.1))
        p = profile(Lk, N, sp, ws)
        prof["$Lk|$N|$sp"] = p
        @printf("profile %s %s %.2f sigma_F %.4f : %s  (%.0fs)\n", Lk, N, sp, p["sigma_F"],
                join([@sprintf("w=%+.3f dchi2=%.6g pred=%.6g", e["w"], e["dchi2"], e["fisher_pred_w"]) for e in p["profile"]], "; "),
                time() - t0)
        flush(stdout)
    end
    res["profile"] = prof

    # 2. linearity of H in tau_R
    lin = Dict{String,Any}()
    for Lk in ("30", "100", "300", "500")
        b = base_for(Lk)
        l0, J0, _ = lnH_jac(b, freqs; rtol = 1e-12)
        H0 = exp.(l0)
        dHdτ = H0 .* J0[:, 1] ./ b.tauR            # dH/dtau_R
        e = Dict{String,Any}()
        for w in (-1.0, -0.5, 1.0, 3.0)
            bw = base_for(Lk; over = Dict{String,Any}("tauR" => b.tauR * (1 + w)))
            Hw = exp.(lnH_val(bw, freqs, zeros(NP); rtol = 1e-12))
            ΔH = Hw .- H0
            linp = w * b.tauR .* dHdτ
            e[@sprintf("%+.1f", w)] = maximum(abs.(ΔH .- linp)) / maximum(abs.(linp))
        end
        lin[Lk] = e
    end
    res["linearity_rel_dev_of_dH"] = lin
    println("linearity ", lin); flush(stdout)

    # 3. frequency grid
    grid = Dict{String,Any}()
    for Lk in ("100", "500")
        b = base_for(Lk)
        rows = Dict{String,Any}()
        for N in (20, 40, 80, 160, 320, 640, 1280)
            for kind in ("endpoints", "cellcentres")
                if kind == "endpoints"
                    f = exp.(range(log(1e4), log(2e8), length = N))
                else
                    e = range(log(1e4), log(2e8), length = N + 1)
                    f = exp.((e[1:end-1] .+ e[2:end]) ./ 2)
                end
                _, J, _ = lnH_jac(b, collect(f); rtol = 1e-11)
                tab = fisher_table(J, collect(f); sigs = (0.01, 0.1))
                rows["$kind|$N"] = Dict(Nm => Dict(sp => tab[Nm][sp]["sigma"] * sqrt(N) for sp in ("0.01", "0.1")) for Nm in ("N0", "N1", "N2", "N3"))
            end
        end
        grid[Lk] = rows
        println("grid ", Lk, " N2|0.1: ", [round(rows["endpoints|$N"]["N2"]["0.1"], digits = 2) for N in (20, 40, 80, 160, 320, 640, 1280)],
                " centres: ", [round(rows["cellcentres|$N"]["N2"]["0.1"], digits = 2) for N in (20, 40, 80, 160, 320, 640, 1280)])
        flush(stdout)
    end
    res["frequency_grid_sigma_sqrtN"] = grid

    # 4. units and information split
    units = Dict{String,Any}()
    for Lk in ("300", "500")
        b = base_for(Lk)
        _, J, _ = lnH_jac(b, freqs; rtol = 1e-12)
        jt = J[:, 1]
        σφ = deg2rad(0.01)
        correct = 1 / sqrt(sum(abs2, real.(jt) ./ 0.01) + sum(abs2, imag.(jt) ./ σφ))
        phase_deg_noise_rad = 1 / sqrt(sum(abs2, real.(jt) ./ 0.01) + sum(abs2, imag.(jt) ./ 0.01))           # sigma_phi 0.01 taken as rad
        J_in_deg = 1 / sqrt(sum(abs2, real.(jt) ./ 0.01) + sum(abs2, rad2deg.(imag.(jt)) ./ σφ))             # phase derivative in deg, noise in rad
        phase_only = 1 / sqrt(sum(abs2, imag.(jt) ./ σφ))
        amp_only = 1 / sqrt(sum(abs2, real.(jt) ./ 0.01))
        units[Lk] = Dict("N0_0.01deg_correct" => correct, "if_sigma_phi_read_as_rad" => phase_deg_noise_rad,
                         "if_J_in_deg_noise_in_rad" => J_in_deg, "phase_only" => phase_only, "amplitude_only_1pct" => amp_only)
    end
    res["units"] = units
    println("units ", units); flush(stdout)

    # 5. scenarios (model-internal sensitivities)
    scen = Dict{String,Any}()
    function sig_for(b, f; sa = 0.01, prior = Dict{String,Float64}())
        _, J, _ = lnH_jac(b, f; rtol = 1e-11)
        out = Dict{String,Any}()
        for Nm in ("N0", "N1", "N2", "N3"), sp in (0.01, 0.1)
            A = whiten(design_cols(J, f, SETS[Nm]), sa, sp)
            # Gaussian priors (calibration) as extra rows
            for (n, sd) in prior
                k = findfirst(==(n), SETS[Nm])
                k === nothing && continue
                row = zeros(1, size(A, 2)); row[1, k] = 1 / sd
                A = vcat(A, row)
            end
            out["$Nm|$sp"] = sigma_proj(A, 1, collect(2:size(A, 2)))
        end
        return out
    end
    for Lk in ("300", "500")
        b = base_for(Lk)
        scen["$Lk|base"] = sig_for(b, freqs)
        scen["$Lk|sigmaA_0.1pct"] = sig_for(b, freqs; sa = 0.001)
        scen["$Lk|R1_3e-9_R2_5e-9"] = sig_for(base_for(Lk; over = Dict{String,Any}("R1" => 3e-9, "R2" => 5e-9)), freqs)
        scen["$Lk|R1_1e-9_R2_1e-9"] = sig_for(base_for(Lk; over = Dict{String,Any}("R1" => 1e-9, "R2" => 1e-9)), freqs)
        scen["$Lk|dAu_20nm"] = sig_for(base_for(Lk; over = Dict{String,Any}("dAu" => 20e-9)), freqs)
        scen["$Lk|beams_3_2.4um"] = sig_for(base_for(Lk; over = Dict{String,Any}("wp" => 3e-6, "wr" => 2.4e-6)), freqs)
        scen["$Lk|beams_30_24um"] = sig_for(base_for(Lk; over = Dict{String,Any}("wp" => 30e-6, "wr" => 24e-6)), freqs)
        scen["$Lk|320pts"] = sig_for(b, collect(exp.(range(log(1e4), log(2e8), length = 320))))
        scen["$Lk|fmax_1GHz_same_density"] = sig_for(b, collect(exp.(range(log(1e4), log(1e9), length = 93))))
        scen["$Lk|tau_x1.044"] = sig_for(base_for(Lk; over = Dict{String,Any}("tauR" => b.tauR * 1.044)), freqs)
        scen["$Lk|prior_phi0_0.01deg_td_1ps"] = sig_for(b, freqs; prior = Dict("phi0" => deg2rad(0.01), "td_ns" => 1e-3))
        scen["$Lk|prior_phi0_0.01deg_td_1ps_R_1pct"] = sig_for(b, freqs; prior = Dict("phi0" => deg2rad(0.01), "td_ns" => 1e-3,
                                                                                      "ln_R1" => 0.01, "ln_R2" => 0.01))
        scen["$Lk|all_favourable"] = sig_for(base_for(Lk; over = Dict{String,Any}("R1" => 3e-9, "R2" => 5e-9, "dAu" => 20e-9)),
                                             collect(exp.(range(log(1e4), log(2e8), length = 320))); sa = 0.001)
        println("scenarios ", Lk, " done ", time() - t0); flush(stdout)
    end
    res["scenarios"] = scen
    writejson(joinpath(OUT, "aud_julia_part2.json"), res)
    println("done ", time() - t0)
end

main2()
