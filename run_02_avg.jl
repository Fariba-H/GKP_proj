import Plots
import Dates


include("utils.jl")
include("2mode_output_state.jl")
include("pdf_sampling.jl")


function run(Δ=0.2, K=30; zr_thr=0, prc=0.95)
    if zr_thr < 0 || zr_thr > 1
        error("zr_thr must be in [0, 1], got $zr_thr")
    end
    a = √(π / 2)
    q_scl = √2a
    rng_q_max = 40 # ceil(estimate_qrng(Δ) / q_scl)
    println("rng_q_max: $rng_q_max")
    n_qs = 200001
    n_bins = 300
    θrange = collect(0.38012pi:0.00001π:0.38248pi)



    date_time = Dates.format(Dates.now(), "yyyymmddHHMMSS")
    output_dir = "outputs/$(date_time)/"
    mkpath(output_dir)
    println("Output directory: $output_dir")
    # Write parameters to text file
    open(joinpath(output_dir, "params.txt"), "a") do f
        write(f, "th_rot_range/π: $(θrange[1]/pi),$(θrange[end]/pi)\nΔ: $Δ\nK: $K\nn_qs: $n_qs\nq_scl: $q_scl\nq_range_max: $rng_q_max\nn_bins: $n_bins\n")
        write(f, "zr_thr: $zr_thr, prc: $prc\n")
        # write(f, "n_samples: $n_samples\n")
    end

    pp, pm, mm = pp_pm_mm(Δ, K)

    pdf_thph_all = zeros(n_bins, n_bins)

    # Define bin edges for θ and φ
    θ_grid = range(0, π, length=n_bins + 1) |> collect
    φ_grid = range(-π, π, length=n_bins + 1) |> collect
    θ_bin = (θ_grid[1:end-1] + θ_grid[2:end]) / 2
    φ_bin = (φ_grid[1:end-1] + φ_grid[2:end]) / 2

    l_prog = -1
    for (i, θr) in enumerate(θrange)

        rng_q = 40 # ceil(estimate_qrng(θr, Δ) / q_scl)

        SpSp = rotated_00(θr, Δ, K)
        SmSm = rotated_11(θr, Δ, K)
        SpSm = rotated_01(θr, Δ, K)

        pdf_q, qs, ths, phs = calc_pdf_q(θr, Δ, SpSp, SpSm, SmSm, pp, pm, mm, n_qs, q_scl, rng_q)
        pdf_thph, _ = calc_pdf_2d(pdf_q, qs, ths, phs, n_bins)

        if zr_thr > 0
            hprc = Statistics.quantile(vec(pdf_thph), prc)
            zr_inds = pdf_thph .< zr_thr * hprc
            pdf_thph[zr_inds] .= 0
            pdf_thph ./= sum(pdf_thph) * (π / n_bins) * (2π / n_bins)  # normalize
        end
        pdf_thph_all .+= pdf_thph
        # sampled_qs, sampled_thetas, sampled_phis = sampling(θr, Δ, pdf_q, qs, n_samples)

        prog = 100 * i / length(θrange)
        if prog >= l_prog + 1
            l_prog = floor(Int, prog)
            print("\rprog: $(lpad(l_prog, 3))% \tθr/π = $(rpad(round(θr/π; digits=5), 7, '0'))   ")
        end
    end
    pdf_thph_all = pdf_thph_all ./ length(θrange)
    save_heatmap(pdf_thph_all, θ_bin, φ_bin, 
    "Average joint probability density of θ and φ\n
     over θr in [$(round(θrange[1]/π; digits=6))...$(round(θrange[end]/π; digits=6))]π, Δ: $Δ\n",
     output_dir, "pdf_thph_hm_all"; prc=0.98, color=Plots.cgrad(:greys, rev=true))
    println()
end


for (Δ, K) in [(sqrt(0.001), 80)] #(0.1, 30), (0.05, 50), (0.075, 30) , (0.3, 30)
    β = round(Δ^2; digits=10)
    println("Δ = $Δ, β = $β, K = $K ----", Dates.format(Dates.now(), " yyyy-mm-dd HH:MM:SS"))
    run(Δ, K;zr_thr=0)
end

## END
