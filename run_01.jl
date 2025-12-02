import Plots
import Dates

include("utils.jl")
include("2mode_output_state.jl")
include("pdf_sampling.jl")

function run(Δ=0.1, K=30, θr=0.38467pi, output_dir="outputs/")
    β = Δ^2
    a = √(π / 2)
    q_scl = √2a
    rng_q_max = 40 
    n_qs = 2000001

    n_bins = 300

    # Magic state in plus-minus basis
    Trg_p, Trg_m = cos(pi / 8) + 0.0im, 0.0 - im * sin(pi / 8) 

    pp, pm, mm = pp_pm_mm(Δ, K)

    rng_q = 40

    SpSp = rotated_00(θr, Δ, K)
    SmSm = rotated_11(θr, Δ, K)
    SpSm = rotated_01(θr, Δ, K)

    pdf_q, qs, ths, phs = calc_pdf_q(θr, Δ, SpSp, SpSm, SmSm, pp, pm, mm, n_qs, q_scl, rng_q)
    pdf_theta, theta_bins = calc_pdf_theta(pdf_q, qs, ths, n_bins)
    pdf_phi, phi_bins = calc_pdf_phi(pdf_q, qs, phs, n_bins)
    pdf_θφ, g_θs, g_φs = calc_pdf_2d(pdf_q, qs, ths, phs, n_bins)
    fidelities = calc_fidelities_q(qs, θr, Δ, Trg_p, Trg_m, pp, pm, mm)
    local_maxima_indices = find_local_maxima(pdf_q, 1)  # n_q_samples ÷ 10)

    # Get q values and probabilities at maxima
    q_maxima = qs[local_maxima_indices]
    # Get corresponding theta values
    theta_maxima = ths[local_maxima_indices]
    phi_maxima = phs[local_maxima_indices]

    cmn_title = "θr: $(round(θr/π; digits=7))π, β: $(round(β; digits=6))"
    cmn_file_name = "rot$(rpad(round(θr/π; digits=7), 7, '0'))pi_beta$(round(β; digits=6))"

    open(joinpath(output_dir, cmn_file_name * "_maxima.txt"), "a") do f
        write(
            f,
            "theta/pi: $(round((θr/pi); digits=7)), sin: $(round(sin(θr); digits=5))
   , cos: $(round(cos(θr); digits=5)), 1/sin: $(round(1 / sin(θr); digits=5))\n"
        )
        write(f, "=================================================\n")
        write(f, "q_max / q_scl, theta/pi, phi/pi\n")
        write(f, "-------------\n")
        for i in eachindex(q_maxima)
            q = q_maxima[i]
            th = theta_maxima[i]
            ph = phi_maxima[i]
            if 0 <= q <= 10 * q_scl
                write(
                    f,
                    "$(round(q/q_scl, digits=5)), $(round(th/pi; digits=3)),$(round(ph/pi; digits=3))\n"
                )
            end
        end
        # New lines: printing combinations of n and m, sorted by their value
        x_rng = 4 * q_scl  # The upper bound
        range_n = -5:5
        range_m = -5:5

        # Collect all valid combinations in an array
        combos = Tuple{Float64,Int,Int}[]  # (value, n, m)

        for n in range_n
            for m in range_m
                val = n * abs(sin(θr)) + m * abs(cos(θr))
                if 0 <= val <= x_rng
                    push!(combos, (val, n, m))
                end
            end
        end

        # Sort by 'val' (the first item in each tuple)
        sort!(combos, by=x -> x[1])

        write(f, "\n--- n|sin(θr)| + m|cos(θr)| (with 0 ≤ value ≤ 4√π), sorted by value ---\n")

        #Now print in ascending order
        for (val, n, m) in combos
            write(f, "$(n)|sin(θr)| + $(m)|cos(θr)| = $(round(val, digits=5))   (n=$(n), m=$(m))\n")
        end
    end
    # Write parameters to text file
    open(joinpath(output_dir, cmn_file_name * "_params.txt"), "a") do f
        write(f, "th_rot/π: $(θr/pi)\nβ: $β\nK: $K\nn_qs: $n_qs\nq_scl: $q_scl\nq_range_max: $rng_q_max\nn_bins: $n_bins\n")
        #  write(f, "n_samples: $n_samples\n")
    end
    a = round(q_scl; digits=5)
    save_scatter_plot(qs, pdf_q, q_scl, 1, "Probability density of q, $cmn_title",
        "q (/√2 a = $a )", "PDF q", output_dir, cmn_file_name * "_pdf_q", x_rng=(0, 2); x_ticks=-6:6)
    save_scatter_plot(theta_bins, pdf_theta, π, 1, "Probability density of θ, " * cmn_title, "θ / π", "PDF θ", output_dir, cmn_file_name * "_pdf_theta")
    # Plot q vs theta at maxima points
    save_scatter_plot(q_maxima, theta_maxima, q_scl, π, "θ vs q at PDF_q maxima, " * cmn_title, "q / $(round(q_scl, digits=4))", "θ / π", output_dir, cmn_file_name * "_max_q_theta", x_rng=(0, 6), ylim=(0, 1))
    save_scatter_plot(phi_bins, pdf_phi, π, 1, "Probability density of ϕ, " * cmn_title, "ϕ / π", "PDF ϕ", output_dir, cmn_file_name * "_pdf_phi")
    # Plot q vs phi at maxima points
    save_scatter_plot(q_maxima, phi_maxima, q_scl, π, "ϕ vs q at PDF_q maxima, " * cmn_title, "q / $(round(q_scl, digits=4))", "ϕ / π", output_dir, cmn_file_name * "_max_q_phi", x_rng=(0, 6), ylim=(-1, 1))
    save_scatter_plot(qs, fidelities, q_scl, 1, "Fidelity vs q, $cmn_title",
    "q (/√2 a = $a )", "Fidelity", output_dir, cmn_file_name * "_fidelity_q", x_rng=(0, 10); x_ticks=-6:6, ylim=(0, 1))

    save_heatmap(pdf_θφ, g_θs, g_φs, "Joint PDF of θ and ϕ, $cmn_title",
        output_dir, cmn_file_name * "_pdf_2d"; prc=0.99, color=Plots.cgrad(:greys, rev=true))

end


## START

date_time = Dates.format(Dates.now(), "yyyymmddHHMMSS")
output_dir = "outputs/$(date_time)/"
mkpath(output_dir)
println("Output directory: $output_dir")

βrange = [0.04, 0.01, 0.001]
θrange = [π / 4, 0.0681π, 0.38467π]
for (β, θr) in zip(βrange, θrange)
    Δ = √β
    println("Δ= $(round(Δ; digits=6)), θr: $(rpad(round(θr/π; digits=7), 7, '0'))π ---", Dates.format(Dates.now(), " yyyy-mm-dd HH:MM:SS"))
    run(Δ, 70, θr, output_dir)
end

## END
