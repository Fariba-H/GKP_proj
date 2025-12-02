
function calc_pdf_q(θ, Δ, SpSp, SpSm, SmSm, pp, pm, mm, n_qs, scl, rng)
    if iseven(n_qs)
        n_qs += 1
    end
    qs = range(-rng * scl, rng * scl, length=n_qs) |> collect
    ths = zeros(n_qs)
    phs = zeros(n_qs)
    pdf_q = zeros(n_qs)
    for (i, q) in enumerate(qs)
        Cp = Cplus(q, θ, Δ)
        Cm = Cminus(q, θ, Δ)
        if isnan(Cp) || isnan(Cm)
            println("in calc_pdf_q:  θ/π= $(round(θ/pi; digits=3)), q= $(round(q; digits=3)), Cp= $Cp, Cm= $Cm")
        end
        # |out> = cos th/2 |+> + e^iph sin th/2 |->
        ths[i] = 2 * atan(abs(Cm), abs(Cp))
        phs[i] = (Cp == 0) ? 0 : angle(Cm / Cp)
        pdf_q[i] = probability(SpSp, SpSm, SmSm, pp, pm, mm, Cp, Cm)
        if isnan(pdf_q[i])
            println("in calc_pdf_q:  θ/π= $(round(θ/pi; digits=3)), q= $(round(q; digits=3)), Cp= $Cp, Cm= $Cm, pdf_q= $pdf_q")
        end
    end
    # Normalize the probability density
    sprq = simpson_integrate(qs, pdf_q)
    # println("Integral of probability density of q: $(round(sprq; digits=3))")
    pdf_q = pdf_q ./ sprq
    return pdf_q, qs, ths, phs
end

#finding probability of the each θ(outputstate)
function calc_pdf_theta(pdf_q, qs, ths, n_bins=1000)
    theta_grid = range(0, π, length=n_bins + 1) |> collect
    pdf_theta = zeros(n_bins)           # will store g(θ) values
    n_qs = length(qs)
    δθ = π / n_bins                      # θ spacing
    for i in 1:n_qs-1
        prbi = (pdf_q[i] + pdf_q[i+1]) * (qs[i+1] - qs[i]) / 2
        thi_min, thi_max = minmax(ths[i], ths[i+1])
        Δθ = thi_max - thi_min
        j_min = Int(floor(thi_min / δθ)) + 1
        j_max = Int(ceil(thi_max / δθ))
        if j_min == j_max || Δθ == 0
            pdf_theta[j_min] += prbi / δθ
            continue
        end
        ΔI = theta_grid[j_min+1] - thi_min
        pdf_theta[j_min] += prbi / δθ * (ΔI / Δθ)
        for j in j_min+1:j_max-1
            pdf_theta[j] += prbi / Δθ
        end
        ΔI = thi_max - theta_grid[j_max]
        pdf_theta[j_max] += prbi / δθ * (ΔI / Δθ)
    end
    return pdf_theta, δθ / 2 .+ theta_grid[1:n_bins]
end

# Finding probability of each φ (output state)
function calc_pdf_phi(pdf_q, qs, phs, n_bins=1000)
    phg_min = -π
    phg_max = π
    phi_grid = range(phg_min, phg_max, length=n_bins + 1) |> collect
    pdf_phi = zeros(n_bins)                # will store g(φ) values
    δφ = (phg_max - phg_min) / n_bins      # φ spacing
    for i in 1:length(qs)-1
        prbi = (pdf_q[i] + pdf_q[i+1]) * (qs[i+1] - qs[i]) / 2
        phi_min, phi_max = minmax(phs[i], phs[i+1])
        Δφ = phi_max - phi_min
        j_min = Int(floor((phi_min - phg_min) / δφ)) + 1
        j_max = Int(ceil((phi_max - phg_min) / δφ))
        if j_min == j_max || Δφ == 0
            pdf_phi[j_min] += prbi / δφ
            continue
        end
        ΔI = phi_grid[j_min+1] - phi_min
        pdf_phi[j_min] += prbi / δφ * (ΔI / Δφ)
        for j in j_min+1:j_max-1
            pdf_phi[j] += prbi / Δφ
        end
        ΔI = phi_max - phi_grid[j_max]
        pdf_phi[j_max] += prbi / δφ * (ΔI / Δφ)
    end
    return pdf_phi, δφ / 2 .+ phi_grid[1:n_bins]
end


function calc_pdf_2d(pdf_q, qs, ths, phs, n_bins=1000)
    # Define the number of bins for θ and φ
    np_bins_θ = n_bins
    np_bins_φ = n_bins
    # Define bin edges for θ and φ
    θ_grid = range(0, π, length=np_bins_θ + 1) |> collect
    φ_grid = range(-π, π, length=np_bins_φ + 1) |> collect
    # Initialize joint probability density array
    prbs_θφ = zeros(np_bins_θ, np_bins_φ)
    δθ = π / np_bins_θ   # θ bin spacing
    δφ = 2π / np_bins_φ  # φ bin spacing
    n_qs = length(qs)  # Number of q samples
    for i in 1:n_qs-1
        prbi = (pdf_q[i] + pdf_q[i+1]) * (qs[i+1] - qs[i]) / 2
        # Find θ bin indices
        θ_min, θ_max = minmax(ths[i], ths[i+1])
        Δθ = θ_max - θ_min
        j_min_θ = Int(floor(θ_min / δθ)) + 1
        j_max_θ = Int(ceil(θ_max / δθ))
        # Find φ bin indices
        φ_min, φ_max = minmax(phs[i], phs[i+1])
        Δφ = φ_max - φ_min
        j_min_φ = Int(floor((φ_min + π) / δφ)) + 1
        j_max_φ = Int(ceil((φ_max + π) / δφ))
        for jθ in j_min_θ:j_max_θ
            ΔI_θ = min(θ_grid[jθ+1], θ_max) - max(θ_grid[jθ], θ_min)
            w_θ = (Δθ == 0) ? 1 : ΔI_θ / Δθ
            w_θ = (w_θ > 1) ? 1 : w_θ
            for jφ in j_min_φ:j_max_φ
                ΔI_φ = min(φ_grid[jφ+1], φ_max) - max(φ_grid[jφ], φ_min)
                w_φ = (Δφ == 0) ? 1 : ΔI_φ / Δφ
                w_φ = (w_φ > 1) ? 1 : w_φ
                prbs_θφ[jθ, jφ] += prbi / (δθ * δφ) * w_θ * w_φ
            end
        end
    end
    # # Print indices where NaN values occur
    # nan_indices = findall(isnan.(prbs_θφ))
    # if !isempty(nan_indices)
    #     println("θr: $θ; 5 NaN vals:")
    #     for idx in nan_indices[1:5]
    #         println("(θ,φ) = ($(round(θ_grid[idx[1]]/π, digits=3))π, $(round(φ_grid[idx[2]]/π, digits=3))π)")
    #     end
    #     continue
    # end
    # Normalize the joint probability density
    # prb_int = sum(prbs_θφ) * δθ * δφ
    # println("Integral of joint probability density: $(round(prb_int; digits=3))")
    # if prb_int != 0
    #     prbs_θφ /= prb_int
    # end
    return prbs_θφ, δθ / 2 .+ θ_grid[1:end-1], δφ / 2 .+ φ_grid[1:end-1]
end


function sampling_q(cdf, qs)
    # Sample a random number from a uniform distribution between 0 and 1
    r = rand()
    # Perform bisection search to find the interval in which r falls
    index = searchsortedfirst(cdf, r)
    if index == 1 || cdf[index] == cdf[index-1]
        return qs[index]
    end
    # Linear interpolation
    q = ((cdf[index] - r) * qs[index] + (r - cdf[index-1]) * qs[index-1]) / (cdf[index] - cdf[index-1])
    return q
end


function sampling(θ, Δ, pdf_q, qs, ns=10001)
    δq = qs[2] - qs[1]
    # Compute the CDF using the cumulative sum. Multiply by the bin width δq to approximate the integral.
    cdf = cumsum(pdf_q) * δq
    sampled_qs = [sampling_q(cdf, qs) for _ in 1:ns]
    sampled_thetas = zeros(ns)
    sampled_phis = zeros(ns)
    for i in 1:ns
        q = sampled_qs[i]
        Cp = Cplus(q, θ, Δ)
        Cm = Cminus(q, θ, Δ)
        sampled_thetas[i] = 2 * atan(abs(Cm), abs(Cp))
        sampled_phis[i] = angle(Cm / Cp)
    end
    return sampled_qs, sampled_thetas, sampled_phis
end


function calc_fidelities_q(qs::Vector{<:Real}, θr::Real, Δ::Real, Trg_p::Complex, Trg_m::Complex, pp::Complex, pm::Complex, mm::Complex; K::Integer=100, cutoff::Real=1e-10)
    out_fids = zeros(length(qs))
    for (i, q) in enumerate(qs)
        Cp = Cplus(q, θr, Δ, K; cutoff=cutoff)
        Cm = Cminus(q, θr, Δ, K; cutoff=cutoff)
        fid = fidelity(Cp, Cm, Trg_p, Trg_m, pp, pm, mm)
        if !isfinite(fid)
            println(stderr, "in calc_fidelities_q:  θ/π= $(round(θr/pi; digits=8)), q= $(round(q; digits=4)), Cp= $Cp, Cm= $Cm, fid= $fid")
            continue
        end
        out_fids[i] = fid
    end
    return out_fids
end
