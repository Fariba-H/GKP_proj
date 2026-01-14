# <q| e^((-Delta^2 + i theta_r) hat{n}) |qp>
# w = e^(-Delta^2 + i theta_r),   Delta^2 = beta
kernel(w::Number, q::Real, qp::Real) =
    exp(-(1 + w^2) / (2 * (1 - w^2)) * (q^2 + qp^2) + (2 * w) / (1 - w^2) * q * qp) / sqrt(π * (1 - w^2))

function estim_K(w::Real, cutoff::Real=1e-10)
    epsilon = cutoff * kernel(w, 0, 0)
    C = -(2 * (1 - w^2)) * log(sqrt(π * (1 - w^2)) * epsilon)
    ceil(Int, sqrt(C / (π * (1 - w)^2)))
end

# |psi_out> = Cplus |~+> + Cminus |~->
function Cplus(q::Real, θr::Real, Δ::Real, K::Integer=1000; cutoff::Real=1e-10, min_j=5)
    gamma = ComplexF64(-Δ^2, θr)
    w = exp(gamma)
    x = kernel(w, q, 0)
    x_m = abs(x)
    s = x
    breaked = false
    for j in 2:2:K
        x = kernel(w, q, j * √π) + kernel(w, q, -j * √π)
        s += x
        if abs(x) > x_m
            x_m = abs(x)
        elseif j >= min_j && abs(x) < cutoff * x_m
            breaked=true
            break
        end
    end
    if !breaked
        println(stderr, "Warning: Cplus did not converge for q=$(q), θr=$(θr/pi)π, Δ=$(Δ). Consider increasing K.")
    end
    s
end

# |psi_out> = Cplus |~+> + Cminus |~->
function Cminus(q::Real, θr::Real, Δ::Real, K::Integer=1000; cutoff::Real=1e-10, min_j=5)
    x_m = 0
    gamma = ComplexF64(-Δ^2, θr)
    w = exp(gamma)
    s = Complex{Float64}(0, 0)
    breaked = false
    for j in 1:2:K
        x = kernel(w, q, j * √π) + kernel(w, q, -j * √π)
        s += x
        if abs(x) > x_m
            x_m = abs(x)
        elseif j >= min_j && abs(x) < cutoff * x_m
            breaked=true
            break
        end
    end
    if !breaked
        println(stderr, "Warning: Cminus did not converge for q=$(q), θr=$(θr/pi)π, Δ=$(Δ). Consider increasing K.")
    end
    s
end

# <~0|~0>
function zero_zero(Δ::Real, K::Integer=300; cutoff::Real=1e-10)
    w = exp(-2 * Δ^2)
    nK = min(K, estim_K(w, cutoff))
    if nK == K
        println(stderr, "Warning: K=$K may be insufficient for Δ=$(round(Δ; digits=4)). Consider increasing K to at least $(estim_K(w, cutoff)).")
    end
    K = nK
    s = 0
    for j in -K÷2:K÷2
        for jp in -K÷2:K÷2  # Using jp instead of j'
            s += kernel(w, 2j * √π, 2jp * √π)
        end
    end
    s
end

# <~1|~1>
function one_one(Δ::Real, K::Integer=300; cutoff::Real=1e-10)
    w = exp(-2 * Δ^2)
    nK = min(K, estim_K(w, cutoff))
    if nK == K
        println(stderr, "Warning: K=$K may be insufficient for Δ=$(round(Δ; digits=4)). Consider increasing K to at least $(estim_K(w, cutoff)).")
    end
    K = nK
    s = 0
    for j in -K÷2:K÷2
        for jp in -K÷2:K÷2  # Using jp instead of j'
            s += kernel(w, (2j + 1) * √π, (2jp + 1) * √π)
        end
    end
    s
end

# <~0|~1>
function zero_one(Δ::Real, K::Integer=300; cutoff::Real=1e-10)
    w = exp(-2 * Δ^2)
    nK = min(K, estim_K(w, cutoff))
    if nK == K
        println(stderr, "Warning: K=$K may be insufficient for Δ=$(round(Δ; digits=4)). Consider increasing K to at least $(estim_K(w, cutoff)).")
    end
    K = nK
    s = 0
    for j in -K÷2:K÷2
        for jp in -K÷2:K÷2  # Using jp instead of j'
            s += kernel(w, 2j * √π, (2jp + 1) * √π)
        end
    end
    s
end

