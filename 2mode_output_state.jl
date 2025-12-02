
function rotated_00(θ::Real, Δ::Real, K::Integer)
    w = exp(2 * (im * θ - Δ^2))
    s = Complex{Float64}(0, 0)
    for j in -K:K
        for jp in -K:K  # Using jp instead of j'  
            exponent = π * ((2 * (2j) * 2jp * w - ((2j)^2 + (2jp)^2) * w^2) / (1 - w^2) - (2j)^2 / 2 - (2jp)^2 / 2)
            x = exp(exponent)
            s += x
        end
    end
    s
end

function rotated_11(θ::Real, Δ::Real, K::Integer)
    w = exp(2 * (im * θ - Δ^2))
    s = Complex{Float64}(0, 0)
    for j in -K:K
        for jp in -K:K  # Using jp instead of j' 
            exponent = π * ((2 * (2j + 1) * (2jp + 1) * w - ((2j + 1)^2 + (2jp + 1)^2) * w^2) / (1 - w^2) - (2j + 1)^2 / 2 - (2jp + 1)^2 / 2)
            x = exp(exponent)
            s += x
        end
    end
    s
end

function rotated_01(θ::Real, Δ::Real, K::Integer)
    w = exp(2 * (im * θ - Δ^2))
    s = Complex{Float64}(0, 0)
    for j in -K:K
        for jp in -K:K  # Using jp instead of j' 
            exponent = π * ((2 * (2j) * (2jp + 1) * w - ((2j)^2 + (2jp + 1)^2) * w^2) / (1 - w^2) - (2j)^2 / 2 - (2jp + 1)^2 / 2)
            x = exp(exponent)
            s += x
        end
    end
    s
end


# prob. measuring q
function probability(SpSp, SpSm, SmSm, pp, pm, mm, Cp, Cm)
    mp = pm'
    SmSp = SpSm'
    abs(Cp' * Cp * pp + Cp' * Cm * pm + Cm' * Cp * mp + Cm' * Cm * mm) / # |<q|⊗I |ψ_2>|^2
    abs(SpSp * pp + SpSm * pm + SmSp * mp + SmSm * mm)  # <ψ_2|ψ_2>
end


# Inner product of fock-damped plus-plus, plus-minus, minus-minus states
function pp_pm_mm(Δ::Real, K::Integer)
    w = exp(-2 * Δ^2)
    s = zeros(Complex{Float64}, 3)
    for j in -K:K
        fj = (j % 2 == 0) ? 1 : -1
        for jp in -K:K
            factor = [1, fj, fj * ((jp % 2 == 0) ? 1 : -1)]
            exponent = π * ((2 * j * jp * w - (j^2 + jp^2) * w^2) / (1 - w^2) - j^2 / 2 - jp^2 / 2)
            s += factor .* exp(exponent)
        end
    end
    s
end


# Fidelity with target state; both written in plus-minus basis
function fidelity(Cp::Complex, Cm::Complex, Trg_p::Complex, Trg_m::Complex, pp::Complex, pm::Complex, mm::Complex)
    # Magic target:
    # Trg_p = cos(pi / 8)
    # Trg_m = -im * sin(pi / 8)
    abs(Trg_p' * Cp * pp + Trg_m' * Cm * mm + Trg_p' * Cm * pm + Trg_m' * Cp * pm')^2 /          # |<Target|Out>|^2
    abs(Trg_p' * Trg_p * pp + Trg_m' * Trg_m * mm + Trg_p' * Trg_m * pm + Trg_m' * Trg_p * pm') /  # <Target|Target>
    abs(Cp' * Cp * pp + Cm' * Cm * mm + Cp' * Cm * pm + Cm' * Cp * pm')                          # <Out|Out>
end


function Cplus(q::Real, θ::Real, Δ::Real, K::Integer=1000; cutoff::Real=1e-5, m=nothing, n=nothing)
    ed = exp(-Δ^2)
    if isnothing(n) || isnothing(m)
        w = exp(im * θ)
        f1 = π / ((w * ed)^(-2) - 1)
        f2 = 2q / ((w * ed) * √π)
        p1 = -q^2 / 2 - q^2 / ((w * ed)^(-2) - 1)
    else
        # w = ((n^2 - m^2) + im(2m * n)) / (n^2 + m^2)
        f1 = π * (m - im * n)^2 / ((m + im * n)^2 / ed^2 - (m - im * n)^2)
        f2 = 2q / √π / ed * (m + im * n) / (im * n - m)
        p1 = -q^2 / 2 - q^2 * (m - im * n)^2 / ((m + im * n)^2 / ed^2 - (m - im * n)^2)
    end

    #f3 = exp(p1) / √(1 - w^2)
    x0 = exp(p1)
    s = x0
    for j in 1:K÷2
        x = exp(f1 * (f2 * (2j) - (2j)^2) - ((2j)^2 * π) / 2 + p1) + exp(f1 * (f2 * (-2j) - (2j)^2) - ((2j)^2 * π) / 2 + p1)
        if isnan(x)
            println("In Cplus: θ= $(round(θ/pi; digits=3)), q= $(round(q; digits=3)), j= $j, x= $x")
        end
        s += x
        if isnan(s)
            println("In Cplus: θ= $(round(θ/pi; digits=3)), q= $(round(q; digits=3)), j= $j, s= $s")
            println("In Cplus: f1= $(round(f1; digits=3)), f2= $(round(f2; digits=3)), p1= $(round(p1; digits=3))")
            return s
        end
        if j >= 10 && abs(x) < cutoff * abs(x0)
            # println("In Cplus: break j= $j, x= $x, x0= $x0, s= $s")
            break
        end
    end
    s
end

function Cminus(q::Real, θ::Real, Δ::Real, K::Integer=1000; cutoff::Real=1e-5, m=nothing, n=nothing)
    ed = exp(-Δ^2)
    if isnothing(n) || isnothing(m)
        w = exp(im * θ)
        f1 = π / ((w * ed)^(-2) - 1)
        f2 = 2q / ((w * ed) * √π)
        p1 = -q^2 / 2 - q^2 / ((w * ed)^(-2) - 1)
    else
        # w = ((n^2 - m^2) + im(2m * n)) / (n^2 + m^2)
        f1 = π * (m - im * n)^2 / ((m + im * n)^2 / ed^2 - (m - im * n)^2)
        f2 = 2q / √π / ed * (m + im * n) / (im * n - m)
        p1 = -q^2 / 2 - q^2 * (m - im * n)^2 / ((m + im * n)^2 / ed^2 - (m - im * n)^2)
    end
    #f3 = exp(p1) / √(1 - w^2)
    x0 = exp(f1 * (f2 - 1) - π / 2 + p1)
    s = x0
    for j in 1:K÷2
        x = exp(f1 * (f2 * (2j + 1) - (2j + 1)^2) - ((2j + 1)^2 * π) / 2 + p1) + exp(f1 * (f2 * (-2j + 1) - (-2j + 1)^2) - ((-2j + 1)^2 * π) / 2 + p1)
        if isnan(x)
            println("In Cminus: θ= $(round(θ/pi; digits=3)), q= $(round(q; digits=3)), j= $j, x= $x")
        end
        s += x
        if isnan(s)
            println("In Cminus: θ= $(round(θ/pi; digits=3)), q= $(round(q; digits=3)), j= $j, s= $s")
            println("In Cminus: f1= $(round(f1; digits=3)), f2= $(round(f2; digits=3)), p1= $(round(p1; digits=3))")
            return s
        end
        if j >= 10 && abs(x) < cutoff * abs(x0)
            break
        end
    end
    s
end