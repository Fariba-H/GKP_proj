
import Statistics as st
import Plots
import DelimitedFiles


"""
    simpson_nonuniform(x, y)

Composite Simpson-like integration for sorted, not-necessarily-uniform x.
- Uses quadratic interpolation on each consecutive triple (x[i],x[i+1],x[i+2])
  and integrates that quadratic exactly on [x[i], x[i+2]].
- If only two points are given, falls back to the trapezoidal rule.
- If the number of points is even, applies the quadratic rule on the first n-1
  points and a trapezoid on the final interval.

Requires x strictly increasing and length(x) == length(y).
"""
function simpson_nonuniform(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    if n != length(y)
        throw(ArgumentError("x and y must have the same length"))
    end
    if n == 0
        return zero(promote_type(eltype(x), eltype(y)))
    elseif n == 1
        return zero(promote_type(eltype(x), eltype(y)))
    elseif n == 2
        return 0.5 * (y[1] + y[2]) * (x[2] - x[1])
    end

    for i in 1:n-1
        if !(x[i+1] > x[i])
            throw(ArgumentError("x must be strictly increasing"))
        end
    end

    T = promote_type(eltype(x), eltype(y))
    total = zero(T)

    m = n
    use_trap_last = false
    if iseven(n)
        m = n - 1
        use_trap_last = true
    end

    i = 1
    while i <= m - 2
        x0 = float(x[i])
        x1 = float(x[i+1])
        x2 = float(x[i+2])
        y0 = float(y[i])
        y1 = float(y[i+1])
        y2 = float(y[i+2])

        # Newton divided differences form:
        f01 = (y1 - y0) / (x1 - x0)
        f12 = (y2 - y1) / (x2 - x1)
        f012 = (f12 - f01) / (x2 - x0)

        # Integrals of basis terms on [x0, x2]
        I0 = x2 - x0                              # ∫ 1 dx
        I1 = (x2 - x0)^2 / 2                      # ∫ (x-x0) dx
        # I2 = ∫ (x-x0)*(x-x1) dx = [x^3/3 - (x0+x1)x^2/2 + x0*x1*x]_{x0}^{x2}
        F(x) = x^3 / 3 - (x0 + x1) * x^2 / 2 + x0 * x1 * x
        I2 = F(x2) - F(x0)

        # Integral of quadratic interpolant
        total += y0 * I0 + f01 * I1 + f012 * I2

        i += 2
    end

    if use_trap_last
        total += 0.5 * (y[end] + y[end-1]) * (x[end] - x[end-1])
    end

    return convert(T, total)
end


# Integrate using Simpson's rule; assumes uniform spacing in x
function simpson_integrate(x, y)
    n = length(x)
    if !isodd(n)
        # throw(ArgumentError("Number of points must be odd for Simpson's rule"))
        println(stderr, "Warning: Number of points is even; dropping last point for Simpson's rule.")
        n -= 1
    end
    # if isodd(n)
    #     yn = y[n]
    # else
    #     n += 1
    #     yn = 4y[n] - 6y[n-1] + 4y[n-2] - y[n-3]
    # end
    h = x[2] - x[1]
    sum = (y[1] + y[n])
    for i in 2:2:n-1
        sum += 4 * y[i]
    end
    for i in 3:2:n-2
        sum += 2 * y[i]
    end
    return sum * h / 3
end

"""
    trapezoid_nonuinform(x, y)

Composite trapezoidal rule for (possibly non-uniform) sorted x.
- Requires x strictly increasing and length(x) == length(y).
- Returns the numeric integral using the trapezoid rule.
"""
function trapezoid_nonuinform(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    if n != length(y)
        throw(ArgumentError("x and y must have the same length"))
    end
    if n < 2
        throw(ArgumentError("At least two points are required for trapezoidal integration"))
    end

    for i in 1:n-1
        if !(x[i+1] > x[i])
            throw(ArgumentError("x must be strictly increasing"))
        end
    end

    T = promote_type(eltype(x), eltype(y))
    total = zero(T)
    for i in 1:n-1
        total += 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])
    end

    return convert(T, total)
end


# Composite trapezoidal rule for uniform spacing in x
function trapezoid_integrate(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    if n != length(y)
        throw(ArgumentError("x and y must have the same length"))
    end
    if n < 2
        throw(ArgumentError("At least two points are required for trapezoidal integration"))
    end

    T = promote_type(eltype(x), eltype(y))
    h = x[2] - x[1]
    if n == 2
        return convert(T, 0.5 * (y[1] + y[2]) * h)
    end
    s = 0.5 * (y[1] + y[n]) + sum(y[2:n-1])
    return convert(T, s * h)
end


function save_scatter_plot(x, y, x_scl, y_scl, title, xlabel, ylabel, output_dir, file_name; n_samples=nothing, x_rng=nothing, ylim=nothing, xlim=nothing, x_ticks=nothing, save_data=false)
    # Create output directory if it doesn't exist
    mkpath(output_dir)
    # Select n_samples linearly spaced points from (x, y) if needed
    if !isnothing(n_samples) && n_samples < length(x)
        inds = round.(Int, range(1, length(x), length=n_samples))
        x = x[inds]
        y = y[inds]
    end
    if x_rng !== nothing
        x_m = x_rng[1] * x_scl
        x_M = x_rng[2] * x_scl
        inds = (x .>= x_m) .& (x .<= x_M)
        x = x[inds] ./ x_scl
        y = y[inds] ./ y_scl
    else
        x = x ./ x_scl
        y = y ./ y_scl
    end

    if isnothing(xlim)
        xmin = minimum(x)
        xmax = maximum(x)
        xlim = (xmin - 0.1 * abs(xmin), xmax + 0.1 * abs(xmax))
    end

    if isnothing(ylim)
        ymin = minimum(y)
        ymax = maximum(y)
        ylim = (ymin - 0.1 * abs(ymin), ymax + 0.1 * abs(ymax))
    elseif ylim isa Integer
        ylim = abs(ylim)
        if ylim == 0
            throw(ArgumentError("ylim must be non-zero."))
        end
        avg_y = st.mean(y)
        std_y = st.std(y)
        ylim_min = avg_y - ylim * std_y
        ylim_max = avg_y + ylim * std_y
        ylim = (ylim_min, ylim_max)
    end

    if x_ticks === nothing
        x_ticks = :auto
    end

    if save_data
        # Save the (x, y) data as a CSV file for later analysis
        csv_file_path = joinpath(output_dir, file_name * ".csv")
        DelimitedFiles.writedlm(csv_file_path, [x y], ',')
    end

    # Add showaxis=false to remove legend from plot
    plt = Plots.scatter(x, y,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        dpi=300,
        color=:black,
        markerstrokewidth=0,
        markersize=.5,
        titlealign=:left,
        xlims=xlim,
        ylims=ylim,
        label=nothing,  # This removes the legend
        xticks=x_ticks
    )
    Plots.plot!(legend=:outertopright, titlefont=8)
    file_path = joinpath(output_dir, file_name * ".png")
    Plots.savefig(plt, file_path)
end


function save_plot_multiline(x, ys, labels, colors, x_scl, y_scl, title, xlabel, ylabel, output_dir, file_name; x_rng=nothing, ylim=nothing)
    # function save_plot_multiline(x::Vector{<:Real}, ys::Vector{Vector{<:Real}}, labels::Vector{String}, colors::Vector,
    #     x_scl::Real, y_scl::Real, title::String, xlabel::String, ylabel::String,
    #     output_dir::String, file_name::String; x_rng=nothing, ylim=nothing)

    @assert length(ys) == length(labels) == length(colors) "Length mismatch between ys, labels, and colors"

    # Create output directory if needed
    mkpath(output_dir)

    # Trim x and ys if x_rng is given
    if x_rng !== nothing
        x_m = x_rng[1] * x_scl
        x_M = x_rng[2] * x_scl
        inds = (x .>= x_m) .& (x .<= x_M)
        x = x[inds] ./ x_scl
        ys = [y[inds] ./ y_scl for y in ys]
    else
        x = x ./ x_scl
        ys = [y ./ y_scl for y in ys]
    end

    # Compute ylim if not provided
    if isnothing(ylim)
        y_all = reduce(vcat, ys)
        ymin = minimum(y_all)
        ymax = maximum(y_all)
        ylim = (ymin - 0.1 * abs(ymin), ymax + 0.1 * abs(ymin))
    elseif ylim isa Integer
        ylim = abs(ylim)
        if ylim == 0
            throw(ArgumentError("ylim must be non-zero."))
        end
        y_all = reduce(vcat, ys)
        avg_y = st.mean(y_all)
        std_y = st.std(y_all)
        ylim = (avg_y - ylim * std_y, avg_y + ylim * std_y)
    end

    # Plot the first line
    plt = Plots.plot(x, ys[1], label=labels[1], color=colors[1],
        dpi=300, markerstrokewidth=0, markersize=1.2,
        title=title, xlabel=xlabel, ylabel=ylabel,
        titlealign=:left, ylims=ylim)

    # Add the rest
    for i in 2:length(ys)
        Plots.plot!(plt, x, ys[i], label=labels[i], color=colors[i])
    end

    Plots.plot!(plt, legend=:outertopright, titlefont=8)
    file_path = joinpath(output_dir, file_name * ".png")
    Plots.savefig(plt, file_path)
end


using Interpolations: interpolate, Gridded, Linear

function save_heatmap(pdf_thph, ths, phis, title, output_dir, file_name;
    prc=0.98, zero_thresh=nothing, color=:viridis, clim=nothing, save_data=false)
    # Create output directory if it doesn't exist
    mkpath(output_dir)
    if !isnothing(clim)
        clim_min, clim_max = clim
    else
        high_percentile = st.quantile(vec(pdf_thph), prc)
        high_values = pdf_thph[pdf_thph.>=high_percentile]
        avg_high = st.mean(high_values)
        clim_min = 0
        clim_max = avg_high
    end
    # println("Color limit max set to: $(round(clim_max; digits=5))")
    if !isnothing(zero_thresh)
        inds_zr = pdf_thph .< zero_thresh * clim_max
        pdf_thph = copy(pdf_thph)
        pdf_thph[inds_zr] .= 0
    end


    if save_data
        # Save the pdf_thph data as a CSV file for later analysis
        csv_file_path = joinpath(output_dir, file_name * ".csv")
        DelimitedFiles.writedlm(csv_file_path, pdf_thph, ',')
    end

    # Create a heatmap plot
    plt = Plots.heatmap(ths / pi, phis / pi, pdf_thph',
        titlefontsize=10,
        title=title,
        xlabel="θ/π",
        ylabel="ϕ/π",
        dpi=300,
        color=color,  # :viridis,
        # transpose=true,
        # aspect_ratio=:equal,
        titlealign=:left,
        # colorbar=false,
        xlims=(0, 1),
        ylims=(-1, 1),
        clims=(clim_min, clim_max)
    )
    file_path = joinpath(output_dir, file_name * ".png")
    Plots.savefig(plt, file_path)
    # Plots.savefig(plt, replace(file_path, ".png" => ".pdf"))
    # Also save as SVG
    #svg_file_path = joinpath(output_dir, file_name * ".svg")
    #Plots.savefig(plt, svg_file_path)
end


function find_local_maxima(arr::AbstractVector{T}, half_w::Integer=1) where {T<:Real}
    # Validate window size
    if half_w <= 0
        throw(ArgumentError("Half window size half_w must be positive."))
    end

    n = length(arr)
    # Handle empty array case immediately
    if n == 0
        return Int[]
    end

    # Initialize vector to store indices of local maxima
    maxima_indices = Int[]


    # Iterate through each element of the array
    for i in 1:n
        # Determine the start and end indices for the window around index i
        # Ensure indices stay within the bounds of the array
        start_idx = max(1, i - half_w)
        end_idx = min(n, i + half_w)

        # Extract the subarray representing the window
        # Note: Slicing creates a view or copy depending on Julia version and context,
        # but `maximum` works correctly on it.
        window = @view arr[start_idx:end_idx] # Use view for efficiency

        # Check if the current element arr[i] is the maximum value in the window
        # `maximum(window)` handles the case of an empty window implicitly if it could occur,
        # but our bounds logic prevents empty windows for n > 0.
        if arr[i] == maximum(window)
            # If it is, add its index to the list
            push!(maxima_indices, i)
        end
    end

    # Return the list of indices
    return maxima_indices
end

function find_local_minima(arr::AbstractVector{T}, half_w::Integer=1) where {T<:Real}
    if half_w <= 0
        throw(ArgumentError("Half window size half_w must be positive."))
    end

    n = length(arr)
    if n == 0
        return Int[]
    end

    minima_indices = Int[]

    for i in 1:n
        start_idx = max(1, i - half_w)
        end_idx = min(n, i + half_w)

        window = @view arr[start_idx:end_idx]

        if arr[i] == minimum(window)
            push!(minima_indices, i)
        end
    end

    return minima_indices
end

function continued_fraction_convergents(x::Real, N::Int; digits=4)
    a = []  # continued fraction coefficients
    p = Int[]  # numerators
    q = Int[]  # denominators

    x0 = x
    for n in 1:N
        an = floor(Int, x0)
        push!(a, an)
        x0 = 1 / (x0 - an)
    end

    # Compute convergents p_n / q_n
    for n in 1:N
        if n == 1
            push!(p, a[1])
            push!(q, 1)
        elseif n == 2
            push!(p, a[2] * a[1] + 1)
            push!(q, a[2])
        else
            pn = a[n] * p[n-1] + p[n-2]
            qn = a[n] * q[n-1] + q[n-2]
            push!(p, pn)
            push!(q, qn)
        end
    end

    # Filter based on number of matching digits
    threshold = 0.5 * 10.0^(-digits)
    uf, vf = Int[], Int[]
    for i in 1:N
        approx = p[i] / q[i]
        if abs(approx - x) < threshold
            push!(uf, p[i])
            push!(vf, q[i])
        end
    end

    return uf, vf
end


if abspath(PROGRAM_FILE) == @__FILE__

    println("Running utils.jl as main script.")

    x = sort(10.0 .* rand(Float64, 1001) .- 5.0)
    y = @. x^3 - 2 * x^2 + x + 1
    integral = simpson_nonuniform(x, y)
    println("Non-uniform Integration Rel. Error: $(abs(3integral/470 + 1))")

    x = collect(range(-5.0, stop=5.0, length=1001))
    y = @. x^3 - 2 * x^2 + x + 1
    integral_nu = simpson_nonuniform(x, y)
    integral = simpson_integrate(x, y)
    println("Non-uniform Integration Rel. Error: $(abs(3integral_nu/470 + 1))")
    println("Integration Rel. Error: $(abs(3integral/470 + 1))")

    x = sqrt(2)
    u, v = continued_fraction_convergents(x, 15; digits=4)
    println("Continued fraction convergents for sqrt(2):")
    println("u = $u")
    println("v = $v")
end
