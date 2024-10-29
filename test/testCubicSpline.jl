#=
Created on Wed 15 May 2024
Updated on Fri 17 May 2024

This test program verifies the cubic spline implemented in CubicSplines.jl by
creating a spline for a sine wave, whose first derivative should produce a
cosine wave, and whose second derivative should produce a negative sine wave.
=#

#------------------------------------------------------------------------------

module testCubicSpline

using
    CairoMakie

import
    CubicSplines:
        CubicSpline,
        Y,
        Y′,
        Y″

#  include("./src/CubicSplines.jl")

export
    run

function run(knots::Int)

    # The vector of independent values.
    xₖ = zeros(Float64, knots)

    dx = 0.5π / (knots - 1)

    # Locations and values associated with the knots.
    for i = 1:knots
        xₖ[i] = (i - 1) * dx
    end

    # The vector of dependent values.
    yₖ = zeros(Float64, knots)

    # Values associated with the knots.
    for i = 1:knots
        yₖ[i] = sin(xₖ[i])
    end

    spline = CubicSpline(xₖ, yₖ)

    # Nodes are at the mid-points between neighboring knots.
    nodes = knots - 1
    xₙ  = zeros(Float64, nodes)
    eₙ  = zeros(Float64, nodes)
    e′ₙ = zeros(Float64, nodes)
    e″ₙ = zeros(Float64, nodes)
    yₙ  = zeros(Float64, nodes)
    y′ₙ = zeros(Float64, nodes)
    y″ₙ = zeros(Float64, nodes)
    zₙ  = zeros(Float64, nodes)
    z′ₙ = zeros(Float64, nodes)
    z″ₙ = zeros(Float64, nodes)

    for n = 1:nodes
        xₙ[n]  = (n - 0.5) * dx
        yₙ[n]  = sin(xₙ[n])
        y′ₙ[n] = cos(xₙ[n])
        y″ₙ[n] = -sin(xₙ[n])
    end

    # Compute spline error at the nodes, i.e., mid points between knots.
    for n = 1:nodes
        zₙ[n]  = Y(spline,  xₙ[n])
        z′ₙ[n] = Y′(spline, xₙ[n])
        z″ₙ[n] = Y″(spline, xₙ[n])
        eₙ[n]  = abs(zₙ[n]  - yₙ[n])
        e′ₙ[n] = abs(z′ₙ[n] - y′ₙ[n])
        e″ₙ[n] = abs(z″ₙ[n] - y″ₙ[n])
    end

    dirPath = string(pwd(), "/figures/")
    if !isdir(dirPath)
        mkdir(dirPath)
    end

    CairoMakie.activate!(type = "png")
    fig = Figure(; size = (1000, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax1 = Axis(fig[1, 1];
        title  = string("Values at Mid Points: ", knots, " knots"),
        xlabel = "x",
        ylabel = "y",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax1, xₙ, zₙ;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "y=sin(x)")
    lines!(ax1, xₙ, z′ₙ;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "y′=cos(x)")
    lines!(ax1, xₙ, z″ₙ;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "y″=-sin(x)")
    axislegend("Locations",
        position = :lb)

    ax2 = Axis(fig[1,2];
        title  = string("Errors at Mid Points: ", knots, " knots"),
        xlabel = "x",
        ylabel = "error",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20,
        yscale = log10)
    lines!(ax2, xₙ, eₙ;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "y=sin(x)")
    lines!(ax2, xₙ, e′ₙ;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "y′=cos(x)")
    lines!(ax2, xₙ, e″ₙ;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "y″=-sin(x)")
    axislegend("Locations",
        position = :ct)

    figName = string("testCubicSpline", knots, "knots.png")
    figPath = string(dirPath, figName)
    save(figPath, fig)
end # run

end # testCubicSpline
