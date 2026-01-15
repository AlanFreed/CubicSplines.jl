#=
Created on Wed 15 May 2024
Updated on Wed 14 Jan 2026   # switched from CairoMakie to Plots

This test program verifies the cubic spline implemented in CubicSplines.jl by
creating a spline for a sine wave, whose first derivative should produce a
cosine wave, and whose second derivative should produce a negative sine wave.
=#

#------------------------------------------------------------------------------

module testCubicSpline

using
    Plots

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
    
    # set the graphics backend to GR
    ENV["QT_QPA_PLATFORM"] = "wayland"
    gr()
    
    p1 = plot(xₙ, [zₙ z′ₙ z″ₙ], linewidth=3, legend=:inside, legendcolumns=3)
    title!(string("y = sin(x); Values at Mid Points: ", knots, " knots"))
    xlabel!("x")
    ylabel!(string("y, y', y", '"'))
    
    p2 = plot(xₙ, [eₙ e′ₙ e″ₙ], linewidth=3)
    plot!(yscale=:log10, minorgrid=true, legend=:top, legendcolumns=3)
    ylims!(1e-10, 1)
    title!(string("Errors at Mid Points: ", knots, " knots"))
    xlabel!("x")
    ylabel!(string("e, e', e", '"'))

    dirPath = string(pwd(), "/figures/")
    if !isdir(dirPath)
        mkdir(dirPath)
    end
    figName = string("testCubicSpline", knots, "knots.png")
    figPath = string(dirPath, figName)
    
    plot(p1, p2, layout=(2,1))
    savefig(figPath)

end # run

end # testCubicSpline
