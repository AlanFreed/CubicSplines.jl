#=
Created on Wed 15 May 2024
Updated on Fri 17 May 2024
Translated from Python (code dated 09/24/2019).
=#

"""
This module was taken from the author's course in numerical methods at TAMU.

Given x and y data, where x represents an ordered ascending array of independent variables and y represents an array of dependent variables of like dimension, this module finds the coefficients (a, b, c, d) for a cubic interpolation that spline these data. Exported are interpolations for

    y(x)  = a + b*x + c*x^2 + d*x^3

    y′(x) = b + 2*c*x + 3*d*x^2

    y″(x) = 2*c + 6*d*x

The four extra equations that associate with a cubic spline are chosen to ensure that the spline goes through its two end points with a slope that is appropriate for the two nodes at each end of the spline. This is often called a clamped spline.

struct CubicSpline
    a::Vector{Float64}  # vector of length n for the constant coefficients
    b::Vector{Float64}  # vector of length n for the linear coefficients
    c::Vector{Float64}  # vector of length n for the quadratic coefficients
    d::Vector{Float64}  # vector of length n for the cubic coefficients
    x::Vector{Float64}  # vector of length n+1 of independent variables

constructors

    spline = CubicSpline(xData, yData)
or
    spline = CubicSpline(a, b, c, d, x)

methods

    y = Y(spline, atX)
        Given
            spline  an instance of type CubicSpline
            atX     independent variable where dependent variable is sought
        returns
            y       interpolated valued for the dependent variable.

    y′ = Y′(spline, atX)
        Given
            spline  an instance of type CubicSpline
            atX     independent variable where dependent variable is sought
        returns
            y′      interpolated valued for the derivative dy(x)/dx.

    y″ = Y″(spline, atX)
        Given
            spline  an instance of type CubicSpline
            atX     independent variable where dependent variable is sought
        returns
            y″      interpolated valued for the derivative d²y(x)/dx².
"""
module CubicSplines

using
    BandedMatrices

import
    Printf:
        @sprintf

export
    CubicSpline,
    Y,
    Y′,
    Y″

struct CubicSpline
    a::Vector{Float64}  # Vector of constant coefficients.
    b::Vector{Float64}  # Vector of linear coefficients.
    c::Vector{Float64}  # Vector of quadratic coefficients.
    d::Vector{Float64}  # Vector of cubic coefficients.
    x::Vector{Float64}  # Vector containing the knots of interpolation.

    # constructors

    function CubicSpline(x::Vector{Float64}, y::Vector{Float64})
        # Verify the input.
        if length(x) == length(y)
            n = length(x) - 1
            if n < 3
                msg = "A cubic spline requires at least three pairs of points."
                throw(ErrorException(msg))
            end
        else
            msg = "The supplied vectors x and y must have the same length."
            throw(ErrorException(msg))
        end
        for i in 1:n
            if x[i+1] ≤ x[i]
                msg = "The independent vector x must be an ascending array."
                throw(ErrorException(msg))
            end
        end

        # Create the arrays holding the cubic coefficients.
        a = zeros(Float64, n)   # the constant coefficients
        b = zeros(Float64, n)   # the linear coefficients
        c = zeros(Float64, n)   # the quadratic coefficients
        d = zeros(Float64, n)   # the cubic coefficients

        # Create the solver arrays.
        rhs = zeros(Float64, 4n)
        mtx = BandedMatrix(Zeros(4n,4n), (4,4))  # 4 bands below and 4 above

        # Assign values to this banded matrix. Do this by using index values
        # for assignment to the full matrix, not the banded matrix.

        # Populate the first two rows in this banded matrix.
        mtx[1,1] = 1.0
        mtx[1,2] = x[1]
        mtx[1,3] = x[1]^2
        mtx[1,4] = x[1]^3
        mtx[2,2] = 1.0
        mtx[2,3] = 2x[1]
        mtx[2,4] = 3x[1]^2

        # Populate the repeating blocks in this banded matrix, row by row.
        for i in 1:n-1
            k = 4i + 1
            mtx[k-2, k-4] = 1.0
            mtx[k-2, k-3] = x[i+1]
            mtx[k-2, k-2] = x[i+1]^2
            mtx[k-2, k-1] = x[i+1]^3
            mtx[k-1, k-2] = -2.0
            mtx[k-1, k-1] = -6x[i+1]
            mtx[k-1, k+2] = 2.0
            mtx[k-1, k+3] = 6x[i+1]
            mtx[k,   k  ] = 1.0
            mtx[k,   k+1] = x[i+1]
            mtx[k,   k+2] = x[i+1]^2
            mtx[k,   k+3] = x[i+1]^3
            mtx[k+1, k-3] = -1.0
            mtx[k+1, k-2] = -2x[i+1]
            mtx[k+1, k-1] = -3x[i+1]^2
            mtx[k+1, k+1] = 1.0
            mtx[k+1, k+2] = 2x[i+1]
            mtx[k+1, k+3] = 3x[i+1]^2
        end

        # Populate the last two rows in this banded matrix.
        k = 4n
        mtx[k-1, k-3] = 1.0
        mtx[k-1, k-2] = x[n+1]
        mtx[k-1, k-1] = x[n+1]^2
        mtx[k-1, k  ] = x[n+1]^3
        mtx[k,   k-2] = 1.0
        mtx[k,   k-1] = 2x[n+1]
        mtx[k,   k  ] = 3x[n+1]^2

        # Assign values to the right-hand side vector.
        # Populate the first two rows.
        rhs[1] = y[1]
        rhs[2] = (-3y[1] + 4y[2] - y[3]) / (-3x[1] + 4x[2] - x[3])

        # Populate the repeating blocks.
        for i in 1:n-1
            rhs[4i-1] = y[i+1]
            rhs[4i+1] = y[i+1]
        end

        # Populate the last two rows.
        rhs[4n-1] = y[n+1]
        rhs[4n]   = (y[n-1] - 4y[n] + 3y[n+1]) / (x[n-1] - 4x[n] + 3x[n+1])

        # Solve the linear system of equations mtx*lhs = rhs for lhs.

        lhs = mtx \ rhs
        for i in 1:n
            j = 4i
            a[i] = lhs[j-3]
            b[i] = lhs[j-2]
            c[i] = lhs[j-1]
            d[i] = lhs[j]
        end

        new(a, b, c, d, x)
    end

    function CubicSpline(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, d::Vector{Float64}, x::Vector{Float64})

        new(a, b, c, d, x)
    end
end # CubicSpline

function _toString(x::Real)::String
    return @sprintf "%.4e" x;
end

function _findSegment(spline::CubicSpline, atX::Float64)::Integer
    # Finds the left index/node for that segment which contains datum atX.
    iLeft  = 1
    iRight = length(spline.a) + 1
    if (atX < spline.x[iLeft]) || (atX > spline.x[iRight])
        msg = string("atX = ", _toString(atX), " lies outside the data range of [", _toString(spline.x[iLeft]), ", ", _toString(spline.x[iRight]), "].")
        throw(ErrorException(msg))
    end

    searching = true
    while searching
        if (iRight - iLeft) ≤ 1
            searching = false
        else
            i = (iLeft + iRight) ÷ 2
            if atX < spline.x[i]
                iRight = i
            else
                iLeft = i
            end
        end
    end

    return iLeft
end  # _findSegment

function Y(spline::CubicSpline, atX::Float64)::Float64
    i = _findSegment(spline, atX)
    y = spline.a[i] + spline.b[i]*atX + spline.c[i]*atX^2 + spline.d[i]*atX^3
    return y
end  # Y

function Y′(spline::CubicSpline, atX::Float64)::Float64
    i  = _findSegment(spline, atX)
    y′ = spline.b[i] + 2spline.c[i]*atX + 3spline.d[i]*atX^2
    return y′
end  # Y′

function Y″(spline::CubicSpline, atX::Float64)::Float64
    i  = _findSegment(spline, atX)
    y″ = 2spline.c[i] + 6spline.d[i]*atX
    return y″
end  # Y″

end # module CubicSplines
