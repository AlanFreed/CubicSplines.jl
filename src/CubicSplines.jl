#=
Created on Wed 15 May 2024
Updated on Tue 29 Oct 2024
Translated from Python (code dated 09/24/2019).
=#

"""
# Module CubicSplines

This module was taken from the author's course in numerical methods at TAMU.

Given x and y data, where x represents an ordered ascending array of independent variables and y represents an array of dependent variables of like dimension, this module finds the coefficients (a, b, c, d) for a cubic interpolation that spline these data. Exported are splined interpolations for

    y(x)  = a + bx + cx² + dx³
    y′(x) = b + 2cx + 3dx²
    y″(x) = 2c + 6dx

The four extra equations needed to construct a cubic spline are chosen to ensure that the spline goes through its two, end, nodal points at slopes that are appropriate for the three nodes at each end of the spline. This is often called a *clamped spline*.
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

"""
# CubicSpline

```julia
struct CubicSpline
    a::Vector{Float64}  # vector of length N for the constant coefficients
    b::Vector{Float64}  # vector of length N for the linear coefficients
    c::Vector{Float64}  # vector of length N for the quadratic coefficients
    d::Vector{Float64}  # vector of length N for the cubic coefficients
    x::Vector{Float64}  # vector of length N+1 of independent variables
end
```

## Constructors

```julia
spline = CubicSpline(xdata, ydata)
```
where

    *xdata* is a vector of independent variables *x* with ascending order.
    *ydata* is a vector of dependent variables *y*.
    
or
```julia
spline = CubicSpline(a, b, c, d, x)
```

## Methods

### Interpolate a spline.

```julia
y = Y(spline, x)
```
given

    spline is an instance of type CubicSpline.
    x      is the independent variable whose dependent variable is sought.
    
returns

    y      is an interpolated valued for the dependent variable.

### Interpolate a spline for its first-order derivative.

```julia
y′ = Y′(spline, x)
```
given

    spline is an instance of type CubicSpline.
    x      is the independent variable whose dependent variable is sought.
    
returns

    y′     is an interpolated valued for the derivative dy(x)/dx.

### Interpolate a spline for its second-order derivative.

```julia
y″ = Y″(spline, x)
```
given

    spline is an instance of type CubicSpline.
    x      is the independent variable whose dependent variable is sought.
    
returns
            
    y″     is an interpolated valued for the derivative d²y(x)/dx².
"""
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
            N = length(x) - 1
            if N < 3
                error("A cubic spline requires at least three pairs of points.")
            end
        else
            error("The supplied vectors x and y must have the same length.")
        end
        for n in 1:N
            if x[n+1] ≤ x[n]
                error("The independent vector x must be an ascending array.")
            end
        end

        # Create the arrays holding the cubic coefficients.
        a = zeros(Float64, N)   # the constant coefficients
        b = zeros(Float64, N)   # the linear coefficients
        c = zeros(Float64, N)   # the quadratic coefficients
        d = zeros(Float64, N)   # the cubic coefficients

        # Create the solver arrays.
        rhs = zeros(Float64, 4N)
        mtx = BandedMatrix(Zeros(4N,4N), (4,4))  # 4 bands below and 4 above

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
        for n in 1:N-1
            k = 4n + 1
            mtx[k-2, k-4] = 1.0
            mtx[k-2, k-3] = x[n+1]
            mtx[k-2, k-2] = x[n+1]^2
            mtx[k-2, k-1] = x[n+1]^3
            mtx[k-1, k-2] = -2.0
            mtx[k-1, k-1] = -6x[n+1]
            mtx[k-1, k+2] = 2.0
            mtx[k-1, k+3] = 6x[n+1]
            mtx[k,   k  ] = 1.0
            mtx[k,   k+1] = x[n+1]
            mtx[k,   k+2] = x[n+1]^2
            mtx[k,   k+3] = x[n+1]^3
            mtx[k+1, k-3] = -1.0
            mtx[k+1, k-2] = -2x[n+1]
            mtx[k+1, k-1] = -3x[n+1]^2
            mtx[k+1, k+1] = 1.0
            mtx[k+1, k+2] = 2x[n+1]
            mtx[k+1, k+3] = 3x[n+1]^2
        end

        # Populate the last two rows in this banded matrix.
        k = 4N
        mtx[k-1, k-3] = 1.0
        mtx[k-1, k-2] = x[N+1]
        mtx[k-1, k-1] = x[N+1]^2
        mtx[k-1, k  ] = x[N+1]^3
        mtx[k,   k-2] = 1.0
        mtx[k,   k-1] = 2x[N+1]
        mtx[k,   k  ] = 3x[N+1]^2

        # Assign values to the right-hand side vector.
        # Populate the first two rows.
        rhs[1] = y[1]
        rhs[2] = (-3y[1] + 4y[2] - y[3]) / (-3x[1] + 4x[2] - x[3])

        # Populate the repeating blocks.
        for n in 1:N-1
            rhs[4n-1] = y[n+1]
            rhs[4n+1] = y[n+1]
        end

        # Populate the last two rows.
        rhs[4N-1] = y[N+1]
        rhs[4N]   = (y[N-1] - 4y[N] + 3y[N+1]) / (x[N-1] - 4x[N] + 3x[N+1])

        # Solve the linear system of equations mtx*lhs = rhs for lhs.

        lhs = mtx \ rhs
        for n in 1:N
            j = 4n
            a[n] = lhs[j-3]
            b[n] = lhs[j-2]
            c[n] = lhs[j-1]
            d[n] = lhs[j]
        end

        new(a, b, c, d, x)::CubicSpline
    end

    function CubicSpline(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, d::Vector{Float64}, x::Vector{Float64})

        new(a, b, c, d, x)::CubicSpline
    end
end # CubicSpline

function _toString(x::Real)::String
    return @sprintf "%.4e" x;
end

function _findSegment(spline::CubicSpline, atX::Float64)::Integer
    # Finds the left index/node for that segment which contains datum atX.
    iLeft  = 1
    iRight = length(spline.a) + 1
    if atX < spline.x[iLeft] || atX > spline.x[iRight]
        msg = string("atX = ", _toString(atX))
        msg = string(msg, " lies outside the data range of [")
        msg = string(msg, _toString(spline.x[iLeft]), ", ")
        msg = string(msg, _toString(spline.x[iRight]), "].")
        error(msg)
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

function Y(spline::CubicSpline, x::Float64)::Float64
    i = _findSegment(spline, x)
    y = spline.a[i] + spline.b[i]*x + spline.c[i]*x^2 + spline.d[i]*x^3
    return y
end  # Y

function Y′(spline::CubicSpline, x::Float64)::Float64
    i  = _findSegment(spline, x)
    y′ = spline.b[i] + 2spline.c[i]*x + 3spline.d[i]*x^2
    return y′
end  # Y′

function Y″(spline::CubicSpline, x::Float64)::Float64
    i  = _findSegment(spline, x)
    y″ = 2spline.c[i] + 6spline.d[i]*x
    return y″
end  # Y″

end # module CubicSplines

