


This module originated from a course in numerical methods that the author (Alan Freed) taught at TAMU.

To use this module, you will need to add the following repository to your project:

```
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/CubicSplines.jl")
```

# CubicSplines

Given x and y data, where x represents an ordered ascending array of independent variables, and y represents an array of dependent variables of like dimension, this module finds the coefficients $(a, b, c, d)$ for a cubic interpolation that spline these data. Exported are interpolations for

$$ y(x)  = a + b x + c x^2 + d x^3 $$

$$ y′(x) = b + 2 c x + 3 d x^2 $$

$$ y″(x) = 2 c + 6 d x $$

where $y′ = \mathrm{d}y/\mathrm{d}x$ and $y″ = \mathrm{d}^2 y/\mathrm{d}x^2$.

The four extra equations that associate with a cubic spline were chosen to ensure that the spline goes through its two, end, nodal points with slopes that are appropriate for the three nodes at each end of the spline. This is often called a *clamped spline*.

The type created for working with cubic splines has a data structure of
```
struct CubicSpline
    a::Vector{Float64}  # Vector of constant coefficients.
    b::Vector{Float64}  # Vector of linear coefficients.
    c::Vector{Float64}  # Vector of quadratic coefficients.
    d::Vector{Float64}  # Vector of cubic coefficients.
    x::Vector{Float64}  # Vector containing the knots of interpolation.
end
```
where the vector of independent values `x` is needed internally to retrieve the appropriate polynomial for interpolating a $y(x)$ value that lies between two neighboring knots, say $x_n$ and $x_{n+1}$ with $x_n \le x < x_{n+1}$.

## Constructors

To create a new cubic spline, one will typically call
```
CubicSpline(x::Vector{Float64}, y::Vector{Float64})
```
where vector `x` contains an ascending array of independent values, while vector `y` contains an array of dependent values that is to be interpolated. These two vectors must have the same dimension, which must be 3 or greater in length.

Or one can use the general constructor
```
function CubicSpline(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, d::Vector{Float64}, x::Vector{Float64})
```

## Methods

To interpolate for a value of the data, call
```
function Y(spline::CubicSpline, x::Float64)::Float64
```
e.g., `y = Y(spline, x),` where `y` is an interpolated value derived from a `spline` whose interpolant is located `x,` which must lie somewhere within the span of `spline.x.`

To interpolate for a first derivative of the data, call
```
function Y′(spline::CubicSpline, x::Float64)::Float64
```
e.g., `y′ = Y′(spline, x),` where `y′` is an interpolated value for the derivative of `y` derived from a `spline` whose interpolant is located `x,` which must lie somewhere within the span of `spline.x.`

To interpolate for a second derivative of the data, call
```
function Y″(spline::CubicSpline, x::Float64)::Float64
```
e.g., `y″ = Y″(spline, x),` where `y″` is an interpolated value for the second derivative of `y` derived from a `spline` whose interpolant is located `x,` which must lie somewhere within the span of `spline.x.`

## Test

In the subdirectory `test` is a file `testCubicSpline.jl` that when compiled will export a function
```
function run(knots::Int)
```
where `knots` is to be a value greater than 2. These knots will span half of a sine wave that is to be interpolated. The output will be a figure saved to a directory `figures` under your current working directory. Interpolations are at the midpoints of the knots, with errors of approximation being computed there.

# Version History

## Version 0.1.3

Switched from CairoMakie to Plots for creating graphs.

## Version 0.1.2

Improved the documentation.

## Version 0.1.1

Used second order (instead of first order) forward and backward difference formula to estimate derivatives at the end points of a spline.

## Version 0.1.0

The original version, dated 17 May 2024. It is a port from the author's python code from a class that he taught in numerical methods to mechanical engineers when he was a professor at Texas A&M University (TAMU). This port was made because of issues that the author had with the BSplineKit.jl package that kept reoccuring with new releases of the compliler. It turned out to be more efficient for me to just code what I needed in my own applications.
