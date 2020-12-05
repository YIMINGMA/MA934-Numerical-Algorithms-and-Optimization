module Starter
# List of dependencies
using Plots
using LaTeXStrings
# List of functions to be exported
export example_plot, new_function
# Definitions of any new types provided

# Function definitions
"""
example_plot(n)

This function evaluates Sin(x) at n points in the interval [0, 2π], creates a
plot and then returns the plot.
"""
function example_plot(n)
    title = "This is an example plot."
    x1 = collect(range(0.0, 2*π, length=n))
    y1 = sin.(x1)

    # Plot the points
    p = plot(x1, y1, seriestype=:scatter,label="Some sample points",
    title=title, xlabel=L"x", ylabel=L"sin(x)", markersize=10, markercolor="red")

    # Now plot the true function
    x2 = collect(range(0.0, 2*π, length=1000))
    y2 = sin.(x2)
    plot!(x2, y2, label="Underlying function.", color="green", linewidth=2)
    return p
end

"""
This function provides scatter log2 plot of f(x) = x^α log(x).

Input: α, n
Output: the log2 plot of f with x = 2^0, x = 2^1,..., x = 2^n
"""
function new_function(α, n)
    title = L"This is a log plot of $f(x) = x^\alpha \log x$"
    x = 0: n
    x1 = [2^i for i in x]
    y1 = x1.^α .* log.(x1)
    y = log.(y1)

    p = plot(x, y, seriestype=:scatter,label=L"$f(x) = x^\alpha \log x$",
    title=title, xlabel=L"$x$", ylabel=L"$f(x)$", markersize=10, markercolor="red")

    return p
end

# End the module definition
end
