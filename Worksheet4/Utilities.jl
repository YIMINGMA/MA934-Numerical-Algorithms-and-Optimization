module Utilities
export plotinequalities
export Dantzig, segment, rearrange, find_basic, swap_back, next_vertex
export gradient_descent, golden_section_search
export gradient_descent_auto_diff
export stochastic_gradient_descent

using Plots
using LaTeXStrings
using LinearAlgebra
using DualNumbers
using Random
Random.seed!(12345)
pyplot()


"""
This function plots the feasible set of a linear programming problem.
Input:
    - fs: A function that defines a line.
    - x_range: The range of x axis. x_range should be a tuple (x_min, x_max)
    - inequalities: If f ≥ 0, then inequalties = "≥".
    - labels: The label of the line.
    - colors: The color of the line
Output:
    - p: The plot of the feasible region.
"""
function plotinequalities(fs, x_range, inequalities, labels, colors)
    x_min = x_range[1]
    x_max = x_range[2]
    x = x_min: 0.1: x_max
    p = plot(xlim = (x_min, x_max), ylim = (x_min, x_max))
    for i ∈ 1: length(fs)
        y = fs[i].(x)
        if inequalities[i] == "≤"
            plot!(x, y, fill=(0, 0.3, colors[i]), label = labels[i]*" < 0")
        elseif inequalities[i] == "≥"
            plot!(x, y, fill=(x_max, 0.3, colors[i]), label = labels[i]*" > 0")
        else
            error("Error")
        end
    end
    return p
end


function swap_back(x, order)
    x_copy = x[:]
    order_copy = order[:]
    for i in 1: (length(order_copy)-1)
        for j in (i+1): length(order_copy)
            if order_copy[i] > order_copy[j]
                order_copy[i], order_copy[j] = order_copy[j], order_copy[i]
                x_copy[i], x_copy[j] = x_copy[j], x_copy[i]
            end
        end
    end
    return x_copy
end


function rearrange(A, c, x, order)
    for p in 1: (length(x) - 1)
        if x[p] == 0
            for q in (p+1): length(x)
                if x[q] != 0
                    A[:, p], A[:, q] = A[:, q], A[:, p]
                    c[p], c[q] = c[q], c[p]
                    x[p], x[q] = x[q], x[p]
                    order[p], order[q] = order[q], order[p]
                    break
                end
            end
        end
    end

    return (A, c, x, order)
end


function segment(A, c, x)
    position = 0
    for j in 1: length(x)
        if x[j] == 0
            position = j - 1
            break
        end
    end
    A_B = A[:, begin: position]
    A_N = A[:, (position+1): end]
    c_B = c[begin: position]
    c_N = c[(position+1): end]
    x_B = x[begin: position]
    x_N = x[(position+1): end]

    return (A_B, A_N, c_B, c_N, x_B, x_N)
end


function next_vertex(A_B, A_N, c_B, c_N, x_B, x_N)
    found_nonbasic = false
    Δf = []
    for k in 1: length(x_N)
        found_basic, i, x_k = find_basic(A_B, A_N, x_B, x_N, k)
        if found_basic
            c_k = c_N[k]
            Δf_k = x_k * (c_k - dot(c_B, A_B\A_N[:, k]))
            push!(Δf, (Δf_k, k, i, x_k, found_basic))
        end
    end

    min_Δf = 0
    nonbasic_ind = 0
    basic_ind = 0
    step = 0
    found_basic = false
    for ele in Δf
        if ele[1] < min_Δf
            min_Δf = ele[1]
            nonbasic_ind = ele[2]
            basic_ind = ele[3]
            step = ele[4]
            found_basic = ele[5]
            found_nonbasic = true
        end
    end

    found = found_nonbasic && found_basic
    return (found, nonbasic_ind, basic_ind, step)
end


function find_basic(A_B, A_N, x_B, x_N, k)
    w = A_B \ A_N[:, k]
    ratio_ind_pairs = []
    i = 0
    for j in 1: length(w)
        if w[j] > 0
            ratio = x_B[j] / w[j]
            push!(ratio_ind_pairs, (ratio, j))
        end
    end

    min_ratio = Inf
    found_basic = false
    for pair in ratio_ind_pairs
        if pair[1] < min_ratio
            min_ratio = pair[1]
            i = pair[2]
            found_basic = true
        end
    end
    return (found_basic, i, min_ratio)
end


function Dantzig(A, b, c, x)
    A_copy = A[:, :]
    c_copy = c[:]
    x_copy = x[:]
    order = zeros(length(x))
    for i in 1: length(x)
        order[i] = i
    end

    found = true
    nonbasic_ind = 1
    basic_ind = 1
    step = 0

    A, c, x, order = rearrange(A_copy, c_copy, x_copy, order)
    A_B, A_N, c_B, c_N, x_B, x_N = segment(A, c, x)

    basic_feasible_vec = []

    while found
        x_B_prime = x_B - step * (A_B\A_N[:, nonbasic_ind])
        x_N_prime = x_N
        x_N_prime[nonbasic_ind] = x_N_prime[nonbasic_ind] + step

        x_copy = vcat(x_B_prime, x_N_prime)

        x_ordered = swap_back(x_copy, order)
        push!(basic_feasible_vec, x_ordered)

        A, c, x, order = rearrange(A_copy, c_copy, x_copy, order)
        A_B, A_N, c_B, c_N, x_B, x_N = segment(A, c, x)
        found, nonbasic_ind, basic_ind, step = next_vertex(A_B, A_N, c_B, c_N, x_B, x_N)
    end

    return basic_feasible_vec
end


"""
This function uses the gradient descent algorithm to solve a optimization problem.
Notice that in this function, the step size is NOT fixed! Instead, a line minimi-
zation is performed to find it.
Input:
    - f: The objective function.
    - ∇f: The gradient of the objective function.
    - step: The step size of updating xₙ.
    - init: The initial value of xₙ, namely, x₀.
    - record: Whether to record the values xₙ during the algorithm.
    - tol: Tolerance. If |∇f(xₙ)| ≤ tol, then the algorithm should stop to avoid
           large marchine errors.
    - max_iter: The maximum number of iterations desired.
Output:
    - x⃰: The solution of the optimization problem.
    - x_vals: The values xₙ took during the algorithm
"""
function gradient_descent(f, ∇f, init, record=false, tol=10*eps(Float64), max_iter=10^6,
                          show_counter=false)
    counter = 0
    xₙ = init
    x_vals = []
    while norm(∇f(xₙ), 2) ≥ tol && counter ≤ max_iter
        vₙ = - (∇f(xₙ) ./ norm(∇f(xₙ), 2))
        g(λ) = f(xₙ + λ.*vₙ)
        a, c, b =  -100, 0, 100
        step = golden_section_search(g, (a, c, b), tol, max_iter)
        xₙ = xₙ + step .* vₙ
        counter = counter + 1
        if record
            push!(x_vals, xₙ)
        end
        if show_counter
            println(counter)
        end
    end
    x⃰ = xₙ
    return (x⃰, counter, x_vals)
end


function golden_section_search(f, triple, tol=10*eps(Float64), max_iter=10^6)
    (a, c, b) = triple
    ϕ = 0.5*(sqrt(5.0)-1.0)
    c = b - (b-a)*ϕ

    counter = 0

    while b - a ≥ tol && counter ≤ max_iter
        x = a + b - c
        if c-a > b-c
            if f(a) ≥ f(x) && f(x) ≤ f(c)
                a, c, b = a, x, c
            else
                a, c, b = x, c, b
            end
        else
            if f(c) ≥ f(x) && f(x) ≤ f(b)
                a, c, b = c, x, b
            else
                a, c, b = a, c, x
            end
        end
        counter = counter + 1
    end
    return c
end


"""
This function implements gradient descent with auto differentiation.
Input:
    - f: The objective function.
    - x₀: The initial value of x.
    - tol: The tolerance of the algorithm, set to eps(Float64)*10 by default.
    - max_iter: The maximum number of itereations, set to 10^6 by default.
Output:
    - xₙ: The desired solution of miminization of f.

Comment: Given f: ℝⁿ → ℝ, how to calculate its derivative by auto differentiation?
D = Matrix{Dual}(Dual(0.0, 1.0)I, n, n) # I is the identity matrix from
                                        # LinearAlgebra.jl.
F = [f(x + D[:, k]) for k in 1: n)]
realpart(F[1]) # The value of f(x).
∇f = epsilon.(F) # The value of
"""
function gradient_descent_auto_diff(f, x₀, tol=eps(Float64)*10, max_iter=10^6,
                                    show_counter=false)
    n = length(x₀)
    D = Matrix{Dual}(Dual(0.0, 1.0)I, n, n)
    Fₙ = [f(x₀ + D[:, k]) for k in 1: n]
    ∇fₙ = epsilon.(Fₙ)
    xₙ = x₀
    counter = 0
    while norm(∇fₙ, 2) ≥ tol && counter ≤ max_iter
        vₙ = - (∇fₙ ./ norm(∇fₙ, 2))
        g(λ) = f(xₙ + λ.*vₙ)
        a, c, b =  -100, 0, 100
        step = golden_section_search(g, (a, c, b), tol, max_iter)
        xₙ = xₙ + step .* vₙ
        Fₙ = [f(xₙ + D[:, k]) for k in 1: n]
        ∇fₙ = epsilon.(Fₙ)
        counter = counter + 1
        if show_counter
            println(counter)
        end
    end
    return xₙ
end



function stochastic_gradient_descent(∇f, x₀, X, Y, ξ₀, ξ₁, tol, max_iter)
    n = 0 # counter
    xₙ = x₀
    m = rand(1: length(Y))
    Xₙ = X[m, :]
    Yₙ = Y[m]
    while norm(∇f(xₙ, Xₙ, Yₙ), 2) ≥ tol && n ≤ max_iter
        # Choose one sample randomly
        m = rand(1: length(Y))
        Xₙ = X[m, :]
        Yₙ = Y[m]
        vₙ = - (∇f(xₙ, Xₙ, Yₙ) ./ norm(∇f(xₙ, Xₙ, Yₙ), 2))
        ξₙ = ξ₀ / (1 + ξ₁ * n)
        xₙ = xₙ + ξₙ .* vₙ
        n = n + 1
    end
    return xₙ
end
end
