module Fibonacci

export fibo_recursion, fibo_iteration, fibo_dictmemo

"""
This function calculates the element of the Fibonacci sequence at n, by using recursion.

Input:
    n: the index of the element.
    F₀: the value of the first Fibonacci element, set to 0 by default.
    F₁: the value of the second Fibonacci element, set to 1 by default.

Output: the value of the element.
"""
function fibo_recursion(n, F₀ = 0, F₁ = 1)
    if n == 0
        return F₀
    elseif n == 1
        return F₁
    else
        return fibo_recursion(n-1, F₀, F₁) + fibo_recursion(n-2, F₀, F₁)
    end
end

"""
This function calculates the element of the Fibonacci sequence at n, by using iteration.

Input:
    n: the index of the element.
    F₀: the value of the first Fibonacci element, set to 0 by default.
    F₁: the value of the second Fibonacci element, set to 1 by default.

Output: the value of the element.
"""
function fibo_iteration(n, F₀ = 0, F₁ = 1)
    fibos = zeros(n+1)
    for i in 0:n
        if i == 0
            fibos[i+1] = F₀
        elseif i == 1
            fibos[i+1] = F₁
        else
            fibos[i+1] = fibos[i] + fibos[i-1]
        end
    end
    return fibos[n+1]
end

"""
This function calculates the element of the Fibonacci sequence at n, by using dictionary memorization.

Input:
    n: the index of the element.
    F₀: the value of the first Fibonacci element, set to 0 by default.
    F₁: the value of the second Fibonacci element, set to 1 by default.

Output: the value of the element.
"""
function fibo_dictmemo(n, F₀ = 0, F₁ = 1)
    fibos = Dict(0 => F₀, 1 => F₁)

    function dict_builder(n)
        if n in keys(fibos)
            return fibos[n]
        else
            result = dict_builder(n-1) + dict_builder(n-2)
            fibos[n] = result
            return result
        end
    end

    dict_builder(n)
    return fibos[n]
end

end
