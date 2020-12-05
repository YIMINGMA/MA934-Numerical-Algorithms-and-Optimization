using Revise
using BenchmarkTools
using LaTeXStrings
using JLD
using Plots
using DataFrames, GLM
using Random
Random.seed!(1234)

##########################################################################################
############################## Task 1 ####################################################
##########################################################################################
println("""
##########################################################################################
############################## Task 1 ####################################################
##########################################################################################
""")


"""
This function takes a list a and an optional parameter reverse,
then return the sorted list according to reverse.

Input:
    - a: An array, whose data should be comparable, such as Char, String, Int64, Float64.
    - reverse: If reverse=true, then sort a in the reverse order (largest -> smallest).
        - The default value for reverse is false.

Output:
    - copy_a: The sorted list of a
"""
function insertion_sort(a, reverse=false)
    copy_a = a[:] # Avoid making any changes to the orginal list a.
    if reverse
        for i in 2: length(copy_a)
            # Extract the i-th element, to compare it with its predecessors.
            ele = copy_a[i]
            idx = i-1 # Find a proper postion for ele.
            while idx >= 1
                if ele < copy_a[idx]
                    break
                else
                    idx = idx - 1
                end
            end
            # After the while loop, ele should be inserted into idx1+1
            idx = idx + 1
            # Move the element between idx1 and (i-1), to (idx1+1) and i.
            for j in i: (-1): (idx+1)
                copy_a[j] = copy_a[j-1]
            end
            # Insert ele.
            copy_a[idx] = ele
        end
    else
        for i in 2: length(copy_a)
            # Extract the i-th element, to compare it with its predecessors.
            ele = copy_a[i]
            idx = i-1 # Find a proper postion for ele.
            while idx >= 1
                if ele > copy_a[idx]
                    break
                else
                    idx = idx - 1
                end
            end
            # After the while loop, ele should be inserted into idx1+1
            idx = idx + 1
            # Move the element between idx1 and (i-1), to (idx1+1) and i.
            for j in i: (-1): (idx+1)
                copy_a[j] = copy_a[j-1]
            end
            # Insert ele.
            copy_a[idx] = ele
        end
    end
    return copy_a
end


println("Try to sort an integer array")
list = [2, 3, 4, 1, 0]
println("Sorting $(list) results in $(insertion_sort(list))")
println("Sorting $(list) reversely results in $(insertion_sort(list, true))")
println("")

println("Try to sort a float array")
list = [π^i for i in list]
println("Sorting $(list) results in $(insertion_sort(list))")
println("Sorting $(list) reversely results in $(insertion_sort(list, true))")
println("")

println("Try to sort a string array")
list = ["hello", "wolrd", "this", "is", "a", "test"]
println("Sorting $(list) results in $(insertion_sort(list))")
println("Sorting $(list) reversely results in $(insertion_sort(list, true))")
println("")

println("Try to sort another string array")
println("The Unicode of 🍌 (banana) is $(codepoint('🍌')).")
println("The Unicode of 🍎 (apple) is $(codepoint('🍎')).")
println("The Unicode of 🍵 (tea) is $(codepoint('🍵')).")
list = ["🍵", "🍌", "🍎"]
println("Sorting $(list) results in $(insertion_sort(list))")
println("Sorting $(list) reversely results in $(insertion_sort(list, true))")
println("")

println("Try to sort a char array")
list = ['g', 'h', 'e', 'f', 'a', 'b', 'c', 'd']
println("Sorting $(list) results in $(insertion_sort(list))")
println("Sorting $(list) reversely results in $(insertion_sort(list, true))")
println("")

println("Try to sort a random-integer array")
list = rand(1:200, 10)
println("Sorting $(list) results in $(insertion_sort(list))")
println("Sorting $(list) reversely results in $(insertion_sort(list, true))")
println("")


##########################################################################################
############################## Task 2 ####################################################
##########################################################################################
println("""
##########################################################################################
############################## Task 2 ####################################################
##########################################################################################
""")


function interlace1(A::Array{Int64,1}, B::Array{Int64,1})
    if length(A) == 0
        return B
    elseif length(B) == 0
        return A
    elseif A[1] < B[1]
        return vcat([A[1]], interlace1(A[2:end], B))
    else
        return vcat([B[1]], interlace1(A, B[2:end]))
    end
end


"""
The following function also interlaces two Int64 arrays, but it does not use recursion.
"""
function interlace2(A::Array{Int64, 1}, B::Array{Int64, 1})
    interlaced = zeros(length(A)+length(B))
    counter_A = 1
    counter_B = 1
    while counter_A <= length(A) || counter_B <= length(B)
        if counter_A <= length(A) && counter_B <= length(B)
            if A[counter_A] < B[counter_B]
                interlaced[counter_A+counter_B-1] = A[counter_A]
                counter_A = counter_A + 1
            else
                interlaced[counter_A+counter_B-1] = B[counter_B]
                counter_B = counter_B + 1
            end
        elseif counter_A > length(A) && counter_B <= length(B)
            interlaced[counter_A+counter_B-1:end] = B[counter_B:end]
            break
        else
            interlaced[counter_A+counter_B-1:end] = A[counter_A:end]
            break
        end
    end
    return Int.(interlaced)
end


"""
The following function implements merge sort, using interlace1( ).

Input:
    a: An Int64 array.
Output:
    sorted: The sorted version of a.
"""
function merge_sort1(a::Array{Int64, 1})
    n = length(a)
    if n == 1
        sorted = a
    else
        m = n ÷ 2
        sorted = interlace1(merge_sort1(a[begin: m]), merge_sort1(a[(m+1): end]))
    end
    return sorted
end


"""
The following function implements merge sort, using interlace2( ).

Input:
    a: An Int64 array.
Output:
    sorted: The sorted version of a.
"""
function merge_sort2(a::Array{Int64, 1})
    n = length(a)
    if n == 1
        sorted = a
    else
        m = n ÷ 2
        sorted = interlace2(merge_sort2(a[begin: m]), merge_sort2(a[(m+1): end]))
    end
    return sorted
end


println("Try to sort an integer array")
list = [2, 3, 4, 1, 0]
println("Use the recursion version of interlace( )")
@btime merge_sort1($list)
println("Use the iteration version of interlace( )")
@btime merge_sort2($list)
sorted1 = merge_sort1(list)
sorted2 = merge_sort2(list)
println("sorted1 == sorted2 ? ", sorted1 == sorted2)
println("Sorting $(list) results in $(sorted1)")
println("")

println("Try to sort a random-integer array")
list = rand(1:200, 10)
println("Use the recursion version of interlace( )")
@btime merge_sort1($list)
println("Use the iteration version of interlace( )")
@btime merge_sort2($list)
sorted1 = merge_sort1(list)
sorted2 = merge_sort2(list)
println("sorted1 == sorted2 ? ", sorted1 == sorted2)
println("Sorting $(list) results in $(sorted1)")
println("")


##########################################################################################
############################## Task 3 ####################################################
##########################################################################################
println("""
##########################################################################################
############################## Task 3 ####################################################
##########################################################################################
""")


# Task 3 has been accomplished in Task 2.
# The functions needed are called interlace2() and merge_sort2()
# so I will copy these two functions and rename them as interlace() and mergesort()


"""
The following function also interlaces two Int64 arrays, but it does not use recursion.
"""
function interlace(A::Array{Int64, 1}, B::Array{Int64, 1})
    interlaced = zeros(length(A)+length(B))
    counter_A = 1
    counter_B = 1
    while counter_A <= length(A) || counter_B <= length(B)
        if counter_A <= length(A) && counter_B <= length(B)
            if A[counter_A] < B[counter_B]
                interlaced[counter_A+counter_B-1] = A[counter_A]
                counter_A = counter_A + 1
            else
                interlaced[counter_A+counter_B-1] = B[counter_B]
                counter_B = counter_B + 1
            end
        elseif counter_A > length(A) && counter_B <= length(B)
            # Copy the rest of B to interlaced
            interlaced[counter_A+counter_B-1:end] = B[counter_B:end]
            break
        else
            # Copy the rest of A to interlaced
            interlaced[counter_A+counter_B-1:end] = A[counter_A:end]
            break
        end
    end
    return Int.(interlaced)
end


"""
The following function implements merge sort, using interlace( ).

Input:
    a: An Int64 array.
Output:
    sorted: The sorted version of a.
"""
function mergesort(a::Array{Int64, 1})
    n = length(a)
    if n == 1
        sorted = a
    else
        m = n ÷ 2
        sorted = interlace(mergesort(a[begin: m]), mergesort(a[(m+1): end]))
    end
    return sorted
end


list = rand(1: 5*(2^15), 2^15)
@btime mergesort($list)


##########################################################################################
############################## Task 4 ####################################################
##########################################################################################
println("""
##########################################################################################
############################## Task 4 ####################################################
##########################################################################################
""")


len = 20
times_insertion = zeros(len)
times_merge = zeros(len)

for i in 1: len
    println(i)
    rand_arr = rand(Int, 2^i)
    bm1 = @benchmark insertion_sort(arr) setup=(arr = copy($rand_arr))
    bm2 = @benchmark mergesort(arr) setup=(arr = copy($rand_arr))
    times_insertion[i] = median(bm1.times)
    times_merge[i] = median(bm2.times)
end

save("benchmarks1.jld", "Insertion Sort", times_insertion, "Merge Sort", times_merge)


##########################################################################################
############################## Task 5 ####################################################
##########################################################################################
println("""
##########################################################################################
############################## Task 5 ####################################################
##########################################################################################
""")


dict = load("benchmarks1.jld")
t_insertion = dict["Insertion Sort"]
t_merge = dict["Merge Sort"]
x = [2^i for i in 1: len]
log_x = log10.(x)
log_t_insertion = log10.(t_insertion)
log_t_merge = log10.(t_merge)

p = plot(log_x, log_t_insertion, seriestype=:scatter,label="Insertion Sort",
    legend=:topleft, title="log - log Benchmark Plot", xlabel=L"$\log_{10}(n)$",
    ylabel= L"$\log_{10} (Benchmarks)$", markersize=5, markercolor="red")
plot!(log_x, log_t_merge, seriestype=:scatter,label="Merge Sort",
    markersize=5, markercolor="blue")

# Now fit the data with Odinary least square
# Fit the data of Insertion Sort first
data1 = DataFrame(X = log_x, Y = log_t_insertion)
ols1 = lm(@formula(Y ~ X), data1)
coef1 = coef(ols1) #coef1[1] is the intercept, coef1[2] is the slope
a1 = coef1[2]
b1 = coef1[1]
f1(x) = a1 * x + b1
plot!(f1, log_x[1], log_x[len], label="Fitted Line of Insertion Sort", color="red")
# Fit the Data of Merge Sort
data2 = DataFrame(X = log_x, Y = log_t_merge)
ols2 = lm(@formula(Y ~ X), data2)
coef2 = coef(ols2) #coef2[1] is the intercept, coef2[2] is the slope
a2 = coef2[2]
b2 = coef2[1]
f2(x) = a2 * x + b2
plot!(f2, log_x[1], log_x[len], label = "Fitted Line of Merge Sort", color="blue")

savefig(p, "benchmarkplot1.pdf")

############################## Conclusion ################################################
# From the plot, one can easily find that when n is small, insertion sort is actually
# faster than merge sort. From my view, I think this is because merge sort uses a recursion,
# which is slower than iteration.
# However, when n is large, the difference of computational complexity becomes more and more
# obvious. The computaional cost of insertion sort is O(n^2), while that of merge sort is
# O(nlog(n)). Thus, it is reasonable to see that the fitted line of merge sort is below
# that of insertion sort.
# Another point to make here is that after taking log transformation, the data of insertion
# sort should be expected to lie on a line ( log(n^2)=2log(n) ) in the log-log plot, but the
# data of merge sort should not ( log(nlog(n))=log(n)+log(log(n)) ). It seems that the blue
# points in the plot fall on a line, but it is actually not. This effect is due to comparatively
# small n, whih makes log(log(n)) almost negligible.
println("""
############################## Conclusion ################################################
From the plot, one can easily find that when n is small, insertion sort is actually
faster than merge sort. From my view, I think this is because merge sort uses a recursion,
which is slower than iteration.

However, when n is large, the difference of computational complexity becomes more and more
obvious. The computaional cost of insertion sort is O(n^2), while that of merge sort is
O(nlog(n)). Thus, it is reasonable to see that the fitted line of merge sort is below
that of insertion sort.

Another point to make here is that after taking log transformation, the data of insertion
sort should be expected to lie on a line ( log(n^2)=2log(n) ) in the log-log plot, but the
data of merge sort should not ( log(nlog(n))=log(n)+log(log(n)) ). It seems that the blue
points in the plot fall on a line, but it is actually not. This effect is due to comparatively
small n, whih makes log(log(n)) almost negligible. Actually, we can print out slopes of the
above two lines.
""")
println("")

println("The slope of Insertion Sort Line is", " ", coef1[2])
println("The intercept of Insertion Sort Line is", " ", coef1[1])
println("The slope of Insertion Sort Line is", " ", coef2[2])
println("The intercept of Insertion Sort Line is", " ", coef2[1])

println("")
println("""
The slope of the fitted line corresponding to insertion sort is about 1.82, which is very near 2. 
The slight difference is probably casued by a relatively small sample space.

As for the slope of the line of merge sort, its slope is greater than 1, due to the existence of 
log(log(𝑛)).
""")

println("")
##########################################################################################
############################## Done! #####################################################
##########################################################################################
println("""
##########################################################################################
############################## Done! #####################################################
##########################################################################################
""")
