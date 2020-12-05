module Utilities

export insertion_sort, interlace_recursion, interlace_iteration, merge_sort1, merge_sort2

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

function interlace_recursion(A::Array{Int64,1}, B::Array{Int64,1})
    if length(A) == 0
        return B
    elseif length(B) == 0
        return A
    elseif A[1] < B[1]
        return vcat([A[1]], interlace_recursion(A[2:end], B))
    else
        return vcat([B[1]], interlace_recursion(A, B[2:end]))
    end
end


"""
The following function also interlaces two Int64 arrays, but it does not use recursion.
"""
function interlace_iteration(A::Array{Int64, 1}, B::Array{Int64, 1})
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
The following function implements merge sort, using interlace_recursion( ).

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
        sorted = interlace_recursion(merge_sort1(a[begin: m]), merge_sort1(a[(m+1): end]))
    end
    return sorted
end


"""
The following function implements merge sort, using interlace_iteration( ).

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
        sorted = interlace_iteration(merge_sort2(a[begin: m]), merge_sort2(a[(m+1): end]))
    end
    return sorted
end

end
