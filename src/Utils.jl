
"""
   gen_exp_vec(n, d, order; vars_reversed=false)

Returns all nonnegative integer lists of length n who entires sum to d

These are the exponent vectors for all the homogeneous monomials of
degree d, in n variables.

TODO: give this function @memoize. perhaps for some big examples the 
storage required to store the result isn't worth it. However, for 
tests where we're running many similar examples of medium size,
I think it'll be nicer.

INPUTS: 
* "n" -- integer
* "d" -- integer
* "order" -- string, monomial ordering, defaulted to lexicographic ordering. Also supports neglex 
"""
function gen_exp_vec(n, d, order::Symbol=:lex; vars_reversed=false)
    result = Vector{Vector{Int64}}(undef, binomial(n+d-1,d))
    for i in 1:binomial(n+d-1,d)
        result[i] = zeros(Int64,n)
    end
    if order == :lex
        for i in 1:n
            dtemp = copy(d)
            k = 0
            while k <= (length(result) - 1)
                if i > 1
                    dtemp = copy(d)
                    for j in 1:n
                        dtemp = dtemp - result[length(result)-k][j]
                    end
                end
                if i == n && dtemp > 0
                    result[length(result)-k][i] = dtemp
                    k = k + 1
                    continue
                end
                if dtemp == 0
                    k = k + 1
                    continue
                end
                if dtemp == 1 && i > 1
                    for j in i:n
                        result[length(result)-(k+j-i)][j] = 1
                    end
                    k = k + n - i + 1
                    continue
                end
                dtemp2 = copy(d - dtemp)
                if i == 1 || dtemp2 == 0
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-1,d-dtemp)
                            result[length(result)-(j+k-1)][i] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-1,d-dtemp)
                        dtemp = dtemp - 1
                    end
                else
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                            result[length(result)-(j+k-1)][i] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                        dtemp = dtemp - 1
                    end
                end
            end
        end
    elseif order == :neglex
        for i in 1:n
            dtemp = copy(d)
            k = 1
            while k <= length(result)
                if i > 1
                    dtemp = copy(d)
                    for j in 1:n
                        dtemp = dtemp - result[k][j]
                    end
                end
                if i == n && dtemp > 0
                    result[k][i] = dtemp
                    k = k + 1
                    continue
                end
                if dtemp == 0
                    k = k + 1
                    continue
                end
                if dtemp == 1 && i > 1
                    for j in i:n
                        result[k+j-i][j] = 1
                    end
                    k = k + n - i + 1
                    continue
                end
                dtemp2 = copy(d - dtemp)
                if i == 1 || dtemp2 == 0
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-1,d-dtemp)
                            result[j+k-1][i] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-1,d-dtemp)
                        dtemp = dtemp - 1
                    end
                else
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                            result[j+k-1][i] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                        dtemp = dtemp - 1
                    end
                end
            end
        end
    elseif order == :invlex
        for i in 1:n
            dtemp = copy(d)
            k = 0
            while k <= (length(result) - 1)
                if i > 1
                    dtemp = copy(d)
                    for j in 1:n
                        dtemp = dtemp - result[length(result)-k][j]
                    end
                end
                if i == n && dtemp > 0
                    result[length(result)-k][n-i+1] = dtemp
                    k = k + 1
                    continue
                end
                if dtemp == 0
                    k = k + 1
                    continue
                end
                if dtemp == 1 && i > 1
                    for j in i:n
                        result[length(result)-(k+j-i)][n-j+1] = 1
                    end
                    k = k + n - i + 1
                    continue
                end
                dtemp2 = copy(d - dtemp)
                if i == 1 || dtemp2 == 0
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-1,d-dtemp)
                            result[length(result)-(j+k-1)][n-i+1] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-1,d-dtemp)
                        dtemp = dtemp - 1
                    end
                else
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                            result[length(result)-(j+k-1)][n-i+1] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                        dtemp = dtemp - 1
                    end
                end
            end
        end
    else
        throw(ArgumentError("Unsupported order '$order'"))
    end
    
    if vars_reversed
        for i in 1:binomial(n+d-1,d)
            reverse!(result[i])
        end
    end
    return result
end

"""
FIXME/DOCUMENTME: exp_vec seems to be a 2d array here?

It will also work if exp_vec is an array of arrays.
"""
function gen_mon(exp_vecs, PR)
    result = zeros(PR,length(exp_vecs))
    R = base_ring(PR)
    B = MPolyBuildCtx(PR)
    for i in axes(exp_vecs,1)
        push_term!(B, one(R), exp_vecs[i])
        result[i] = finish(B)
    end
    result
end

"""
Computes the monomials in n variables, of degree d, in the
variable order order, as Oscar polynonmials.
"""
function compute_monomials(n,d,PR,order=:lex,cache=nothing; vars_reversed=false)
    if n < 0 || d < 0
        return []
    end
    if cache != nothing
        exp_vecs = cache[d]
    else
        exp_vecs = gen_exp_vec(n,d,order)
    end

    if vars_reversed
        reverse!.(exp_vecs)
        res = gen_mon(exp_vecs,PR)
        reverse!.(exp_vecs)
        res
    else
        gen_mon(exp_vecs,PR)
    end
end

"""
    polynomial_to_vector!(v, f, n, order=:lex, cache=nothing; vars_reversed=false)

Convert (homogeneous) polynomial to vector form with specified order. Default is lexicographic.
Update v with this vector. Assumes v has the correct length.

f - an oscar polynomial
n - the number of variables (minus 1???)
order - a symbol which denotes the term order
cache - a GradedExpCache which gives us the variables
vars_reversed - should we reverse the variables?
"""
function polynomial_to_vector!(v, f, n, order=:lex,cache=nothing; vars_reversed=false)

    #TODO: remove vars_reversed and set it based on the cache's property

    #TODO: if f turns out to be zero, we don't know what the degree should be.
    #
    #How best to fix this?
    R = base_ring(parent(f))
    d = total_degree(f)

    if cache != nothing
        mon_vec = cache[d]
        #println(cache[d])
    else
        mon_vec = gen_exp_vec(n,d,order,vars_reversed=vars_reversed)
    end
    res = v
    for i in eachindex(mon_vec)
        if vars_reversed
            reverse!(mon_vec[i])
            res[i] = coeff(f, mon_vec[i])
            reverse!(mon_vec[i])
        else
            res[i] = coeff(f, mon_vec[i])
        end
    end

    res
end

"""
    polynomial_to_vector(f, n, order=:lex, cache=nothing; vars_reversed=false)

Convert (homogeneous) polynomial to vector form with specified order. Default is lexicographic.

f - an oscar polynomial
n - the number of variables (minus 1???)
order - a symbol which denotes the term order
cache - a GradedExpCache which gives us the variables
vars_reversed - should we reverse the variables?
"""
function polynomial_to_vector(f, n, order=:lex,cache=nothing; vars_reversed=false)

    #TODO: remove vars_reversed and set it based on the cache's property

    #TODO: if f turns out to be zero, we don't know what the degree should be.
    #
    #How best to fix this?
    R = base_ring(parent(f))
    d = total_degree(f)

    if cache != nothing
        mon_vec = cache[d]
        #println(cache[d])
    else
        mon_vec = gen_exp_vec(n,d,order,vars_reversed=vars_reversed)
    end
    res = fill(R(0), length(mon_vec))
    for i in eachindex(mon_vec)
        if vars_reversed
            reverse!(mon_vec[i])
            res[i] = coeff(f, mon_vec[i])
            reverse!(mon_vec[i])
        else
            res[i] = coeff(f, mon_vec[i])
        end
    end

    res
end

# 
"""
vector_to_polynomial(vect, n, d, PR, order=:lex)

Convert vector to polynomial with specified order. Default is lexicographic.

vect - vector of coefficients
n - number of variables-1
d - homogeneous degree
PR - polynomial ring to be the parent of the return value
order - a symbol which denotes the term order
"""
function vector_to_polynomial(vect, n, d, PR, order=:lex, vars_reversed=false)
    
    C = MPolyBuildCtx(PR)
    R = base_ring(PR)
    exp_vecs = gen_exp_vec(n+1,d,order)
    @assert length(vect) == length(exp_vecs) "vector has incorrect length for the specified degree"
    for i in eachindex(vect)
        #res += PR(vect[i]) * mon[i]
        v = vect[i]
        if v isa AbstractFloat
            v = Integer(vect[i])
        end

        if vars_reversed
            push_term!(C, R(v), reverse(exp_vecs[i]))
        else 
            push_term!(C, R(v), exp_vecs[i])
        end
    end
    res = finish(C)

    res
end

# Converts vector of homogeneous polynomials to a matrix of their coefficents
function convert_p_to_m(polys, expvecs; vars_reversed=false)
    R = coefficient_ring(parent(polys[1]))
    MS = matrix_space(R, length(polys), length(expvecs))
    result = MS()
    for i in axes(polys,1)
        for j in axes(expvecs,1)
            exp_vec = expvecs[j]

            if vars_reversed
                reverse!(exp_vec)
            end
                
            result[i,j] = coeff(polys[i], exp_vec)

            if vars_reversed
                reverse!(exp_vec)
            end
        end
    end
    return result
    #=
    result = []
    for i in axes(polys,1)
        temp = []
        for j in axes(expvec,1)
            push!(temp, coeff(polys[i], expvec[j]))
        end
        push!(result, transpose(temp))
    end
    vcat(result...)
    =#
end

# Converts Matrix of coefficents to vector of polynomials, each row is one polynomial
function convert_m_to_p(mat, expvec, R, PR)
    result = []
    B = MPolyBuildCtx(PR)
    for i in axes(mat,1)
        for j in axes(expvec,1)
            push_term!(B, mat[i,j], expvec[j])
        end
        push!(result,finish(B))
    end
    result
end

