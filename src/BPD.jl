##
using LightGraphs
using DataStructures
using Random
using GraphPlot
using Statistics
using StatsBase
##

function initialize_message(g, active=fill(true,nv(g)), deg=degree(g))
    mes = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    for i in vertices(g)
        if active[i]
            for j in neighbors(g,i)
                if active[j]
                    mes[(i,j)] = (1/(deg[i]+2),1/(deg[i]+2))
                end
            end
        end
    end
    mes
end

function update_message!(mes, wa, wb, wc, w0, i, g, active, maxdiff, dampingfactor)
    for (index, j) in enumerate(neighbors(g,i))
        if active[j]
            wa[index] = 1.0
            wb[index] = 0.0
            wc[index] = w0
        end
    end
    for j in neighbors(g,i)
        if active[j]
            q0, qroot = mes[(j,i)]
            q0root = q0 + qroot
            for (index, k) in enumerate(neighbors(g,i))
                if active[k] && k != j
                    wb[index] = wa[index]*(1-q0) + wb[index]*q0root
                    wa[index] *= q0root
                    maxval = max(wa[index], wb[index], wc[index])
                    wa[index] /= maxval
                    wb[index] /= maxval
                    wc[index] /= maxval
                end
            end
        end
    end
    for (index, j) in enumerate(neighbors(g,i))
        if active[j]
            norm = wa[index] + wb[index] + wc[index]
            mes_diff = (wc[index]/norm, wa[index]/norm) .- mes[(i,j)]
            maxdiff = max(sum(abs.(mes_diff)), maxdiff)
            mes[(i,j)] = mes[(i,j)] .+ mes_diff .* dampingfactor
        end
    end
    return maxdiff
end

function belief_propagation!(mes, wa, wb, wc, w0, g, active, ϵ, dampingfactor, T)
    for t in 1:T
        maxdiff = 0
        for i in shuffle(vertices(g))
            if active[i]
                maxdiff = update_message!(mes, wa, wb, wc, w0, i, g, active, maxdiff, dampingfactor)
            end
        end
        maxdiff < ϵ && break
    end
end

function empty_probability!(q, mes, w0, g, active)
    for i in vertices(g)
        if active[i]
            wa = 1.0
            wb = 0.0
            wc = w0
            for j in neighbors(g,i)
                if active[j]
                    q0, qroot = mes[(j,i)]
                    q0root = q0 + qroot
                    wb = wa * (1-q0) + wb * q0root
                    wa *= q0root
                    maxval = max(wa, wb, wc)
                    wa /= maxval
                    wb /= maxval
                    wc /= maxval
                end
            end
            norm = wa + wb + wc
            q[i] = wc/norm
        end
    end
end

function simplify!(leafset, deg, active, q, g)
    k, i = top_with_handle(deg)
    while k < 2
        active[i] = false
        push!(leafset, i)
        deg[i] = nv(g)
        q[i] = -1.0
        for j in neighbors(g,i)
            if active[j]
                deg[j] -= 1
            end
        end
        k, i = top_with_handle(deg)
    end
end

function BPD(g; fix_fraction=0.01, T=500, x=10, damping_factor=0.95)
    n = nv(g)
    w0 = exp(-x)
    deg = MutableBinaryMinHeap(degree(g))
    q = MutableBinaryMaxHeap(fill(-1.0,n))
    active = fill(true, n)
    fvs = Int[]
    leaf = Int[]
    maxdeg = maximum(degree(g))
    wa = zeros(maxdeg)
    wb = zeros(maxdeg)
    wc = zeros(maxdeg)

    simplify!(leaf, deg, active, q, g)
    mes = initialize_message(g, active, deg)
    belief_propagation!(mes, wa, wb, wc, w0, g, active, 1.0e-10, damping_factor, T)
    while length(leaf) + length(fvs) < n
        simplify!(leaf, deg, active, q, g)
        belief_propagation!(mes, wa, wb, wc, w0, g, active, 1.0e-7, damping_factor, 10)
        empty_probability!(q, mes, w0, g, active)

        max_fix_num = max(1,round(fix_fraction*(n-length(leaf)-length(fvs))))
        cnt = 0
        q0, i = top_with_handle(q)
        while cnt < max_fix_num && q0 >= 0
            cnt += 1
            active[i] = false
            deg[i] = n
            for j in neighbors(g,i)
                if active[j]
                    deg[j] -= 1
                end
            end
            q[i] = -1.0
            push!(fvs, i)
            q0, i = top_with_handle(q)
        end
    end   
    fvs, leaf
end

function BPD1(g; fix_fraction=0.01, T=500, x=10, damping_factor=0.95)
    n = nv(g)
    w0 = exp(-x)
    deg = MutableBinaryMinHeap(degree(g))
    q = MutableBinaryMaxHeap(fill(-1.0,n))
    active = fill(true, n)
    fvs = Int[]
    leaf = Int[]
    maxdeg = maximum(degree(g))
    wa = zeros(maxdeg)
    wb = zeros(maxdeg)
    wc = zeros(maxdeg)

    simplify!(leaf, deg, active, q, g)
    mes = initialize_message(g, active, deg)
    belief_propagation!(mes, wa, wb, wc, w0, g, active, 1.0e-10, damping_factor, T)
    while length(leaf) + length(fvs) < n
        simplify!(leaf, deg, active, q, g)
        empty_probability!(q, mes, w0, g, active)

        max_fix_num = max(1,round(fix_fraction*(n-length(leaf)-length(fvs))))
        cnt = 0
        q0, i = top_with_handle(q)
        while cnt < max_fix_num && q0 >= 0
            cnt += 1
            active[i] = false
            deg[i] = n
            for j in neighbors(g,i)
                if active[j]
                    deg[j] -= 1
                end
            end
            q[i] = -1.0
            push!(fvs, i)
            q0, i = top_with_handle(q)
        end
    end   
    fvs, leaf
end
##