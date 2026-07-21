# Utilities to work with contact forces
#
# ------------------------------------------------
# Author: Nico Romeo
# Date: Dec 2025

using LinearAlgebra



function make_neighbour_points(points::Matrix{Float64}, radii::Vector{Float64}; cutoff=1.5, fact=1.0)
    # Create a cell list for efficient neighbor searching.
    #
    # We must distinguish two ways points can be close:
    #   (a) they are adjacent *along the rod* (the tube is thick, so consecutive
    #       points are naturally within a radius of each other) -> NOT a collision
    #   (b) a distant part of the rod has folded back and touches/penetrates -> a
    #       genuine self-contact.
    # Euclidean distance alone cannot tell these apart, so we exclude case (a) by
    # ARCLENGTH separation along the rod: points whose arclength separation is
    # smaller than (rᵢ+rⱼ)*fact are "expected" neighbours and are skipped. Crucially
    # we do NOT gate on a lower Euclidean distance -- otherwise a deeply penetrating
    # segment (dist -> 0) would drop out of the neighbour list and tunnel straight
    # through instead of being pushed apart.
    n = size(points, 1)
    # cumulative arclength position of each point
    arclen = zeros(n)
    for i = 2:n
        arclen[i] = arclen[i-1] + norm(points[i,:] - points[i-1,:])
    end
    neighbour_list = Vector{Vector{Int}}(undef, n)
    for i = 1:n
        neighbour_points = Vector{Int}()
        for j = 1:n
            i == j && continue
            dist = norm(points[i,:] - points[j,:])
            if abs(arclen[i] - arclen[j]) > fact*(radii[i]+radii[j]) # skip along-rod neighbours
                if dist <= cutoff*(radii[i]+radii[j])
                    push!(neighbour_points, j)
                end
            end
        end
        neighbour_list[i] = neighbour_points
    end
    return neighbour_list# vector of vector of indices of points
end



function get_contact(point::Vector{Float64}, points::Matrix{Float64}, radius_a::Float64, radii::Vector{Float64}, neighbours::Vector{Int})
    dist0 = 1e6
    jc = 0
    for j =  neighbours
        point_j = points[j,:]
        dist = norm(point - point_j)
        if dist < dist0
            dist0 = dist
            jc = j
        end
    end
    overlap = max(radius_a + radii[jc] - dist0, 0.)
    orient = normalize(points[jc,:] - point)
    return overlap, orient, jc
end

function compute_contact_force_clamped(points::Matrix{Float64}, neighbour_list::Vector{Vector{Int}}, radii::Vector{Float64}, alpha::Float64)
    f_cont = zeros(size(points,1),3)
    N = size(points,1)-1 # exclude the first and last point as they are clamped here...
    for i = 2:N
        if !isempty(neighbour_list[i])
            overlap, orient, j = get_contact(points[i,:], points, radii[i], radii, neighbour_list[i])
            f_cont[i,:] += -alpha * overlap * orient
            f_cont[j,:] += +alpha * overlap * orient
        end
    end
    f_cont[1,:] .= 0. # make sure clamped
    f_cont[N+1,:] .= 0. # make sure clamped
    return f_cont
end



