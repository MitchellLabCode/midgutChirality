# Utilities to work with the boundary defined as an intersection of halfspaces
# Walls are defined by their normal vectors and offsets from the origin
# Implicitly assumes that the walls define a convex Polyhedron, and that the given normals 
# point outwards, i.e. that for all i, dot(n_i, x) <= d_i for all points x inside the polyhedron
# ------------------------------------------------
# Author: Nico Romeo
# Date: Dec 2025

using LinearAlgebra


struct Plane
    normal::Vector{Float64}
    d::Float64
end


function make_walls(normals::AbstractArray{Float64}, ds::AbstractArray{Float64}) # ds really is a vector of floats
    planes = Vector{Plane}()
    for i in 1:size(normals, 1)
        normal = normalize(normals[i, :])
        d = ds[i]
        push!(planes, Plane(normal, d))
    end
    return planes
end

function point_in_polyhedron(point::AbstractArray{Float64}, poly::Vector{Plane})
    for face in poly
        if  face.d-dot(face.normal, point) < -1e-8
            return false
        end
    end
    return true
end

function distance_to_polyhedron(point::Vector{Float64}, poly::Vector{Plane})
    min_dist = Inf
    for face in poly
        dist = face.d - dot(face.normal, point)
        if dist < min_dist
            min_dist = dist
        end
    end
    return min_dist
end

function distance_to_face(point::AbstractArray{Float64}, face::Plane)
    dist = face.d - dot(face.normal, point)
    return dist
end

function distance_to_faces(point::AbstractArray{Float64}, poly::Vector{Plane})
    dists = Vector{Float64}(undef, size(faces, 1))
    for (i, face) in enumerate(poly)
        dists[i] = distance_to_face(point, face)
    end
    return dists
end


function project_point_to_polyhedron(point::AbstractArray{Float64}, poly::Vector{Plane})
    if isempty(poly) || point_in_polyhedron(point, poly)
        return point, 0.0
    end
    min_dist = Inf
    closest_point = point
    for face in poly
        dist = face.d- dot(face.normal, point)
        if dist < min_dist
            min_dist = dist
            closest_point = point - dist * face.normal
        end
    end
    return closest_point, min_dist
end


function get_force_components(point::AbstractArray{Float64}, poly::Vector{Plane})
    #if point_in_polyhedron(point, poly)
    #    return zeros(3), 0.0
    #end
    min_dist = 0.0
    force_orient = zeros(3)
    for face in poly
        dist = face.d- dot(face.normal, point)
        if dist < min_dist
            min_dist = dist
            force_orient = - face.normal
        end
    end
    return force_orient, -min_dist
end

function make_neighbour_walls(points::Matrix{Float64}, walls::Vector{Plane}, cutoff::Float64)
    # Create a cell list for efficient neighbor searching
    cell_size = cutoff  
    neighbour_walls = Vector{Vector{Plane}}(undef, size(points, 1))
    for (i, point) in enumerate(eachrow(points))
        neighbour_faces = Vector{Plane}()
        for face in walls
            dist = distance_to_face(point, face) # positive if inside
            if dist <= cutoff
                push!(neighbour_faces, face)
            end
        end
        neighbour_walls[i] = neighbour_faces
    end
    return neighbour_walls# vector of vector of face
end