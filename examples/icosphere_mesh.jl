using LinearAlgebra

"""
    icosphere_mesh(subdivision_level)

Build a triangulated unit sphere by recursively subdividing a regular
icosahedron `subdivision_level` times and projecting new vertices onto
the unit sphere.

Level 0 → 20 triangles, 1 → 80, 2 → 320, 3 → 1280.
"""
function icosphere_mesh(subdivision_level::Int)
    φ = (1 + sqrt(5)) / 2

    ico_points_raw = [
         0.0   1.0    φ;    0.0   1.0   -φ;
         0.0  -1.0    φ;    0.0  -1.0   -φ;
         1.0    φ   0.0;    1.0   -φ   0.0;
        -1.0    φ   0.0;   -1.0   -φ   0.0;
           φ  0.0   1.0;      φ  0.0  -1.0;
          -φ  0.0   1.0;     -φ  0.0  -1.0
    ]

    points = similar(ico_points_raw)
    for i in 1:size(ico_points_raw, 1)
        v = ico_points_raw[i, :]
        points[i, :] = v / norm(v)
    end

    faces = [
         1  3  9;   1  9  5;   1  5  7;   1  7 11;   1 11  3;
         4  2 10;   4 10  6;   4  6  8;   4  8 12;   4 12  2;
         3  6  9;   9  6 10;   9 10  5;   5 10  2;   5  2  7;
         7  2 12;   7 12 11;  11 12  8;  11  8  3;   3  8  6
    ]

    for _ in 1:subdivision_level
        midpoint_cache = Dict{Tuple{Int,Int},Int}()
        new_points = [points[i, :] for i in 1:size(points, 1)]

        function get_midpoint(i::Int, j::Int)
            key = (min(i, j), max(i, j))
            haskey(midpoint_cache, key) && return midpoint_cache[key]
            m = (points[i, :] + points[j, :]) / 2
            m = m / norm(m)
            push!(new_points, m)
            midpoint_cache[key] = length(new_points)
            return length(new_points)
        end

        n_faces = size(faces, 1)
        new_faces = Matrix{Int}(undef, 4 * n_faces, 3)
        for k in 1:n_faces
            a, b, c = faces[k, 1], faces[k, 2], faces[k, 3]
            ab = get_midpoint(a, b)
            bc = get_midpoint(b, c)
            ca = get_midpoint(c, a)
            new_faces[4k - 3, :] = [a,  ab, ca]
            new_faces[4k - 2, :] = [ab, b,  bc]
            new_faces[4k - 1, :] = [ca, bc, c ]
            new_faces[4k,     :] = [ab, bc, ca]
        end

        points = reduce(vcat, (p' for p in new_points))
        faces  = new_faces
    end

    return points, faces
end