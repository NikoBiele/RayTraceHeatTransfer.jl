# Updated constructor with spectral support
function PolyFace3D(vertices::Vector{Point3{G}}, solidFace::Bool, domain_midpoint::Point3{G}, 
                epsilon::Union{G, Vector{G}}, q_in_w::G, T_in_w::G;
                n_spectral_bins::Int=1) where {G}
    
    midPoint = sum(vertices)/length(vertices)
    
    if length(vertices) == 3
        area = norm(cross(vertices[2] - vertices[1], vertices[3] - vertices[1]))/2
    else
        area = norm(cross(vertices[2] - vertices[1], vertices[4] - vertices[1]))
    end
    
    inwardNormal = calculateInwardNormal(vertices[1], vertices[2], vertices[3], domain_midpoint)
    
    # Initialize spectral properties based on epsilon type
    is_spectral = isa(epsilon, Vector)
    
    if is_spectral
        n_bins = length(epsilon)
        # Spectral mode - radiative quantities are vectors
        j_w = zeros(G, n_bins)
        g_a_w = zeros(G, n_bins)
        e_w = zeros(G, n_bins)
        r_w = zeros(G, n_bins)
        g_w = zeros(G, n_bins)
        i_w = zeros(G, n_bins)
        # Physical quantities remain scalar
        q_w = nothing
        T_w = nothing
    else
        # Grey mode - all scalar values
        j_w = nothing
        g_a_w = nothing
        e_w = nothing
        r_w = nothing
        g_w = nothing
        i_w = nothing
        q_w = nothing
        T_w = nothing
    end
    
    return PolyFace3D(vertices, solidFace, midPoint, inwardNormal, area, nothing, epsilon, 
                  j_w, g_a_w, e_w, r_w, g_w, i_w,
                  q_in_w, q_w, T_in_w, T_w)
end