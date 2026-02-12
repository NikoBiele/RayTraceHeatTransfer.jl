# Quadrilateral constructor with spectral support
function PolyVolume2D{G}(p::SVector{4,Point2{G}}, b::SVector{4,Bool}, 
                         n_spectral_bins::Int=1,
                         kappa_default::G=zero(G), sigma_s_default::G=zero(G)) where {G}
    
    # Calculate geometric properties
    vertices = [p[1], p[2], p[3], p[4]]
    solidWalls = [b[1], b[2], b[3], b[4]]
    midPoint = convert(Point2{G}, (p[1]+p[2]+p[3]+p[4])/4)
    wallMidPoints = [convert(Point2{G}, (p[1]+p[2])/2), 
                     convert(Point2{G}, (p[2]+p[3])/2), 
                     convert(Point2{G}, (p[3]+p[4])/2), 
                     convert(Point2{G}, (p[4]+p[1])/2)]
    
    inwardNormals = [calculateInwardNormal(p[1], p[2], midPoint),
                      calculateInwardNormal(p[2], p[3], midPoint),
                      calculateInwardNormal(p[3], p[4], midPoint),
                      calculateInwardNormal(p[4], p[1], midPoint)]
    
    volume = convert(G, 0.5)*(p[1][1]*(p[2][2]-p[3][2])+p[2][1]*(p[3][2]-p[1][2])+p[3][1]*(p[1][2]-p[2][2])) + 
             convert(G, 0.5)*(p[3][1]*(p[4][2]-p[1][2])+p[4][1]*(p[1][2]-p[3][2])+p[1][1]*(p[3][2]-p[4][2]))
    
    area = [convert(G, norm(p[i]-p[mod(i, 4)+1])) for i = 1:4]
    subVolumes = PolyVolume2D{G}[]
    
    # Initialize extinction properties based on spectral mode
    if n_spectral_bins == 1
        kappa_g = kappa_default
        sigma_s_g = sigma_s_default
        epsilon = zeros(G, 4)
        
        # Grey mode - scalar values for volume properties
        j_g = zero(G)
        g_a_g = zero(G)
        e_g = zero(G)
        r_g = zero(G)
        g_g = zero(G)
        i_g = zero(G)
        q_in_g = zero(G)
        q_g = zero(G)
        T_in_g = zero(G)
        T_g = zero(G)
        
        # Grey mode - scalar values for wall properties  
        j_w = [zero(G), zero(G), zero(G), zero(G)]
        g_a_w = [zero(G), zero(G), zero(G), zero(G)]
        e_w = [zero(G), zero(G), zero(G), zero(G)]
        r_w = [zero(G), zero(G), zero(G), zero(G)]
        g_w = [zero(G), zero(G), zero(G), zero(G)]
        i_w = [zero(G), zero(G), zero(G), zero(G)]
        q_in_w = [zero(G), zero(G), zero(G), zero(G)]
        q_w = [zero(G), zero(G), zero(G), zero(G)]
        T_in_w = [zero(G), zero(G), zero(G), zero(G)]
        T_w = [zero(G), zero(G), zero(G), zero(G)]
    else
        kappa_g = fill(kappa_default, n_spectral_bins)
        sigma_s_g = fill(sigma_s_default, n_spectral_bins)
        epsilon = [zeros(G, n_spectral_bins) for _ in 1:4]
    
        # Spectral mode - vector values for volume properties
        j_g = zeros(G, n_spectral_bins)
        g_a_g = zeros(G, n_spectral_bins)
        e_g = zeros(G, n_spectral_bins)
        r_g = zeros(G, n_spectral_bins)
        g_g = zeros(G, n_spectral_bins)
        i_g = zeros(G, n_spectral_bins)
        q_in_g = zero(G)
        q_g = zero(G)
        T_in_g = zero(G)
        T_g = zero(G)
        
        # Spectral mode - vector values for wall properties
        j_w = [zeros(G, n_spectral_bins) for _ in 1:4]
        g_a_w = [zeros(G, n_spectral_bins) for _ in 1:4]
        e_w = [zeros(G, n_spectral_bins) for _ in 1:4]
        r_w = [zeros(G, n_spectral_bins) for _ in 1:4]
        g_w = [zeros(G, n_spectral_bins) for _ in 1:4]
        i_w = [zeros(G, n_spectral_bins) for _ in 1:4]
        q_in_w = [zero(G) for _ in 1:4]
        q_w = [zero(G) for _ in 1:4]
        T_in_w = [zero(G) for _ in 1:4]
        T_w = [zero(G) for _ in 1:4]
    end
    
    # Call the default constructor with all fields
    return PolyVolume2D{G}(
        vertices, solidWalls, midPoint, wallMidPoints, inwardNormals,
        volume, area, subVolumes,
        epsilon, kappa_g, sigma_s_g,
        j_g, g_a_g, e_g, r_g, g_g, i_g, q_in_g, q_g, T_in_g, T_g,
        j_w, g_a_w, e_w, r_w, g_w, i_w, q_in_w, q_w, T_in_w, T_w
    )
end

# Triangle constructor with spectral support
function PolyVolume2D{G}(p::SVector{3, Point2{G}}, b::SVector{3,Bool}, 
                         n_spectral_bins::Int=1,
                         kappa_default::G=zero(G), sigma_s_default::G=zero(G)) where {G}
    
    # Calculate geometric properties
    vertices = [p[1], p[2], p[3]]
    solidWalls = [b[1], b[2], b[3]]
    midPoint = convert(Point2{G}, (p[1]+p[2]+p[3])./3)
    wallMidPoints = [convert(Point2{G}, (p[1]+p[2])/2), 
                     convert(Point2{G}, (p[2]+p[3])/2), 
                     convert(Point2{G}, (p[3]+p[1])/2)]
    
    inwardNormals = [calculateInwardNormal(p[1], p[2], midPoint),
                      calculateInwardNormal(p[2], p[3], midPoint),
                      calculateInwardNormal(p[3], p[1], midPoint)]
    
    volume = convert(G, 0.5)*(p[1][1]*(p[2][2]-p[3][2])+p[2][1]*(p[3][2]-p[1][2])+p[3][1]*(p[1][2]-p[2][2]))
    
    area = [convert(G, norm(p[i]-p[mod(i, 3)+1])) for i = 1:3]
    subVolumes = PolyVolume2D{G}[]
    
    # Initialize extinction properties based on spectral mode
    if n_spectral_bins == 1
        kappa_g = kappa_default
        sigma_s_g = sigma_s_default
        epsilon = zeros(G, 3)
        
        # Grey mode - scalar values for volume properties
        j_g = zero(G)
        g_a_g = zero(G)
        e_g = zero(G)
        r_g = zero(G)
        g_g = zero(G)
        i_g = zero(G)
        q_in_g = zero(G)
        q_g = zero(G)
        T_in_g = zero(G)
        T_g = zero(G)
        
        # Grey mode - scalar values for wall properties  
        j_w = [zero(G), zero(G), zero(G)]
        g_a_w = [zero(G), zero(G), zero(G)]
        e_w = [zero(G), zero(G), zero(G)]
        r_w = [zero(G), zero(G), zero(G)]
        g_w = [zero(G), zero(G), zero(G)]
        i_w = [zero(G), zero(G), zero(G)]
        q_in_w = [zero(G), zero(G), zero(G)]
        q_w = [zero(G), zero(G), zero(G)]
        T_in_w = [zero(G), zero(G), zero(G)]
        T_w = [zero(G), zero(G), zero(G)]
    else
        kappa_g = fill(kappa_default, n_spectral_bins)
        sigma_s_g = fill(sigma_s_default, n_spectral_bins)
        epsilon = [zeros(G, n_spectral_bins) for _ in 1:3]
    
        # Spectral mode - vector values for volume properties
        j_g = zeros(G, n_spectral_bins)
        g_a_g = zeros(G, n_spectral_bins)
        e_g = zeros(G, n_spectral_bins)
        r_g = zeros(G, n_spectral_bins)
        g_g = zeros(G, n_spectral_bins)
        i_g = zeros(G, n_spectral_bins)
        q_in_g = zero(G)
        q_g = zero(G)
        T_in_g = zero(G)
        T_g = zero(G)
        
        # Spectral mode - vector values for wall properties
        j_w = [zeros(G, n_spectral_bins) for _ in 1:3]
        g_a_w = [zeros(G, n_spectral_bins) for _ in 1:3]
        e_w = [zeros(G, n_spectral_bins) for _ in 1:3]
        r_w = [zeros(G, n_spectral_bins) for _ in 1:3]
        g_w = [zeros(G, n_spectral_bins) for _ in 1:3]
        i_w = [zeros(G, n_spectral_bins) for _ in 1:3]
        q_in_w = [zero(G) for _ in 1:3]
        q_w = [zero(G) for _ in 1:3]
        T_in_w = [zero(G) for _ in 1:3]
        T_w = [zero(G) for _ in 1:3]
    end
    
    # Call the default constructor with all fields
    return PolyVolume2D{G}(
        vertices, solidWalls, midPoint, wallMidPoints, inwardNormals,
        volume, area, subVolumes,
        epsilon, kappa_g, sigma_s_g,
        j_g, g_a_g, e_g, r_g, g_g, i_g, q_in_g, q_g, T_in_g, T_g,
        j_w, g_a_w, e_w, r_w, g_w, i_w, q_in_w, q_w, T_in_w, T_w
    )
end
