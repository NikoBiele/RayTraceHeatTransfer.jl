function plotTrapezoids(Nx_fine,Ny_fine,N_subs,Tg,point1,point2,point3,point4)
    # this function enables contour-like plots of custom geometries

    Tg_matrix = Array{Float64}(undef, Nx_fine, Ny_fine*N_subs)
    Tg_count = 0
    for i = 1:Nx_fine
        for j = 1:Ny_fine*N_subs
            Tg_count += 1
            Tg_matrix[i,j] = Tg[Tg_count]
        end
    end

    function trapezoid(point1, point2, point3, point4)
        Shape([point1[1], point2[1], point3[1], point4[1]],[point1[2], point2[2], point3[2], point4[2]])
    end
    
    display(plot([trapezoid(point1[i,j], point2[i,j], point3[i,j], point4[i,j]) for i in 1:Nx_fine for j in 1:Ny_fine*N_subs],
        legend=false,
        fill_z = [Tg_matrix[i,j] for i in 1:Nx_fine for j in 1:Ny_fine*N_subs]',
        color = :thermal,
        aspect_ratio = 1.0,
        title = "Temperature distribution"
    ))

    return Tg_matrix

end