function plotTemperatureField(mesh::TracingMesh,Tg::Vector{Float64},Tw=nothing)
    # this function enables contour-like plots of custom geometries

    Tg_matrix = Array{Float64}(undef, mesh.Nx, mesh.Ny*mesh.N_subs)
    Tg_count = 0
    for i = 1:mesh.Nx
        for j = 1:mesh.Ny*mesh.N_subs
            Tg_count += 1
            Tg_matrix[i,j] = Tg[Tg_count]
        end
    end

    function trapezoid(point1, point2, point3, point4)
        Shape([point1[1], point2[1], point3[1], point4[1]],[point1[2], point2[2], point3[2], point4[2]])
    end
    
    display(plot([trapezoid(mesh.point1[i,j], mesh.point2[i,j],
                            mesh.point3[i,j], mesh.point4[i,j]) for i in 1:mesh.Nx for j in 1:mesh.Ny*mesh.N_subs],
        legend=false,
        fill_z = permutedims([Tg_matrix[i,j] for i in 1:mesh.Nx for j in 1:mesh.Ny*mesh.N_subs]),
        color = :thermal,
        aspect_ratio = 1.0,
        cbar = true,
        title = "Temperature distribution",
        colorbar_title = "Temperature / K",
        xlabel = "Position / m",
        ylabel = "Position / m"
    ))

    if Tw !== nothing
        # plot the left wall
        display(plot!([trapezoid(mesh.point1[1,j].-(mesh.point1[2,j].-mesh.point1[1,j])*0.5, mesh.point1[1,j],
                            mesh.point4[1,j],mesh.point4[1,j].-(mesh.point4[2,j].-mesh.point4[1,j])*0.5) for j = 1:mesh.Ny*mesh.N_subs],
                            fill_z = permutedims(Tw[1:mesh.N_subs*mesh.Ny]),color = :thermal))
        # plot the right wall
        Nx = mesh.Nx
        display(plot!([trapezoid(mesh.point2[Nx,j],mesh.point2[Nx,j].+(mesh.point2[Nx,j].-mesh.point2[Nx-1,j])*0.5,
                                mesh.point3[Nx,j].+(mesh.point3[Nx,j].-mesh.point3[Nx-1,j])*0.5,mesh.point3[Nx,j]) for j = 1:mesh.Ny*mesh.N_subs],
                                fill_z = permutedims(Tw[mesh.N_subs*mesh.Ny+1:2*mesh.N_subs*mesh.Ny]),color = :thermal))
        # plot the bottom wall
        display(plot!([trapezoid(mesh.point1[i,1].-(mesh.point1[i,2].-mesh.point1[i,1])*0.5,
                                mesh.point2[i,1].-(mesh.point2[i,2].-mesh.point2[i,1])*0.5,
                                mesh.point2[i,1],mesh.point1[i,1]) for i = 1:mesh.Nx],
                                fill_z = permutedims(Tw[2*mesh.N_subs*mesh.Ny+1:2*mesh.N_subs*mesh.Ny+mesh.Nx]),color = :thermal))

        # plot top wall
        Ny = mesh.Ny*mesh.N_subs
        display(plot!([trapezoid(mesh.point4[i,Ny],mesh.point3[i,Ny],
                                mesh.point3[i,Ny]+(mesh.point3[i,Ny]-mesh.point3[i,Ny-1])*0.5,
                                mesh.point4[i,Ny]+(mesh.point4[i,Ny]-mesh.point4[i,Ny-1])*0.5) for i = 1:mesh.Nx],
                                fill_z = permutedims(Tw[2*mesh.N_subs*mesh.Ny+mesh.Nx+1:end]), color = :thermal))
    end

    return Tg_matrix

end