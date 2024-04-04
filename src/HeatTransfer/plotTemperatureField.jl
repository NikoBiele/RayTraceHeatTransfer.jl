"""
    plotTemperatureField(mesh::RayTracingMesh,Tg::Vector{Float64},Tw=nothing)

This function enables contour-like plots of custom geometries.
It plots the gas temperature distribution and optionally also wall temperatures.
The two temperature input vectors are returned by the function 'steadyState'.
The wall temperature vector 'Tw' is optional.
"""
function plotTemperatureField(mesh::RayTracingMesh,Tg::Vector{Float64},Tw=nothing)

    # organize the input Tg vector for plotting
    Tg_array = Array{Float64}(undef, mesh.N_subs, mesh.Nx, mesh.Ny)
    Tg_count = 0
    for k = 1:mesh.N_subs
        for i = 1:mesh.Nx
            for j = 1:mesh.Ny
                Tg_count += 1
                Tg_array[k,i,j] = Tg[Tg_count]
            end
        end
    end

    # function for plotting each cell
    function trapezoid(point1, point2, point3, point4)
        Shape([point1[1], point2[1], point3[1], point4[1]],[point1[2], point2[2], point3[2], point4[2]])
    end
    
    plot()
    Tg_count = 1

    for k = 1:mesh.N_subs # loop over each SubEnclosure

        # plot gas temperatures
        display(plot!([trapezoid(mesh.point1[i,j+(k-1)*mesh.Ny], mesh.point2[i,j+(k-1)*mesh.Ny],
                                mesh.point3[i,j+(k-1)*mesh.Ny], mesh.point4[i,j+(k-1)*mesh.Ny]) for i in 1:mesh.Nx for j in 1:mesh.Ny],
            legend=false,
            fill_z = permutedims([Tg_array[k,i,j] for i in 1:mesh.Nx for j in 1:mesh.Ny]),
            color = :thermal,
            aspect_ratio = 1.0,
            cbar = true,
            title = "Temperature distribution",
            colorbar_title = "Temperature / K",
            xlabel = "Position / m",
            ylabel = "Position / m")
        )

        if Tw !== nothing

            # plot wall temperatures

            if mesh.solidWalls[k][1]
                # plot the bottom wall
                display(plot!(
                    [trapezoid(mesh.point1[i,1+(k-1)*mesh.Ny].-(mesh.point1[i,2+(k-1)*mesh.Ny].-mesh.point1[i,1+(k-1)*mesh.Ny])*0.5,
                    mesh.point2[i,1+(k-1)*mesh.Ny].-(mesh.point2[i,2+(k-1)*mesh.Ny].-mesh.point2[i,1+(k-1)*mesh.Ny])*0.5,
                    mesh.point2[i,1+(k-1)*mesh.Ny],
                    mesh.point1[i,1+(k-1)*mesh.Ny])
                        for i = 1:mesh.Nx],
                            fill_z = permutedims(Tw[Tg_count:Tg_count+mesh.Ny-1]),color = :thermal))
                Tg_count += mesh.Nx
            end

            if mesh.solidWalls[k][2]
                # plot the right wall
                Nx = mesh.Nx
                display(plot!(
                    [trapezoid(mesh.point2[Nx,j+(k-1)*mesh.Ny],
                    mesh.point2[Nx,j+(k-1)*mesh.Ny].+(mesh.point2[Nx,j+(k-1)*mesh.Ny].-mesh.point2[Nx-1,j+(k-1)*mesh.Ny])*0.5,
                    mesh.point3[Nx,j+(k-1)*mesh.Ny].+(mesh.point3[Nx,j+(k-1)*mesh.Ny].-mesh.point3[Nx-1,j+(k-1)*mesh.Ny])*0.5,
                    mesh.point3[Nx,j+(k-1)*mesh.Ny])
                        for j = 1:mesh.Ny],
                            fill_z = permutedims(Tw[Tg_count:Tg_count+mesh.Ny-1]),color = :thermal))
                Tg_count += mesh.Nx
            end

            if mesh.solidWalls[k][3]
                # plot top wall
                Ny = mesh.Ny
                display(plot!(
                    [trapezoid(mesh.point4[i,Ny+(k-1)*mesh.Ny],mesh.point3[i,Ny+(k-1)*mesh.Ny],
                    mesh.point3[i,Ny+(k-1)*mesh.Ny]+(mesh.point3[i,Ny+(k-1)*mesh.Ny]-mesh.point3[i,Ny-1+(k-1)*mesh.Ny])*0.5,
                    mesh.point4[i,Ny+(k-1)*mesh.Ny]+(mesh.point4[i,Ny+(k-1)*mesh.Ny]-mesh.point4[i,Ny-1+(k-1)*mesh.Ny])*0.5)
                        for i = 1:mesh.Nx],
                            fill_z = permutedims(Tw[Tg_count:Tg_count+mesh.Ny-1]),color = :thermal))
                Tg_count += mesh.Nx
            end

            if mesh.solidWalls[k][4]
                # plot the left wall
                display(plot!(
                    [trapezoid(mesh.point1[1,j+(k-1)*mesh.Ny].-(mesh.point1[2,j+(k-1)*mesh.Ny].-mesh.point1[1,j+(k-1)*mesh.Ny])*0.5,
                    mesh.point1[1,j+(k-1)*mesh.Ny],
                    mesh.point4[1,j+(k-1)*mesh.Ny],
                    mesh.point4[1,j+(k-1)*mesh.Ny].-(mesh.point4[2,j+(k-1)*mesh.Ny].-mesh.point4[1,j+(k-1)*mesh.Ny])*0.5)
                        for j = 1:mesh.Ny],
                            fill_z = permutedims(Tw[Tg_count:Tg_count+mesh.Ny-1]),color = :thermal))
                Tg_count += mesh.Nx
            end

        end

    end

    return Tg_array

end