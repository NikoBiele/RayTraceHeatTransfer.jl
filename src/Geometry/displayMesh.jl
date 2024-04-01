"""
    displayMesh(mesh::RayTracingMesh)

This function plots the mesh stored in a RayTracingMesh struct.
"""
function displayMesh(mesh::RayTracingMesh)

    # function to plot one cell
    function trapezoid(pointOne, pointTwo, pointThree, pointFour)
        Shape([pointOne[1], pointTwo[1], pointThree[1], pointFour[1]],[pointOne[2], pointTwo[2], pointThree[2], pointFour[2]])
    end

    Nx = mesh.Nx
    Ny = mesh.Ny
    plot()

    for k = 1:mesh.N_subs
        # display entire mesh
        display(Plots.plot!([trapezoid(mesh.point1[i,j+Ny*(k-1)],
                                        mesh.point2[i,j+Ny*(k-1)],
                                        mesh.point3[i,j+Ny*(k-1)],
                                        mesh.point4[i,j+Ny*(k-1)]) for i in 1:Nx for j in 1:Ny],
                    legend=false,color=:white,aspect_ratio=1.0,
                    xlabel="Position / m",ylabel="Position / m",title="Meshed ray tracing geometry"))
    end
end