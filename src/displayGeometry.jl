function displayGeometry(mesh::TracingMesh)

    # lookup mesh in struct
    N_subs = mesh.N_subs
    point1 = mesh.point1
    point2 = mesh.point2
    point3 = mesh.point3
    point4 = mesh.point4
    Nx = mesh.Nx
    Ny = mesh.Ny

    # function to plot one cell
    function trapezoid(pointOne, pointTwo, pointThree, pointFour)
        Shape([pointOne[1], pointTwo[1], pointThree[1], pointFour[1]],[pointOne[2], pointTwo[2], pointThree[2], pointFour[2]])
    end

    # display entire mesh
    plot([trapezoid(point1[i,j], point2[i,j],point3[i,j],point4[i,j]) for i in 1:Nx for j in 1:Ny*N_subs],
                legend=false,color=:white,aspect_ratio=1.0,
                xlabel="Position / m",ylabel="Position / m",title="Meshed ray tracing geometry")
    
end