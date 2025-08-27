function plotMesh2D(ax, faces::RayTracingMeshOptim; superfaces=false, wallNumbers=[], volumeNumbers=[], color=:lightblue, strokecolor=:black, strokewidth=1)

    function plot_volume(ax, face::PolyFace2D, index=nothing)
        Makie.poly!(ax, Point2f[face.vertices...], color=(color, 0.5), strokecolor=strokecolor, strokewidth=strokewidth)
        if index !== nothing
            Makie.text!(ax, "g$(index)", position=face.midPoint, space = :data)
        end
        # arrows!(ax, Point2f[face.wallMidPoints...], Point2f[face.outwardNormals...], color=:black, lengthscale=0.1)
        # text!(ax, String[[face.solidWalls[i] ? "solid" : "gas" for i in 1:length(face.wallMidPoints)]...], position=face.wallMidPoints, align=(:center, :center), space = :data)
        # text!(ax, "$(face.T_g)", position=face.midPoint, align=(:center, :center), space = :data)
        # text!(ax, String[["$(face.T_w[i])" for i in 1:length(face.wallMidPoints)]...], position=[face.wallMidPoints[i] for i = 1:length(face.wallMidPoints)], align=(:center, :center), space = :data)

    end

    function plot_surface(ax, face::PolyFace2D, index=nothing, wallNumber=nothing)
        # poly!(ax, Point2f[face.vertices...], color=(color, 0.5), strokecolor=strokecolor, strokewidth=strokewidth)
        if index !== nothing && wallNumber !== nothing
            Makie.text!(ax, "w$(index)", position=face.wallMidPoints[wallNumber], space = :data)
        end
        # arrows!(ax, Point2f[face.wallMidPoints...], Point2f[face.outwardNormals...], color=:black, lengthscale=0.1)
        # text!(ax, String[[face.solidWalls[i] ? "solid" : "gas" for i in 1:length(face.wallMidPoints)]...], position=face.wallMidPoints, align=(:center, :center), space = :data)
        # text!(ax, "$(face.T_g)", position=face.midPoint, align=(:center, :center), space = :data)
        # text!(ax, String[["$(face.T_w[i])" for i in 1:length(face.wallMidPoints)]...], position=[face.wallMidPoints[i] for i = 1:length(face.wallMidPoints)], align=(:center, :center), space = :data)

    end

    # plot entire current mesh (including specific numberings)
    volume_count = 0
    surface_count = 0
    for i in eachindex(faces.coarse_mesh)
        if superfaces == true
            plot_volume(ax, faces.coarse_mesh[i], i)
        else
            for face in faces.fine_mesh[i] # face in faces.fine_mesh[i]
                volume_count += 1
                if volumeNumbers != [] && volume_count in volumeNumbers
                    plot_volume(ax, face, volume_count)
                else
                    plot_volume(ax, face)
                end
                for (i, solidWall) in enumerate(face.solidWalls)
                    if solidWall
                        surface_count += 1
                        if wallNumbers != [] && surface_count in wallNumbers
                            plot_surface(ax, face, surface_count, i)
                        end
                    end
                end
            end
        end
    end

    hidespines!(ax)
end