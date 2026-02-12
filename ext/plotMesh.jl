# multiple dispatch based on domain type

# 2d plotMesh
function RayTraceHeatTransfer.plotMesh(ax, domain::RayTracingDomain2D; supervolumes=false,
                wallNumbers=[], volumeNumbers=[], color=:lightblue, strokecolor=:black, strokewidth=1)

    function plot_volume(ax, volume::PolyVolume2D, index=nothing)
        Makie.poly!(ax, Point2f[volume.vertices...], color=(color, 0.5), strokecolor=strokecolor, strokewidth=strokewidth)
        if index !== nothing
            Makie.text!(ax, "g$(index)", position=volume.midPoint, space = :data)
        end
        # arrows!(ax, Point2f[volume.wallMidPoints...], Point2f[volume.outwardNormals...], color=:black, lengthscale=0.1)
        # text!(ax, String[[volume.solidWalls[i] ? "solid" : "gas" for i in 1:length(volume.wallMidPoints)]...], position=volume.wallMidPoints, align=(:center, :center), space = :data)
        # text!(ax, "$(volume.T_g)", position=volume.midPoint, align=(:center, :center), space = :data)
        # text!(ax, String[["$(volume.T_w[i])" for i in 1:length(volume.wallMidPoints)]...], position=[volume.wallMidPoints[i] for i = 1:length(volume.wallMidPoints)], align=(:center, :center), space = :data)

    end

    function plot_surface(ax, volume::PolyVolume2D, index=nothing, wallNumber=nothing)
        # poly!(ax, Point2f[volume.vertices...], color=(color, 0.5), strokecolor=strokecolor, strokewidth=strokewidth)
        if index !== nothing && wallNumber !== nothing
            Makie.text!(ax, "w$(index)", position=volume.wallMidPoints[wallNumber], space = :data)
        end
        # arrows!(ax, Point2f[volume.wallMidPoints...], Point2f[volume.outwardNormals...], color=:black, lengthscale=0.1)
        # text!(ax, String[[volume.solidWalls[i] ? "solid" : "gas" for i in 1:length(volume.wallMidPoints)]...], position=volume.wallMidPoints, align=(:center, :center), space = :data)
        # text!(ax, "$(volume.T_g)", position=volume.midPoint, align=(:center, :center), space = :data)
        # text!(ax, String[["$(volume.T_w[i])" for i in 1:length(volume.wallMidPoints)]...], position=[volume.wallMidPoints[i] for i = 1:length(volume.wallMidPoints)], align=(:center, :center), space = :data)

    end

    # plot entire current mesh (including specific numberings)
    volume_count = 0
    surface_count = 0
    for i in eachindex(domain.coarse_mesh)
        if supervolumes == true
            plot_volume(ax, domain.coarse_mesh[i], i)
        else
            for subvolume in domain.fine_mesh[i] # volume in domain.fine_mesh[i]
                volume_count += 1
                if volumeNumbers != [] && volume_count in volumeNumbers
                    plot_volume(ax, subvolume, volume_count)
                else
                    plot_volume(ax, subvolume)
                end
                for (i, solidWall) in enumerate(subvolume.solidWalls)
                    if solidWall
                        surface_count += 1
                        if wallNumbers != [] && surface_count in wallNumbers
                            plot_surface(ax, subvolume, surface_count, i)
                        end
                    end
                end
            end
        end
    end

    hidespines!(ax)
end

# 3d plotMesh
function RayTraceHeatTransfer.plotMesh(ax1, domain::Union{ViewFactorDomain3D,Vector{PolyFace3D{G}}}; superfaces::Bool=false) where {G}

    # get vector of superfaces
    if typeof(domain) == Vector{PolyFace3D{G}}
        facesVector = domain
    elseif typeof(domain) == ViewFactorDomain3D
        facesVector = domain.facesMesh
    end

    for i in 1:length(facesVector) # vector of superfaces

        # plot only superface if subfaces are empty
        if superfaces && facesVector[i].vertices !== nothing
            # colors of faces
            facecolors = rand(3)
            # Create triangles using the quad logic, but with linear indexing
            tri1_x = [facesVector[i].vertices[1][1], facesVector[i].vertices[2][1], facesVector[i].vertices[3][1]]
            tri1_y = [facesVector[i].vertices[1][2], facesVector[i].vertices[2][2], facesVector[i].vertices[3][2]]
            tri1_z = [facesVector[i].vertices[1][3], facesVector[i].vertices[2][3], facesVector[i].vertices[3][3]]
            Makie.mesh!(ax1, tri1_x, tri1_y, tri1_z, color = RGBf(facecolors...))
            if length(facesVector[i].vertices) == 4
                tri2_x = [facesVector[i].vertices[1][1], facesVector[i].vertices[3][1], facesVector[i].vertices[4][1]]
                tri2_y = [facesVector[i].vertices[1][2], facesVector[i].vertices[3][2], facesVector[i].vertices[4][2]]
                tri2_z = [facesVector[i].vertices[1][3], facesVector[i].vertices[3][3], facesVector[i].vertices[4][3]]
                Makie.mesh!(ax1, tri2_x, tri2_y, tri2_z, color = RGBf(facecolors...))
            end
        else
            for subface in facesVector[i].subFaces               
                facecolors = rand(3) 
                # Create triangles using the quad logic, but with linear indexing
                tri1_x = [subface.vertices[1][1], subface.vertices[2][1], subface.vertices[3][1]]
                tri1_y = [subface.vertices[1][2], subface.vertices[2][2], subface.vertices[3][2]]
                tri1_z = [subface.vertices[1][3], subface.vertices[2][3], subface.vertices[3][3]]
                Makie.mesh!(ax1, tri1_x, tri1_y, tri1_z, color = RGBf(facecolors...))
                if length(subface.vertices) == 4
                    tri2_x = [subface.vertices[1][1], subface.vertices[3][1], subface.vertices[4][1]]
                    tri2_y = [subface.vertices[1][2], subface.vertices[3][2], subface.vertices[4][2]]
                    tri2_z = [subface.vertices[1][3], subface.vertices[3][3], subface.vertices[4][3]]
                    Makie.mesh!(ax1, tri2_x, tri2_y, tri2_z, color = RGBf(facecolors...))
                end
            end
        end
    end
end