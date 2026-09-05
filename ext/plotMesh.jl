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