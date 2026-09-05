# multiple dispatch based on domain type

# 2d plotField
function RayTraceHeatTransfer.plotField(domain::RayTracingDomain2D; field::Symbol=:T,
                            include_walls::Bool=false, transparent_interfaces::Bool=false, 
                            minmax::Union{Nothing, Tuple{Float64, Float64}}=nothing,
                            xlims::Union{Nothing, Tuple{Float64, Float64}}=nothing,
                            ylims::Union{Nothing, Tuple{Float64, Float64}}=nothing,
                            title::Union{Nothing, String}=nothing)
    
    # Function to create a shape for plotting
    function shape(vertices)
        Plots.Shape([vert[1] for vert in vertices], [vert[2] for vert in vertices])
    end

    # Function to get the field value
    function get_field_value(face, field, is_wall=false, index=nothing)
        gas_field = Symbol(string(field) * "_g")
        wall_field = Symbol(string(field) * "_w")
        if is_wall
            return index !== nothing ? getfield(face, wall_field)[index] : NaN
        else
            return getfield(face, gas_field)
        end
    end

    # Collect all field values to determine color range
    all_values = Vector()
    for i in 1:length(domain.coarse_mesh)
        for subvolume in domain.fine_mesh[i] # superface.subvolumes
            push!(all_values, get_field_value(subvolume, field))
            if include_walls
                append!(all_values, filter(!isnan, [get_field_value(subvolume, field, true, i) for i in 1:length(subvolume.solidWalls) if subvolume.solidWalls[i]]))
            end
        end
    end
    vmin, vmax = extrema(all_values)

    # Initialize the plot
    p = Plots.plot(legend=false)

    # Plot volumes (gas properties)
    for i in 1:length(domain.coarse_mesh)
        for subvolume in domain.fine_mesh[i] # face.subvolumes
            color_value = get_field_value(subvolume, field)
            line_color = transparent_interfaces ? :transparent : :black
            Plots.plot!(shape(subvolume.vertices), fill_z=color_value, colorbar=false, linecolor=line_color)
        end
    end

    # Plot surfaces (wall properties) if included
    if include_walls
        for i in 1:length(domain.coarse_mesh)
            for subvolume in domain.fine_mesh[i] #face.subvolume
                for (i, is_solid) in enumerate(subvolume.solidWalls)
                    if is_solid
                        v1 = subvolume.vertices[i]
                        v2 = subvolume.vertices[mod1(i+1, length(subvolume.vertices))]
                        midpoint = (v1 + v2) / 2
                        normal = subvolume.inwardNormals[i]
                        
                        # Create a small rectangle for the wall
                        wall_width = 0.1 * norm(v2 - v1)
                        p1 = midpoint - normal * wall_width / 2
                        p2 = midpoint + normal * wall_width / 2
                        p3 = p2 + (v2 - v1) / 2
                        p4 = p1 + (v2 - v1) / 2
                        
                        color_value = get_field_value(subvolume, field, true, i)
                        Plots.plot!(p, shape([p1, p2, p3, p4]), fill_z=color_value, colorbar=false)
                    end
                end
            end
        end
    end

    # Add colorbar
    if minmax !== nothing
        vmin, vmax = minmax
    end
    Plots.plot!(p, colorbar=true, clims=(vmin, vmax), color=:thermal, dpi=1000)

    Plots.plot!(p, xlabel="Position (m)", 
                    ylabel="Position (m)",
                    guidefontsize=14,
                    tickfontsize=14,
                    colorbar_titlefontsize=14,
                    colorbar_tickfontsize=14,
                    right_margin=15Plots.mm,
                    xtickfont=font(14, "Computer Modern"),
      ytickfont=font(14, "Computer Modern"),
      colorbar_tickfont=font(14, "Computer Modern"),
      colorbar_title_margin=10)

    if xlims !== nothing
        Plots.xlims!(p, xlims)
    end
    if ylims !== nothing
        Plots.ylims!(p, ylims)
    end
    if title !== nothing
        Plots.title!(p, title)
    else
        Plots.title!(p,"Distribution of $field")
    end

    display(p)
    return p
end