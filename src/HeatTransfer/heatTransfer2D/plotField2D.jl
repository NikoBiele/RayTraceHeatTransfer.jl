function plotField2D(mesh::RayTracingMeshOptim; field::Symbol=:T, include_walls::Bool=false, minmax::Union{Nothing, Tuple{Float64, Float64}}=nothing) # where {G<:AbstractFloat, T} # Vector{PolyFace2D{G,T}}
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
    for i in 1:length(mesh.coarse_mesh)
        for subface in mesh.fine_mesh[i] # superface.subFaces
            push!(all_values, get_field_value(subface, field))
            if include_walls
                append!(all_values, filter(!isnan, [get_field_value(subface, field, true, i) for i in 1:length(subface.solidWalls) if subface.solidWalls[i]]))
            end
        end
    end
    vmin, vmax = extrema(all_values)

    # Initialize the plot
    p = Plots.plot(aspect_ratio=1.0, legend=false)

    # Plot volumes (gas properties)
    for i in 1:length(mesh.coarse_mesh)
        for subface in mesh.fine_mesh[i] # face.subFaces
            color_value = get_field_value(subface, field)
            Plots.plot!(shape(subface.vertices), fill_z=color_value, colorbar=false)
        end
    end

    # Plot surfaces (wall properties) if included
    if include_walls
        for i in 1:length(mesh.coarse_mesh)
            for subface in mesh.fine_mesh[i] #face.subFaces
                for (i, is_solid) in enumerate(subface.solidWalls)
                    if is_solid
                        v1 = subface.vertices[i]
                        v2 = subface.vertices[mod1(i+1, length(subface.vertices))]
                        midpoint = (v1 + v2) / 2
                        normal = subface.outwardNormals[i]
                        
                        # Create a small rectangle for the wall
                        wall_width = 0.1 * norm(v2 - v1)
                        p1 = midpoint - normal * wall_width / 2
                        p2 = midpoint + normal * wall_width / 2
                        p3 = p2 + (v2 - v1) / 2
                        p4 = p1 + (v2 - v1) / 2
                        
                        color_value = get_field_value(subface, field, true, i)
                        Plots.plot!(p, shape([p1, p2, p3, p4]), fill_z=color_value, colorbar=false)
                    end
                end
            end
        end
    end

    # Add colorbar
    if minmax !== nothing
        # vmin, vmax = minmax
        vmin, vmax = (500.0, 900.0)
    end
    Plots.plot!(p, colorbar=true, clims=(vmin, vmax), color=:thermal, dpi=1000) #, colorbar_title="Temperature / K")

    # Set labels and title
    # Plots.xlabel!(p, "position (m)")
    # Plots.ylabel!(p, "position (m)")
    # Plots.title!(p, "Distribution of $field" * (include_walls ? " (Gas and Walls)" : " (Gas)"))
    Plots.plot!(p, xlabel=L"\mathrm{position\;(m)}", 
                    ylabel=L"\mathrm{position\;(m)}",
                    # colorbar_title=L"Temperature\;(K)",
                    guidefontsize=20,
                    tickfontsize=18,
                    colorbar_titlefontsize=20,
                    colorbar_tickfontsize=18,
                    right_margin=15Plots.mm,
                    xtickfont=font(18, "Computer Modern"),
      ytickfont=font(18, "Computer Modern"),
      colorbar_tickfont=font(18, "Computer Modern"),
      colorbar_title_margin=10)    # More space between colorbar and title)

      Plots.xlims!(p, 0.0, 1.0)
      Plots.ylims!(p, 0.0, 1.0)
    #   # Get plot dimensions for annotation positioning
    #     xrange = [0.0, 1.0]  # Your position range
    #     yrange = [0.0, 1.0]  # Your position range
    #     plot_width = 1.0
    #     plot_height = 1.0

# # Add custom colorbar title annotation
# annotate!(p, [(xrange[2] + plot_width*0.6, yrange[1] + plot_height*0.5, 
#                Plots.text(L"\mathrm{Temperature\;(K)}", 20, :center, rotation=90))])

    display(p)
    return p
end