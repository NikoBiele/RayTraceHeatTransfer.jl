"""
    displayGeometry(geometry::Vector{SubEnclosure})

This function plots one or more SubEnclosures.
It is used for building and validating the geometry.
It automatically shows which walls are solid and empty, and shows the number of each SubEnclosure.
This information is important when defining boundary conditions when solving heat transfer problems.
"""
function displayGeometry(geometry::Vector{SubEnclosure})

    # function to plot one SubEnclosure
    function trapezoid(pointOne, pointTwo, pointThree, pointFour)
        Shape([pointOne[1], pointTwo[1], pointThree[1], pointFour[1]],[pointOne[2], pointTwo[2], pointThree[2], pointFour[2]])
    end

    # plot all SubEnclosures stored in the geometry vector
    Plots.plot()
    for i = 1:length(geometry)
        # display entire mesh
        display(Plots.plot!(trapezoid(geometry[i].point1, geometry[i].point2, geometry[i].point3, geometry[i].point4),
                                legend=false,aspect_ratio=1.0,color=:white,
                                xlabel="Position / m",ylabel="Position / m",title="A collection of SubEnclosures"))
        display(annotate!([geometry[i].point1[1].+geometry[i].point2[1]]./2, [geometry[i].point1[2].+geometry[i].point2[2]]./2, Plots.text(geometry[i].solidWall1 ? "Solid" : "Empty", geometry[i].solidWall1 ? :red : :blue, :mid, 20)))
        display(annotate!([geometry[i].point2[1]+geometry[i].point3[1]]/2, [geometry[i].point2[2]+geometry[i].point3[2]]/2, Plots.text(geometry[i].solidWall2 ? "Solid" : "Empty", geometry[i].solidWall2 ? :red : :blue, :mid, 20)))
        display(annotate!([geometry[i].point3[1]+geometry[i].point4[1]]/2, [geometry[i].point3[2]+geometry[i].point4[2]]/2, Plots.text(geometry[i].solidWall3 ? "Solid" : "Empty", geometry[i].solidWall3 ? :red : :blue, :mid, 20)))
        display(annotate!([geometry[i].point4[1]+geometry[i].point1[1]]/2, [geometry[i].point4[2]+geometry[i].point1[2]]/2, Plots.text(geometry[i].solidWall4 ? "Solid" : "Empty", geometry[i].solidWall4 ? :red : :blue, :mid, 20)))
        display(annotate!([geometry[i].midpoint[1]], [geometry[i].midpoint[2]], Plots.text("Sub #$i", :green, :mid, 20)))
    end
end