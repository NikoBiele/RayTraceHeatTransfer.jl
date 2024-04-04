using Documenter, Documenter.HTMLWriter
using StaticArrays
using LinearAlgebra
using Plots
using CSV
using DataFrames
include("../src/RayTraceHeatTransfer.jl")
using .RayTraceHeatTransfer

# Ensure the Documenter knows where to find the source files for documentation. hej
docroot = joinpath(@__DIR__, "src")

makedocs(
    sitename="Documentation",
    # Specify the module for which documentation is being generated.
    # This tells Documenter to generate documentation for the RayTraceHeatTransfer module,
    # including its submodules and exported members.
    modules=[RayTraceHeatTransfer],
    # Specify the format and any format-specific options.
    format=HTMLWriter.HTML(sidebar_sitename=false),
    # Ensure that the source directory is correctly specified.
    # This option tells Documenter where to find the markdown files and other documentation sources.
    source=docroot,
    pages = [
        "Home" => "index.md",
        "Getting started" => "gettingstarted.md",
        "Usage" => [
            "Generate and mesh geometry" => "geometry.md",
            "Ray trace the geometry" => "raytracing.md",
            "Calculate heat transfer" => "heattransfer.md",
            "View Factors in 3D" => "viewfactor.md"
        ],
        "Code Validation" => "validation.md",
        "Code" => "code.md"
        ]
)

deploydocs(
    repo = "github.com/NikoBiele/RayTraceHeatTransfer.jl.git",
    target = "site",
    branch = "gh-pages",
    julia = "1.x",
    deploy_config = Dict("GITHUB_TOKEN" => ENV["GITHUB_TOKEN"]),
)