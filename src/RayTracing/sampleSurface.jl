"""
    sampleSurface(mesh::RayTracingMesh, subNumber::Int64, wallNumber::Int64,
                            xCountSample::Int64, yCountSample::Int64,
                            sampleLeftRight::Bool, sampleTopBottom::Bool)

This function samples an emission position on a 2D surface (a line) as well as an emission direction,
given the endpoints of the line and position in the mesh.
The function samples from a diffuse (Lambertian) distribution and rotates the sample
from local to global coordinates using a rotation matrix.
"""
function sampleSurface(mesh::RayTracingMesh, subNumber::Int64, wallNumber::Int64,
                            xCountSample::Int64, yCountSample::Int64,
                            sampleLeftRight::Bool, sampleTopBottom::Bool)

    # Unit vectors in global coordinate system.
    xVecGlobal = SVector(1.0, 0.0)
    yVecGlobal = SVector(0.0, 1.0) # Used for rotating from local to global with rotation matrix.

    # vector based
    if sampleLeftRight
        if wallNumber == 2
            firstPoint = mesh.point2[mesh.Nx,yCountSample+(subNumber-1)*mesh.Ny]
            secondPoint = mesh.point3[mesh.Nx,yCountSample+(subNumber-1)*mesh.Ny]
        else # wallNumber == 4
            firstPoint = mesh.point4[1,yCountSample+(subNumber-1)*mesh.Ny]
            secondPoint = mesh.point1[1,yCountSample+(subNumber-1)*mesh.Ny]
        end
    elseif sampleTopBottom
        if wallNumber == 1
            firstPoint = mesh.point1[xCountSample,1+(subNumber-1)*mesh.Ny] # remember to change mesh1 to mesh
            secondPoint = mesh.point2[xCountSample,1+(subNumber-1)*mesh.Ny]
        else # if wallNumber == 3
            firstPoint = mesh.point3[xCountSample,mesh.Ny+(subNumber-1)*mesh.Ny]
            secondPoint = mesh.point4[xCountSample,mesh.Ny+(subNumber-1)*mesh.Ny]
        end
    end

    # sample position to find emission point
    point = firstPoint + (secondPoint - firstPoint) * rand()
    # nudge the point a tiny bit towards the midpoint, to ensure we are inside cell
    point = point + (mesh.subMidpoints[subNumber] - point) * 1e-9

    # sample direction of ray (local coordinate system)
    i1_loc = lambertSample3D()
    # unit vectors in local coordinate system
    xVecLocal = normalize(secondPoint-firstPoint)
    yVecLocal = SVector{2}([-xVecLocal[2], xVecLocal[1]])
    # rotate according to rotation matrix
    RotationMatrix = SMatrix{2,2}([dot(xVecGlobal, xVecLocal) dot(xVecGlobal, yVecLocal);
                                    dot(yVecGlobal, xVecLocal) dot(yVecGlobal, yVecLocal)])
    i1 = SVector{2}(RotationMatrix*i1_loc) # emission direction (global coordinate system)

    return point, i1
    
end