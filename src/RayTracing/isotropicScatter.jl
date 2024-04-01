"""
    isotropicScatter()

This function generates a uniformly distributed random sample
of a direction in 3D and projects it onto the 2D plane.
"""
function isotropicScatter()

    # Find isotropic scatter direction (3D spherical projected onto 2D)
    theta = acos(2.0*rand()-1.0) # cone angle
    phi = 2*pi*rand() # circumferential angle (turn the cone angle around the plane)
    # direction vector according to conventions
    r = 1
    xdir = r*sin(theta)*cos(phi)
    # ydir = r*sin(theta)*sin(phi) # this one is projected onto plane (set to zero)
    zdir = r*cos(theta) 
    i1 = SVector(xdir, zdir)
    
    return i1
end