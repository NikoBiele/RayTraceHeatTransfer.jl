"""
    imagLi_2(mag::Float64, angle::Float64,
                    almostZero::Float64)

Eq.(24) from the paper.
"""
function imagLi_2(mag, angle, almostZero)

    if mag > almostZero
        omega = atan(mag*sin(angle),(1 - mag*cos(angle)))
        imaginaryPart = ( (0.5*Cl(2.0*angle, almostZero)) + (0.5*Cl(2.0*omega, almostZero))
                            - (0.5*Cl(2.0*omega + 2.0*angle, almostZero)) + (log(mag)*omega)  )
    else
        imaginaryPart = mag*sin(angle)
    end

    return Float32(imaginaryPart)
end