# =============== IJK Sync Rad Reaction ================== #

function IFluxFunction(force::SyncRadReact,space_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p,u0,u1,phi0,phi1,name::String)

    # Evaluate the flux through the I surface i.e. surface of constant p

    mode = force.mode
    B = force.B

    # convert p to Float64 to ensure no floating point issues
    p64 = Float64(p)

    μ0 = getfield(DC,Symbol("μ0"))
    m = getfield(DC,Symbol("mu"*name))
    q = getfield(DC,Symbol("q"))
    Z = getfield(DC,Symbol("z"*name))
    mEle = getfield(DC,Symbol("mEle"))
    c = getfield(DC,Symbol("c"))

    fluxScale = (Z^4*B^2)/(μ0*m^2*mEle*c^2) / (4pi^2) # normalised by σT*c, factor of 4pi^2 is to make sure timescale is correct, no idea why this is needed 

    if m == 0 || Z == 0
        return flux = 0f0
    else
        if typeof(mode) == Ani
            flux = 1/(6*m^2)
            flux *= p64*sqrt(m^2 + p64^2) * (-3*u0 + u0^3 + 3*u1 - u1^3) * (phi0 - phi1)
            flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
            flux *= fluxScale
        elseif typeof(mode) == Axi # average over azimuthal angles
            flux = 1/(6*m^2)
            flux *= p64*sqrt(m^2 + p64^2) * (-3*u0 + u0^3 + 3*u1 - u1^3) * (phi0 - phi1)
            flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
            flux *= fluxScale
        elseif typeof(mode) == Iso # averaged force over all angles
            flux = -1/(3*m^2)
            flux *= p64*sqrt(m^2 + p64^2) * (phi0 - phi1) * (u0 - u1)
            flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
            flux *= fluxScale
        else
            error("Synchrotron mode not recognised.")
        end
    end

    return flux

end

function JFluxFunction(force::SyncRadReact,space_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u,phi0,phi1,name::String)

    # Evaluate the flux through the J surface i.e. surface of constant u

    mode = force.mode
    B = force.B

    # convert p to Float64 to ensure no floating point issues
    p064 = Float64(p0)
    p164 = Float64(p1)

    μ0 = getfield(DC,Symbol("μ0"))
    m = getfield(DC,Symbol("mu"*name))
    q = getfield(DC,Symbol("q"))
    Z = getfield(DC,Symbol("z"*name))
    mEle = getfield(DC,Symbol("mEle"))
    c = getfield(DC,Symbol("c"))

    fluxScale = (Z^4*B^2)/(μ0*m^2*mEle*c^2) / (2pi) # normalised by σT*c, factor of 4pi^2 is to make sure timescale is correct, no idea why this is needed  

    if m == 0 || Z == 0
        flux = 0f0
    else
        if typeof(mode) == Ani
            flux = 1/4
            flux *= u * (-1 + u^2) * (phi0 - phi1) * log(((-p064 + sqrt(m^2 + p064^2)) * (p164 + sqrt(m^2 + p164^2))^2) / (m^2 * (p064 + sqrt(m^2 + p064^2))))
            flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
            flux *= fluxScale
        elseif typeof(mode) == Axi
            flux = 1/4
            flux *= u * (-1 + u^2) * (phi0 - phi1) * log(((-p064 + sqrt(m^2 + p064^2)) * (p164 + sqrt(m^2 + p164^2))^2) / (m^2 * (p064 + sqrt(m^2 + p064^2))))
            flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
            flux *= fluxScale
        elseif typeof(mode) == Iso
            flux = 0f0
        end
    end

    return flux

end

function KFluxFunction(force::SyncRadReact,space_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u0,u1,phi,name::String)

    # Evaluate the flux through the K surface i.e. surface of constant phi

    mode = force.mode
    B = force.B

    # convert p to Float64 to ensure no floating point issues
    p064 = Float64(p0)
    p164 = Float64(p1)

    μ0 = getfield(DC,Symbol("μ0"))
    m = getfield(DC,Symbol("mu"*name))
    q = getfield(DC,Symbol("q"))
    Z = getfield(DC,Symbol("z"*name))
    mEle = getfield(DC,Symbol("mEle"))
    c = getfield(DC,Symbol("c"))

    fluxScale = (Z^4*B^2)/(μ0*m^2*mEle*c^2) / (2pi) # normalised by σT*c, factor of 4pi^2 is to make sure timescale is correct, no idea why this is needed 

    # convert p to Float64 to ensure no floating point issues
    p064 = Float64(p0)
    p164 = Float64(p1)

    if m == 0 || Z == 0
        flux = 0f0
    else
        flux = 0f0
        flux *= fluxScale
    end

    return flux

end

# ======================================================== #