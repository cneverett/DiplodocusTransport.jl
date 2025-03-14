# =============== IJK Sync Rad Reaction ================== #

function IFluxFunction(force::SyncRadReact,space_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p,u0,u1,phi0,phi1,name::String)

    # Evaluate the flux through the I surface i.e. surface of constant p

    # convert p to Float64 to ensure no floating point issues
    p64 = Float64(p)

    μ0 = getfield(BCI,Symbol("μ0"))
    m = getfield(BCI,Symbol("mu"*name))
    q = getfield(BCI,Symbol("q"))
    Z = getfield(BCI,Symbol("z"*name))
    mEle = getfield(BCI,Symbol("mEle"))
    c = getfield(BCI,Symbol("c"))

    # TOFIX Variable B
    B = 1f-4

    fluxScale = (Z^4*B^2)/(μ0*m^3*mEle*c^2) # normalised by σT*c 

    if m == 0 || Z == 0
        return flux = 0f0
    else
        #flux = (Z^4*q^4*μ0*B^2)/(72*pi^2*m^2*mEle^2) * 1/sqrt(m^2 + p64^2) * (p64 + p64^3) * (t0 - t1) * (-3*u0 + u0^3 + 3*u1 - u1^3) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1) *(phi0 - phi1)

        flux = 1/(12*m^2*pi)
        flux *= p64*sqrt(m^2 + p64^2) * (-3*u0 + u0^3 + 3*u1 - u1^3)  * (phi0 - phi1)
        flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
        flux *= fluxScale
    end

    return Float32(flux)

end

function JFluxFunction(force::SyncRadReact,space_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u,phi0,phi1,name::String)

    # Evaluate the flux through the J surface i.e. surface of constant u

    # convert p to Float64 to ensure no floating point issues
    p064 = Float64(p0)
    p164 = Float64(p1)

    μ0 = getfield(BCI,Symbol("μ0"))
    m = getfield(BCI,Symbol("mu"*name))
    q = getfield(BCI,Symbol("q"))
    Z = getfield(BCI,Symbol("z"*name))
    mEle = getfield(BCI,Symbol("mEle"))
    c = getfield(BCI,Symbol("c"))

    # TOFIX Variable B
    B = 1f-4

    fluxScale = (Z^4*B^2)/(μ0*m^3*mEle*c^2) # normalised by σT*c 

    if m == 0 || Z == 0
        flux = 0f0
    else
#=         flux = -(Z^4*q^4*μ0*B^2)/(48*pi^2*m^2*mEle^2)
        flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1) 
        flux *= (-u + u^3) * (phi0 - phi1) * log(((-p064 + sqrt(m^2 + p064^2)) * (p164 + sqrt(m^2 + p164^2))^2)/(m^2 * (p064 + sqrt(m^2 + p064^2))))
        flux /= getfield(BCI,Symbol("σT"))#*getfield(BCI,Symbol("c")) =#

        flux = -1/(8*pi)
        flux *= u * (-1 + u^2) * (phi0 - phi1) * log(((-p064 + sqrt(m^2 + p064^2)) * (p164 + sqrt(m^2 + p164^2))^2) / (m^2 * (p064 + sqrt(m^2 + p064^2))))
        flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
        flux *= fluxScale
    end

    return Float32(flux)

end

function KFluxFunction(force::SyncRadReact,space_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u0,u1,phi,name::String)

    # Evaluate the flux through the K surface i.e. surface of constant phi

    # convert p to Float64 to ensure no floating point issues
    p064 = Float64(p0)
    p164 = Float64(p1)

    μ0 = getfield(BCI,Symbol("μ0"))
    m = getfield(BCI,Symbol("mu"*name))
    q = getfield(BCI,Symbol("q"))
    Z = getfield(BCI,Symbol("z"*name))
    mEle = getfield(BCI,Symbol("mEle"))
    c = getfield(BCI,Symbol("c"))

    # TOFIX Variable B
    B = 1f-4

    fluxScale = (Z^4*B^2)/(μ0*m^3*mEle*c^2) # normalised by σT*c

    # convert p to Float64 to ensure no floating point issues
    p064 = Float64(p0)
    p164 = Float64(p1)

    if m == 0 || Z == 0
        flux = 0f0
    else
        flux = 0f0
        flux *= fluxScale
    end

    return Float32(flux)

end

# ======================================================== #