# =============== ABCD Vol Cartesian ===================== #

    function AFluxFunction(space_coords::Cartesian,momentum_coords::Spherical,t,x0,x1,y0,y1,z0,z1,p0,p1,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the A surface i.e. surface of constant t

        m = getfield(DC,Symbol("mu"*name))

        flux::Float64 = (-x0 + x1) * (-y0 + y1) * (-z0 + z1) 
        flux *= (-p0 + p1) * (-u0 + u1) * (-phi0 + phi1)

        return flux

    end

    function BFluxFunction(space_coords::Cartesian,momentum_coords::Spherical,t0,t1,x,y0,y1,z0,z1,p0,p1,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the B surface i.e. surface of constant x

        m = getfield(DC,Symbol("mu"*name))

        flux::Float64 =  (t0 - t1) * (y0 - y1) * (z0 - z1)
        flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1)) * (sinpi(phi0/pi) - sinpi(phi1/pi)) / 2

        return flux

    end

    function CFluxFunction(space_coords::Cartesian,momentum_coords::Spherical,t0,t1,x0,x1,y,z0,z1,p0,p1,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the C surface i.e. surface of constant y

        m = getfield(DC,Symbol("mu"*name))

        flux::Float64 =  (t0 - t1) * (x0 - x1) * (z0 - z1) 
        flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (-u0*sqrt(1 - u0^2) + u1*sqrt(1 - u1^2) + 2acot_mod(u0) - 2acot_mod(u1)) * (cospi(phi0/pi) - cospi(phi1/pi)) / 2 

        return flux

    end

    function DFluxFunction(space_coords::Cartesian,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z,p0,p1,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the D surface i.e. surface of constant z

        m = getfield(DC,Symbol("mu"*name))

        flux::Float64 = (t0 - t1) * (x0 - x1) * (y0 - y1) 
        flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (u0^2 - u1^2) * (phi0 - phi1) / 2

        return flux

    end

    function VolFunction(space_coords::Cartesian,t0,t1,x0,x1,y0,y1,z0,z1)

        # Evaluate the spacetime volume element

        vol::Float64 = (-t0 + t1) * (-x0 + x1) * (-y0 + y1) * (-z0 + z1)

        return vol

    end

# ======================================================== #

# =============== IJK Cartesian Ricci ==================== #

    function IFluxFunction(force::CoordinateForce,space_coords::Cartesian,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the I surface i.e. surface of constant spherical p

        flux::Float64 = 0e0

        return flux

    end

    function JFluxFunction(force::CoordinateForce,space_coords::Cartesian,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u,phi0,phi1,name::String)

        # Evaluate the flux through the J surface i.e. surface of constant spherical u

        flux::Float64 = 0e0

        return flux

    end

    function KFluxFunction(force::CoordinateForce,space_coords::Cartesian,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u0,u1,phi,name::String)

        # Evaluate the flux through the K surface i.e. surface of constant spherical phi

        flux::Float64 = 0e0

        return flux

    end

# ======================================================== #


# =============== ABCD Vol Cylindrical ===================== #

    function AFluxFunction(space_coords::Cylindrical,momentum_coords::Spherical,t,x0,x1,y0,y1,z0,z1,p0,p1,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the A surface i.e. surface of constant t

        # convert p to Float64 to ensure no floating point issues
        p064 = Float64(p0)
        p164 = Float64(p1)

        m = getfield(DC,Symbol("mu"*name))

        flux = (-x0^2 + x1^2) / 2 * (-y0 + y1) * (-z0 + z1) 
        flux *= (-p064 + p164) * (-u0 + u1) * (-phi0 + phi1)

        return flux

    end

    function BFluxFunction(space_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x,y0,y1,z0,z1,p0,p1,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the B surface i.e. surface of constant x, spherical r, cylindrical r

        # convert p to Float64 to ensure no floating point issues
        p064 = Float64(p0)
        p164 = Float64(p1)

        m = getfield(DC,Symbol("mu"*name))

        flux = (t0 - t1) * x * (y0 - y1) * (z0 - z1) 
        flux *= (1/2)*(sqrt(m^2 + p064^2) - sqrt(m^2 + p164^2)) * (u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - 2*acot_mod(u0) + 2*acot_mod(u1))*(sinpi(phi0/pi) - sinpi(phi1/pi))

        return flux

    end

    function CFluxFunction(space_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y,z0,z1,p0,p1,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the C surface i.e. surface of constant y, spherical θ, cylindrical ϕ

        # convert p to Float64 to ensure no floating point issues
        p064 = Float64(p0)
        p164 = Float64(p1)

        m = getfield(DC,Symbol("mu"*name))

        flux = (t0 - t1) * (x0 - x1) * (z0 - z1) 
        flux *= (1/2)*(-sqrt(m^2 + p064^2) + sqrt(m^2 + p164^2)) * (u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - 2*acot_mod(u0) + 2*acot_mod(u1))*(cospi(phi0/pi) - cospi(phi1/pi))

        return flux

    end

    function DFluxFunction(space_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z,p0,p1,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the D surface i.e. surface of constant z, spherical ϕ, cylindrical z

        # convert p to Float64 to ensure no floating point issues
        p064 = Float64(p0)
        p164 = Float64(p1)

        m = getfield(DC,Symbol("mu"*name))

        flux = (t0 - t1) * (x0^2 - x1^2) / 2 * (y0 - y1) 
        flux *= (1/2) * (sqrt(m^2 + p064^2) - sqrt(m^2 + p164^2)) * (u0^2 - u1^2) * (phi0 - phi1)

        return flux

    end

    function VolFunction(space_coords::Cylindrical,t0,t1,x0,x1,y0,y1,z0,z1)

        # Evaluate the spacetime volume element

        return (-t0 + t1) * (-x0^2 + x1^2) * (-y0 + y1) * (-z0 + z1) / 2

    end

# ======================================================== #

# =============== IJK Cylindrical Ricci ================== #

    function IFluxFunction(force::CoordinateForce,space_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the I surface i.e. surface of constant spherical p

        α::Float64 = space_coords.α
        β::Float64 = space_coords.β
        γ::Float64 = space_coords.γ

        flux::Float64 = 0e0

        return flux

    end

    function JFluxFunction(force::CoordinateForce,space_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u,phi0,phi1,name::String)

        # Evaluate the flux through the J surface i.e. surface of constant spherical u

        α::Float64 = space_coords.α
        β::Float64 = space_coords.β
        γ::Float64 = space_coords.γ

        if α != 0.0 || γ != 0.0
            error("JFluxFunction not implemented for non-zero α or γ in Cylindrical coordinates.")
        end

        # convert p to Float64 to ensure no floating point issues
        p064 = Float64(p0)
        p164 = Float64(p1)
        phi0pi::Float64 = Float64(phi0/pi) # using sinpi rather than sin to avoid float issues 
        phi1pi::Float64 = Float64(phi1/pi) # using sinpi rather than sin to avoid float issues

        m = getfield(DC,Symbol("mu"*name))

        flux::Float64 = -(1/2)*(-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2))*(-t0 + t1)*(-x0 + x1)*(-y0 + y1)*(-z0 + z1)
        flux *=  sinpi(β)*(sinpi(phi0pi) - sinpi(phi1pi))*(-2u*sqrt(1 - u^2)*sinpi(β) + (-1 + u^2)*cospi(β)*(sinpi(phi0pi) + sinpi(phi1pi)))

        return flux

    end

    function KFluxFunction(force::CoordinateForce,space_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u0,u1,phi,name::String)

        # Evaluate the flux through the K surface i.e. surface of constant spherical phi

        α::Float64 = space_coords.α
        β::Float64 = space_coords.β
        γ::Float64 = space_coords.γ

        if α != 0.0 || γ != 0.0
            error("KFluxFunction not implemented for non-zero α or γ in Cylindrical coordinates.")
        end

        # convert p to Float64 to ensure no floating point issues
        p064 = Float64(p0)
        p164 = Float64(p1)
        phipi::Float64 = Float64(phi/pi) # using sinpi rather than sin to avoid float issues 

        m = getfield(DC,Symbol("mu"*name))

        #flux = (1/2)*(-sqrt(m^2 + p064^2) + sqrt(m^2 + p164^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2*acot_mod(u0) + 2*acot_mod(u1))*sinpi(phipi)
        #flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)

        flux::Float64 = (1/4)*(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*(t0 - t1)*(x0 - x1)*(y0 - y1)*(z0 - z1)
        flux *= (-3u0^2*cospi(β)*sinpi(β) + 3u1^2*cospi(β)*sinpi(β) + u0^2*cospi(β)*cospi(phi)^2*sinpi(β) - u1^2*cospi(β)*cospi(phi)^2*sinpi(β) + 4acot_mod(u1)*sinpi(phi) - 4acot_mod(u0)*sinpi(phi) - 2*u0*sqrt(1 - u0^2)*cospi(β)^2*sinpi(phi) + 2u1*sqrt(1 - u1^2)*cospi(β)^2*sinpi(phi) + 2u0*sqrt(1 - u0^2)*sinpi(β)^2*sinpi(phi) - 2u1*sqrt(1 - u1^2)*sinpi(β)^2*sinpi(phi) - u0^2*cospi(β)*sinpi(β)*sinpi(phi)^2 + u1^2*cospi(β)*sinpi(β)*sinpi(phi)^2)

        return flux

    end

# ======================================================== #

# =============== ABCD Vol Spherical ===================== #

    function AFluxFunction(space_coords::Spherical,momentum_coords::Spherical,t,x0,x1,y0,y1,z0,z1,p0,p1,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the A surface i.e. surface of constant t

        m = getfield(DC,Symbol("mu"*name))

        flux::Float64 = -(1/3) * (x0^3 - x1^3) * (z0 - z1) * (cospi(y0/pi) - cospi(y1/pi))
        flux *= (p0 - p1) * (u0 - u1) * (phi0 - phi1)

        return flux

    end

    function BFluxFunction(space_coords::Spherical,momentum_coords::Spherical,t0,t1,x,y0,y1,z0,z1,p0,p1,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the B surface i.e. surface of constant spherical r

        m = getfield(DC,Symbol("mu"*name))

        flux::Float64 = (t0 - t1) * x^2 * (z0 - z1) * (cospi(y0/pi) - cospi(y1/pi)) * (sinpi(phi0/pi) - sinpi(phi1/pi)) / 2
        flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (-u0*sqrt(1 - u0^2) + u1*sqrt(1 - u1^2) + 2acot_mod(u0) - 2acot_mod(u1)) * (sinpi(phi0/pi) - sinpi(phi1/pi))
        return flux

    end

    function CFluxFunction(space_coords::Spherical,momentum_coords::Spherical,t0,t1,x0,x1,y,z0,z1,p0,p1,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the C surface i.e. surface of constant spherical θ 

        m = getfield(DC,Symbol("mu"*name))

        flux::Float64 = (t0 - t1) * (x0^2 - x1^2) * (z0 - z1) * sin(y) / 2
        flux *= (-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1)) * (cospi(phi0/pi) - cospi(phi1/pi)) / 2

        return flux

    end

    function DFluxFunction(space_coords::Spherical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z,p0,p1,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the D surface i.e. surface of constant spherical ϕ

        m = getfield(DC,Symbol("mu"*name))

        flux::Float64 = (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) / 2
        flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (u0^2 - u1^2) * (phi0 - phi1) / 2

        return flux

    end

    function VolFunction(space_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1)

        # Evaluate the spacetime volume element

        vol::Float64 = -(1/3) * (t0 - t1) * (x0^3 - x1^3) * (z0 - z1) * (cospi(y0/pi) - cospi(y1/pi))

        return vol

    end

# ======================================================== #

# =============== IJK Spherical Ricci ==================== #

    function IFluxFunction(force::CoordinateForce,space_coords::Spherical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the I surface i.e. surface of constant spherical p

        flux::Float64 = 0e0

        return flux

    end

    function JFluxFunction(force::CoordinateForce,space_coords::Spherical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u,phi0,phi1,name::String)

        # Evaluate the flux through the J surface i.e. surface of constant spherical u

        flux::Float64 = -2(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (t0 - t1) * u*sqrt(1 - u^2) * (x0^2 - x1^2) * (z0 - z1) * sinpi((y0 - y1)/(2pi)) * sinpi((phi0 - phi1)/(2pi)) * sinpi((y0 + y1 + phi0 + phi1)/(2pi))

        return flux

    end

    function KFluxFunction(force::CoordinateForce,space_coords::Spherical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u0,u1,phi,name::String)

        # Evaluate the flux through the K surface i.e. surface of constant spherical phi

        flux::Float64 =  (t0 - t1) * (x0^2 - x1^2) * (z0 - z1) / 2
        flux *= (-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * ((u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2)) * cospi(phi/pi) * (sinpi(y0/pi) - sinpi(y1/pi)) + arcsin(u1) * (cospi(phi/pi) * (sinpi(y0/pi) - sinpi(y1/pi)) + 2(cospi(y0/pi) - cospi(y1/pi)) * sinpi(phi/pi)) + arcsin(u0) * (cospi(phi/pi) * (-sinpi(y0/pi) + sinpi(y1/pi)) + 2(-cospi(y0/pi) + cospi(y1/pi))*sinpi(phi/pi))) / 2

        return flux

    end

# ======================================================== #

# ======= IJK Cylindrical Ricci aligned to B field ======= #

    #=
    B field is formed of a axial field and azimuthal field, such that B^α = (0,0,b1/r^2,b2), local coordinates are aligned to this field such that B^a = (0,0,0,B) where B= sqrt(b1^2/r^2+b2^2)
    =#

    function IFluxFunction(force::CoordinateForce,space_coords::CylindricalMag,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p,u0,u1,phi0,phi1,name::String)

        # Evaluate the flux through the I surface i.e. surface of constant spherical p

        b1::Float64 = space_coords.b1
        b2::Float64 = space_coords.b2

        B::Float64 = b1/b2

        if b1 == 0 # no azimuthal field
            flux = 0e0
        elseif b2 == 0 # no axial field
            flux = (1/(48*sqrt(m^2 + p^2)))*p^2*(t0 - t1)*(x0 - x1)*(y0 - y1)*(z0 - z1)*(sinpi(phi0/pi) - sinpi(phi1/pi))*(-2* u0*sqrt(1 - u0^2)*(-11 + x0^2 + x0*x1 + x1^2) + 4*u0^3*sqrt(1 - u0^2)*(-7 + x0^2 + x0*x1 + x1^2) + 2* u1*sqrt(1 - u1^2)*(-11 + x0^2 + x1*(x0 + x1) - 2*u1^2*(-7 + x0^2 + x0*x1 + x1^2)) + (u0*sqrt(1 - u0^2)*(-5 + 2*u0^2) + u1*(5 - 2*u1^2)*sqrt(1 - u1^2))*(cos(2*phi0) + cos(2*phi1) - 2*sinpi(phi0/pi)*sinpi(phi1/pi)) - 2*acot_mod(u0)*(-3*cos(2*phi0) - 3*cos(2*phi1) + 2*(-3 + x0^2 + x0*x1 + x1^2 + 3*sinpi(phi0/pi)*sinpi(phi1/pi))) + 2*acot_mod(u1)*(-3*cos(2*phi0) - 3*cos(2*phi1) + 2*(-3 + x0^2 + x0*x1 + x1^2 + 3*sinpi(phi0/pi)*sinpi(phi1/pi))))
        else
            #flux = 
        end

        return flux

    end

    function JFluxFunction(force::CoordinateForce,space_coords::CylindricalMag,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u,phi0,phi1,name::String)

        # Evaluate the flux through the J surface i.e. surface of constant spherical u

        return 0f0

    end

    function KFluxFunction(force::CoordinateForce,space_coords::CylindricalMag,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u0,u1,phi,name::String)

        # Evaluate the flux through the K surface i.e. surface of constant spherical phi

        # convert p to Float64 to ensure no floating point issues
        p064 = Float64(p0)
        p164 = Float64(p1)
        phipi = Float64(phi/pi) # using sinpi rather than sin to avoid float issues 

        m = getfield(DC,Symbol("mu"*name))

        flux = (1/2)*(-sqrt(m^2 + p064^2) + sqrt(m^2 + p164^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2*acot_mod(u0) + 2*acot_mod(u1))*sinpi(phipi)
        flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)

        return flux

    end

# ======================================================== #

function acot_mod(x)

    # ArcCot[(1+x)/Sqrt[1-x^2]]
    # this is actually just acos(x)/2

    if x == -1
        return pi/2
    elseif x == 1
        return 0.0
    else
        return acot((1 + x)/sqrt(1 - x^2))
    end

end