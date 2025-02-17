# =============== ABCD Cylindrical ===================== #

    function AFluxFunction(spacetime_coords::Cylindrical,momentum_coords::Spherical,t,x0,x1,y0,y1,z0,z1,p0,p1,u0,u1,phi0,phi1,mass::Float32)

        # Evaluate the flux through the A surface i.e. surface of constant t

        # convert p to Float64 to ensure no floating point issues
        p064 = Float64(p0)
        p164 = Float64(p1)

        flux = ((-p064 + p164)*(-u0 + u1)*(-(x0^2/2) + x1^2/2)*(-y0 + y1)*(-z0 + z1)*(-phi0 + phi1))/(2*pi)

        return Float32(flux)

    end

    function BFluxFunction(spacetime_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x,y0,y1,z0,z1,p0,p1,u0,u1,phi0,phi1,mass::Float32)

        # Evaluate the flux through the B surface i.e. surface of constant x

        # convert p to Float64 to ensure no floating point issues
        p064 = Float64(p0)
        p164 = Float64(p1)

        flux = (1/(4*pi))*(sqrt(m^2 + p064^2) - sqrt(m^2 + p164^2))*(t0 - t1) * x * (y0 - y1) * (z0 - z1) * (u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - 2*acot_mod(u0) + 2*acot_mod(u1))*(sin(phi0) - sin(phi1))

        return Float32(flux)

    end

    function CFluxFunction(spacetime_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y,z0,z1,p0,p1,u0,u1,phi0,phi1,mass::Float32)

        # Evaluate the flux through the I surface i.e. surface of constant phi

        # convert p to Float64 to ensure no floating point issues
        p064 = Float64(p0)
        p164 = Float64(p1)

        flux = (1/(4*pi))*(-sqrt(m^2 + p064^2) + sqrt(m^2 + p164^2))*(t0 - t1) * (x0 - x1) * (z0 - z1) * (u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - 2*acot_mod(u0) + 2*acot_mod(u1))*(cos(phi0) - cos(phi1))

        return Float32(flux)

    end

    function DFluxFunction(spacetime_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z,p0,p1,u0,u1,phi0,phi1,mass::Float32)

        # Evaluate the flux through the D surface i.e. surface of constant z

        # convert p to Float64 to ensure no floating point issues
        p064 = Float64(p0)
        p164 = Float64(p1)

        flux = ((sqrt(m^2 + p064^2) - sqrt(m^2 + p164^2))*(t0 - t1)*(u0^2 - u1^2)*(x0^2 - x1^2)*(y0 - y1)*(phi0 - phi1))/(8*pi)

        return Float32(flux)

    end

# ======================================================== #

# =============== IJK Cylindrical Ricci ================== #

    function IFluxFunction(spacetime_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p,u0,u1,phi0,phi1,mass::Float32)

        # Evaluate the flux through the I surface i.e. surface of constant p

        return 0f0

    end

    function JFluxFunction(spacetime_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u,phi0,phi1,mass::Float32)

        # Evaluate the flux through the J surface i.e. surface of constant u

        return 0f0

    end

    function KFluxFunction(spacetime_coords::Cylindrical,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u0,u1,phi,m::Float32)

        # Evaluate the flux through the K surface i.e. surface of constant phi

        # convert p to Float64 to ensure no floating point issues
        p064 = Float64(p0)
        p164 = Float64(p1)

        flux = (1/(4*pi))*(-sqrt(m^2 + p064^2) + sqrt(m^2 + p164^2)) * (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2*acot_mod(u0) + 2*acot_mod(u1))*sin(phi)

        return Float32(flux)

    end

# ======================================================== #

function acot_mod(x)

    if x == -1
        return pi/2
    elseif x == 1
        return 0
    else
        return acot((1 + x)/sqrt(1 - x^2))
    end

end