#= A Note of Dimensions 

    If spatial coordinates have dimension L and time coordinates have dimensions of cT then:
        - for A fluxes the dimensions of `flux` are L^3
        - for B, C, D fluxes the dimensions of `flux` are L^2*cT
        - for I, J, and K fluxes the dimensions of `flux` are cT*L^2
        - for volume elements the dimensions are L^3*cT

=#

# =============== ABCD Vol  ===================== #

    function AFluxFunction(PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        m = Grids.mass_list[species_idx]

        if plus_minus == "plus"
            t = Grids.tr[t_idx+1]
        elseif plus_minus == "minus"
            t = Grids.tr[t_idx]
        else
            error("plus_minus string not recognised.")
        end
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]
        p0 = Grids.pxr_list[species_idx][px_idx]
        p1 = Grids.pxr_list[species_idx][px_idx+1]
        u0 = Grids.pyr_list[species_idx][py_idx]
        u1 = Grids.pyr_list[species_idx][py_idx+1]
        h0 = Grids.pzr_list[species_idx][pz_idx]
        h1 = Grids.pzr_list[species_idx][pz_idx+1]

        # Evaluate the flux through the A surface i.e. surface of constant t

        flux::Float64 = 1.0

        if space_coords isa Cartesian && momentum_coords isa Spherical

            flux *= (-x0 + x1) * (-y0 + y1) * (-z0 + z1) 
            flux *= (-p0 + p1) * (-u0 + u1) * (-h0 + h1)

        elseif space_coords isa Cylindrical && momentum_coords isa Spherical

            flux *= (-x0^2 + x1^2) / 2 * (-y0 + y1) * (-z0 + z1) 
            flux *= (-p0 + p1) * (-u0 + u1) * (-h0 + h1)

        elseif space_coords isa Spherical && momentum_coords isa Spherical

            flux *= -(1/3) * (x0^3 - x1^3) * (z0 - z1) * (cospi(y0/pi) - cospi(y1/pi))
            flux *= (p0 - p1) * (u0 - u1) * (h0 - h1)

        else
            error("A Flux function not implemented for this combination of co-ordinate systems.")
        end

        return flux

    end

    function BFluxFunction(PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        m = Grids.mass_list[species_idx]

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        if plus_minus == "plus"
            x = Grids.xr[x_idx+1]
        elseif plus_minus == "minus"
            x = Grids.xr[x_idx]
        else
            error("plus_minus string not recognised.")
        end
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]
        p0 = Grids.pxr_list[species_idx][px_idx]
        p1 = Grids.pxr_list[species_idx][px_idx+1]
        u0 = Grids.pyr_list[species_idx][py_idx]
        u1 = Grids.pyr_list[species_idx][py_idx+1]
        h0 = Grids.pzr_list[species_idx][pz_idx]
        h1 = Grids.pzr_list[species_idx][pz_idx+1]

        # Evaluate the flux through the B surface i.e. surface of constant x coordinate: cartesian x, spherical r, cylindrical r

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        if space_coords isa Cartesian && momentum_coords isa Spherical

            h1 = h1/pi
            h0 = h0/pi
        
            flux *=  (t0 - t1) * (y0 - y1) * (z0 - z1)
            flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1)) * (sinpi(h0) - sinpi(h1)) / 2

        elseif space_coords isa Cylindrical && momentum_coords isa Spherical

            α::Float64 = Rotation[1]
            β::Float64 = Rotation[2]
            γ::Float64 = Rotation[3]

            if α != 0.0
                error("JFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            h1 = h1/pi
            h0 = h0/pi

            flux *= 1/2 * (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (t0 - t1) * x * (y0 - y1) * (z0 - z1)
            flux *= (-u0 * sqrt(1 - u0^2) + u1 * sqrt(1 - u1^2) + 2*acot_mod(u0) - 2*acot_mod(u1)) * (sinpi(γ - h0) - sinpi(γ - h1))



        elseif space_coords isa Spherical && momentum_coords isa Spherical 

            h1 = h1/pi
            h0 = h0/pi
            y1 = y1/pi
            y0 = y0/pi

            flux *= (t0 - t1) * x^2 * (z0 - z1) * (cospi(y0) - cospi(y1)) / 2
            flux *= (-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * (u0^2 - u1^2) * pi*(h0 - h1)

        else
            error("B Flux function not implemented for this combination of co-ordinate systems.")
        end

        return flux

    end

    function CFluxFunction(PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        m = Grids.mass_list[species_idx]

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        if plus_minus == "plus"
            y = Grids.yr[y_idx+1]
        elseif plus_minus == "minus"
            y = Grids.yr[y_idx]
        else
            error("plus_minus string not recognised.")
        end
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]
        p0 = Grids.pxr_list[species_idx][px_idx]
        p1 = Grids.pxr_list[species_idx][px_idx+1]
        u0 = Grids.pyr_list[species_idx][py_idx]
        u1 = Grids.pyr_list[species_idx][py_idx+1]
        h0 = Grids.pzr_list[species_idx][pz_idx]
        h1 = Grids.pzr_list[species_idx][pz_idx+1]

        # Evaluate the flux through the C surface i.e. surface of constant y coordinate: cartesian y, spherical θ, cylindrical ϕ

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        if space_coords isa Cartesian && momentum_coords isa Spherical

            h1 = h1/pi
            h0 = h0/pi

            flux *=  (t0 - t1) * (x0 - x1) * (z0 - z1) 
            flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (-u0*sqrt(1 - u0^2) + u1*sqrt(1 - u1^2) + 2acot_mod(u0) - 2acot_mod(u1)) * (cospi(h0) - cospi(h1)) / 2 

        elseif space_coords isa Cylindrical && momentum_coords isa Spherical

            α::Float64 = Rotation[1]
            β::Float64 = Rotation[2]
            γ::Float64 = Rotation[3]

            if α != 0.0
                error("JFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            h1 = h1/pi
            h0 = h0/pi

            flux *= 1/2*(-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * (-t0 + t1) * (-x0 + x1) * (-z0 + z1) 
            flux *= ((-u0 * sqrt(1 - u0^2) + u1 * sqrt(1 - u1^2) + 2*acot_mod(u0) - 2*acot_mod(u1)) * cospi(β) * (cospi(γ - h0) - cospi(γ - h1)) + u0^2 * h0 * sinpi(β) - u1^2 * h0 * sinpi(β) - u0^2 * h1 * sinpi(β) + u1^2 * h1 * sinpi(β))

        elseif space_coords isa Spherical && momentum_coords isa Spherical

            h1 = h1/pi
            h0 = h0/pi

            flux *= (t0 - t1) * (x0^2 - x1^2) * (z0 - z1) * sin(y) / 2
            flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1)) * (sinpi(h0) - sinpi(h1)) / 2

        else 
            error("C Flux function not implemented for this combination of co-ordinate systems.")
        end

        return flux

    end

    function DFluxFunction(PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        if plus_minus == "plus"
            z = Grids.zr[z_idx+1]
        elseif plus_minus == "minus"
            z = Grids.zr[z_idx]
        else
            error("plus_minus string not recognised.")
        end
        p0 = Grids.pxr_list[species_idx][px_idx]
        p1 = Grids.pxr_list[species_idx][px_idx+1]
        u0 = Grids.pyr_list[species_idx][py_idx]
        u1 = Grids.pyr_list[species_idx][py_idx+1]
        h0 = Grids.pzr_list[species_idx][pz_idx]
        h1 = Grids.pzr_list[species_idx][pz_idx+1]

        # Evaluate the flux through the D surface i.e. surface of constant z coordinate: cartesian z, spherical ϕ, cylindrical z

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        if space_coords isa Cartesian && momentum_coords isa Spherical

            flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) 
            flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (u0^2 - u1^2) * (h0 - h1) / 2

        elseif space_coords isa Cylindrical && momentum_coords isa Spherical

            α::Float64 = Rotation[1]
            β::Float64 = Rotation[2]
            γ::Float64 = Rotation[3]

            if α != 0.0
                error("JFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            flux *= (1/4) * (-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) 
            flux *= ((-u0^2 + u1^2) * (h0 - h1) * cospi(β) + (-u0 * sqrt(1 - u0^2) + u1 * sqrt(1 - u1^2) + 2*acot_mod(u0) - 2*acot_mod(u1)) * (cospi(γ - h0) - cospi(γ - h1)) * sinpi(β))

        elseif space_coords isa Spherical && momentum_coords isa Spherical

            h1 = h1/pi
            h0 = h0/pi

            flux *= (t0 - t1) * (x0^2 - x1^2) * (z0 - z1) * (y0 - y1) / 2
            flux *= (-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1)) * (cospi(h0) - cospi(h1)) / 2

        else
            error("D Flux function not implemented for this combination of co-ordinate systems.")
        end

        return flux

    end

    function VolFunction(PhaseSpace::PhaseSpaceStruct,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64)

        Grids = PhaseSpace.Grids
        space_coords = PhaseSpace.Space.space_coordinates

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]

        # Evaluate the spacetime volume element

        vol::Float64 = 1.0

        if space_coords isa Cartesian

            vol *= (-t0 + t1) * (-x0 + x1) * (-y0 + y1) * (-z0 + z1)

        elseif space_coords isa Cylindrical

            vol *=(-t0 + t1) * (-x0^2 + x1^2) * (-y0 + y1) * (-z0 + z1) / 2

        elseif space_coords isa Spherical

            y1 = y1/pi
            y0 = y0/pi

            vol *= -(1/3) * (t0 - t1) * (x0^3 - x1^3) * (z0 - z1) * (cospi(y0) - cospi(y1))

        else
            error("Volume function not implemented for this co-ordinate system.")
        end

        return vol

    end

# ======================================================== #

# =============== IJK Cartesian Ricci ==================== #

    function IFluxFunction(force::CoordinateForce,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]
        if plus_minus == "plus"
            p = Grids.pxr_list[species_idx][px_idx+1]
        elseif plus_minus == "minus"
            p = Grids.pxr_list[species_idx][px_idx]
        else
            error("plus_minus string not recognised.")
        end
        u0 = Grids.pyr_list[species_idx][py_idx]
        u1 = Grids.pyr_list[species_idx][py_idx+1]
        h0 = Grids.pzr_list[species_idx][pz_idx]
        h1 = Grids.pzr_list[species_idx][pz_idx+1]

        # Evaluate the flux through the I surface i.e. surface of constant spherical p

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        if space_coords isa Cartesian && momentum_coords isa Spherical

            flux *= 0e0

        elseif space_coords isa Cylindrical && momentum_coords isa Spherical

            flux *= 0e0

        elseif space_coords isa Spherical && momentum_coords isa Spherical

            flux *= 0e0

        else
            error("IFluxFunction not implemented for this combination of co-ordinate systems.")
        end

        return flux

    end

    function JFluxFunction(force::CoordinateForce,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        m = Grids.mass_list[species_idx]

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]
        p0 = Grids.pxr_list[species_idx][px_idx]
        p1 = Grids.pxr_list[species_idx][px_idx+1]
        if plus_minus == "plus"
            u = Grids.pyr_list[species_idx][py_idx+1]
        elseif plus_minus == "minus"
            u = Grids.pyr_list[species_idx][py_idx]
        else
            error("plus_minus string not recognised.")
        end
        h0 = Grids.pzr_list[species_idx][pz_idx]
        h1 = Grids.pzr_list[species_idx][pz_idx+1]

        # Evaluate the flux through the J surface i.e. surface of constant spherical u

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        if space_coords isa Cartesian && momentum_coords isa Spherical

            flux *= 0e0

        elseif space_coords isa Cylindrical && momentum_coords isa Spherical

            α::Float64 = Rotation[1]
            β::Float64 = Rotation[2]
            γ::Float64 = Rotation[3]

            if α != 0.0
                error("JFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            h1 = h1/pi
            h0 = h0/pi

            flux *= 1/2 * (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1) 
            flux *= cospi(γ - h0/2 - h1/2) * sinpi(β) * sinpi((h0 - h1)/2) * (-4u*sqrt(1 - u^2)*sinpi(β) -4(-1 + u^2)*cospi(β)*cospi((h0 - h1)/2)*sinpi(γ - h0/2 - h1/2))


        elseif space_coords isa Spherical && momentum_coords isa Spherical

            h1 = h1/pi
            h0 = h0/pi
            y1 = y1/pi
            y0 = y0/pi

            flux *= (t0 - t1) * (x0^2 - x1^2) * (z0 - z1) * (cospi(y0) - cospi(y1)) / 2
            flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (-1+u^2) * pi*(h0 - h1)

        else
            error("JFluxFunction not implemented for this combination of co-ordinate systems.")
        end

        return flux

    end

    function KFluxFunction(force::CoordinateForce,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        Rotation = PhaseSpace.GlobalToLocalRotation
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        m = Grids.mass_list[species_idx]

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]
        p0 = Grids.pxr_list[species_idx][px_idx]
        p1 = Grids.pxr_list[species_idx][px_idx+1]
        u0 = Grids.pyr_list[species_idx][py_idx]
        u1 = Grids.pyr_list[species_idx][py_idx+1]
        if plus_minus == "plus"
            h = Grids.pzr_list[species_idx][pz_idx+1]
        elseif plus_minus == "minus"
            h = Grids.pzr_list[species_idx][pz_idx]
        else
            error("plus_minus string not recognised.")
        end

        # Evaluate the flux through the K surface i.e. surface of constant spherical phi

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        if space_coords isa Cartesian && momentum_coords isa Spherical

            flux *= 0e0

        elseif space_coords isa Cylindrical && momentum_coords isa Spherical

            α::Float64 = Rotation[1]
            β::Float64 = Rotation[2]
            γ::Float64 = Rotation[3]

            if α != 0.0
                error("KFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            h = h/pi

            flux *= -(1/16) * (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)
            flux *= (6u0^2*sinpi(2*β) - 6u1^2*sinpi(2*β) - 8*asin(u0)*sinpi(γ - h) + 8*asin(u1)*sinpi(γ - h) - u0^2*sinpi(2*(β + γ - h)) + u1^2*sinpi(2*(β + γ - h)) - 4*u0*sqrt(1 - u0^2)*sinpi(2*β + γ - h) + 4*u1*sqrt(1 - u1^2)*sinpi(2*β + γ - h) - u0^2*sinpi(2*(β - γ + h)) + u1^2*sinpi(2*(β - γ + h)) + 4*u0*sqrt(1 - u0^2)*sinpi(2*β - γ + h) - 4*u1*sqrt(1 - u1^2)*sinpi(2*β - γ + h))

        elseif space_coords isa Spherical && momentum_coords isa Spherical

            h = h/pi
            y1 = y1/pi
            y0 = y0/pi

            flux *=  (t0 - t1) * (x0^2 - x1^2) * (sinpi(y0) - sinpi(y1)) * (z0 - z1) / 2
            flux *= (-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * ((u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2*acot_mod(u0) + 2*acot_mod(u1)) * sinpi(h)) / 2

        else
            error("KFluxFunction not implemented for this combination of co-ordinate systems.")
        end

        return flux

    end

# ======================================================== #

function acot_mod(x)

    # ArcCot[(1+x)/Sqrt[1-x^2]]
    # this is actually just acos(x)/2

    #=if x == -1
        return pi/2
    elseif x == 1
        return 0.0
    else
        return acot((1 + x)/sqrt(1 - x^2))
    end=#

    return acos(x)/2

end



# OUTDATED BELOW HERE
#= # ======= IJK Cylindrical Ricci aligned to B field ======= #

    #=
    B field is formed of a axial field and azimuthal field, such that B^α = (0,0,b1/r^2,b2), local coordinates are aligned to this field such that B^a = (0,0,0,B) where B= sqrt(b1^2/r^2+b2^2)
    =#

    function IFluxFunction(force::CoordinateForce,space_coords::CylindricalMag,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p,u0,u1,h0,h1,name::String)

        # Evaluate the flux through the I surface i.e. surface of constant spherical p

        b1::Float64 = space_coords.b1
        b2::Float64 = space_coords.b2

        B::Float64 = b1/b2

        if b1 == 0 # no azimuthal field
            flux = 0e0
        elseif b2 == 0 # no axial field
            flux = (1/(48*sqrt(m^2 + p^2)))*p^2*(t0 - t1)*(x0 - x1)*(y0 - y1)*(z0 - z1)*(sinpi(h0/pi) - sinpi(h1/pi))*(-2* u0*sqrt(1 - u0^2)*(-11 + x0^2 + x0*x1 + x1^2) + 4*u0^3*sqrt(1 - u0^2)*(-7 + x0^2 + x0*x1 + x1^2) + 2* u1*sqrt(1 - u1^2)*(-11 + x0^2 + x1*(x0 + x1) - 2*u1^2*(-7 + x0^2 + x0*x1 + x1^2)) + (u0*sqrt(1 - u0^2)*(-5 + 2*u0^2) + u1*(5 - 2*u1^2)*sqrt(1 - u1^2))*(cos(2*h0) + cos(2*h1) - 2*sinpi(h0/pi)*sinpi(h1/pi)) - 2*acot_mod(u0)*(-3*cos(2*h0) - 3*cos(2*h1) + 2*(-3 + x0^2 + x0*x1 + x1^2 + 3*sinpi(h0/pi)*sinpi(h1/pi))) + 2*acot_mod(u1)*(-3*cos(2*h0) - 3*cos(2*h1) + 2*(-3 + x0^2 + x0*x1 + x1^2 + 3*sinpi(h0/pi)*sinpi(h1/pi))))
        else
            #flux = 
        end

        return flux

    end

    function JFluxFunction(force::CoordinateForce,space_coords::CylindricalMag,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u,h0,h1,name::String)

        # Evaluate the flux through the J surface i.e. surface of constant spherical u

        return 0f0

    end

    function KFluxFunction(force::CoordinateForce,space_coords::CylindricalMag,momentum_coords::Spherical,t0,t1,x0,x1,y0,y1,z0,z1,p0,p1,u0,u1,phi,name::String)

        # Evaluate the flux through the K surface i.e. surface of constant spherical phi

        # convert p to Float64 to ensure no floating point issues
        p064 = Float64(p0)
        p164 = Float64(p1)
        phipi = Float64(phi/pi) # using sinpi rather than sin to avoid float issues 

        m = eval(Symbol("CONST_mu"*name))

        flux = (1/2)*(-sqrt(m^2 + p064^2) + sqrt(m^2 + p164^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2*acot_mod(u0) + 2*acot_mod(u1))*sinpi(phipi)
        flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)

        return flux

    end

# ======================================================== # =#