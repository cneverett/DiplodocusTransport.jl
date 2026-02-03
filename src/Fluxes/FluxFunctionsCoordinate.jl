#= A Note of Dimensions 

    If spatial coordinates have dimension L and time coordinates have dimensions of cT then:
        - for A fluxes the dimensions of `flux` are L^3
        - for B, C, D fluxes the dimensions of `flux` are L^2*cT
        - for I, J, and K fluxes the dimensions of `flux` are cT*L^2
        - for volume elements the dimensions are L^3*cT

=#

# =============== ABCD Vol Cartesian ===================== #

    function AFluxFunction(space_coords::Cartesian,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        flux::Float64 = (-x0 + x1) * (-y0 + y1) * (-z0 + z1) 
        flux *= (-p0 + p1) * (-u0 + u1) * (-h0 + h1)

        return flux

    end

    function BFluxFunction(space_coords::Cartesian,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        # Evaluate the flux through the B surface i.e. surface of constant x

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        flux *=  (t0 - t1) * (y0 - y1) * (z0 - z1)
        flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1)) * (sinpi(h0/pi) - sinpi(h1/pi)) / 2

        return flux

    end

    function CFluxFunction(space_coords::Cartesian,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        # Evaluate the flux through the C surface i.e. surface of constant y

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        flux *=  (t0 - t1) * (x0 - x1) * (z0 - z1) 
        flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (-u0*sqrt(1 - u0^2) + u1*sqrt(1 - u1^2) + 2acot_mod(u0) - 2acot_mod(u1)) * (cospi(h0/pi) - cospi(h1/pi)) / 2 

        return flux

    end

    function DFluxFunction(space_coords::Cartesian,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        # Evaluate the flux through the D surface i.e. surface of constant z

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) 
        flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (u0^2 - u1^2) * (h0 - h1) / 2

        return flux

    end

    function VolFunction(space_coords::Cartesian,Grids::GridsStruct,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64)

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]

        # Evaluate the spacetime volume element

        vol::Float64 = (-t0 + t1) * (-x0 + x1) * (-y0 + y1) * (-z0 + z1)

        return vol

    end

# ======================================================== #

# =============== IJK Cartesian Ricci ==================== #

    function IFluxFunction(force::CoordinateForce,space_coords::Cartesian,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        flux *= 0e0

        return flux

    end

    function JFluxFunction(force::CoordinateForce,space_coords::Cartesian,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        flux *= 0e0

        return flux

    end

    function KFluxFunction(force::CoordinateForce,space_coords::Cartesian,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        flux *= 0e0

        return flux

    end

# ======================================================== #


# =============== ABCD Vol Cylindrical ===================== #

    function AFluxFunction(space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        flux::Float64 = (-x0^2 + x1^2) / 2 * (-y0 + y1) * (-z0 + z1) 
        flux *= (-p0 + p1) * (-u0 + u1) * (-h0 + h1)

        return flux

    end

    function BFluxFunction(space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        # Evaluate the flux through the B surface i.e. surface of constant x, spherical r, cylindrical r

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        flux *= (t0 - t1) * x * (y0 - y1) * (z0 - z1) 
        flux *= (1/2)*(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - 2*acot_mod(u0) + 2*acot_mod(u1))*(sinpi(h0/pi) - sinpi(h1/pi))

        return flux

    end

    function CFluxFunction(space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        # Evaluate the flux through the C surface i.e. surface of constant y, spherical θ, cylindrical ϕ
         
        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        flux *= (t0 - t1) * (x0 - x1) * (z0 - z1) 
        flux *= (1/2)*(-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * (u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - 2*acot_mod(u0) + 2*acot_mod(u1))*(cospi(h0/pi) - cospi(h1/pi))

        return flux

    end

    function DFluxFunction(space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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
        
        # Evaluate the flux through the D surface i.e. surface of constant z, spherical ϕ, cylindrical z

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        flux *= (t0 - t1) * (x0^2 - x1^2) / 2 * (y0 - y1) 
        flux *= (1/2) * (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (u0^2 - u1^2) * (h0 - h1)

        return flux

    end

    function VolFunction(space_coords::Cylindrical,Grids::GridsStruct,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64)

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]

        # Evaluate the spacetime volume element

        return (-t0 + t1) * (-x0^2 + x1^2) * (-y0 + y1) * (-z0 + z1) / 2

    end

# ======================================================== #

# =============== IJK Cylindrical Ricci ================== #

    function IFluxFunction(force::CoordinateForce,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        α::Float64 = space_coords.α
        β::Float64 = space_coords.β
        γ::Float64 = space_coords.γ

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        flux *= 0e0

        return flux

    end

    function JFluxFunction(force::CoordinateForce,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        α::Float64 = space_coords.α
        β::Float64 = space_coords.β
        γ::Float64 = space_coords.γ

        if α != 0.0 || γ != 0.0
            error("JFluxFunction not implemented for non-zero α or γ in Cylindrical coordinates.")
        end

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        flux *= -(1/2)*(-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2))*(-t0 + t1)*(-x0 + x1)*(-y0 + y1)*(-z0 + z1)
        flux *=  sinpi(β)*(sinpi(h0/pi) - sinpi(h1/pi))*(-2u*sqrt(1 - u^2)*sinpi(β) + (-1 + u^2)*cospi(β)*(sinpi(h0/pi) + sinpi(h1/pi)))

        return flux

    end

    function KFluxFunction(force::CoordinateForce,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        α::Float64 = space_coords.α
        β::Float64 = space_coords.β
        γ::Float64 = space_coords.γ

        if α != 0.0 || γ != 0.0
            error("KFluxFunction not implemented for non-zero α or γ in Cylindrical coordinates.")
        end

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        flux *= (1/4)*(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))*(t0 - t1)*(x0 - x1)*(y0 - y1)*(z0 - z1)
        flux *= (-3u0^2*cospi(β)*sinpi(β) + 3u1^2*cospi(β)*sinpi(β) + u0^2*cospi(β)*cospi(phi/pi)^2*sinpi(β) - u1^2*cospi(β)*cospi(phi/pi)^2*sinpi(β) + 4acot_mod(u1)*sinpi(phi/pi) - 4acot_mod(u0)*sinpi(phi/pi) - 2*u0*sqrt(1 - u0^2)*cospi(β)^2*sinpi(phi/pi) + 2u1*sqrt(1 - u1^2)*cospi(β)^2*sinpi(phi/pi) + 2u0*sqrt(1 - u0^2)*sinpi(β)^2*sinpi(phi/pi) - 2u1*sqrt(1 - u1^2)*sinpi(β)^2*sinpi(phi/pi) - u0^2*cospi(β)*sinpi(β)*sinpi(phi/pi)^2 + u1^2*cospi(β)*sinpi(β)*sinpi(phi/pi)^2)

        return flux

    end

# ======================================================== #

# =============== ABCD Vol Spherical ===================== #

    function AFluxFunction(space_coords::Spherical,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        flux::Float64 = -(1/3) * (x0^3 - x1^3) * (z0 - z1) * (cospi(y0/pi) - cospi(y1/pi))
        flux *= (p0 - p1) * (u0 - u1) * (h0 - h1)

        return flux

    end

    function BFluxFunction(space_coords::Spherical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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
        # Evaluate the flux through the B surface i.e. surface of constant spherical r

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        flux *= (t0 - t1) * x^2 * (z0 - z1) * (cospi(y0/pi) - cospi(y1/pi)) * (sinpi(h0/pi) - sinpi(h1/pi)) / 2
        flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (-u0*sqrt(1 - u0^2) + u1*sqrt(1 - u1^2) + 2acot_mod(u0) - 2acot_mod(u1)) * (sinpi(h0/pi) - sinpi(h1/pi))
        return flux

    end

    function CFluxFunction(space_coords::Spherical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        # Evaluate the flux through the C surface i.e. surface of constant spherical θ 

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        flux *= (t0 - t1) * (x0^2 - x1^2) * (z0 - z1) * sin(y) / 2
        flux *= (-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1)) * (cospi(h0/pi) - cospi(h1/pi)) / 2

        return flux

    end

    function DFluxFunction(space_coords::Spherical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        # Evaluate the flux through the D surface i.e. surface of constant spherical ϕ

        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        flux::Float64 = c * T / L

        flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) / 2
        flux *= (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (u0^2 - u1^2) * (h0 - h1) / 2

        return flux

    end

    function VolFunction(space_coords::Spherical,Grids::GridsStruct,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64)

        t0 = Grids.tr[t_idx]
        t1 = Grids.tr[t_idx+1]
        x0 = Grids.xr[x_idx]
        x1 = Grids.xr[x_idx+1]
        y0 = Grids.yr[y_idx]
        y1 = Grids.yr[y_idx+1]
        z0 = Grids.zr[z_idx]
        z1 = Grids.zr[z_idx+1]

        # Evaluate the spacetime volume element

        vol::Float64 = -(1/3) * (t0 - t1) * (x0^3 - x1^3) * (z0 - z1) * (cospi(y0/pi) - cospi(y1/pi))

        return vol

    end

# ======================================================== #

# =============== IJK Spherical Ricci ==================== #

    function IFluxFunction(force::CoordinateForce,space_coords::Spherical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        flux *= 0e0

        return flux

    end

    function JFluxFunction(force::CoordinateForce,space_coords::Spherical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        flux *= -2(sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (t0 - t1) * u*sqrt(1 - u^2) * (x0^2 - x1^2) * (z0 - z1) * sinpi((y0 - y1)/(2pi)) * sinpi((h0 - h1)/(2pi)) * sinpi((y0 + y1 + h0 + h1)/(2pi))

        return flux

    end

    function KFluxFunction(force::CoordinateForce,space_coords::Spherical,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        flux *=  (t0 - t1) * (x0^2 - x1^2) * (z0 - z1) / 2
        flux *= (-sqrt(m^2 + p0^2) + sqrt(m^2 + p1^2)) * ((u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2)) * cospi(phi/pi) * (sinpi(y0/pi) - sinpi(y1/pi)) + asin(u1) * (cospi(phi/pi) * (sinpi(y0/pi) - sinpi(y1/pi)) + 2(cospi(y0/pi) - cospi(y1/pi)) * sinpi(phi/pi)) + asin(u0) * (cospi(phi/pi) * (-sinpi(y0/pi) + sinpi(y1/pi)) + 2(-cospi(y0/pi) + cospi(y1/pi))*sinpi(phi/pi))) / 2

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