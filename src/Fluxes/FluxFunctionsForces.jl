#= A Note of Dimensions 

    If spatial coordinates have dimension L and time coordinates have dimensions of cT.
    All force fluxes i.e. `flux` calculated in the functions below have dimensions of: either c*T*L^2 or L^3T/T_c depending on if the flux depends on external parameters that denote a characteristic time (e.g. magnetic field in synchrotron radiation reaction). 

=#

# ======= IJK GradB Inv Z Decay Cylindrical ============== #

    function IFluxFunction(force::GradBInvZDecay,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        EMField = PhaseSpace.ElectroMagneticField
        @assert EMField isa ElectroMagneticField_InvZDecay "ElectroMagneticField must be of type ElectroMagneticField_InvZDecay for GradBInvZDecay force."
        @assert space_coords isa Cylindrical "Space coordinates must be cylindrical for GradBInvZDecay force."

        B0 = EMField.parameters[1]
        L0 = EMField.parameters[2]

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

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

        # Evaluate the flux through the I surface i.e. surface of constant p

        # fluxScale must have no dimensions
        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        #invTc =  c * (log(z0) - log(z1)) / (z0-z1) / L  # c removed by normalising mF by mEle^2*c in the definition of the force

        #flux::Float64 = T*invTc
        flux::Float64 = L0 * c * T / L # L0 is the characteristic length in terms of L 

        if m == 0 || Z == 0
            flux *= 0e0
        else
            
            α::Float64 = space_coords.α
            β::Float64 = space_coords.β
            γ::Float64 = space_coords.γ
            if α != 0.0
                error("IFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            flux *= -(1/(32(m^2 + p^2)^(1/2))) * p^2 * (2(u0 - u1)*(u0 + u1)*(-2 + u0^2 + u1^2)*(h0 - h1)*cospi(β) + 2(u0*sqrt(1 - u0^2)*(-5 + 2u0^2) + u1*(5 - 2u1^2)*sqrt(1 - u1^2) - 3asin(u0) + 3asin(u1))*sinpi(β)*sinpi((h0 - h1)/(2pi))*sinpi(γ - h0/(2pi) - h1/(2pi))) 

            # space part (independent of momentum)
            flux *= (log(z0) - log(z1))
            flux *= (t0 - t1) * (x0 - x1) * (x0 + x1) * (y0 - y1)

        end

        return flux

    end

    function JFluxFunction(force::GradBInvZDecay,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)
        
        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        EMField = PhaseSpace.ElectroMagneticField
        @assert EMField isa ElectroMagneticField_InvZDecay "ElectroMagneticField must be of type ElectroMagneticField_InvZDecay for GradBInvZDecay force."
        @assert space_coords isa Cylindrical "Space coordinates must be cylindrical for GradBInvZDecay force."

        B0 = EMField.parameters[1]
        L0 = EMField.parameters[2]

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

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

        # Evaluate the flux through the J surface i.e. surface of constant u

        # fluxScale must have no dimensions
        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        #invTc =  c * (log(z0) - log(z1)) / (z0-z1) / L  # c removed by normalising mF by mEle^2*c in the definition of the force

        flux::Float64 = L0 * c * T / L # L0 is the characteristic length in terms of L 
        #flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0e0
        else
            
            α::Float64 = space_coords.α
            β::Float64 = space_coords.β
            γ::Float64 = space_coords.γ
            if α != 0.0
                error("IFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            # momentum part (independent of space)
            #flux *= m^2 * ((sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (-1 + u^2) * ((-1 + u^2) * (h0 - h1) * cospi(β) + u * sqrt(1 - u^2) * (cospi(γ - h0/pi) - cospi(γ - h1/pi)) * sinpi(β)))/(2 * sqrt((m^2 + p0^2) * (m^2 + p1^2))) 

            flux *= ((sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2)) * (-1 + u^2) * ((-1 + u^2) * (h0 - h1) * cospi(β) + u * sqrt(1 - u^2) * (cospi(γ - h0/pi) - cospi(γ - h1/pi)) * sinpi(β)))/(4) 


            # space part (independent of momentum)
            flux *= (log(z0) - log(z1)) 
            flux *= (t0 - t1) * (x0 - x1) * (x0 + x1) * (y0 - y1)

        end

        return flux

    end

    function KFluxFunction(force::GradBInvZDecay,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
        space_coords = PhaseSpace.Space.space_coordinates
        momentum_coords = PhaseSpace.Momentum.momentum_coordinates

        EMField = PhaseSpace.ElectroMagneticField
        @assert EMField isa ElectroMagneticField_InvZDecay "ElectroMagneticField must be of type ElectroMagneticField_InvZDecay for GradBInvZDecay force."
        @assert space_coords isa Cylindrical "Space coordinates must be cylindrical for GradBInvZDecay force."

        B0 = EMField.parameters[1]
        L0 = EMField.parameters[2]        

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

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

        # Evaluate the flux through the K surface i.e. surface of constant phi

        # fluxScale must have no dimensions
        T = Characteristic.CHAR_time
        L = Characteristic.CHAR_length
        c = CONST_c

        #invTc =  c * (log(z0) - log(z1)) / (z0-z1) / L  # c removed by normalising mF by mEle^2*c in the definition of the force

        #flux::Float64 = T*invTc
        flux::Float64 = L0 * c * T / L # L0 is the characteristic length in terms of L 

        if m == 0 || Z == 0
            flux *= 0e0
        else
            
            α::Float64 = space_coords.α
            β::Float64 = space_coords.β
            γ::Float64 = space_coords.γ
            if α != 0.0
                error("IFluxFunction not implemented for non-zero α in Cylindrical coordinates.")
            end

            # momentum part (independent of space)
            #flux *= -(m^2 * (sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))  * (u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - acos(u0) + acos(u1)) * cospi(γ - h/pi) * sinpi(β))/(4 * sqrt((m^2 + p0^2) * (m^2 + p1^2)))

            flux *= -((sqrt(m^2 + p0^2) - sqrt(m^2 + p1^2))  * (u0 * sqrt(1 - u0^2) - u1 * sqrt(1 - u1^2) - acos(u0) + acos(u1)) * cospi(γ - h/pi) * sinpi(β))/(8)

            # space part (independent of momentum)
            flux *= (log(z0) - log(z1))
            flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1)
        end

        return flux

    end


# ================= IJK ExB Cartesian ==================== #

    function IFluxFunction(force::ExB,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
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

        # Evaluate the flux through the I surface i.e. surface of constant p

        E0 = force.E0
        B0 = force.B0

        invTc = (Z*CONST_q) / (m*CONST_mEle) # old normalisation / (CONST_σT*CONST_c)
        T = Characteristic.CHAR_time

        flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0e0
        else
            # momentum part (independent of space)
            flux *= -(1/2) * E0* (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1)) * (cospi(h0/pi) - cospi(h1/pi))

            # space part (independent of momentum)
            if space_coords isa Cartesian
                flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)
            else
                error("Space co-ordinate system not recognised.")
            end

        end

        return flux

    end

    function JFluxFunction(force::ExB,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
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

        # Evaluate the flux through the J surface i.e. surface of constant u

        E0 = force.E0
        B0 = force.B0

        invTc = (Z*CONST_q) / (m*CONST_mEle) # old normalisation / (CONST_σT*CONST_c)
        T = Characteristic.CHAR_time    

        flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0e0
        else

            # momentum part (independent of space)
            flux *= E0 * log(p1/p0) * u*sqrt(1 - u^2) * (cospi(h0/pi) - cospi(h1/pi))

            # space part (independent of momentum)
            if space_coords isa Cartesian
                flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)
            else
                error("Space co-ordinate system not recognised.")
            end

        end

        return flux

    end

    function KFluxFunction(force::ExB,PhaseSpace::PhaseSpaceStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        Grids = PhaseSpace.Grids
        Characteristic = PhaseSpace.Characteristic
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

        # Evaluate the flux through the K surface i.e. surface of constant phi

        E0 = force.E0
        B0 = force.B0

        invTc = (Z*CONST_q) / (m*CONST_mEle) # old normalisation / (CONST_σT*CONST_c)
        T = Characteristic.CHAR_time    

        flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0e0
        else

            # momentum part (independent of space)
            flux *= -E0 * (-2acot_mod(u0)+2acot_mod(u1)) * cospi(phi/pi) * log(p1/p0) + B0 * (u0 - u1) * (-asinh(p0/m) + asinh(p1/m)) # sign of E0 flipped to make sure direction is correct ?????

            # space part (independent of momentum)
            if space_coords isa Cartesian
                flux *= (t0 - t1) * (x0 - x1) * (y0 - y1) * (z0 - z1)
            else
                error("Space co-ordinate system not recognised.")
            end

        end

        return flux

    end

# ======================================================== #


# TODO: Update this function for aximuthal B and radial E to meet standards of functions above
# =============== IJK ExB Cylindrical ================== #

    function IFluxFunction(force::ExB,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

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

        # Evaluate the flux through the I surface i.e. surface of constant p

        E0 = force.E0
        B0 = force.B0

        xi = 1.0 # value of x at which E0 and B0 are defined

        invTc = (Z*CONST_q) / (m*CONST_mEle) # old normalisation / (CONST_σT*CONST_c) 
        T = Characteristic.CHAR_time

        flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0.0
        else
            flux *= 1/2 * E0 * (-t0 + t1) * (-x0 + x1) * (-y0 + y1) * (-z0 + z1) * xi
            flux *= (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1))
            flux *= fluxScale
        end

        return flux

    end

    function JFluxFunction(force::ExB,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

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
        # Evaluate the flux through the J surface i.e. surface of constant u

        E0 = force.E0
        B0 = force.B0

        xi = 1.0 # value of x at which E0 and B0 are defined

        invTc = (Z*CONST_q) / (m*CONST_mEle) # old normalisation / (CONST_σT*CONST_c)
        T = Characteristic.CHAR_time
        
        flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0.0
        else
            flux *= (-t0 + t1)*(-x0 + x1)*(-y0 + y1)*(-z0 + z1) * xi
            flux *= E0 * u * sqrt(1 - u^2) * log(p1/p0)*(-sinpi(h0/pi) + sinpi(h1/pi))
            flux *= fluxScale
        end

        return flux

    end

    function KFluxFunction(force::ExB,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,Characteristic::CharacteristicStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

        m = Grids.mass_list[species_idx]
        Z = Grids.charge_list[species_idx]

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

        # Evaluate the flux through the K surface i.e. surface of constant phi

        E0 = force.E0
        B0 = force.B0

        xi = 1.0 # value of x at which E0 and B0 are defined

        invTc = (Z*CONST_q) / (m*CONST_mEle) # old normalisation / (CONST_σT*CONST_c)
        T = Characteristic.CHAR_time

        flux::Float64 = T*invTc

        if m == 0 || Z == 0
            flux *= 0.0
        else
            flux *= (-t0 + t1)*(-x0 + x1)*(-y0 + y1)*(-z0 + z1) * xi
            flux *= (1/2 * B0 * (-u0 + u1) * (-2asinh(p0/m) + 2asinh(p1/m)) + 2*E0 * acot_mod(u0) * log(p1/p0) * sinpi(phi/pi) - 2*E0 * acot_mod(u1) * log(p1/p0) * sinpi(phi/pi))
            flux *= fluxScale
        end

        return flux

    end

# ======================================================== 
