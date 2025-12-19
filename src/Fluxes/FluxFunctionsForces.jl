# ======= IJK Sync Rad Reaction Cylindrical ============== #

    function IFluxFunction(force::SyncRadReact,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        mode = force.mode
        B = force.B

        fluxScale = (Z^4*B^2)/(CONST_μ0*m^2*CONST_mEle*CONST_c^2) # normalised by σT*c
        
        flux::Float64 = fluxScale

        if m == 0 || Z == 0
            flux = 0e0
        else
            if typeof(mode) == Ani
                flux *= 1/(6*m^2)
                flux *= p*sqrt(m^2 + p^2) * (-3*u0 + u0^3 + 3*u1 - u1^3) * (h0 - h1)
                flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
            elseif typeof(mode) == Axi # average over azimuthal angles
                flux *= 1/(6*m^2)
                flux *= p*sqrt(m^2 + p^2) * (-3*u0 + u0^3 + 3*u1 - u1^3) * (h0 - h1)
                flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
            elseif typeof(mode) == Iso # averaged force over all angles
                flux *= -1/(3*m^2)
                flux *= p*sqrt(m^2 + p^2) * (h0 - h1) * (u0 - u1)
                flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
            else
                error("Synchrotron mode not recognised.")
            end
        end

        return flux

    end

    function JFluxFunction(force::SyncRadReact,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        mode = force.mode
        B = force.B

        fluxScale = (Z^4*B^2)/(CONST_μ0*m^2*CONST_mEle*CONST_c^2) # normalised by σT*c

        flux::Float64 = fluxScale

        if m == 0 || Z == 0
            flux = 0e0
        else
            if typeof(mode) == Ani
                flux *= 1/4
                flux *= u * (-1 + u^2) * (h0 - h1) * (-2asinh(p0/m)+2asinh(p1/m))
                flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
            elseif typeof(mode) == Axi
                flux *= 1/4
                flux *= u * (-1 + u^2) * (h0 - h1) * (-2asinh(p0/m)+2asinh(p1/m))
                flux *= (t0 - t1) * (x0^2 - x1^2) * (y0 - y1) * (z0 - z1)
            elseif typeof(mode) == Iso
                flux *= 0e0
            end
        end

        return flux

    end

    function KFluxFunction(force::SyncRadReact,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        mode = force.mode
        B = force.B

        fluxScale = (Z^4*B^2)/(CONST_μ0*m^2*CONST_mEle*CONST_c^2) # normalised by σT*c

        flux::Float64 = fluxScale

        if m == 0 || Z == 0
            flux = 0e0
        else
            flux *= 0e0
        end

        return flux

    end

# ======================================================== #

# ======= IJK Sync Rad Reaction Spherical ============== #

    function IFluxFunction(force::SyncRadReact,space_coords::Spherical,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)


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

        mode = force.mode
        B = force.B

        fluxScale = (Z^4*B^2)/(CONST_μ0*m^2*CONST_mEle*CONST_c^2) # normalised by σT*c

        flux::Float64 = fluxScale

        if m == 0 || Z == 0
            flux = 0e0
        else
            if typeof(mode) == Ani
                flux *= -1/(9*m^2)
                flux *= p*sqrt(m^2 + p^2) * (-3*u0 + u0^3 + 3*u1 - u1^3) * (h0 - h1)
                flux *= (t0 - t1) * (x0^3 - x1^3) * (cospi(y0/pi) - cospi(y1/pi)) * (z0 - z1)
            elseif typeof(mode) == Axi # average over azimuthal angles
                flux *= -1/(9*m^2)
                flux *= p*sqrt(m^2 + p^2) * (-3*u0 + u0^3 + 3*u1 - u1^3) * (h0 - h1)
                flux *= (t0 - t1) * (x0^3 - x1^3) * (cospi(y0/pi) - cospi(y1/pi)) * (z0 - z1)
            elseif typeof(mode) == Iso # averaged force over all angles
                flux *= 2/(9*m^2)
                flux *= p*sqrt(m^2 + p^2) * (h0 - h1) * (u0 - u1)
                flux *= (t0 - t1) * (x0^3 - x1^3) * (cospi(y0/pi) - cospi(y1/pi)) * (z0 - z1)
            else
                error("Synchrotron mode not recognised.")
            end
        end

        return flux

    end

    function JFluxFunction(force::SyncRadReact,space_coords::Spherical,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        mode = force.mode
        B = force.B

        fluxScale = (Z^4*B^2)/(CONST_μ0*m^2*CONST_mEle*CONST_c^2) # normalised by σT*c

        flux::Float64 = fluxScale

        if m == 0 || Z == 0
            flux = 0e0
        else
            if typeof(mode) == Ani
                flux *= -1/6
                flux *= u * (-1 + u^2) * (h0 - h1) * (-2asinh(p0/m)+2asinh(p1/m))
                flux *= (t0 - t1) * (x0^3 - x1^3) * (cospi(y0/pi) - cospi(y1/pi)) * (z0 - z1)
            elseif typeof(mode) == Axi
                flux *= -1/6
                flux *= u * (-1 + u^2) * (h0 - h1) * (-2asinh(p0/m)+2asinh(p1/m))
                flux *= (t0 - t1) * (x0^3 - x1^3) * (cospi(y0/pi) - cospi(y1/pi)) * (z0 - z1)
            elseif typeof(mode) == Iso
                flux *= 0e0
            end
        end

        return flux

    end

    function KFluxFunction(force::SyncRadReact,space_coords::Spherical,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        mode = force.mode
        B = force.B

        fluxScale = (Z^4*B^2)/(CONST_μ0*m^2*CONST_mEle*CONST_c^2) # normalised by σT*c

        flux::Float64 = fluxScale

        if m == 0 || Z == 0
            flux = 0e0
        else
            flux *= 0e0
        end

        return flux

    end

# ======================================================== #

# =============== IJK ExB Cylindrical ================== #

    function IFluxFunction(force::ExB,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        fluxScale = (Z*CONST_q) / (m*CONST_mEle) / (CONST_σT*CONST_c) 

        if m == 0 || Z == 0
            return flux = 0f0
        else
            flux = 1/2 * E0 * (-t0 + t1) * (-x0 + x1) * (-y0 + y1) * (-z0 + z1)
            flux *= (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1))
            flux *= fluxScale
        end

        return flux

    end

    function JFluxFunction(force::ExB,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        fluxScale = (Z*CONST_q) / (m*CONST_mEle) / (CONST_σT*CONST_c) 

        if m == 0 || Z == 0
            return flux = 0f0
        else
            flux = (-t0 + t1)*(-x0 + x1)*(-y0 + y1)*(-z0 + z1) 
            flux *= E0 * u * sqrt(1 - u^2) * log(p1/p0)*(-sinpi(h0/pi) + sinpi(h1/pi))
            flux *= fluxScale
        end

        return flux

    end

    function KFluxFunction(force::ExB,space_coords::Cylindrical,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        fluxScale = (Z*CONST_q) / (m*CONST_mEle) / (CONST_σT*CONST_c) 

        if m == 0 || Z == 0
            return flux = 0f0
        else
            flux = (-t0 + t1)*(-x0 + x1)*(-y0 + y1)*(-z0 + z1) 
            flux *= (1/2 * B0 * (-u0 + u1) * (-2asinh(p0/m) + 2asinh(p1/m)) + 2*E0 * acot_mod(u0) * log(p1/p0) * sinpi(phi/pi) - 2*E0 * acot_mod(u1) * log(p1/p0) * sinpi(phi/pi))
            flux *= fluxScale
        end

        return flux

    end

# ======================================================== #

# ================= IJK ExB Cartesian ==================== #

    function IFluxFunction(force::ExB,space_coords::Cartesian,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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
            p = pxr_list[species_idx][px_idx+1]
        elseif plus_minus == "minus"
            p = pxr_list[species_idx][px_idx]
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

        fluxScale = (Z*CONST_q) / (m*CONST_mEle) / (CONST_σT*CONST_c) 

        flux::Float64 = fluxScale

        if m == 0 || Z == 0
            flux *= 0e0
        else
            flux *= (-t0 + t1) * (-x0 + x1) * (-y0 + y1) * (-z0 + z1)  
            flux *= -(1/2) * E0* (u0*sqrt(1 - u0^2) - u1*sqrt(1 - u1^2) - 2acot_mod(u0) + 2acot_mod(u1)) * (cospi(h0/pi) - cospi(h1/pi))
        end

        return flux

    end

    function JFluxFunction(force::ExB,space_coords::Cartesian,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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
        h0 = pzr_list[species_idx][pz_idx]
        h1 = pzr_list[species_idx][pz_idx+1]

        # Evaluate the flux through the J surface i.e. surface of constant u

        E0 = force.E0
        B0 = force.B0

        fluxScale = (Z*CONST_q) / (m*CONST_mEle) / (CONST_σT*CONST_c) 

        flux::Float64 = fluxScale

        if m == 0 || Z == 0
            flux *= 0e0
        else
            flux *= (-t0 + t1) * (-x0 + x1) * (-y0 + y1) * (-z0 + z1) 
            flux *= E0 * log(p1/p0) * u*sqrt(1 - u^2) * (cospi(h0/pi) - cospi(h1/pi))
        end

        return flux

    end

    function KFluxFunction(force::ExB,space_coords::Cartesian,momentum_coords::Spherical,Grids::GridsStruct,species_idx::Int64,plus_minus::String,t_idx::Int64,x_idx::Int64,y_idx::Int64,z_idx::Int64,px_idx::Int64,py_idx::Int64,pz_idx::Int64)

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

        fluxScale = (Z*CONST_q) / (m*CONST_mEle) / (CONST_σT*CONST_c) 

        flux::Float64 = fluxScale

        if m == 0 || Z == 0
            flux *= 0e0
        else
            flux *= (-t0 + t1) * (-x0 + x1) * (-y0 + y1) *  (-z0 + z1) 
            flux *= -E0 * (-2acot_mod(u0)+2acot_mod(u1)) * cospi(phi/pi) * log(p1/p0) + B0 * (u0 - u1) * (-asinh(p0/m) + asinh(p1/m)) # sign of E0 flipped to make sure direction is correct ?????
        end

        return flux

    end

# ======================================================== #