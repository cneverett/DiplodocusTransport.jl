# Functions for the force terms due to Synchrotron Radiation Reaction with a constant magnetic field along the axis of symmetry

function Force_p_Sync(p,u,m)
    #three force normalised by σT*c and (Z^4 B^2)/(m^3*mu0)

    Fp = -(p^3/m^2+p) * (1-u^2) 

    return Fp
end

function Force_u_Sync(p,u)
    #three force normalised by σT*c and (Z^4 B^2)/(m^3*mu0)

    Fu = -p * u * sqrt(1-u^2) 

    return Fu
end

# integrated functions

    function Force_p_Sync_du(p,u,m)
        #three force normalised by σT*c  and (Z^4 B^2)/(m^3*mu0) - integrated over u

        Fpdu = -(1/3) * (p^3/m^2+p) * (3*u-u^3) 

        return Fpdu
    end

    function Force_u_Sync_dp(p,u)
        # three-force normalised by σT*c and (Z^4 B^2)/(m^3*mu0) times by sqrt(1-u^2)/p integrated over p

        Fudp =  p * (-u+u^3) 

        return Fudp

    end


# Matrix of the forces due to synchrotron radiation reaction
function Force_Sync_Matrix(pr,ur,B,m,Z)

    scale = Z^4*B^2/(m^3*getfield(BCI,Symbol(μ0)))

    Fp_array = [Force_p_Sync(p,u,m) for p in pr, u in ur] 
    Fu_array = [Force_u_Sync(p,u) for p in pr, u in ur]

    Force_array = zeros(Float64,length(pr)-1,length(ur)-1,length(pr)-1,length(ur)-1)

    for i in axes(Force_array,1), j in axes(Force_array,2), k in axes(Force_array,3), l in axes(Force_array,4)
        if k == i && l == j
            Force_array[i,j,k,l] += Fu_array[i+1,j+1] - Fu_array[i,j+1] - Fu_array[i+1,j] + Fu_array[i,j] 
            Force_array[i,j,k,l] += Fp_array[i+1,j+1] - Fp_array[i+1,j] - Fp_array[i,j+1] + Fp_array[i,j] 
        elseif k == i && l+1 == j
            Force_array[i,j,k,l] += Fu_array[i+1,j+1] - Fu_array[i,j+1]
        elseif k == i && l-1 == j
            Force_array[i,j,k,l] -= Fu_array[i+1,j] - Fu_array[i,j] 
        elseif k+1 == i && l == j
            Force_array[i,j,k,l] += Fp_array[i+1,j+1] - Fp_array[i+1,j] 
        elseif k-1 == i && l == j
            Force_array[i,j,k,l] -= Fp_array[i,j+1] - Fp_array[i,j]
        end
        # normalise for multiplication by gij rather than fij
        Force_array[i,j,k,l] /= (pr[k+1]-pr[k])*(ur[l+1]-ur[l]) 
    end

    Force_array .*= scale/2

    Force_array_2D = reshape(Force_array,(size(Force_array,1)*size(Force_array,2),size(Force_array,3)*size(Force_array,4)))

    return Force_array_2D  

end
