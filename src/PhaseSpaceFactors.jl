function PhaseSpaceFactors_Sync_Undo!(SMatrix::Array{Float64,4},p1val::Vector{Float64},t1val::Vector{Float64},p2val::Vector{Float64},t2val::Vector{Float64})

    for ii in axes(SMatrix,4), jj in axes(SMatrix,3), kk in axes(SMatrix,2), ll in axes(SMatrix,1)
        SMatrix[ll,kk,jj,ii] /= (t2val[ii+1]-t2val[ii])*(p2val[jj+1]-p2val[jj]) #dp2dmu2
        SMatrix[ll,kk,jj,ii] *= (t1val[kk+1]-t2val[kk])*(p1val[ll+1]-p1val[ll])
    end

end # function

function PhaseSpaceFactors_Binary_Undo!(p3val::Vector{Float64},t3val::Vector{Float64},p4val::Vector{Float64},t4val::Vector{Float64},p1val::Vector{Float64},t1val::Vector{Float64},p2val::Vector{Float64},t2val::Vector{Float64};SMatrix3=0f0,SMatrix4=0f0,TMatrix1=0f0,TMatrix2=0f0)

    # Function that applies the correct phase space factors to SMatrix and TMatrix derived from Stotal and Ttotal arrays such that the symmetries can be applied.

    if SMatrix3 == 0f0 && SMatrix4 == 0f0 && TMatrix1 == 0f0 && TMatrix2 == 0f0
        error("No Matricies to apply phase space factors to.")
    end

    # === SMatrix3 === #
    # Momentum space volume elements
    if SMatrix3 != 0f0
        for ii in axes(SMatrix3,6), jj in axes(SMatrix3,5), kk in axes(SMatrix3,4), ll in axes(SMatrix3,3)
            for mm in axes(SMatrix3,2), nn in 1:size(SMatrix3,1)-1
                SMatrix3[nn,mm,ll,kk,jj,ii] *= (t3val[mm+1]-t3val[mm])*(p3val[nn+1]-p3val[nn]) # dmu3
                SMatrix3[nn,mm,ll,kk,jj,ii] /= (t1val[kk+1]-t1val[kk])*(p1val[ll+1]-p1val[ll]) # dp1dmu1
                SMatrix3[nn,mm,ll,kk,jj,ii] /= (t2val[ii+1]-t2val[ii])*(p2val[jj+1]-p2val[jj]) # dp2dmu2
            end
            if TMatrix1 != 0f0
                TMatrix1[ll,kk,jj,ii] /= (t2val[ii+1]-t2val[ii])*(p2val[jj+1]-p2val[jj]) # dp2dmu2
            end
            if TMatrix2 != 0f0
                TMatrix2[jj,ii,ll,kk] /= (t1val[kk+1]-t1val[kk])*(p1val[ll+1]-p1val[ll]) # dp2dmu2 
            end
        end
    end

    # === SMatrix4 === #
    # Momentum space volume elements
    if SMatrix4 != 0f0
        for ii in axes(SMatrix4,6), jj in axes(SMatrix4,5), kk in axes(SMatrix4,4), ll in axes(SMatrix4,3)
            for mm in axes(SMatrix4,2), nn in 1:size(SMatrix4,1)-1
                SMatrix4[nn,mm,ll,kk,jj,ii] *= (t4val[mm+1]-t4val[mm])*(p4val[nn+1]-p4val[nn]) # dmu3
                SMatrix4[nn,mm,ll,kk,jj,ii] /= (t1val[kk+1]-t1val[kk])*(p1val[ll+1]-p1val[ll]) # dp1dmu1
                SMatrix4[nn,mm,ll,kk,jj,ii] /= (t2val[ii+1]-t2val[ii])*(p2val[jj+1]-p2val[jj]) # dp2dmu2
            end
        end
    end

    return nothing

end
