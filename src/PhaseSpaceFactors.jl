function PhaseSpaceFactors_Emi_Undo!(SMatrix::Array{Float64,4},pxr2::Vector{Float64},pyr2::Vector{Float64},pxr1::Vector{Float64},pyr1::Vector{Float64})

    for ii in axes(SMatrix,4), jj in axes(SMatrix,3), kk in axes(SMatrix,2), ll in axes(SMatrix,1)
        SMatrix[ll,kk,jj,ii] /= (pyr1[ii+1]-pyr1[ii])*(pxr1[jj+1]-pxr1[jj]) #dp2dmu2
        SMatrix[ll,kk,jj,ii] *= (pyr2[kk+1]-pyr2[kk])*(pxr2[ll+1]-pxr2[ll])
    end

end # function

function PhaseSpaceFactors_Binary_Undo!(pxr3::Vector{Float64},pyr3::Vector{Float64},pxr4::Vector{Float64},pyr4::Vector{Float64},pxr1::Vector{Float64},pyr1::Vector{Float64},pxr2::Vector{Float64},pyr2::Vector{Float64};SMatrix3=0f0,SMatrix4=0f0,TMatrix1=0f0,TMatrix2=0f0)

    # Function that applies the correct phase space factors to SMatrix and TMatrix derived from Stotal and Ttotal arrays such that the symmetries can be applied.

    if SMatrix3 == 0f0 && SMatrix4 == 0f0 && TMatrix1 == 0f0 && TMatrix2 == 0f0
        error("No Matrices to apply phase space factors to.")
    end

    # === SMatrix3 === #
    # Momentum space volume elements
    if SMatrix3 != 0f0
        for ii in axes(SMatrix3,6), jj in axes(SMatrix3,5), kk in axes(SMatrix3,4), ll in axes(SMatrix3,3)
            for mm in axes(SMatrix3,2), nn in 1:size(SMatrix3,1)-1
                SMatrix3[nn,mm,ll,kk,jj,ii] *= (pyr3[mm+1]-pyr3[mm])*(pxr3[nn+1]-pxr3[nn]) # dmu3
                SMatrix3[nn,mm,ll,kk,jj,ii] /= (pyr1[kk+1]-pyr1[kk])*(pxr1[ll+1]-pxr1[ll]) # dp1dmu1
                SMatrix3[nn,mm,ll,kk,jj,ii] /= (pyr2[ii+1]-pyr2[ii])*(pxr2[jj+1]-pxr2[jj]) # dp2dmu2
            end
            if TMatrix1 != 0f0
                TMatrix1[ll,kk,jj,ii] /= (pyr2[ii+1]-pyr2[ii])*(pxr2[jj+1]-pxr2[jj]) # dp2dmu2
            end
            if TMatrix2 != 0f0
                TMatrix2[jj,ii,ll,kk] /= (pyr1[kk+1]-pyr1[kk])*(pxr1[ll+1]-pxr1[ll]) # dp2dmu2 
            end
        end
    end

    # === SMatrix4 === #
    # Momentum space volume elements
    if SMatrix4 != 0f0
        for ii in axes(SMatrix4,6), jj in axes(SMatrix4,5), kk in axes(SMatrix4,4), ll in axes(SMatrix4,3)
            for mm in axes(SMatrix4,2), nn in 1:size(SMatrix4,1)-1
                SMatrix4[nn,mm,ll,kk,jj,ii] *= (pyr4[mm+1]-pyr4[mm])*(pxr4[nn+1]-pxr4[nn]) # dmu3
                SMatrix4[nn,mm,ll,kk,jj,ii] /= (pyr1[kk+1]-pyr1[kk])*(pxr1[ll+1]-pxr1[ll]) # dp1dmu1
                SMatrix4[nn,mm,ll,kk,jj,ii] /= (pyr2[ii+1]-pyr2[ii])*(pxr2[jj+1]-pxr2[jj]) # dp2dmu2
            end
        end
    end

    return nothing

end
