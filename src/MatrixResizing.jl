function SixDtoThreeD(matrix::Array{Float32,6})

    # Flattens a 6D matrix to a 3D matrix where S[i,j,k,l,m,n] -> S[(j-1)*nump+i,(l-1)*nump+k,(n-1)*nump+m]

    for n in axes(Matrix,6), m in axes(Matrix,5), l in axes(Matrix,4), k in axes(Matrix,3), j in axes(Matrix,2), i in axes(Matrix,1)
        matrix3D[(j-1)*nump+i,(l-1)*nump+k,(n-1)*nump+m] = matrix[i,j,k,l,m,n]
    end

    return matrix3D

end

function FourDtoTwoD(matrix::Array{Float32,4})

    # Flattens a 4D matrix to a 2D matrix where T[i,j,k,l] -> T[(j-1)*nump+i,(l-1)*nump+k]

    for l in axes(Matrix,4), k in axes(Matrix,3), j in axes(Matrix,2), i in axes(Matrix,1)
        matrix2D[(j-1)*nump+i,(l-1)*nump+k] = matrix[i,j,k,l]
    end

    return matrix2D
    
end