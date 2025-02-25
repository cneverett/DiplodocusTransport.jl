function Allocate_A_Binary(Lists::ListStruct)

    # allocates space for the BIG matrix A for all binary interactions

    p_num_list = Lists.p_num_list
    u_num_list = Lists.u_num_list

    n = sum(p_num_list.*u_num_list)
    m = n*n

    A_Binary::AbstractArray{Float32,2} = zeros(Float32,m,n)

    size_A_Binary = sizeof(A_Binary)

    if size_A_Binary > 1e9
        println("A_Binary is approx. $(size_A_Binary/1e9) GB in memory")
    elseif size_A_Binary > 1e6
        println("A_Binary is approx. $(size_A_Binary/1e6) MB in memory")
    elseif size_A_Binary > 1e3
        println("A_Binary is approx. $(size_A_Binary/1e3) KB in memory")
    else
        println("A_Binary is approx. $size_A_Binary bytes in memory")
    end

    return A_Binary

end

function Fill_A_Binary!(A_Binary::AbstractArray{Float32,2},interaction::Vector{String},Lists::ListStruct;SMatrix3=nothing,SMatrix4=nothing,TMatrix1=nothing,TMatrix2=nothing)

    name_list = Lists.name_list
    p_num_list = Lists.p_num_list
    u_num_list = Lists.u_num_list

    # first find offsets of the different species in the big matrix A

    # assign elements of each matrix to the big matrix A, converting to Float32

    offset = zeros(Int64,size(name_list)[1])
    for i in eachindex(offset)
        if i == 1
            offset[i] = 0
        else
            offset[i] = offset[i-1]+p_num_list[i-1]*u_num_list[i-1]
        end
    end

    # interaction = ["name1","name2","name3","name4"]
    name1_loc = findfirst(==(interaction[1]),name_list)
    name2_loc = findfirst(==(interaction[2]),name_list)
    name3_loc = findfirst(==(interaction[3]),name_list)
    name4_loc = findfirst(==(interaction[4]),name_list)

    if interaction[1] == interaction[2] && interaction[3] == interaction[4] # name1=name2 and name3=name4

        SMatrix_to_A_Binary!(A_Binary,SMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc])
        TMatrix_to_A_Binary!(A_Binary,TMatrix1,offset[name1_loc],offset[name2_loc])

        SMatrix3 = nothing
        TMatrix1 = nothing

        GC.gc()

    end

    if interaction[1] == interaction[2] && interaction[3] != interaction[4] # name1=name2 and name3=name4

        SMatrix_to_A_Binary!(A_Binary,SMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc])
        SMatrix_to_A_Binary!(A_Binary,SMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc])
        TMatrix_to_A_Binary!(A_Binary,TMatrix1,offset[name1_loc],offset[name2_loc])

        SMatrix3 = nothing
        SMatrix4 = nothing
        TMatrix1 = nothing

        GC.gc()

    end

    if interaction[1] != interaction[2] && interaction[3] == interaction[4] # name1=name2 and name3=name4

        SMatrix_to_A_Binary!(A_Binary,SMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc])
        TMatrix_to_A_Binary!(A_Binary,TMatrix1,offset[name1_loc],offset[name2_loc])
        TMatrix_to_A_Binary!(A_Binary,TMatrix2,offset[name2_loc],offset[name1_loc])

        SMatrix3 = nothing
        TMatrix1 = nothing
        TMatrix2 = nothing

        GC.gc()

    end

    if interaction[1] != interaction[2] && interaction[3] != interaction[4] # name1=name2 and name3=name4

        SMatrix_to_A_Binary!(A_Binary,SMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc])
        SMatrix_to_A_Binary!(A_Binary,SMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc])
        TMatrix_to_A_Binary!(A_Binary,TMatrix1,offset[name1_loc],offset[name2_loc])
        TMatrix_to_A_Binary!(A_Binary,TMatrix2,offset[name2_loc],offset[name1_loc])

        SMatrix3 = nothing
        SMatrix4 = nothing
        TMatrix1 = nothing
        TMatrix2 = nothing

        GC.gc()

    end

end

function Allocate_A_Emi(Lists::ListStruct)

    # allocates space for the BIG matrix A for all binary interactions

    p_num_list = Lists.p_num_list
    u_num_list = Lists.u_num_list

    n = sum(p_num_list.*u_num_list)

    A_Emi::AbstractArray{Float32,2} = zeros(Float32,n,n)

    size_A_Emi = sizeof(A_Emi)

    if size_A_Emi > 1e9
        println("A_Emi is approx. $(size_A_Emi/1e9) GB in memory")
    elseif size_A_Emi > 1e6
        println("A_Emi is approx. $(size_A_Emi/1e6) MB in memory")
    elseif size_A_Emi > 1e3
        println("A_Emi is approx. $(size_A_Emi/1e3) KB in memory")
    else
        println("A_Emi is approx. $size_A_Emi bytes in memory")
    end

    return A_Emi

end

function Fill_A_Emi!(A_Emi::AbstractArray{Float32,2},interaction::Vector{String},Lists::ListStruct;SMatrix2=nothing,SMatrix3=nothing,TMatrix1=nothing)

    name_list = Lists.name_list
    p_num_list = Lists.p_num_list
    u_num_list = Lists.u_num_list

    # first find offsets of the different species in the big matrix A

    # assign elements of each matrix to the big matrix A, converting to Float32

    offset = zeros(Int64,size(name_list)[1])
    for i in eachindex(offset)
        if i == 1
            offset[i] = 0
        else
            offset[i] = offset[i-1]+p_num_list[i-1]*u_num_list[i-1]
        end
    end

    # interaction = ["name1","name2","name3","interaction"]
    name1_loc = findfirst(==(interaction[1]),name_list)
    name2_loc = findfirst(==(interaction[2]),name_list)
    name3_loc = findfirst(==(interaction[3]),name_list)

    # interaction is name1 -> name2 + name3
    # absorption of name1 and emission of name2 not implemented with name1 == name2

    #SMatrix_to_A_Emi!(A_Emi,SMatrix2,offset[name2_loc],offset[name1_loc])
    SMatrix_to_A_Emi!(A_Emi,SMatrix3,offset[name3_loc],offset[name1_loc])
    #TMatrix_to_A_Emi!(A_Emi,TMatrix1,offset[name1_loc],offset[name1_loc])

    #SMatrix2 = nothing
    SMatrix3 = nothing
    #TMatrix1 = nothing

    GC.gc()

end

function SMatrix_to_A_Binary!(A_Binary::AbstractArray{Float32,2},SMatrix::Array{Float64,6},offset3::Int64,offset1::Int64,offset2::Int64)

    p3_num = size(SMatrix,1)-1 # ignore overflow bin
    u3_num = size(SMatrix,2)
    p1_num = size(SMatrix,3)  
    u1_num = size(SMatrix,4)
    p2_num = size(SMatrix,5)
    u2_num = size(SMatrix,6)

    # sanity check 
    if p1_num != p2_num || u1_num != u2_num
        error("Error: SMatrix dimensions not as expected")
    end

    N = size(A_Binary,2)

    for n in axes(SMatrix,6), m in axes(SMatrix,5), l in axes(SMatrix,4), k in axes(SMatrix,3), j in axes(SMatrix,2), i in 1:p3_num

        a = (j-1)*p3_num+i+offset3
        b = (l-1)*p1_num+k+offset1
        c = (n-1)*p2_num+m+offset2

        #if (a-1)*N+(b-1)+1 > size(A_Binary,1)
        #    error(i,j,k,l,m,n,a,b,c)
        #end
        
        A_Binary[(b-1)*N+(a-1)+1,c] += SMatrix[i,j,k,l,m,n]/2
        A_Binary[(c-1)*N+(a-1)+1,b] += SMatrix[i,j,k,l,m,n]/2 # symmetrising 
        #A_Binary[a,b,c] += SMatrix[i,j,k,l,m,n]

    end

end

function TMatrix_to_A_Binary!(A_Binary::Array{Float32},TMatrix::Array{Float64,4},offset1::Int64,offset2::Int64)

    p1_num = size(TMatrix,1) # ignore overflow bin
    u1_num = size(TMatrix,2)
    p2_num = size(TMatrix,3)  
    u2_num = size(TMatrix,4)

    # sanity check 
    if p1_num != p2_num || u1_num != u2_num
        error("Error: SMatrix dimensions not as expected")
    end

    N = size(A_Binary,2)

    for l in axes(TMatrix,4), k in axes(TMatrix,3), j in axes(TMatrix,2), i in axes(TMatrix,1)

        a = (j-1)*p1_num+i+offset1
        b = a
        c = (l-1)*p2_num+k+offset2
    
        A_Binary[(b-1)*N+(a-1)+1,c] -= TMatrix[i,j,k,l]/2
        A_Binary[(c-1)*N+(a-1)+1,b] -= TMatrix[i,j,k,l]/2


        #A_Binary[a,b,c] -= TMatrix[i,j,k,l]
        
    end

end

function SMatrix_to_A_Emi!(A_Emi::Array{Float32},SMatrix::Array{Float32,4},offset2::Int64,offset1::Int64)

    p2_num = size(SMatrix,1)#-1 # ignore overflow bin
    u2_num = size(SMatrix,2)
    p1_num = size(SMatrix,3)  
    u1_num = size(SMatrix,4)

    for l in axes(SMatrix,4), k in axes(SMatrix,3), j in axes(SMatrix,2), i in 1:p2_num

        a = (j-1)*p2_num+i+offset2
        b = (l-1)*p1_num+k+offset1

        A_Emi[a,b] += SMatrix[i,j,k,l]

    end

end

function TMatrix_to_A_Emi!(A_Emi::Array{Float32},TMatrix::Array{Float64,2},offset1::Int64)

    p1_num = size(TMatrix,1)
    u1_num = size(TMatrix,2)

    for j in axes(TMatrix,2), i in axes(TMatrix,1)

        a = (j-1)*p1_num+i+offset1
        b = a

        A_Emi[a,b] -= TMatrix[i,j]
        
    end

end