"""
    Allocate_M_Bin(PhaseSpace::PhaseSpaceStruct)

Allocates a big matrix `M_Bin` which stores the interaction rates for binary interactions between all particles in the simulation. If `n` is the size of the momentum domain, then `M_Bin` is an `n^2 x n` matrix. The size of `M_Bin` in memory is printed to the console upon allocation.
"""
function Allocate_M_Bin(PhaseSpace::PhaseSpaceStruct,loading_check::Bool)

    Momentum = PhaseSpace.Momentum

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n = sum(px_num_list.*py_num_list.*pz_num_list)
    m = n*n

    size = m*n*sizeof(Float32)

    if loading_check

        if size > 1e9
            println("M_Bin will be approx. $(size/1e9) GB in memory")
        elseif size > 1e6
            println("M_Bin will be approx. $(size/1e6) MB in memory")
        elseif size > 1e3
            println("M_Bin will be approx. $(size/1e3) KB in memory")
        else
            println("M_Bin will be approx. $size bytes in memory")
        end

        println("Do you want to proceed? (y/n)")
        proceed = readline()
        if proceed != "y"
            error("Error: User aborted binary interaction matrix loading")
        end

    else

        println("Building binary interaction matrix")
        if size > 1e9
            println("M_Bin is approx. $(size/1e9) GB in memory")
        elseif size > 1e6
            println("M_Bin is approx. $(size/1e6) MB in memory")
        elseif size > 1e3
            println("M_Bin is approx. $(size/1e3) KB in memory")
        else
            println("M_Bin is approx. $size bytes in memory")
        end

    end

    M_Bin::AbstractArray{Float32,2} = zeros(Float32,m,n)

    return M_Bin

end


"""
    Fill_M_Bin!(M_Bin::AbstractArray{Float32,2},interaction::Vector{String},Lists::ListStruct;SMatrix3=nothing,SMatrix4=nothing,TMatrix1=nothing,TMatrix2=nothing)

Fills the big matrix `M_Bin` with the interaction rates for binary interactions between all particles in the simulation.
"""
function Fill_M_Bin!(M_Bin::AbstractArray{Float32,2},interaction::BinaryStruct,PhaseSpace::PhaseSpaceStruct;SMatrix3=nothing,SMatrix4=nothing,TMatrix1=nothing,TMatrix2=nothing)

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Mode = Momentum.momentum_coordinates.Mode

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    name1 = interaction.name1
    name2 = interaction.name2
    name3 = interaction.name3
    name4 = interaction.name4

    # first find offsets of the different species in the big matrix A

    offset = zeros(Int64,size(name_list)[1])
    for i in eachindex(offset)
        if i == 1
            offset[i] = 0
        else
            offset[i] = offset[i-1]+px_num_list[i-1]*py_num_list[i-1]*pz_num_list[i-1]
        end
    end

    # interaction = ["name1","name2","name3","name4"]
    name1_loc = findfirst(==(name1),name_list)
    name2_loc = findfirst(==(name2),name_list)
    name3_loc = findfirst(==(name3),name_list)
    name4_loc = findfirst(==(name4),name_list)

    if name1 == name2 && name3 == name4

        if typeof(Mode) == Axi

            SMatrix_to_M_Bin_Axi!(M_Bin,SMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc])
            TMatrix_to_M_Bin_Axi!(M_Bin,TMatrix1,offset[name1_loc],offset[name2_loc])

        elseif typeof(Mode) == Iso

            Grids = PhaseSpace.Grids
            dpy1 = Grids.dpy[name1_loc]
            dpy2 = Grids.dpy[name2_loc]
            dpy3 = Grids.dpy[name3_loc]

            SMatrix_to_M_Bin_Iso!(M_Bin,SMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],dpy3,dpy1,dpy2)
            TMatrix_to_M_Bin_Iso!(M_Bin,TMatrix1,offset[name1_loc],offset[name2_loc],dpy1,dpy2)

        else
            error("Error: Momentum mode not recognised")
        end

        SMatrix3 = nothing
        TMatrix1 = nothing

        GC.gc()

    end

    if name1 == name2 && name3 != name4

        if typeof(Mode) == AxiT

            SMatrix_to_M_Bin_Axi!(M_Bin,SMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc])
            SMatrix_to_M_Bin_Axi!(M_Bin,SMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc])
            TMatrix_to_M_Bin_Axi!(M_Bin,TMatrix1,offset[name1_loc],offset[name2_loc])

        elseif typeof(Mode) == Iso

            Grids = PhaseSpace.Grids
            dpy1 = Grids.dpy[name1_loc]
            dpy2 = Grids.dpy[name2_loc]
            dpy3 = Grids.dpy[name3_loc]
            dpy4 = Grids.dpy[name4_loc]

            SMatrix_to_M_Bin_Iso!(M_Bin,SMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],dpy1,dpy2,dpy3)
            SMatrix_to_M_Bin_Iso!(M_Bin,SMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc],dpy4,dpy1,dpy2)
            TMatrix_to_M_Bin_Iso!(M_Bin,TMatrix1,offset[name1_loc],offset[name2_loc],dpy1,dpy2)

        end

        SMatrix3 = nothing
        SMatrix4 = nothing
        TMatrix1 = nothing

        GC.gc()

    end

    if name1 != name2 && name3 == name4

        if typeof(Mode) == Axi

            SMatrix_to_M_Bin_Axi!(M_Bin,SMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc])
            TMatrix_to_M_Bin_Axi!(M_Bin,TMatrix1,offset[name1_loc],offset[name2_loc])
            TMatrix_to_M_Bin_Axi!(M_Bin,TMatrix2,offset[name2_loc],offset[name1_loc])

        elseif typeof(Mode) == Iso

            Grids = PhaseSpace.Grids
            dpy1 = Grids.dpy[name1_loc]
            dpy2 = Grids.dpy[name2_loc]
            dpy3 = Grids.dpy[name3_loc]

            SMatrix_to_M_Bin_Iso!(M_Bin,SMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],dpy3,dpy1,dpy2)
            TMatrix_to_M_Bin_Iso!(M_Bin,TMatrix1,offset[name1_loc],offset[name2_loc],dpy1,dpy2)
            TMatrix_to_M_Bin_Iso!(M_Bin,TMatrix2,offset[name2_loc],offset[name1_loc],dpy2,dpy1)

        end

        SMatrix3 = nothing
        TMatrix1 = nothing
        TMatrix2 = nothing

        GC.gc()

    end

    if name1 != name2 && name3 != name4

        if typeof(Mode) == Axi

            SMatrix_to_M_Bin_Axi!(M_Bin,SMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc])
            SMatrix_to_M_Bin_Axi!(M_Bin,SMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc])
            TMatrix_to_M_Bin_Axi!(M_Bin,TMatrix1,offset[name1_loc],offset[name2_loc])
            TMatrix_to_M_Bin_Axi!(M_Bin,TMatrix2,offset[name2_loc],offset[name1_loc])

        elseif typeof(Mode) == Iso

            Grids = PhaseSpace.Grids
            dpy1 = Grids.dpy[name1_loc]
            dpy2 = Grids.dpy[name2_loc]
            dpy3 = Grids.dpy[name3_loc]
            dpy4 = Grids.dpy[name4_loc]

            SMatrix_to_M_Iso!(M_Bin,SMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],dpy3,dpy1,dpy2)
            SMatrix_to_M_Bin!(M_Bin,SMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc],dpy4,dpy1,dpy2)
            TMatrix_to_M_Bin!(M_Bin,TMatrix1,offset[name1_loc],offset[name2_loc],dpy1,dpy2)
            TMatrix_to_M_Bin!(M_Bin,TMatrix2,offset[name2_loc],offset[name1_loc],dpy2,dpy1)

        end

        SMatrix3 = nothing
        SMatrix4 = nothing
        TMatrix1 = nothing
        TMatrix2 = nothing

        GC.gc()

    end

end

"""
    Allocate_M_Emi(Momentum::MomentumStruct)

Allocates a big matrix `M_Emi` which stores the interaction rates for emission interactions between all particles in the simulation.
"""
function Allocate_M_Emi(PhaseSpace::PhaseSpaceStruct,loading_check::Bool)

    Momentum = PhaseSpace.Momentum

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n = sum(px_num_list.*py_num_list.*pz_num_list)

    size = n*n*sizeof(Float32)

    if loading_check

        if size > 1e9
            println("M_Emi will be approx. $(size/1e9) GB in memory")
        elseif size > 1e6
            println("M_Emi will be approx. $(size/1e6) MB in memory")
        elseif size > 1e3
            println("M_Emi will be approx. $(size/1e3) KB in memory")
        else
            println("M_Emi will be approx. $size bytes in memory")
        end

        println("Do you want to proceed? (y/n)")
        #proceed = readline()
        proceed = "y" ## FIX LATER
        if proceed != "y"
            error("Error: User aborted emission interaction matrix loading")
        end

    else

        println("Building emission interaction matrix")
        if size > 1e9
            println("M_Emi is approx. $(size/1e9) GB in memory")
        elseif size > 1e6
            println("M_Emi is approx. $(size/1e6) MB in memory")
        elseif size > 1e3
            println("M_Emi is approx. $(size/1e3) KB in memory")
        else
            println("M_Emi is approx. $size bytes in memory")
        end

    end

    M_Emi::AbstractArray{Float32,2} = zeros(Float32,n,n)

    return M_Emi

end

function Fill_M_Emi!(M_Emi::AbstractArray{Float32,2},interaction::EmiStruct,PhaseSpace::PhaseSpaceStruct;SMatrix2=nothing,SMatrix3=nothing,TMatrix1=nothing)

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    mode = interaction.mode
  
    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    name1 = interaction.name1
    name2 = interaction.name2
    name3 = interaction.name3

    # first find offsets of the different species in the big matrix A

    offset = zeros(Int64,size(name_list)[1])
    for i in eachindex(offset)
        if i == 1
            offset[i] = 0
        else
            offset[i] = offset[i-1]+px_num_list[i-1]*py_num_list[i-1]*pz_num_list[i-1]
        end
    end

    name1_loc = findfirst(==(name1),name_list)
    name2_loc = findfirst(==(name2),name_list)
    name3_loc = findfirst(==(name3),name_list)

    # interaction is name1 -> name2 + name3
    # absorption of name1 and emission of name2 not implemented with name1 == name2

    if typeof(mode) == Axi

        #SMatrix_to_M_Emi_Axi!(M_Emi,SMatrix2,offset[name2_loc],offset[name1_loc])
        SMatrix_to_M_Emi_Axi!(M_Emi,SMatrix3,offset[name3_loc],offset[name1_loc])
        #TMatrix_to_M_Emi_Axi!(M_Emi,TMatrix1,offset[name1_loc])

    elseif typeof(mode) == Iso

        Grids = PhaseSpace.Grids
        dpy1 = Grids.dpy_list[name1_loc]
        dpy2 = Grids.dpy_list[name2_loc]
        dpy3 = Grids.dpy_list[name3_loc]
        
        #SMatrix_to_M_Emi_Iso!(M_Emi,SMatrix2,offset[name2_loc],offset[name1_loc],dpy2,dpy1)
        SMatrix_to_M_Emi_Iso!(M_Emi,SMatrix3,offset[name3_loc],offset[name1_loc],dpy3,dpy1)
        #TMatrix_to_M_Emi_Iso!(M_Emi,TMatrix1,offset[name1_loc],dpy1)

    end
    #SMatrix2 = nothing
    SMatrix3 = nothing
    #TMatrix1 = nothing

    GC.gc()

end

function SMatrix_to_M_Bin_Axi!(M_Bin::AbstractArray{Float32,2},SMatrix::Array{Float64,6},offset3::Int64,offset1::Int64,offset2::Int64)

    px3_num = size(SMatrix,1)-1 # ignore overflow bin
    py3_num = size(SMatrix,2)
    px1_num = size(SMatrix,3)  
    py1_num = size(SMatrix,4)
    px2_num = size(SMatrix,5)
    py2_num = size(SMatrix,6)

    # sanity check 
    if px1_num != px2_num || py1_num != py2_num
        error("Error: SMatrix dimensions not as expected")
    end

    N = size(M_Bin,2)

    for n in axes(SMatrix,6), m in axes(SMatrix,5), l in axes(SMatrix,4), k in axes(SMatrix,3), j in axes(SMatrix,2), i in 1:px3_num

        a = (j-1)*px3_num+i+offset3
        b = (l-1)*px1_num+k+offset1
        c = (n-1)*px2_num+m+offset2
        
        # M_Bin terms allocated symmetrically
        M_Bin[(b-1)*N+(a-1)+1,c] += SMatrix[i,j,k,l,m,n]/2
        M_Bin[(c-1)*N+(a-1)+1,b] += SMatrix[i,j,k,l,m,n]/2 

    end

end

function SMatrix_to_M_Bin_Iso!(M_Bin::AbstractArray{Float32,2},SMatrix::Array{Float64,6},offset3::Int64,offset1::Int64,offset2::Int64,dpy3::Vector{Float64},dpy1::Vector{Float64},dpy2::Vector{Float64})

    px3_num = size(SMatrix,1)-1 # ignore overflow bin
    py3_num = size(SMatrix,2)
    px1_num = size(SMatrix,3)  
    py1_num = size(SMatrix,4)
    px2_num = size(SMatrix,5)
    py2_num = size(SMatrix,6)

    pxr1_grid

    # sanity check 
    if px1_num != px2_num || py1_num != py2_num
        error("Error: SMatrix dimensions not as expected")
    end

    N = size(M_Bin,2)

    for n in axes(SMatrix,6), m in axes(SMatrix,5), l in axes(SMatrix,4), k in axes(SMatrix,3), j in axes(SMatrix,2), i in 1:px3_num

        val = 0.0 # py (pz) averaged SMatrix term 
        for p in axes(SMatrix,6), q in axes(SMatrix,4), r in axes(SMatrix,2)
            val += SMatrix[i,r,k,q,m,p] * dpy3[r] * dpy1[q] * dpy2[p] # check order
        end
        val /= sum(dpy3)*sum(dpy1)*sum(dpy2) 

        a = (j-1)*px3_num+i+offset3
        b = (l-1)*px1_num+k+offset1
        c = (n-1)*px2_num+m+offset2
        
        # M_Bin terms allocated symmetrically
        M_Bin[(b-1)*N+(a-1)+1,c] += val/2
        M_Bin[(c-1)*N+(a-1)+1,b] += val/2 

    end

end

function TMatrix_to_M_Bin_Axi!(M_Bin::Array{Float32},TMatrix::Array{Float64,4},offset1::Int64,offset2::Int64)

    p1_num = size(TMatrix,1) # ignore overflow bin
    u1_num = size(TMatrix,2)
    p2_num = size(TMatrix,3)  
    u2_num = size(TMatrix,4)

    # sanity check 
    if p1_num != p2_num || u1_num != u2_num
        error("Error: SMatrix dimensions not as expected")
    end

    N = size(M_Bin,2)

    for l in axes(TMatrix,4), k in axes(TMatrix,3), j in axes(TMatrix,2), i in axes(TMatrix,1)

        a = (j-1)*p1_num+i+offset1
        b = a
        c = (l-1)*p2_num+k+offset2
    
        # M_Bin terms allocated symmetrically
        M_Bin[(b-1)*N+(a-1)+1,c] -= TMatrix[i,j,k,l]/2
        M_Bin[(c-1)*N+(a-1)+1,b] -= TMatrix[i,j,k,l]/2

    end

end

function TMatrix_to_M_Bin_Iso!(M_Bin::Array{Float32},TMatrix::Array{Float64,4},offset1::Int64,offset2::Int64,dpy1,dpy2)

    p1_num = size(TMatrix,1) # ignore overflow bin
    u1_num = size(TMatrix,2)
    p2_num = size(TMatrix,3)  
    u2_num = size(TMatrix,4)

    # sanity check 
    if p1_num != p2_num || u1_num != u2_num
        error("Error: SMatrix dimensions not as expected")
    end

    N = size(M_Bin,2)

    for l in axes(TMatrix,4), k in axes(TMatrix,3), j in axes(TMatrix,2), i in axes(TMatrix,1)

        val = 0.0 # py (pz) averaged TMatrix term 
        for p in axes(TMatrix,4), q in axes(TMatrix,2)
            val += TMatrix[i,q,k,p] * dpy1[q] * dpy2[p] # check order
        end
        val /= sum(dpy1)*sum(dpy2) 

        a = (j-1)*p1_num+i+offset1
        b = a
        c = (l-1)*p2_num+k+offset2
    
        # M_Bin terms allocated symmetrically
        M_Bin[(b-1)*N+(a-1)+1,c] -= val/2
        M_Bin[(c-1)*N+(a-1)+1,b] -= val/2

    end

end

function SMatrix_to_M_Emi_Axi!(M_Emi::Array{Float32},SMatrix::Array{Float64,4},offset2::Int64,offset1::Int64)

    px2_num = size(SMatrix,1)#-1 # ignore overflow bin
    py2_num = size(SMatrix,2)
    px1_num = size(SMatrix,3)  
    py1_num = size(SMatrix,4)

    for l in axes(SMatrix,4), k in axes(SMatrix,3), j in axes(SMatrix,2), i in 1:px2_num

        a = (j-1)*px2_num+i+offset2
        b = (l-1)*px1_num+k+offset1

        M_Emi[a,b] += SMatrix[i,j,k,l]

    end

end

function SMatrix_to_M_Emi_Iso!(M_Emi::Array{Float32},SMatrix::Array{Float64,4},offset2::Int64,offset1::Int64,dpy2,dpy1)

    px2_num = size(SMatrix,1)#-1 # ignore overflow bin
    py2_num = size(SMatrix,2)
    px1_num = size(SMatrix,3)  
    py1_num = size(SMatrix,4)

    for l in axes(SMatrix,4), k in axes(SMatrix,3), j in axes(SMatrix,2), i in 1:px2_num

        val = 0.0 # py (pz) averaged SMatrix term
        for p in axes(SMatrix,4), q in axes(SMatrix,2)
            val += SMatrix[i,q,k,p] * dpy2[q] * dpy1[p] # check order
        end
        val /= sum(dpy2)*sum(dpy1)

        a = (j-1)*px2_num+i+offset2
        b = (l-1)*px1_num+k+offset1

        M_Emi[a,b] += val

    end

end

function TMatrix_to_M_Emi_Axi!(M_Emi::Array{Float32},TMatrix::Array{Float64,2},offset1::Int64)

    px1_num = size(TMatrix,1)
    py1_num = size(TMatrix,2)

    for j in axes(TMatrix,2), i in axes(TMatrix,1)

        a = (j-1)*px1_num+i+offset1
        b = a

        M_Emi[a,b] -= TMatrix[i,j]
        
    end

end

function TMatrix_to_M_Emi_Iso!(M_Emi::Array{Float32},TMatrix::Array{Float64,2},offset1::Int64)

    px1_num = size(TMatrix,1)
    py1_num = size(TMatrix,2)

    for j in axes(TMatrix,2), i in axes(TMatrix,1)

        val = 0.0 # py (pz) averaged TMatrix term
        for p in axes(TMatrix,2)
            val += TMatrix[i,p] * dpy1[p] # check order
        end
        val /= sum(dpy1)

        a = (j-1)*px1_num+i+offset1
        b = a

        M_Emi[a,b] -= val
        
    end

end