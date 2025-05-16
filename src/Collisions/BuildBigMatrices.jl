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
        #proceed = readline()
        #if proceed != "y"
        #    error("Error: User aborted binary interaction matrix loading")
        #end

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
    Fill_M_Bin!(M_Bin::AbstractArray{Float32,2},interaction::Vector{String},Lists::ListStruct;GainMatrix3=nothing,GainMatrix4=nothing,LossMatrix1=nothing,LossMatrix2=nothing)

Fills the big matrix `M_Bin` with the interaction rates for binary interactions between all particles in the simulation.
"""
function Fill_M_Bin!(M_Bin::AbstractArray{Float32,2},interaction::BinaryStruct,PhaseSpace::PhaseSpaceStruct,GainMatrix3::Array{Float64,9},GainMatrix4::Array{Float64,9},LossMatrix1::Array{Float64,6},LossMatrix2::Array{Float64,6})

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    #Mode = Momentum.momentum_coordinates.Mode

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

            GainMatrix_to_M_Bin_Axi!(M_Bin,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc])
            LossMatrix_to_M_Bin_Axi!(M_Bin,LossMatrix1,offset[name1_loc],offset[name2_loc])

        elseif typeof(Mode) == Iso

            Grids = PhaseSpace.Grids
            dpy1 = Grids.dpy[name1_loc]
            dpy2 = Grids.dpy[name2_loc]
            dpy3 = Grids.dpy[name3_loc]

            GainMatrix_to_M_Bin_Iso!(M_Bin,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],dpy3,dpy1,dpy2)
            LossMatrix_to_M_Bin_Iso!(M_Bin,LossMatrix1,offset[name1_loc],offset[name2_loc],dpy1,dpy2)

        else
            error("Error: Momentum mode not recognised")
        end

        GainMatrix3 = nothing
        LossMatrix1 = nothing

        GC.gc()

    end

    if name1 == name2 && name3 != name4

        if typeof(Mode) == Axi

            GainMatrix_to_M_Bin_Axi!(M_Bin,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc])
            GainMatrix_to_M_Bin_Axi!(M_Bin,GainMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc])
            LossMatrix_to_M_Bin_Axi!(M_Bin,LossMatrix1,offset[name1_loc],offset[name2_loc])

        elseif typeof(Mode) == Iso

            Grids = PhaseSpace.Grids
            dpy1 = Grids.dpy[name1_loc]
            dpy2 = Grids.dpy[name2_loc]
            dpy3 = Grids.dpy[name3_loc]
            dpy4 = Grids.dpy[name4_loc]

            GainMatrix_to_M_Bin_Iso!(M_Bin,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],dpy1,dpy2,dpy3)
            GainMatrix_to_M_Bin_Iso!(M_Bin,GainMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc],dpy4,dpy1,dpy2)
            LossMatrix_to_M_Bin_Iso!(M_Bin,LossMatrix1,offset[name1_loc],offset[name2_loc],dpy1,dpy2)

        end

        GainMatrix3 = nothing
        GainMatrix4 = nothing
        LossMatrix1 = nothing

        GC.gc()

    end

    if name1 != name2 && name3 == name4

        if typeof(Mode) == Axi

            GainMatrix_to_M_Bin_Axi!(M_Bin,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc])
            LossMatrix_to_M_Bin_Axi!(M_Bin,LossMatrix1,offset[name1_loc],offset[name2_loc])
            LossMatrix_to_M_Bin_Axi!(M_Bin,LossMatrix2,offset[name2_loc],offset[name1_loc])

        elseif typeof(Mode) == Iso

            Grids = PhaseSpace.Grids
            dpy1 = Grids.dpy[name1_loc]
            dpy2 = Grids.dpy[name2_loc]
            dpy3 = Grids.dpy[name3_loc]

            GainMatrix_to_M_Bin_Iso!(M_Bin,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],dpy3,dpy1,dpy2)
            LossMatrix_to_M_Bin_Iso!(M_Bin,LossMatrix1,offset[name1_loc],offset[name2_loc],dpy1,dpy2)
            LossMatrix_to_M_Bin_Iso!(M_Bin,LossMatrix2,offset[name2_loc],offset[name1_loc],dpy2,dpy1)

        end

        GainMatrix3 = nothing
        LossMatrix1 = nothing
        LossMatrix2 = nothing

        GC.gc()

    end

    if name1 != name2 && name3 != name4

        #=if typeof(Mode) == Axi

            GainMatrix_to_M_Bin_Axi!(M_Bin,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc])
            GainMatrix_to_M_Bin_Axi!(M_Bin,GainMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc])
            LossMatrix_to_M_Bin_Axi!(M_Bin,LossMatrix1,offset[name1_loc],offset[name2_loc])
            LossMatrix_to_M_Bin_Axi!(M_Bin,LossMatrix2,offset[name2_loc],offset[name1_loc])

        elseif typeof(Mode) == Iso

            Grids = PhaseSpace.Grids
            dpy1 = Grids.dpy[name1_loc]
            dpy2 = Grids.dpy[name2_loc]
            dpy3 = Grids.dpy[name3_loc]
            dpy4 = Grids.dpy[name4_loc]

            GainMatrix_to_M_Iso!(M_Bin,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],dpy3,dpy1,dpy2)
            GainMatrix_to_M_Bin!(M_Bin,GainMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc],dpy4,dpy1,dpy2)
            LossMatrix_to_M_Bin!(M_Bin,LossMatrix1,offset[name1_loc],offset[name2_loc],dpy1,dpy2)
            LossMatrix_to_M_Bin!(M_Bin,LossMatrix2,offset[name2_loc],offset[name1_loc],dpy2,dpy1)

        end=#

        GainMatrix_to_M_Bin!(M_Bin,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc])
        GainMatrix_to_M_Bin!(M_Bin,GainMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc])
        LossMatrix_to_M_Bin!(M_Bin,LossMatrix1,offset[name1_loc],offset[name2_loc])
        LossMatrix_to_M_Bin!(M_Bin,LossMatrix2,offset[name2_loc],offset[name1_loc])

        GainMatrix3 = nothing
        GainMatrix4 = nothing
        LossMatrix1 = nothing
        LossMatrix2 = nothing

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

function Fill_M_Emi!(M_Emi::AbstractArray{Float32,2},interaction::EmiStruct,PhaseSpace::PhaseSpaceStruct;GainMatrix2=nothing,GainMatrix3=nothing,LossMatrix1=nothing)

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

        #GainMatrix_to_M_Emi_Axi!(M_Emi,GainMatrix2,offset[name2_loc],offset[name1_loc])
        GainMatrix_to_M_Emi_Axi!(M_Emi,GainMatrix3,offset[name3_loc],offset[name1_loc])
        #LossMatrix_to_M_Emi_Axi!(M_Emi,LossMatrix1,offset[name1_loc])

    elseif typeof(mode) == Iso

        Grids = PhaseSpace.Grids
        dpy1 = Grids.dpy_list[name1_loc]
        dpy2 = Grids.dpy_list[name2_loc]
        dpy3 = Grids.dpy_list[name3_loc]
        
        #GainMatrix_to_M_Emi_Iso!(M_Emi,GainMatrix2,offset[name2_loc],offset[name1_loc],dpy2,dpy1)
        GainMatrix_to_M_Emi_Iso!(M_Emi,GainMatrix3,offset[name3_loc],offset[name1_loc],dpy3,dpy1)
        #LossMatrix_to_M_Emi_Iso!(M_Emi,LossMatrix1,offset[name1_loc],dpy1)

    end
    #GainMatrix2 = nothing
    GainMatrix3 = nothing
    #LossMatrix1 = nothing

    GC.gc()

end

function GainMatrix_to_M_Bin!(M_Bin::AbstractArray{Float32,2},GainMatrix::Array{Float64,9},offset3::Int64,offset1::Int64,offset2::Int64)

    px3_num = size(GainMatrix,1)-1 # ignore overflow bin
    py3_num = size(GainMatrix,2)
    pz3_num = size(GainMatrix,3)
    px1_num = size(GainMatrix,4)  
    py1_num = size(GainMatrix,5)
    pz1_num = size(GainMatrix,6)
    px2_num = size(GainMatrix,7)
    py2_num = size(GainMatrix,8)
    pz2_num = size(GainMatrix,9)

    N = size(M_Bin,2)
    #println("N = $N")
    #println("$offset2")

    for pz2 in 1:pz2_num, py2 in 1:py2_num, px2 in 1:px2_num, pz1 in 1:pz1_num, py1 in 1:py1_num, px1 in 1:px1_num, pz3 in 1:pz3_num, py3 in 1:py3_num, px3 in 1:px3_num

        a = (pz3-1)*px3_num*py3_num+(py3-1)*px3_num+px3+offset3
        b = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
        c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

        # M_Bin terms allocated symmetrically
        M_Bin[(b-1)*N+(a-1)+1,c] += GainMatrix[px3,py3,pz3,px1,py1,pz1,px2,py2,pz2]/2
        M_Bin[(c-1)*N+(a-1)+1,b] += GainMatrix[px3,py3,pz3,px1,py1,pz1,px2,py2,pz2]/2

    end

end

function LossMatrix_to_M_Bin!(M_Bin::Array{Float32},LossMatrix::Array{Float64,6},offset1::Int64,offset2::Int64)

    px1_num = size(LossMatrix,1)  
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)
    px2_num = size(LossMatrix,4)
    py2_num = size(LossMatrix,5)
    pz2_num = size(LossMatrix,6)

    N = size(M_Bin,2)

    for pz2 in 1:pz2_num, py2 in 1:py2_num, px2 in 1:px2_num, pz1 in 1:pz1_num, py1 in 1:py1_num, px1 in 1:px1_num

        a = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
        b = a
        c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

        # M_Bin terms allocated symmetrically
        M_Bin[(b-1)*N+(a-1)+1,c] -= LossMatrix[px1,py1,pz1,px2,py2,pz2]/2
        M_Bin[(c-1)*N+(a-1)+1,b] -= LossMatrix[px1,py1,pz1,px2,py2,pz2]/2

    end

end



function GainMatrix_to_M_Emi_Axi!(M_Emi::Array{Float32},GainMatrix::Array{Float64,4},offset2::Int64,offset1::Int64)

    px2_num = size(GainMatrix,1)#-1 # ignore overflow bin
    py2_num = size(GainMatrix,2)
    px1_num = size(GainMatrix,3)  
    py1_num = size(GainMatrix,4)

    for l in axes(GainMatrix,4), k in axes(GainMatrix,3), j in axes(GainMatrix,2), i in 1:px2_num

        a = (j-1)*px2_num+i+offset2
        b = (l-1)*px1_num+k+offset1

        M_Emi[a,b] += GainMatrix[i,j,k,l]

    end

end

function GainMatrix_to_M_Emi_Iso!(M_Emi::Array{Float32},GainMatrix::Array{Float64,4},offset2::Int64,offset1::Int64,dpy2,dpy1)

    px2_num = size(GainMatrix,1)#-1 # ignore overflow bin
    py2_num = size(GainMatrix,2)
    px1_num = size(GainMatrix,3)  
    py1_num = size(GainMatrix,4)

    for l in axes(GainMatrix,4), k in axes(GainMatrix,3), j in axes(GainMatrix,2), i in 1:px2_num

        val = 0.0 # py (pz) averaged GainMatrix term
        for p in axes(GainMatrix,4), q in axes(GainMatrix,2)
            val += GainMatrix[i,q,k,p] * dpy2[q] * dpy1[p] # check order
        end
        val /= sum(dpy2)*sum(dpy1)

        a = (j-1)*px2_num+i+offset2
        b = (l-1)*px1_num+k+offset1

        M_Emi[a,b] += val

    end

end

function LossMatrix_to_M_Emi_Axi!(M_Emi::Array{Float32},LossMatrix::Array{Float64,2},offset1::Int64)

    px1_num = size(LossMatrix,1)
    py1_num = size(LossMatrix,2)

    for j in axes(LossMatrix,2), i in axes(LossMatrix,1)

        a = (j-1)*px1_num+i+offset1
        b = a

        M_Emi[a,b] -= LossMatrix[i,j]
        
    end

end

function LossMatrix_to_M_Emi_Iso!(M_Emi::Array{Float32},LossMatrix::Array{Float64,2},offset1::Int64)

    px1_num = size(LossMatrix,1)
    py1_num = size(LossMatrix,2)

    for j in axes(LossMatrix,2), i in axes(LossMatrix,1)

        val = 0.0 # py (pz) averaged LossMatrix term
        for p in axes(LossMatrix,2)
            val += LossMatrix[i,p] * dpy1[p] # check order
        end
        val /= sum(dpy1)

        a = (j-1)*px1_num+i+offset1
        b = a

        M_Emi[a,b] -= val
        
    end

end

##  TO BE REMOVED LATER

function GainMatrix_to_M_Bin_Axi!(M_Bin::AbstractArray{Float32,2},GainMatrix::Array{Float64,6},offset3::Int64,offset1::Int64,offset2::Int64)

    px3_num = size(GainMatrix,1)-1 # ignore overflow bin
    py3_num = size(GainMatrix,2)
    px1_num = size(GainMatrix,3)  
    py1_num = size(GainMatrix,4)
    px2_num = size(GainMatrix,5)
    py2_num = size(GainMatrix,6)

    # sanity check 
    if px1_num != px2_num || py1_num != py2_num
        error("Error: GainMatrix dimensions not as expected")
    end

    N = size(M_Bin,2)

    for n in axes(GainMatrix,6), m in axes(GainMatrix,5), l in axes(GainMatrix,4), k in axes(GainMatrix,3), j in axes(GainMatrix,2), i in 1:px3_num

        a = (j-1)*px3_num+i+offset3
        b = (l-1)*px1_num+k+offset1
        c = (n-1)*px2_num+m+offset2
        
        # M_Bin terms allocated symmetrically
        M_Bin[(b-1)*N+(a-1)+1,c] += GainMatrix[i,j,k,l,m,n]/2
        M_Bin[(c-1)*N+(a-1)+1,b] += GainMatrix[i,j,k,l,m,n]/2 

    end

end

function GainMatrix_to_M_Bin_Iso!(M_Bin::AbstractArray{Float32,2},GainMatrix::Array{Float64,6},offset3::Int64,offset1::Int64,offset2::Int64,dpy3::Vector{Float64},dpy1::Vector{Float64},dpy2::Vector{Float64})

    px3_num = size(GainMatrix,1)-1 # ignore overflow bin
    py3_num = size(GainMatrix,2)
    px1_num = size(GainMatrix,3)  
    py1_num = size(GainMatrix,4)
    px2_num = size(GainMatrix,5)
    py2_num = size(GainMatrix,6)

    pxr1_grid

    # sanity check 
    if px1_num != px2_num || py1_num != py2_num
        error("Error: GainMatrix dimensions not as expected")
    end

    N = size(M_Bin,2)

    for n in axes(GainMatrix,6), m in axes(GainMatrix,5), l in axes(GainMatrix,4), k in axes(GainMatrix,3), j in axes(GainMatrix,2), i in 1:px3_num

        val = 0.0 # py (pz) averaged GainMatrix term 
        for p in axes(GainMatrix,6), q in axes(GainMatrix,4), r in axes(GainMatrix,2)
            val += GainMatrix[i,r,k,q,m,p] * dpy3[r] * dpy1[q] * dpy2[p] # check order
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

function LossMatrix_to_M_Bin_Axi!(M_Bin::Array{Float32},LossMatrix::Array{Float64,4},offset1::Int64,offset2::Int64)

    p1_num = size(LossMatrix,1) # ignore overflow bin
    u1_num = size(LossMatrix,2)
    p2_num = size(LossMatrix,3)  
    u2_num = size(LossMatrix,4)

    # sanity check 
    if p1_num != p2_num || u1_num != u2_num
        error("Error: GainMatrix dimensions not as expected")
    end

    N = size(M_Bin,2)

    for l in axes(LossMatrix,4), k in axes(LossMatrix,3), j in axes(LossMatrix,2), i in axes(LossMatrix,1)

        a = (j-1)*p1_num+i+offset1
        b = a
        c = (l-1)*p2_num+k+offset2
    
        # M_Bin terms allocated symmetrically
        M_Bin[(b-1)*N+(a-1)+1,c] -= LossMatrix[i,j,k,l]/2
        M_Bin[(c-1)*N+(a-1)+1,b] -= LossMatrix[i,j,k,l]/2

    end

end

function LossMatrix_to_M_Bin_Iso!(M_Bin::Array{Float32},LossMatrix::Array{Float64,4},offset1::Int64,offset2::Int64,dpy1,dpy2)

    p1_num = size(LossMatrix,1) # ignore overflow bin
    u1_num = size(LossMatrix,2)
    p2_num = size(LossMatrix,3)  
    u2_num = size(LossMatrix,4)

    # sanity check 
    if p1_num != p2_num || u1_num != u2_num
        error("Error: GainMatrix dimensions not as expected")
    end

    N = size(M_Bin,2)

    for l in axes(LossMatrix,4), k in axes(LossMatrix,3), j in axes(LossMatrix,2), i in axes(LossMatrix,1)

        val = 0.0 # py (pz) averaged LossMatrix term 
        for p in axes(LossMatrix,4), q in axes(LossMatrix,2)
            val += LossMatrix[i,q,k,p] * dpy1[q] * dpy2[p] # check order
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