"""
    BuildBigMatrices(PhaseSpace,DataDirectory;loading_check)

Function that builds the big matrices associated with binary and emissive interactions. If there are such interactions, first space is allocated for the arrays, then data is loaded into these arrays from the desired `DataDirectory` location and finally the big matrices are returned as an immutable `BigMatricesStruct`.
"""
function BuildBigMatrices(PhaseSpace::PhaseSpaceStruct,DataDirectory::String;loading_check::Bool=false,Precision::DataType=Float32,Bin_Mode::ModeType=Ani(),Bin_corrected::Bool=true,Emi_corrected::Bool=true,Bin_sparse::Bool=false)

    if Bin_sparse == true && Bin_corrected == false
        error("Bin_sparse=true requires Bin_corrected=true")
    end

    if isempty(PhaseSpace.Binary_list)
        if Bin_sparse
            M_Bin = spzeros(Float32,0,0)
        else
            M_Bin = zeros(Float32,0,0)
        end 
    else
        M_Bin = Allocate_M_Bin(PhaseSpace,loading_check,Precision,Bin_sparse)
        LoadMatrices_Binary(M_Bin,DataDirectory,PhaseSpace,Bin_Mode,Bin_corrected,Bin_sparse)
    end

    if isempty(PhaseSpace.Emi_list)
        M_Emi = Matrix{Precision}(undef,0,0)
    else
        M_Emi = Allocate_M_Emi(PhaseSpace,loading_check,Precision)
        LoadMatrices_Emi(M_Emi,DataDirectory,PhaseSpace,corrected=Emi_corrected)
    end

    return BigMatricesStruct{Precision}(M_Bin,M_Emi)

end


"""
    Allocate_M_Bin(PhaseSpace::PhaseSpaceStruct)

Allocates a big matrix `M_Bin` which stores the interaction rates for binary interactions between all particles in the simulation. If `n` is the size of the momentum domain, then `M_Bin` is an `n^2 x n` matrix. The size of `M_Bin` in memory is printed to the console upon allocation.
"""
function Allocate_M_Bin(PhaseSpace::PhaseSpaceStruct,loading_check::Bool,Precision::DataType,Bin_sparse::Bool=false)

    Momentum = PhaseSpace.Momentum

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n = sum(px_num_list.*py_num_list.*pz_num_list)
    m = n*n

    size = m*n*sizeof(Precision)

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

        #println("Do you want to proceed? (y/n)")
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

    if Bin_sparse
        M_Bin = spzeros(Precision,m,n)
    else 
        M_Bin = zeros(Precision,m,n)
    end

    return M_Bin

end


"""
    Fill_M_Bin!(M_Bin::Array{Float32,2},interaction::Vector{String},Lists::ListStruct;GainMatrix3=nothing,GainMatrix4=nothing,LossMatrix1=nothing,LossMatrix2=nothing)

Fills the big matrix `M_Bin` with the interaction rates for binary interactions between all particles in the simulation.
"""
function Fill_M_Bin!(M_Bin::AbstractMatrix{F},name_locs::Tuple{Int64,Int64,Int64,Int64},PhaseSpace::PhaseSpaceStruct,GainMatrix3::Array{Float64,9},GainMatrix4::Array{Float64,9},LossMatrix1::Array{Float64,6},LossMatrix2::Array{Float64,6};mode::ModeType=Ani(),Bin_sparse::Bool=false) where F<:Union{Float32,Float64}

    Grids = PhaseSpace.Grids
    offset = Grids.momentum_species_offset

    (name1_loc,name2_loc,name3_loc,name4_loc) = name_locs

    # first find offsets of the different species in the big matrix A

    if mode isa Ani

        GainMatrix_to_M_Bin!(M_Bin,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],Bin_sparse=Bin_sparse)
        GainMatrix_to_M_Bin!(M_Bin,GainMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc],Bin_sparse=Bin_sparse)
        LossMatrix_to_M_Bin!(M_Bin,LossMatrix1,offset[name1_loc],offset[name2_loc])
        LossMatrix_to_M_Bin!(M_Bin,LossMatrix2,offset[name2_loc],offset[name1_loc])

    elseif mode isa Axi

        dpz1 = Grids.dpz_list[name1_loc]
        dpz2 = Grids.dpz_list[name2_loc]
        dpz3 = Grids.dpz_list[name3_loc]
        dpz4 = Grids.dpz_list[name4_loc]

        GainMatrix_to_M_Bin_Axi!(M_Bin,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],dpz1,dpz2,dpz3)
        GainMatrix_to_M_Bin_Axi!(M_Bin,GainMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc],dpz1,dpz2,dpz4)
        LossMatrix_to_M_Bin_Axi!(M_Bin,LossMatrix1,offset[name1_loc],offset[name2_loc],dpz1,dpz2)
        LossMatrix_to_M_Bin_Axi!(M_Bin,LossMatrix2,offset[name2_loc],offset[name1_loc],dpz1,dpz2)

    elseif mode isa Iso

        dpy1 = Grids.dpy_list[name1_loc]
        dpy2 = Grids.dpy_list[name2_loc]
        dpy3 = Grids.dpy_list[name3_loc]
        dpy4 = Grids.dpy_list[name4_loc]
        dpz1 = Grids.dpz_list[name1_loc]
        dpz2 = Grids.dpz_list[name2_loc]
        dpz3 = Grids.dpz_list[name3_loc]
        dpz4 = Grids.dpz_list[name4_loc]

        GainMatrix_to_M_Bin_Iso!(M_Bin,GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],dpy1,dpz1,dpy2,dpz2,dpy3,dpz3)
        GainMatrix_to_M_Bin_Iso!(M_Bin,GainMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc],dpy1,dpz1,dpy2,dpz2,dpy4,dpz4)
        LossMatrix_to_M_Bin_Iso!(M_Bin,LossMatrix1,offset[name1_loc],offset[name2_loc],dpy1,dpz1,dpy2,dpz2)
        LossMatrix_to_M_Bin_Iso!(M_Bin,LossMatrix2,offset[name2_loc],offset[name1_loc],dpy1,dpz1,dpy2,dpz2)

    end

    GainMatrix3 = nothing
    GainMatrix4 = nothing
    LossMatrix1 = nothing
    LossMatrix2 = nothing

    GC.gc()

    #end

end

"""
    Allocate_M_Emi(Momentum::MomentumStruct)

Allocates a big matrix `M_Emi` which stores the interaction rates for emission interactions between all particles in the simulation.
"""
function Allocate_M_Emi(PhaseSpace::PhaseSpaceStruct,loading_check::Bool,Precision::DataType)

    Momentum = PhaseSpace.Momentum

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n = sum(px_num_list.*py_num_list.*pz_num_list)

    size = n*n*sizeof(Precision)

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

    M_Emi::Matrix{Precision} = Matrix{Precision}(undef,n,n)
    fill!(M_Emi,zero(Precision))

    return M_Emi

end

function Fill_M_Emi!(M_Emi::Matrix{F},name_locs::Tuple{Int64,Int64,Int64},PhaseSpace::PhaseSpaceStruct;GainMatrix2=nothing,GainMatrix3=nothing,LossMatrix1=nothing,mode::ModeType=Ani()) where F<:Union{Float32,Float64}

    Grids = PhaseSpace.Grids
    offset = Grids.momentum_species_offset

    (name1_loc,name2_loc,name3_loc) = name_locs

    # interaction is name1 -> name2 + name3
    # absorption of name1 and emission of name2 not implemented with name1 == name2
    if mode isa Ani

        #GainMatrix_to_M_Emi_Ani!(M_Emi,GainMatrix2,offset[name2_loc],offset[name1_loc])
        GainMatrix_to_M_Emi_Ani!(M_Emi,GainMatrix3,offset[name3_loc],offset[name1_loc])
        #LossMatrix_to_M_Emi_Ani!(M_Emi,LossMatrix1,offset[name1_loc])

    elseif mode isa Axi

        dpz1 = Grids.dpz_list[name1_loc]
        dpz2 = Grids.dpz_list[name2_loc]
        dpz3 = Grids.dpz_list[name3_loc]

        #GainMatrix_to_M_Emi_Axi!(M_Emi,GainMatrix2,offset[name2_loc],offset[name1_loc])
        GainMatrix_to_M_Emi_Axi!(M_Emi,GainMatrix3,offset[name3_loc],offset[name1_loc],dpz3,dpz1)
        #LossMatrix_to_M_Emi_Axi!(M_Emi,LossMatrix1,offset[name1_loc],dpz1)

    elseif mode isa Iso

        dpy1 = Grids.dpy_list[name1_loc]
        dpy2 = Grids.dpy_list[name2_loc]
        dpy3 = Grids.dpy_list[name3_loc]
        dpz1 = Grids.dpz_list[name1_loc]
        dpz2 = Grids.dpz_list[name2_loc]
        dpz3 = Grids.dpz_list[name3_loc]
        
        #GainMatrix_to_M_Emi_Iso!(M_Emi,GainMatrix2,offset[name2_loc],offset[name1_loc],dpy2,dpy1)
        GainMatrix_to_M_Emi_Iso!(M_Emi,GainMatrix3,offset[name3_loc],offset[name1_loc],dpy3,dpy1,dpz3,dpz1)
        #LossMatrix_to_M_Emi_Iso!(M_Emi,LossMatrix1,offset[name1_loc],dpy1,dpz1)

    end
    #GainMatrix2 = nothing
    GainMatrix3 = nothing
    #LossMatrix1 = nothing

    GC.gc()

end

function GainMatrix_to_M_Bin!(M_Bin::AbstractMatrix{F},GainMatrix::Array{Float64,9},offset3::Int64,offset1::Int64,offset2::Int64;Bin_sparse::Bool=false) where F<:Union{Float32,Float64}

    px3_num = size(GainMatrix,1)-2 # ignore underflow and overflow bins
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

    for pz2 in 1:pz2_num, py2 in 1:py2_num, px2 in 1:px2_num, pz1 in 1:pz1_num, py1 in 1:py1_num, px1 in 1:px1_num, pz3 in 1:pz3_num, py3 in 1:py3_num

        if Bin_sparse
            GainMax = maximum(@view(GainMatrix[:,py3,pz3,px1,py1,pz1,px2,py2,pz2]))
        end
        
        for px3 in 1:px3_num

            # First index of GainMatrix is [underflow,1,2,3,4,...,px_num,overflow], so need to offset by 1
            val = GainMatrix[px3+1,py3,pz3,px1,py1,pz1,px2,py2,pz2]

            if Bin_sparse && val < GainMax * 1e-6 # skips values smaller than this value to reduce memory usage
                continue
            end

            a = (pz3-1)*px3_num*py3_num+(py3-1)*px3_num+px3+offset3
            b = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
            c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

            # M_Bin terms allocated symmetrically
            M_Bin[(b-1)*N+(a-1)+1,c] += convert(F,val/2)
            M_Bin[(c-1)*N+(a-1)+1,b] += convert(F,val/2)

        end

    end

end

function LossMatrix_to_M_Bin!(M_Bin::Matrix{F},LossMatrix::Array{Float64,6},offset1::Int64,offset2::Int64) where F<:Union{Float32,Float64}

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
        M_Bin[(b-1)*N+(a-1)+1,c] -= convert(F,LossMatrix[px1,py1,pz1,px2,py2,pz2]/2)
        M_Bin[(c-1)*N+(a-1)+1,b] -= convert(F,LossMatrix[px1,py1,pz1,px2,py2,pz2]/2)

    end

end

function GainMatrix_to_M_Emi_Ani!(M_Emi::Matrix{F},GainMatrix::Array{Float64,6},offset2::Int64,offset1::Int64) where F<:Union{Float32,Float64}

    # 1 is incident particle, 2 is emitted particle

    px2_num = size(GainMatrix,1)
    py2_num = size(GainMatrix,2)
    pz2_num = size(GainMatrix,3)
    px1_num = size(GainMatrix,4)  
    py1_num = size(GainMatrix,5)
    pz1_num = size(GainMatrix,6)

    for pz1 in axes(GainMatrix,6), py1 in axes(GainMatrix,5), px1 in axes(GainMatrix,4), pz2 in axes(GainMatrix,3), py2 in axes(GainMatrix,2), px2 in axes(GainMatrix,1)

        a = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2
        b = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1

        M_Emi[a,b] += convert(F,GainMatrix[px2,py2,pz2,px1,py1,pz1])

    end

end

function GainMatrix_to_M_Emi_Axi!(M_Emi::Matrix{F},GainMatrix::Array{Float64,6},offset2::Int64,offset1::Int64,dpz2,dpz1) where F<:Union{Float32,Float64}

    # 1 is incident particle, 2 is emitted particle

    px2_num = size(GainMatrix,1)
    py2_num = size(GainMatrix,2)
    pz2_num = size(GainMatrix,3)
    px1_num = size(GainMatrix,4)  
    py1_num = size(GainMatrix,5)
    pz1_num = size(GainMatrix,6)

    for px1 in axes(GainMatrix,4), px2 in axes(GainMatrix,1), py1 in axes(GainMatrix,5), py2 in axes(GainMatrix,2)

        val = 0.0 # pz averaged GainMatrix term
        # average over incoming and outgoing phi angles
        for pz1 in axes(GainMatrix,6), pz2 in axes(GainMatrix,3) 
            val += GainMatrix[px2,py2,pz2,px1,py1,pz1] * dpz2[pz2] * dpz1[pz1]
        end

        for pz1 in axes(GainMatrix,6), pz2 in axes(GainMatrix,3) 

            a = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2
            b = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1

            # weight for the size of the incoming and outgoing phase space
            w = 1 / sum(dpz1) / sum(dpz2)

            M_Emi[a,b] += convert(F,val*w)

        end

    end

end

function GainMatrix_to_M_Emi_Iso!(M_Emi::Matrix{F},GainMatrix::Array{Float64,6},offset2::Int64,offset1::Int64,dpy2,dpy1,dpz2,dpz1) where F<:Union{Float32,Float64}

    px2_num = size(GainMatrix,1)
    py2_num = size(GainMatrix,2)
    pz2_num = size(GainMatrix,3)
    px1_num = size(GainMatrix,4)  
    py1_num = size(GainMatrix,5)
    pz1_num = size(GainMatrix,6)

    for px1 in axes(GainMatrix,4), px2 in axes(GainMatrix,1)

        val = 0.0 # py pz averaged GainMatrix term
        # average over incoming and outgoing u and phi angles
        for py1 in axes(GainMatrix,5), pz1 in axes(GainMatrix,6), py2 in axes(GainMatrix,2), pz2 in axes(GainMatrix,3)
            val += GainMatrix[px2,py2,pz2,px1,py1,pz1] * dpy2[py2] * dpz2[pz2] * dpy1[py1] * dpz1[pz1]
        end

        for py1 in axes(GainMatrix,5), pz1 in axes(GainMatrix,6), py2 in axes(GainMatrix,2), pz2 in axes(GainMatrix,3)

            a = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2
            b = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1

            # weight for the size of the incoming and outgoing phase space
            w = 1 / (sum(dpz1) * sum(dpy1)) / (sum(dpz2) * sum(dpy2))

            M_Emi[a,b] += convert(F,val*w) 

        end

    end

end

function LossMatrix_to_M_Emi_Axi!(M_Emi::Matrix{F},LossMatrix::Array{Float64,2},offset1::Int64,dpz1) where F<:Union{Float32,Float64}

    px1_num = size(LossMatrix,1)
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)

    for px1 in axes(LossMatrix,1), py1 in axes(LossMatrix,2)

        val = 0.0 # pz averaged LossMatrix term
        for pz1 in axes(LossMatrix,3)
            val += LossMatrix[px1,py1,pz1] * dpz1[pz1] # check order
        end

        for py1 in axes(LossMatrix,2), pz1 in axes(LossMatrix,3)

            a = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
            b = a

            w = 1 / sum(dpz1)

            M_Emi[a,b] -= convert(F,val*w)

        end
        
    end

end

function LossMatrix_to_M_Emi_Iso!(M_Emi::Matrix{F},LossMatrix::Array{Float64,2},offset1::Int64,dpy1,dpz1) where F<:Union{Float32,Float64}

    px1_num = size(LossMatrix,1)
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)

    for px1 in axes(LossMatrix,1)

        val = 0.0 # py (pz) averaged LossMatrix term
        for py1 in axes(LossMatrix,2), pz1 in axes(LossMatrix,3)
            val += LossMatrix[px1,py1,pz1] * dpy1[py1] * dpz1[pz1] # check order
        end

        for py1 in axes(LossMatrix,2), pz1 in axes(LossMatrix,3)

            a = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
            b = a

            w = 1 / (sum(dpz1) * sum(dpy1))

            M_Emi[a,b] -= convert(F,val*w)

        end
        
    end

end

function GainMatrix_to_M_Bin_Axi!(M_Bin::Matrix{F},GainMatrix::Array{Float64,9},offset3::Int64,offset1::Int64,offset2::Int64,dpz1,dpz2,dpz3) where F<:Union{Float32,Float64}

    px3_num = size(GainMatrix,1)-2 # ignore underflow and overflow bins
    py3_num = size(GainMatrix,2)
    pz3_num = size(GainMatrix,3)
    px1_num = size(GainMatrix,4)  
    py1_num = size(GainMatrix,5)
    pz1_num = size(GainMatrix,6)
    px2_num = size(GainMatrix,7)
    py2_num = size(GainMatrix,8)
    pz2_num = size(GainMatrix,9)

    N = size(M_Bin,2)

    for py2 in axes(GainMatrix,8), px2 in axes(GainMatrix,7), py1 in axes(GainMatrix,5), px1 in axes(GainMatrix,4), py3 in axes(GainMatrix,2), px3 in 1:px3_num

        val = 0.0 # pz averaged GainMatrix term
        for pz2 in axes(GainMatrix,9), pz1 in axes(GainMatrix,6), pz3 in axes(GainMatrix,3)
            val += GainMatrix[px3+1,py3,pz3,px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpz3[pz3]
        end

        for pz2 in axes(GainMatrix,9), pz1 in axes(GainMatrix,6), pz3 in axes(GainMatrix,3)

            a = (pz3-1)*px3_num*py3_num+(py3-1)*px3_num+px3+offset3
            b = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
            c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

            # weight for the size of the incoming and outgoing phase space
            w = 1 / sum(dpz1) / sum(dpz2) / sum(dpz3)
        
            # M_Bin terms allocated symmetrically
            M_Bin[(b-1)*N+(a-1)+1,c] += convert(F,val*w/2)
            M_Bin[(c-1)*N+(a-1)+1,b] += convert(F,val*w/2) 

        end

    end

end

function GainMatrix_to_M_Bin_Iso!(M_Bin::Matrix{F},GainMatrix::Array{Float64,9},offset3::Int64,offset1::Int64,offset2::Int64,dpy1,dpz1,dpy2,dpz2,dpy3,dpz3) where F<:Union{Float32,Float64}

    px3_num = size(GainMatrix,1)-2 # ignore underflow and overflow bins
    py3_num = size(GainMatrix,2)
    pz3_num = size(GainMatrix,3)
    px1_num = size(GainMatrix,4)  
    py1_num = size(GainMatrix,5)
    pz1_num = size(GainMatrix,6)
    px2_num = size(GainMatrix,7)
    py2_num = size(GainMatrix,8)
    pz2_num = size(GainMatrix,9)

    N = size(M_Bin,2)

    for px2 in axes(GainMatrix,7), px1 in axes(GainMatrix,4), px3 in 1:px3_num

        val = 0.0 # py pz averaged GainMatrix term
        for py2 in axes(GainMatrix,8), pz2 in axes(GainMatrix,9), pz1 in axes(GainMatrix,6),py1 in axes(GainMatrix,5), pz3 in axes(GainMatrix,3), py3 in axes(GainMatrix,2)
            val += GainMatrix[px3+1,py3,pz3,px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpy1[py1] * dpz2[pz2] * dpy2[py2] * dpz3[pz3] * dpy3[py3]
        end

        for py2 in axes(GainMatrix,8), pz2 in axes(GainMatrix,9), pz1 in axes(GainMatrix,6),py1 in axes(GainMatrix,5), pz3 in axes(GainMatrix,3), py3 in axes(GainMatrix,2)

            a = (pz3-1)*px3_num*py3_num+(py3-1)*px3_num+px3+offset3
            b = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
            c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

            # weight for the size of the incoming and outgoing phase space
            w = 1 / sum(dpz1) / sum(dpz2) / sum(dpz3) / sum(dpy1) / sum(dpy2) / sum(dpy3)
        
            # M_Bin terms allocated symmetrically
            M_Bin[(b-1)*N+(a-1)+1,c] += convert(F,val*w/2)
            M_Bin[(c-1)*N+(a-1)+1,b] += convert(F,val*w/2) 

        end

    end

end

function LossMatrix_to_M_Bin_Axi!(M_Bin::Matrix{F},LossMatrix::Array{Float64,6},offset1::Int64,offset2::Int64,dpz1,dpz2) where F <:Union{Float32,Float64}

    px1_num = size(LossMatrix,1)  
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)
    px2_num = size(LossMatrix,4)
    py2_num = size(LossMatrix,5)
    pz2_num = size(LossMatrix,6)

    N = size(M_Bin,2)

    for py2 in axes(LossMatrix,5), px2 in axes(LossMatrix,4), py1 in axes(LossMatrix,2), px1 in axes(GainMatrix,1)

        val = 0.0
        for pz2 in axes(GainMatrix,6), pz1 in axes(GainMatrix,3)
            val += LossMatrix[px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2]
        end

        for pz2 in axes(GainMatrix,6), pz1 in axes(GainMatrix,3)

            a = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
            b = a
            c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

            w = 1 / sum(dpz1) / sum(dpz2)
        
            # M_Bin terms allocated symmetrically
            M_Bin[(b-1)*N+(a-1)+1,c] -= convert(F,val*w/2)
            M_Bin[(c-1)*N+(a-1)+1,b] -= convert(F,val*w/2)

        end

    end

end

function LossMatrix_to_M_Bin_Iso!(M_Bin::Matrix{F},LossMatrix::Array{Float64,6},offset1::Int64,offset2::Int64,dpy1,dpz1,dpy2,dpz2) where F<:Union{Float32,Float64}

    px1_num = size(LossMatrix,1)  
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)
    px2_num = size(LossMatrix,4)
    py2_num = size(LossMatrix,5)
    pz2_num = size(LossMatrix,6)

    N = size(M_Bin,2)

    for px2 in axes(LossMatrix,4), px1 in axes(LossMatrix,1)

        val = 0.0
        for pz2 in axes(LossMatrix,6), py2 in axes(LossMatrix,5), pz1 in axes(LossMatrix,3),py1 in axes(LossMatrix,2)
            val += LossMatrix[px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpy1[py1] * dpy2[py2]
        end

        for pz2 in axes(LossMatrix,6), py2 in axes(LossMatrix,5), pz1 in axes(LossMatrix,3),py1 in axes(LossMatrix,2)

            a = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
            b = a
            c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

            w = 1 / sum(dpz1) / sum(dpz2) / sum(dpy1) / sum(dpy2)
        
            # M_Bin terms allocated symmetrically
            M_Bin[(b-1)*N+(a-1)+1,c] -= convert(F,val*w/2)
            M_Bin[(c-1)*N+(a-1)+1,b] -= convert(F,val*w/2)

        end

    end

end