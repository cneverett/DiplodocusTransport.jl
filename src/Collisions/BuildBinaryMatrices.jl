"""
    BuildBigMatrices(PhaseSpace,DataDirectory;loading_check)

Function that builds the big matrices associated with binary and emissive interactions. If there are such interactions, first space is allocated for the arrays, then data is loaded into these arrays from the desired `DataDirectory` location and finally the big matrices are returned as an immutable `BinaryMatricesStruct`.
"""
function BuildBinaryMatrices(PhaseSpace::PhaseSpaceStruct,Binary_list::Vector{BinaryInteraction},Domain::Union{Vector{Int64},Nothing},DataDirectory::String;loading_check::Bool=false,Bin_Mode::AbstractMode=Ani(),Bin_corrected::Bool=true,Bin_sparse::Bool=false)

    Precision::DataType = getfield(Main,Symbol("Precision"))

    @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

    Momentum = PhaseSpace.Momentum

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n = sum(px_num_list.*py_num_list.*pz_num_list)
    m = n*n

    size = m*n*sizeof(Precision)

    if size > 1e9
        println("M_Bin will be approx. $(size/1e9) GB in memory if dense")
    elseif size > 1e6
        println("M_Bin will be approx. $(size/1e6) MB in memory if dense")
    elseif size > 1e3
        println("M_Bin will be approx. $(size/1e3) KB in memory if dense")
    else
        println("M_Bin will be approx. $size bytes in memory if dense")
    end

    if isempty(Binary_list)
        if Bin_sparse
            M_Bin = spzeros(Precision,0,0)
        else
            M_Bin = zeros(Precision,0,0)
        end 
    else
        if Bin_sparse
            M_Bin_I::Vector{Int64} = Int64[]
            M_Bin_J::Vector{Int64} = Int64[]
            M_Bin_V::Vector{Precision} = Precision[]
            LoadMatrices_Binary(Binary_list,DataDirectory,PhaseSpace,Bin_Mode,Bin_corrected;M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)
            println("Building sparse M_Bin")
            M_Bin = sparse(M_Bin_I,M_Bin_J,M_Bin_V,m,n)::SparseMatrixCSC{Precision,Int64}
   
            GC.gc()
        else
            M_Bin = zeros(Precision,m,n)::Matrix{Precision}
            LoadMatrices_Binary(Binary_list,DataDirectory,PhaseSpace,Bin_Mode,Bin_corrected;M_Bin=M_Bin)
        end
    end

    size = Base.summarysize(M_Bin)

    if size > 1e9
        println("M_Bin is approx. $(size/1e9) GB in memory")
    elseif size > 1e6
        println("M_Bin is approx. $(size/1e6) MB in memory")
    elseif size > 1e3
        println("M_Bin is approx. $(size/1e3) KB in memory")
    else
        println("M_Bin is approx. $size bytes in memory")
    end

    return BinaryMatricesStruct{Precision}(M_Bin,Binary_list,Domain)

end

function BuildBinaryMatricesPatankar(PhaseSpace::PhaseSpaceStruct,Binary_list::Vector{BinaryInteraction},Domain::Union{Vector{Int64},Nothing},DataDirectory::String;loading_check::Bool=false,Bin_Mode::AbstractMode=Ani(),Bin_corrected::Bool=true,Bin_sparse::Bool=false)

    Precision::DataType = getfield(Main,Symbol("Precision"))

    @assert Precision == Float32 || Precision == Float64 "Precision must be either Float32 or Float64"

    Momentum = PhaseSpace.Momentum

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n = sum(px_num_list.*py_num_list.*pz_num_list)
    m = n*n

    size = m*n*sizeof(Precision)

    if size > 1e9
        println("M_Bin will be approx. $(size/1e9) GB in memory if dense")
    elseif size > 1e6
        println("M_Bin will be approx. $(size/1e6) MB in memory if dense")
    elseif size > 1e3
        println("M_Bin will be approx. $(size/1e3) KB in memory if dense")
    else
        println("M_Bin will be approx. $size bytes in memory if dense")
    end

    if isempty(Binary_list)
        if Bin_sparse
            Gijk = spzeros(Precision,0,0)
            Liij = spzeros(Precision,0,0)
        else
            Gijk = zeros(Precision,0,0)
            Liij = zeros(Precision,0,0)
        end 
    else
        if Bin_sparse
            Gijk_I::Vector{Int64} = Int64[]
            Gijk_J::Vector{Int64} = Int64[]
            Gijk_V::Vector{Precision} = Precision[]
            
            Liij_I::Vector{Int64} = Int64[]
            Liij_J::Vector{Int64} = Int64[]
            Liij_V::Vector{Precision} = Precision[]
            LoadMatrices_BinaryPatankar(Binary_list,DataDirectory,PhaseSpace,Bin_Mode,Bin_corrected;Gijk_I=Gijk_I,Gijk_J=Gijk_J,Gijk_V=Gijk_V,Liij_I=Liij_I,Liij_J=Liij_J,Liij_V=Liij_V)
            println("Building sparse Gijk and Liij")
            Gijk = sparse(Gijk_I,Gijk_J,Gijk_V,m,n)::SparseMatrixCSC{Precision,Int64}
            Liij = sparse(Liij_I,Liij_J,Liij_V,m,n)::SparseMatrixCSC{Precision,Int64}
   
            GC.gc()
        else
            Gijk = zeros(Precision,m,n)::Matrix{Precision}
            Liij = zeros(Precision,m,n)::Matrix{Precision}
            LoadMatrices_BinaryPatankar(Binary_list,DataDirectory,PhaseSpace,Bin_Mode,Bin_corrected;Gijk=Gijk,Liij=Liij)
        end
    end

    size = Base.summarysize(Gijk)

    if size > 1e9
        println("Gijk is approx. $(size/1e9) GB in memory")
    elseif size > 1e6
        println("Gijk is approx. $(size/1e6) MB in memory")
    elseif size > 1e3
        println("Gijk is approx. $(size/1e3) KB in memory")
    else
        println("Gijk is approx. $size bytes in memory")
    end

    return BinaryMatricesStructPatankar{Precision}(Gijk,Liij,Binary_list,Domain)

end


"""
    Fill_M_Bin!(...,name_locs,PhaseSpace;GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2)

Fills the big matrix `M_Bin` directly if dense or the vectors of rows, columns and values `M_Bin_I`, `M_Bin_J``, `M_Bin_V` if sparse, with the interaction rates for a specific binary interactions given by `name_locs` and the collision arrays `GainMatrix3`, `GainMatrix4`, `LossMatrix1`, `LossMatrix2`.
"""
function Fill_M_Bin!(name_locs::Tuple{Int64,Int64,Int64,Int64},PhaseSpace::PhaseSpaceStruct,GainMatrix3::Array{Float64,9},GainMatrix4::Array{Float64,9},LossMatrix1::Array{Float64,6},LossMatrix2::Array{Float64,6},n_momentum::Int64;mode::AbstractMode=Ani(),M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    Grids = PhaseSpace.Grids
    offset = Grids.momentum_species_offset

    (name1_loc,name2_loc,name3_loc,name4_loc) = name_locs

    dpy1 = Grids.dpy_list[name1_loc]
    dpy2 = Grids.dpy_list[name2_loc]
    dpy3 = Grids.dpy_list[name3_loc]
    dpy4 = Grids.dpy_list[name4_loc]
    dpz1 = Grids.dpz_list[name1_loc]
    dpz2 = Grids.dpz_list[name2_loc]
    dpz3 = Grids.dpz_list[name3_loc]
    dpz4 = Grids.dpz_list[name4_loc]

    GainMatrix_to_M_Bin!(GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,dpy3,dpz3,n_momentum;M_Bin=M_Bin,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)
    GainMatrix_to_M_Bin!(GainMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,dpy4,dpz4,n_momentum;M_Bin=M_Bin,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)
    LossMatrix_to_M_Bin!(LossMatrix1,offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,n_momentum;M_Bin=M_Bin,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)
    LossMatrix_to_M_Bin!(LossMatrix2,offset[name2_loc],offset[name1_loc],mode,dpy2,dpz2,dpy1,dpz1,n_momentum;M_Bin=M_Bin,M_Bin_I=M_Bin_I,M_Bin_J=M_Bin_J,M_Bin_V=M_Bin_V)

    GainMatrix3 = nothing
    GainMatrix4 = nothing
    LossMatrix1 = nothing
    LossMatrix2 = nothing

    GC.gc()

    return nothing

end

function Fill_M_BinPatankar!(name_locs::Tuple{Int64,Int64,Int64,Int64},PhaseSpace::PhaseSpaceStruct,GainMatrix3::Array{Float64,9},GainMatrix4::Array{Float64,9},LossMatrix1::Array{Float64,6},LossMatrix2::Array{Float64,6},n_momentum::Int64;mode::AbstractMode=Ani(),Gijk::Union{Nothing,Matrix{F}}=nothing,Gijk_I::Union{Nothing,Vector{Int64}}=nothing,Gijk_J::Union{Nothing,Vector{Int64}}=nothing,Gijk_V::Union{Nothing,Vector{F}}=nothing,Liij::Union{Nothing,Matrix{F}}=nothing,Liij_I::Union{Nothing,Vector{Int64}}=nothing,Liij_J::Union{Nothing,Vector{Int64}}=nothing,Liij_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    Grids = PhaseSpace.Grids
    offset = Grids.momentum_species_offset

    (name1_loc,name2_loc,name3_loc,name4_loc) = name_locs

    dpy1 = Grids.dpy_list[name1_loc]
    dpy2 = Grids.dpy_list[name2_loc]
    dpy3 = Grids.dpy_list[name3_loc]
    dpy4 = Grids.dpy_list[name4_loc]
    dpz1 = Grids.dpz_list[name1_loc]
    dpz2 = Grids.dpz_list[name2_loc]
    dpz3 = Grids.dpz_list[name3_loc]
    dpz4 = Grids.dpz_list[name4_loc]

    GainMatrix_to_M_Bin!(GainMatrix3,offset[name3_loc],offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,dpy3,dpz3,n_momentum;M_Bin=Gijk,M_Bin_I=Gijk_I,M_Bin_J=Gijk_J,M_Bin_V=Gijk_V)
    GainMatrix_to_M_Bin!(GainMatrix4,offset[name4_loc],offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,dpy4,dpz4,n_momentum;M_Bin=Gijk,M_Bin_I=Gijk_I,M_Bin_J=Gijk_J,M_Bin_V=Gijk_V)
    LossMatrix_to_M_Bin!(LossMatrix1,offset[name1_loc],offset[name2_loc],mode,dpy1,dpz1,dpy2,dpz2,n_momentum;M_Bin=Liij,M_Bin_I=Liij_I,M_Bin_J=Liij_J,M_Bin_V=Liij_V)
    LossMatrix_to_M_Bin!(LossMatrix2,offset[name2_loc],offset[name1_loc],mode,dpy2,dpz2,dpy1,dpz1,n_momentum;M_Bin=Liij,M_Bin_I=Liij_I,M_Bin_J=Liij_J,M_Bin_V=Liij_V)

    GainMatrix3 = nothing
    GainMatrix4 = nothing
    LossMatrix1 = nothing
    LossMatrix2 = nothing

    GC.gc()

    return nothing

end



function GainMatrix_to_M_Bin!(GainMatrix::Array{Float64,9},offset3::Int64,offset1::Int64,offset2::Int64,mode::AbstractMode,dpy1::Vector{Float64},dpz1::Vector{Float64},dpy2::Vector{Float64},dpz2::Vector{Float64},dpy3::Vector{Float64},dpz3::Vector{Float64},n_momentum::Int64;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    px3_num = size(GainMatrix,1)-2 # ignore underflow and overflow bins
    py3_num = size(GainMatrix,2)
    pz3_num = size(GainMatrix,3)
    px1_num = size(GainMatrix,4)  
    py1_num = size(GainMatrix,5)
    pz1_num = size(GainMatrix,6)
    px2_num = size(GainMatrix,7)
    py2_num = size(GainMatrix,8)
    pz2_num = size(GainMatrix,9)

    N = n_momentum
    #println("N = $N")
    #println("$offset2")

    is_sparse = isnothing(M_Bin)

    for px2 in 1:px2_num, px1 in 1:px1_num, px3 in 1:px3_num

        #=if px1 == 1 || px2 == 1
            continue # skip first bin as the occupation of this bin can become very large causing time stepping issues.
        end=#

        if mode isa Iso

            val = 0.0 
            w = 1.0 / (sum(dpz1) * sum(dpz2) * sum(dpz3) * sum(dpy1) * sum(dpy2) * sum(dpy3))

            # average over incoming and outgoing u and phi angles (py,pz)
            for py1 in 1:py1_num, pz1 in 1:pz1_num, py2 in 1:py2_num, pz2 in 1:pz2_num, py3 in 1:py3_num, pz3 in 1:pz3_num
                val += GainMatrix[px3+1,py3,pz3,px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpz3[pz3] * dpy1[py1] * dpy2[py2] * dpy3[py3]
            end

        end

        for py1 in 1:py1_num, py2 in 1:py2_num, py3 in 1:py3_num

            if mode isa Axi

                val = 0.0 
                w = 1.0 / (sum(dpz1) * sum(dpz2) * sum(dpz3))

                # average over incoming and outgoing phi angles (pz)
                for pz1 in 1:pz1_num, pz2 in 1:pz2_num, pz3 in 1:pz3_num
                    val += GainMatrix[px3+1,py3,pz3,px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpz3[pz3]
                end

            end

            for pz1 in 1:pz1_num, pz2 in 1:pz2_num, pz3 in 1:pz3_num
                
                if mode isa Ani

                    if is_sparse
                        #GainMax = maximum(@view(GainMatrix[:,py3,pz3,px1,py1,pz1,px2,py2,pz2]))
                        GainMax = maximum(@view(GainMatrix[:,:,:,px1,py1,pz1,px2,py2,pz2]))
                    end

                    val = GainMatrix[px3+1,py3,pz3,px1,py1,pz1,px2,py2,pz2]
                    w = 1.0

                    #if is_sparse && val*w < eps(GainMax) # skips values smaller than this value to reduce memory usage
                    #    continue
                    #end

                end

                if val == 0.0
                    continue
                end

                a = (pz3-1)*px3_num*py3_num+(py3-1)*px3_num+px3+offset3
                b = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
                c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

                if is_sparse
                    # M_Bin terms allocated symmetrically
                    #push!(M_Bin_I,(b-1)*N+(a-1)+1)
                    #push!(M_Bin_J,c)
                    #push!(M_Bin_V,convert(F,val*w/2))
                    #push!(M_Bin_I,(c-1)*N+(a-1)+1)
                    #push!(M_Bin_J,b)
                    #push!(M_Bin_V,convert(F,val*w/2))
                    # M_Bin terms allocated asymmetrically (save memory and avoids implicit coupling)
                    push!(M_Bin_I,(b-1)*N+(a-1)+1)
                    push!(M_Bin_J,c)
                    push!(M_Bin_V,convert(F,val*w))
                else
                    # M_Bin terms allocated symmetrically
                    #M_Bin[(b-1)*N+(a-1)+1,c] += convert(F,val*w/2)
                    #M_Bin[(c-1)*N+(a-1)+1,b] += convert(F,val*w/2)
                    # M_Bin terms allocated asymmetrically
                    M_Bin[(b-1)*N+(a-1)+1,c] += convert(F,val*w)
                end

            end # pz loop

        end # py loop

    end # px loop

end

function LossMatrix_to_M_Bin!(LossMatrix::Array{Float64,6},offset1::Int64,offset2::Int64,mode::AbstractMode,dpy1::Vector{Float64},dpz1::Vector{Float64},dpy2::Vector{Float64},dpz2::Vector{Float64},n_momentum::Int64;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

    px1_num = size(LossMatrix,1)  
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)
    px2_num = size(LossMatrix,4)
    py2_num = size(LossMatrix,5)
    pz2_num = size(LossMatrix,6)

    N = n_momentum

    is_sparse = isnothing(M_Bin)

    #for pz2 in 1:pz2_num, py2 in 1:py2_num, px2 in 1:px2_num, pz1 in 1:pz1_num, py1 in 1:py1_num, px1 in 1:px1_num
    for px2 in 1:px2_num, px1 in 1:px1_num

        #=if px1 == 1 || px2 == 1
            continue # skip first bin as the occupation of this bin can become very large causing time stepping issues.
        end=#

        if mode isa Iso 

            val = 0.0 
            w = 1.0 / (sum(dpz1) * sum(dpz2) * sum(dpy1) * sum(dpy2))

            # average over incoming and outgoing u and phi angles (py,pz)
            for py2 in 1:py2_num, py1 in 1:py1_num, pz2 in 1:pz2_num, pz1 in 1:pz1_num
                val += LossMatrix[px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2] * dpy1[py1] * dpy2[py2]
            end

        end

        for py2 in 1:py2_num, py1 in 1:py1_num

            if mode isa Axi 

                val = 0.0 
                w = 1.0 / (sum(dpz1) * sum(dpz2))

                # average over incoming and outgoing phi angles (pz)
                for pz2 in 1:pz2_num, pz1 in 1:pz1_num
                    val += LossMatrix[px1,py1,pz1,px2,py2,pz2] * dpz1[pz1] * dpz2[pz2]
                end

            end

            for pz2 in 1:pz2_num, pz1 in 1:pz1_num

                if mode isa Ani

                    val = LossMatrix[px1,py1,pz1,px2,py2,pz2]
                    w = 1.0

                end

                if val == 0.0
                    continue
                end

                # Asymmetric: Labc = Laac δab  
                # Symmetric: (Laac δab + Laab δac) / 2
                # (Laac δab + Laab δac) fb fc / 2
                a = (pz1-1)*px1_num*py1_num+(py1-1)*px1_num+px1+offset1
                b = a
                c = (pz2-1)*px2_num*py2_num+(py2-1)*px2_num+px2+offset2

                # M_Bin terms allocated symmetrically
                if is_sparse
                    # symmetric in b and c
                    #push!(M_Bin_I,(b-1)*N+(a-1)+1)
                    #push!(M_Bin_J,c)
                    #push!(M_Bin_V,-convert(F,val*w/2))
                    #push!(M_Bin_I,(c-1)*N+(a-1)+1)
                    #push!(M_Bin_J,b)
                    #push!(M_Bin_V,-convert(F,val*w/2))
                    # symmetric in b and c
                    #push!(M_Bin_I,(b-1)*N+(a-1)+1)
                    #push!(M_Bin_J,c)
                    #push!(M_Bin_V,-convert(F,val*w/2 * (a==b)))
                    #push!(M_Bin_I,(c-1)*N+(a-1)+1)
                    #push!(M_Bin_J,b)
                    #push!(M_Bin_V,-convert(F,val*w/2 * (a==c)))
                    # asymmetric in b and c (save memory and works better with MPE solver)
                    push!(M_Bin_I,(b-1)*N+(a-1)+1)
                    push!(M_Bin_J,c)
                    push!(M_Bin_V,-convert(F,val*w))
                    #push!(M_Bin_I,(c-1)*N+(c-1)+1)
                    #push!(M_Bin_J,a)
                    #push!(M_Bin_V,-convert(F,val*w/2))
                else
                    # symmetric in b and c
                    #M_Bin[(b-1)*N+(a-1)+1,c] -= convert(F,val*w/2)
                    #M_Bin[(c-1)*N+(a-1)+1,b] -= convert(F,val*w/2)
                    # asymmetric in b and c
                    M_Bin[(b-1)*N+(a-1)+1,c] -= convert(F,val*w)
                end

            end # pz loop

        end # py loop

    end # px loop

end


#===== OLD CODE =======# 
#=
"""
    Allocate_M_Bin(PhaseSpace::PhaseSpaceStruct)

Allocates a big matrix `M_Bin` which stores the interaction rates for binary interactions between all particles in the simulation. If `n` is the size of the momentum domain, then `M_Bin` is an `n^2 x n` matrix. The size of `M_Bin` in memory is printed to the console upon allocation.
"""
function Allocate_M_Bin(PhaseSpace::PhaseSpaceStruct,loading_check::Bool,Precision::DataType,Bin_sparse::Bool)

    Momentum = PhaseSpace.Momentum

    px_num_list = Momentum.px_num_list
    py_num_list = Momentum.py_num_list
    pz_num_list = Momentum.pz_num_list

    n = sum(px_num_list.*py_num_list.*pz_num_list)
    m = n*n

    size = m*n*sizeof(Precision)

    if loading_check

        if size > 1e9
            println("M_Bin will be approx. $(size/1e9) GB in memory if dense")
        elseif size > 1e6
            println("M_Bin will be approx. $(size/1e6) MB in memory if dense")
        elseif size > 1e3
            println("M_Bin will be approx. $(size/1e3) KB in memory if dense")
        else
            println("M_Bin will be approx. $size bytes in memory if dense")
        end

        #println("Do you want to proceed? (y/n)")
        #proceed = readline()
        #if proceed != "y"
        #    error("Error: User aborted binary interaction matrix loading")
        #end

    else

        println("Building binary interaction matrix")
        if size > 1e9
            println("M_Bin will be approx. $(size/1e9) GB in memory if dense")
        elseif size > 1e6
            println("M_Bin will be approx. $(size/1e6) MB in memory if dense")
        elseif size > 1e3
            println("M_Bin will be approx. $(size/1e3) KB in memory if dense")
        else
            println("M_Bin will be approx. $size bytes in memory if dense")
        end

    end

    if Bin_sparse
        M_Bin = spzeros(Precision,m,n)
    else 
        M_Bin = zeros(Precision,m,n)
    end

    return M_Bin

end

function GainMatrix_to_M_Bin_Axi!(GainMatrix::Array{Float64,9},offset3::Int64,offset1::Int64,offset2::Int64,dpz1,dpz2,dpz3;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

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

    is_sparse = isnothing(M_Bin)

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
            if is_sparse
                push!(M_Bin_I,(b-1)*N+(a-1)+1)
                push!(M_Bin_J,c)
                push!(M_Bin_V,convert(F,val*w/2))
                push!(M_Bin_I,(c-1)*N+(a-1)+1)
                push!(M_Bin_J,b)
                push!(M_Bin_V,convert(F,val*w/2))
            else
                M_Bin[(b-1)*N+(a-1)+1,c] += convert(F,val*w/2)
                M_Bin[(c-1)*N+(a-1)+1,b] += convert(F,val*w/2) 
            end

        end

    end

end

function GainMatrix_to_M_Bin_Iso!(GainMatrix::Array{Float64,9},offset3::Int64,offset1::Int64,offset2::Int64,dpy1,dpz1,dpy2,dpz2,dpy3,dpz3;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F<:Union{Float32,Float64}

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
            if is_sparse
                push!(M_Bin_I,(b-1)*N+(a-1)+1)
                push!(M_Bin_J,c)
                push!(M_Bin_V,convert(F,val*w/2))
                push!(M_Bin_I,(c-1)*N+(a-1)+1)
                push!(M_Bin_J,b)
                push!(M_Bin_V,convert(F,val*w/2))
            else
                M_Bin[(b-1)*N+(a-1)+1,c] += convert(F,val*w/2)
                M_Bin[(c-1)*N+(a-1)+1,b] += convert(F,val*w/2) 
            end

        end

    end

end

function LossMatrix_to_M_Bin_Axi!(LossMatrix::Array{Float64,6},offset1::Int64,offset2::Int64,dpz1,dpz2;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F <:Union{Float32,Float64}

    px1_num = size(LossMatrix,1)  
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)
    px2_num = size(LossMatrix,4)
    py2_num = size(LossMatrix,5)
    pz2_num = size(LossMatrix,6)

    N = size(M_Bin,2)

    is_sparse = isnothing(M_Bin)

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
            if is_sparse
                push!(M_Bin_I,(b-1)*N+(a-1)+1)
                push!(M_Bin_J,c)
                push!(M_Bin_V,-convert(F,val*w/2))
                push!(M_Bin_I,(c-1)*N+(a-1)+1)
                push!(M_Bin_J,b)
                push!(M_Bin_V,-convert(F,val*w/2))
            else
                M_Bin[(b-1)*N+(a-1)+1,c] -= convert(F,val*w/2)
                M_Bin[(c-1)*N+(a-1)+1,b] -= convert(F,val*w/2)
            end

        end

    end

end

function LossMatrix_to_M_Bin_Iso!(LossMatrix::Array{Float64,6},offset1::Int64,offset2::Int64,dpy1,dpz1,dpy2,dpz2;M_Bin::Union{Nothing,Matrix{F}}=nothing,M_Bin_I::Union{Nothing,Vector{Int64}}=nothing,M_Bin_J::Union{Nothing,Vector{Int64}}=nothing,M_Bin_V::Union{Nothing,Vector{F}}=nothing) where F <:Union{Float32,Float64}

    px1_num = size(LossMatrix,1)  
    py1_num = size(LossMatrix,2)
    pz1_num = size(LossMatrix,3)
    px2_num = size(LossMatrix,4)
    py2_num = size(LossMatrix,5)
    pz2_num = size(LossMatrix,6)

    N = size(M_Bin,2)

    is_sparse = isnothing(M_Bin)

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
            if is_sparse
                push!(M_Bin_I,(b-1)*N+(a-1)+1)
                push!(M_Bin_J,c)
                push!(M_Bin_V,-convert(F,val*w/2))
                push!(M_Bin_I,(c-1)*N+(a-1)+1)
                push!(M_Bin_J,b)
                push!(M_Bin_V,-convert(F,val*w/2))
            else
                M_Bin[(b-1)*N+(a-1)+1,c] -= convert(F,val*w/2)
                M_Bin[(c-1)*N+(a-1)+1,b] -= convert(F,val*w/2)
            end

        end

    end

end

=#
