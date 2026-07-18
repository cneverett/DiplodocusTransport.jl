mutable struct PadeExpWorkspaceStruct{T}
    dtA::Matrix{T}
    A2::Matrix{T}
    A4::Matrix{T}
    A6::Matrix{T}
    A8::Matrix{T}
    Umat::Matrix{T}
    Vmat::Matrix{T}
    lhs::Matrix{T}
    rhs::Matrix{T}
    tmp1::Matrix{T}
    tmp2::Matrix{T}
    Id::Matrix{T}
end

function PadeExpWorkspaceStruct{T}(A::AbstractMatrix{T}) where {T}
    m = size(A,1)

    mat(r, c) = begin
        X = similar(A, T, r, c)
        fill!(X, zero(T))
        X
    end

    dtA = mat(m, m)
    A2 = mat(m, m)
    A4 = mat(m, m)
    A6 = mat(m, m)
    A8 = mat(m, m)
    Umat = mat(m, m)
    Vmat = mat(m, m)
    lhs = mat(m, m)
    rhs = mat(m, m)
    tmp1 = mat(m, m)
    tmp2 = mat(m, m)
    Id = mat(m, m)
    #Id .= I

    return PadeExpWorkspaceStruct{T}(dtA,A2,A4,A6,A8,Umat,Vmat,lhs,rhs,tmp1,tmp2,Id)
end


# -----------------------------------------------------------------------------
# In-place Padé scaling-and-squaring exponential on the active k x k block only.
# This keeps the original full-matrix exp(dt*H) path, not the e1-action variant.
# -----------------------------------------------------------------------------
function expm_pade!(dt::T, F::M, A::M,ws::PadeExpWorkspaceStruct{T}) where {T,M<:AbstractMatrix{T}}
    
    #nrm = _matrix_one_norm!(ws, A)

    dtA = ws.dtA 
    @. dtA = dt * A

    nrm = norm(dtA,1)
    #println("nrm = $nrm, 1-norm: $(norm(dtA,1))")

    if T === Float32
        theta3 = 4.258730016922831f-1
        theta5 = 1.880152677804762f0
        theta7 = 3.925724783138660f0
        theta9 = 6.42276f0 # more digits?
        if nrm <= theta3
            return _pade3!(F, dtA, ws)
        elseif nrm <= theta5
            return _pade5!(F, dtA, ws)
        else
            s = max(0, ceil(Int, log2(nrm / theta7)))
            #println("nrm = $nrm, theta7 = $theta7, s = $s")
            if isodd(s) # force s to be even to allow mul step below to always give F without need for mem copies
                s += 1
            end
            @. dtA /= T(2)^s
            #println("nrmH: $(norm(ws.Hscaled)), maxH: $(maximum(ws.Hscaled)), minH: $(minimum(ws.Hscaled))")
            _pade7!(F, dtA, ws)
            T1 = ws.tmp1
            for sv in 1:s
                if isodd(sv)
                    mul!(T1, F, F)
                else
                    mul!(F, T1, T1)
                end
                #println("nrmF: $(norm(F)), maxF: $(maximum(F)), minF: $(minimum(F))")
                #mul!(T1, F, F)
                #copyto!(F,T1)
            end
            return F
        end
    else
        theta3  = 1.495585217958292e-2
        theta5  = 2.539398330063230e-1
        theta7  = 9.504178996162932e-1
        theta9  = 2.097847961257068e0
        theta13 = 5.371920351148152e0
        if nrm <= theta3
            return _pade3!(F, dtA, ws)
        elseif nrm <= theta5
            return _pade5!(F, dtA, ws)
        elseif nrm <= theta7
            return _pade7!(F, dtA, ws)
        elseif nrm <= theta9
            return _pade9!(F, dtA, ws)
        else
            s = max(0, ceil(Int, log2(nrm / theta13)))
            println("nrm = $nrm, theta13 = $theta13, s = $s")
            @. dtA /= T(2)^s
            _pade13!(F, dtA, ws)
            for _ in 1:s
                mul!(ws.tmp1, F, F)
                @. F = ws.tmp1
            end
            return F
        end
    end
end

function _finish_pade!(F::M, ws::PadeExpWorkspaceStruct{T}) where {T,M<:AbstractMatrix{T}}
    @. ws.lhs = ws.Vmat - ws.Umat
    @. ws.rhs = T(2) * ws.Umat
    L = ws.lhs
    R = ws.rhs
    #ldiv!(R, lu!(L),R)
    luL = lu(L)
    ldiv!(R,luL,R)
    F .= R + I
    return F
end

function _pade3!(F::M, A::M, ws::PadeExpWorkspaceStruct{T}) where {T,M<:AbstractMatrix{T}}
    A2 = ws.A2
    U = ws.Umat
    V = ws.Vmat
    mul!(A2, A, A)
    U .= T(1) .* A2 + T(60) * I
    mul!(ws.tmp1, A, U)
    @. U = ws.tmp1
    V .= T(12) .* A2 + T(120) * I
    return _finish_pade!(F, ws)
end

function _pade5!(F::M, A::M, ws::PadeExpWorkspaceStruct{T}) where {T,M<:AbstractMatrix{T}}
    A2 = ws.A2
    A4 = ws.A4
    U = ws.Umat
    V = ws.Vmat
    mul!(A2, A, A)
    mul!(A4, A2, A2)
    U .= T(1) .* A4 + T(420) * A2 + T(15120) * I
    mul!(ws.tmp1, A, U)
    @. U = ws.tmp1
    V .= T(30) .* A4 + T(3360) * A2 + T(30240) * I
    return _finish_pade!(F, ws)
end

function _pade7!(F::M, A::M, ws::PadeExpWorkspaceStruct{T}) where {T,M<:AbstractMatrix{T}}
    A2 = ws.A2
    A4 = ws.A4
    A6 = ws.A6
    U = ws.Umat 
    V = ws.Vmat
    T1 = ws.tmp1
    #I1 = ws.Id

    mul!(A2, A, A)
    mul!(A4, A2, A2)
    mul!(A6, A2, A4)
    U .= T(1) .* A6 + T(1512) * A4 + T(277200) * A2 + T(8648640) * I
    mul!(T1, A, U)
    U .= T1
    V .= T(56) .* A6 + T(25200) * A4 + T(1995840) * A2 + T(17297280) * I

    return _finish_pade!(F, ws)
end

function _pade9!(F::M, A::M, ws::PadeExpWorkspaceStruct{T}) where {T,M<:AbstractMatrix{T}}
    A2 = ws.A2
    A4 = ws.A4
    A6 = ws.A6
    A8 = ws.A8
    U = ws.Umat
    V = ws.Vmat
    mul!(A2, A, A)
    mul!(A4, A2, A2)
    mul!(A6, A2, A4)
    mul!(A8, A4, A4)
    U .= T(1) .* A8 + T(3960) * A6 + T(2162160) * A4 + T(302702400) * A2 + T(8821612800) * I
    mul!(ws.tmp1, A, U)
    U .= ws.tmp1
    V .= T(90) .* A8 + T(110880) * A6 + T(30270240) * A4 + T(2075673600) * A2 + T(17643225600) * I
    return _finish_pade!(F, ws)
end

function _pade13!(F::M, A::M, ws::PadeExpWorkspaceStruct{T}) where {T,M<:AbstractMatrix{T}}

        A2 = ws.A2
        A4 = ws.A4
        A6 = ws.A6
        U = ws.Umat
        V = ws.Vmat
        tmp1 = ws.tmp1
        tmp2 = ws.tmp2
        Ik = ws.Id

        mul!(A2, A, A)
        mul!(A4, A2, A2)
        mul!(A6, A2, A4)

        @. tmp1 = T(1) * A6 + T(16380) * A4 + T(40840800) * A2
        mul!(tmp2, A6, tmp1)
        tmp1 .= tmp2 .+ T(33522128640) .* A6 .+ T(10559470521600) .* A4 + T(1187353796428800) .* A2 + T(32382376266240000) * I
        mul!(U, A, tmp1)

        @. tmp1 = T(182) * A6 + T(960960) * A4 + T(1323241920) * A2
        mul!(tmp2, A6, tmp1)
        @. V = tmp2 .+ T(670442572800) .* A6 .+ T(129060195264000) .* A4 .+ T(7771770303897600) .* A2 .+ T(64764752532480000) * I

    return _finish_pade!(F, ws)
end