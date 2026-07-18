using LinearAlgebra


# A is the augmented matrix A = [J u; 0 0] 
function arnoldi_scalar!(A, u, V, H, w; reorth::Bool = true)
    T = eltype(V)
    n = length(u)
    m = size(H, 2) - 1

    @assert size(A, 1) == n && size(A, 2) == n
    @assert size(V, 1) == n && size(V, 2) >= m + 1
    @assert size(H, 1) >= m + 1
    @assert length(w) == n

    β0 = norm(u)
    iszero(β0) && throw(ArgumentError("Initial vector u has zero norm"))

    copyto!(@view(V[:, 1]), u)
    rmul!(@view(V[:, 1]), inv(β0)) # TODO: not needed if u is [0,0,0,0 ... 1]

    hkp1k = zero(T)

    @inbounds for k in 1:m
        Vk  = @view V[:, 1:k]
        vk  = @view V[:, k]
        vkp = @view V[:, k+1]
        hk  = @view H[1:k, k]

        mul!(w, A, vk)
        mul!(hk, adjoint(Vk), w)
        mul!(w, Vk, hk, -one(T), one(T))

        if reorth
            mul!(hk, adjoint(Vk), w)
            mul!(w, Vk, hk, -one(T), one(T))
        end

        hkp1k = norm(w)

        if iszero(hkp1k)
            println("happy")
            fill!(vkp, zero(T))
            # fill last entry for error estimate 
            fill!(@view(H[1, m+1]), one(T))
            break
        end

        fill!(@view(H[k+1, k]), hkp1k)

        copyto!(vkp, w)
        rmul!(vkp, inv(hkp1k))

    end

    # fill last entry for error estimate 
    fill!(@view(H[1, m+1]), one(T))
    # remove last residual for error estimate form 
    fill!(@view(H[m+1, m]), zero(T))

    return m, hkp1k
end

function arnoldi_block4!(A, v, V, H, W, R; reorth::Bool = true)
    T = eltype(V)
    n = length(v)
    @assert size(W, 2) == 4
    @assert size(R, 1) == 4 && size(R, 2) == 4
    @assert size(V, 1) == n
    @assert size(H, 2) % 4 == 0

    nblocks = size(H, 2) ÷ 4
    @assert size(V, 2) >= 4 * (nblocks + 1)
    @assert size(H, 1) >= 4 * (nblocks + 1)

    fill!(H, zero(T))

    @views begin
        copyto!(W[:, 1], v)
        mul!(W[:, 2], A, W[:, 1])
        mul!(W[:, 3], A, W[:, 2])
        mul!(W[:, 4], A, W[:, 3])
    end

    F0 = qr!(W)
    copyto!(W, Matrix(F0.Q))
    copyto!(R, Matrix(F0.R))
    copyto!(@view(H[1:4, 1:4]), R)
    copyto!(@view(V[:, 1:4]), W)

    @inbounds for s in 1:nblocks
        c1 = 4 * (s - 1) + 1
        c2 = 4 * s
        n1 = c2 + 1
        n2 = 4 * (s + 1)

        @views begin
            mul!(W[:, 1], A, V[:, c2])
            mul!(W[:, 2], A, W[:, 1])
            mul!(W[:, 3], A, W[:, 2])
            mul!(W[:, 4], A, W[:, 3])

            Vprev = V[:, 1:c2]
            Hproj = H[1:c2, c1:c2]

            mul!(Hproj, adjoint(Vprev), W)
            mul!(W, Vprev, Hproj, -one(T), one(T))

            if reorth
                mul!(Hproj, adjoint(Vprev), W)
                mul!(W, Vprev, Hproj, -one(T), one(T))
            end
        end

        Fs = qr!(W)
        copyto!(W, Matrix(Fs.Q))
        copyto!(R, Matrix(Fs.R))

        copyto!(@view(V[:, n1:n2]), W)

        if s < nblocks # write residual block to H
            copyto!(@view(H[n1:n2, c1:c2]), R)
        end
    end

    # last residual block returned in R which is used for error estimate 

    return nothing
end