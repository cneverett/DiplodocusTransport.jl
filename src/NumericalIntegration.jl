"""
    DistributionToDIPIntegration(f,pxr,pyr,pzr,method,samples,DistributionFunction,Parameters...;kwargs...)

Numerically integrates a `DistributionFunction`` over momentum space domains to generate the discrete form of the distribution function used in DiplodocusTransport.jl. By default this `HCubature` from the HCubature.jl package, also by default the maximum number of function evaluations is 16^3 but can be increased by increasing `samples` (default is 16). `Parameters...` and `kwargs...` are any additional parameters and kwargs required by the `DistributionFunction`.
"""
function DistributionToDIPIntegration!(f::Array{Float64,3},pxr::Vector{Float64},pyr::Vector{Float64},pzr::Vector{Float64},method::String,samples::Int,DistributionFunction::Function,Parameters...;kwargs...)

    # Integrates a distribution function over momentum space using Simpson's rule
    # f is the 3D array to be filled with the distribution function values
    # DistributionFunction is a function that takes Parameters... and returns the distribution function f(vec{p})

    integrand = x->DistributionFunction(x[1],x[2],x[3],Parameters...;kwargs...)*x[1]^2

    for px in 1:length(pxr)-1, py in 1:length(pyr)-1, pz in 1:length(pzr)-1

        a = [pxr[px], pyr[py], pzr[pz]]
        b = [pxr[px+1], pyr[py+1], pzr[pz+1]]
        n = [samples, samples, samples]

        if method == "simpson"
            f[px,py,pz] += simpson3d(integrand,a,b,n)
        elseif method == "trapezium"
            f[px,py,pz] += trapz3d(integrand,a,b,n)
        elseif method == "hcubature"
            f[px,py,pz] += hcubature(integrand,a,b,rtol=1e-3,atol=1e-10,maxevals=samples^3,initdiv=8)[1]
        else
            error("Invalid integration method specified. Choose from 'simpson', 'trapezium', or 'hcubature'.")
        end

    end

    if sum(f) == 0.0
        error("Distribution function integration revealed no non-zero values. Check parameters.")
    end

    return nothing

end



# Optimized 3D Trapezoidal Rule Integration in Julia
#
# Integrates f(x, y, z) over the box:
#   x in [ax, bx], y in [ay, by], z in [az, bz]
#
# nx, ny, nz are the number of intervals in each dimension.

function trapz3d(f, a, b, n::Vector{Int})
    nx = n[1]
    ny = n[2]
    nz = n[3]

    ax = a[1]
    ay = a[2]
    az = a[3]
    bx = b[1]
    by = b[2]
    bz = b[3]

    nx < 1 && error("nx must be at least 1")
    ny < 1 && error("ny must be at least 1")
    nz < 1 && error("nz must be at least 1")

    hx = (bx - ax) / nx
    hy = (by - ay) / ny
    hz = (bz - az) / nz

    total = 0.0
    xyz= zeros(Float64,3)

    @inbounds for i in 0:nx
        xyz[1] = ax + i * hx
        wx = (i == 0 || i == nx) ? 1.0 : 2.0

        for j in 0:ny
            xyz[2] = ay + j * hy
            wy = (j == 0 || j == ny) ? 1.0 : 2.0

            for k in 0:nz
                xyz[3] = az + k * hz
                wz = (k == 0 || k == nz) ? 1.0 : 2.0

                total += wx * wy * wz * f(xyz)
            end
        end
    end

    return (hx * hy * hz / 8.0) * total
end

# Optimized 3D Simpson's Rule Integration in Julia
#
# Integrates f(x, y, z) over the box:
#   x in [ax, bx], y in [ay, by], z in [az, bz]
#
# nx, ny, nz must be even.

function simpson3d(f, a, b, n::Vector{Int})
    nx = n[1]
    ny = n[2]
    nz = n[3]

    ax = a[1]
    ay = a[2]
    az = a[3]
    bx = b[1]
    by = b[2]
    bz = b[3]

    nx < 2 && error("nx must be at least 2")
    ny < 2 && error("ny must be at least 2")
    nz < 2 && error("nz must be at least 2")

    iseven(nx) || error("nx must be even")
    iseven(ny) || error("ny must be even")
    iseven(nz) || error("nz must be even")

    hx = (bx - ax) / nx
    hy = (by - ay) / ny
    hz = (bz - az) / nz

    total = 0.0
    xyz = zeros(Float64,3)

    @inbounds for i in 0:nx
        xyz[1] = ax + i * hx
        wx = (i == 0 || i == nx) ? 1.0 : (isodd(i) ? 4.0 : 2.0)

        for j in 0:ny
            xyz[2] = ay + j * hy
            wy = (j == 0 || j == ny) ? 1.0 : (isodd(j) ? 4.0 : 2.0)

            for k in 0:nz
                xyz[3] = az + k * hz
                wz = (k == 0 || k == nz) ? 1.0 : (isodd(k) ? 4.0 : 2.0)

                total += wx * wy * wz * f(xyz)
            end
        end
    end

    return (hx * hy * hz / 27.0) * total
end