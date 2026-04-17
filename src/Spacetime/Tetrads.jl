abstract type AbstractTetrad end

""" 
    ZeroRicciTetrad 

A `ZeroRicciTetrad` is a tetrad that maps a coordinate system to a fixed set of basis vectors in Minkowski spacetime. By this all Ricci rotation coefficients are zero across the whole spacetime, simplifying the calculation of IJK fluxes.
"""
abstract type ElectromagneticTetrad <: AbstractTetrad end
abstract type ForceFreeTetrad <: AbstractTetrad end

"""
    StationaryObserverTetrad()

The simplest possible tetrad, well-defined for all stationary spacetimes, and is associated with the local frame of a stationary (Eularian) observer n=-Adt. Given a metric `g`:
```math
g = (-A^2+E^2/B^2+F^2/C^2+G^2/D^2)dt^2 + B^2dx^2 + C^2dy^2 + D^2dz^2 + 2E dtdx + 2F dtdy + 2G dtdz
```
the tetrad is given by 
```math
\\begin{align}
e_a^{~α} &= (n^α,X^α,Y^α,Z^α), \\
n^α &= (1/A, -E/(A*B^2), -F/(A*C^2), -G/(A*D^2)), \\
X^α &= (0, 1/B, 0, 0), \\
Y^α &= (0, 0, 1/C, 0), \\
Z^α &= (0, 0, 0, 1/D)
\\end{align}
"""
struct StationaryObserverTetrad <: AbstractTetrad end

"""
    UniformElectromagneticFieldTetrad(B0::Float64=1.0, E0::Float64=0.0)

A simple uniform force-free electromagnetic field configuration where the locally measured electric and magnetic fields are orthogonal such that this configuration if force-free. The magnetic field is taken to be in the Cartesian z direction and the Electric field in the Cartesian y direction, such that the ExB drift is in the Cartesian x direction. The strength of the magnetic field as measured by a **stationary observer n=-Adt** is given by the parameter B0 (in Tesla) and the strength of the electric field is given by the parameter E0 (in units of c*Tesla). By default, B0 is set to 1.0 Tesla and E0 is set to 0.0, which corresponds to a purely magnetic field with no electric field.
"""
@kwdef struct UniformElectromagneticFieldTetrad{T} <: ForceFreeTetrad
    B0::T = 1.0 # strength of the uniform magnetic field (Tesla)
    E0::T = 0.0 # strength of the uniform electric field (c*Tesla)
end

"""
    ParabolicForceFreeFieldTetrad(B0::Float64=1.0, E0::Float64=0.0)

The force-free field configuration as described by Blandford1976. This is an analytic field configuration for flat spacetime and an arbitrary field line rotation rate `Ω(ρ,z=0)=Ω(ρ0)` where ρ and z are cylindrical coordinates. 
"""
@kwdef struct ParabolicForceFreeFieldTetrad{T} <: ForceFreeTetrad
    B0::T = 1.0 # strength of the uniform magnetic field (Tesla)
    Ω = ΩZero # field line rotation rate as a function of cylindrical radius ρ0 in the equatorial plane (z=0)
end
@inline function ΩZero(ρ0) 
        return 0.0
end



#============= Tetrad Components =============#
#=============================================#
"""
    TetradComponents!(pos,e,::AbstractMetric,::AbstractCoordinates,::AbstractTetrad)

Returns the tetrad components ``e=e_a^{~α}`` at a given point `pos` in spacetime defined by the coordinates `pos=(t,x,y,z)`. The tetrad components are stored in the static matrix e, which are modified in-place. The specific form of the tetrad components depends on the type of metric, coordinate basis, and tetrad being used.

For force-free electromagnetic fields, the tetrad is given by e=(T^α,X^α,Y^α,Z^α), with ``Z^α=B^α/B`` and ``Y^α=E^α/E`` being the magnetic and electric field directions as measured by a static observer ``n_α``. The timelike vector ``T`` is the vector of an observer moving with the "ExB" drift velocity ``T=γ(n-U_⟂)`` where ``U_⟂=*(n∧E∧B)/B^2`` and ``γ=\\frac{B^2}{B^2-E^2}``. The final spacelike vector ``X`` is given by the orthonormality condition ``X=*(T∧Y∧Z)``.

The inverse tetrad is given by inve=transpose(-T_α,X_α,Y_α,Z_α).
"""
TetradComponents!(pos,e,metric::AbstractMetric,coordinates::AbstractCoordinates,tetrad::AbstractTetrad) = error("Tetrad function not defined for this metric $(typeof(metric)), coordinate $(typeof(coordinates)), and tetrad $(typeof(tetrad)).")

"""
    InverseTetradComponents!(pos,inve,::AbstractMetric,::AbstractCoordinates,::AbstractTetrad)

Returns the inverse tetrad components ``inve=inve^α_{~a}`` at a given point `pos` in spacetime defined by the coordinates `pos=(t,x,y,z)`. The tetrad components are stored in the static matrix inve, which are modified in-place. The specific form of the tetrad components depends on the type of metric, coordinate basis, and tetrad being used.

For force-free electromagnetic fields, the tetrad is given by e=(T^α,X^α,Y^α,Z^α), with ``Z^α=B^α/B`` and ``Y^α=E^α/E`` being the magnetic and electric field directions as measured by a static observer ``n_α``. The timelike vector ``T`` is the vector of an observer moving with the "ExB" drift velocity ``T=γ(n-U_⟂)`` where ``U_⟂=*(n∧E∧B)/B^2`` and ``γ=\\frac{B^2}{B^2-E^2}``. The final spacelike vector ``X`` is given by the orthonormality condition ``X=*(T∧Y∧Z)``.

The inverse tetrad is given by inve=transpose(-T_α,X_α,Y_α,Z_α).
"""
InverseTetradComponents!(pos,inve,metric::AbstractMetric,coordinates::AbstractCoordinates,tetrad::AbstractTetrad) = error("Tetrad function not defined for this metric $(typeof(metric)), coordinate $(typeof(coordinates)), and tetrad $(typeof(tetrad)).")

"""
    CoordinateFluxSpaceAIntegrand!(xyzt,A,::AbstractMetric,::AbstractCoordinates,::AbstractTetrad)

Returns the coordinate flux space integrand ``A=A_a=e_a^{~0}χ``, where ``e`` is the tetrad and ``χ`` is the volume element, at a given point in spacetime defined by the coordinates `xyzt=(x,y,z,t)` where the first three are to be integrated over and the last is fixed. The components are stored in the static vector A, which are modified in-place.
"""
CoordinateFluxSpaceAIntegrand!(xyzt,A,metric::AbstractMetric,coordinates::AbstractCoordinates,tetrad::AbstractTetrad) = error("CoordinateFluxSpaceAIntegrand function not defined for this metric $(typeof(metric)), coordinate $(typeof(coordinates)), and tetrad $(typeof(tetrad)).")

"""
    CoordinateFluxSpaceBIntegrand!(yztx,B,::AbstractMetric,::AbstractCoordinates,::AbstractTetrad)

Returns the coordinate flux space integrand ``B=B_a=e_a^{~1}χ``, where ``e`` is the tetrad and ``χ`` is the volume element, at a given point in spacetime defined by the coordinates `pos=(y,z,t,x)` where the first three are to be integrated over and the last is fixed. The components are stored in the static vector B, which are modified in-place.
"""
CoordinateFluxSpaceBIntegrand!(yztx,B,spacetime::AbstractMetric,coordinates::AbstractCoordinates,tetrad::AbstractTetrad) = error("CoordinateFluxSpaceBIntegrand function not defined for this spacetime $(typeof(spacetime)), coordinate $(typeof(coordinates)), and tetrad $(typeof(tetrad)).")

"""
    CoordinateFluxSpaceCIntegrand!(ztxy,C,::AbstractMetric,::AbstractCoordinates,::AbstractTetrad)

Returns the coordinate flux space integrand ``C=C_a=e_a^{~2}χ``, where ``e`` is the tetrad and ``χ`` is the volume element, at a given point in spacetime defined by the coordinates `ztxy=(z,t,x,y)` where the first three are to be integrated over and the last is the fixed coordinate. The components are stored in the static vector C, which are modified in-place.
"""
CoordinateFluxSpaceCIntegrand!(ztxy,C,metric::AbstractMetric,coordinates::AbstractCoordinates,tetrad::AbstractTetrad) = error("CoordinateFluxSpaceCIntegrand function not defined for this metric $(typeof(metric)), coordinate $(typeof(coordinates)), and tetrad $(typeof(tetrad)).")

"""
    CoordinateFluxSpaceDIntegrand!(txyz,D,::AbstractMetric,::AbstractCoordinates,::AbstractTetrad)

Returns the coordinate flux space integrand ``D=D_a=e_a^{~3}χ``, where ``e`` is the tetrad and ``χ`` is the volume element, at a given point in spacetime defined by the coordinates to be integrated over `pos=(t,x,y,z)`  where the first three are to be integrated over and the last is fixed. The components are stored in the static vector D, which are modified in-place.
"""
CoordinateFluxSpaceDIntegrand!(txyz,D,metric::AbstractMetric,coordinates::AbstractCoordinates,tetrad::AbstractTetrad) = error("CoordinateFluxSpaceDIntegrand function not defined for this metric $(typeof(metric)), coordinate $(typeof(coordinates)), and tetrad $(typeof(tetrad)).")

#================= Stationary Observer ==================#
#========================================================#
    # ============ Minkowski Cartesian ============ #
        @inline function TetradComponents!(pos::MVector{4,T},e::MMatrix{4,4,T,16},::Minkowski,::Cartesian,::StationaryObserverTetrad) where T

            #=     n = -dt     X = dx      Y = dy     Z = dz    =#
            # T components T^α = (1, 0, 0, 0)
            e[1,1] = 1.0
            # X components X^α = (0, 1, 0, 0)
            e[2,2] = 1.0
            # Y components Y^α = (0, 0, 1, 0)
            e[3,3] = 1.0
            # Z components Z^α = (0, 0, 0, 1)
            e[4,4] = 1.0

            return nothing

        end
        @inline function InverseTetradComponents!(pos::MVector{4,T},inve::MMatrix{4,4,T,16},::Minkowski,::Cartesian,::StationaryObserverTetrad) where T
            
            #=     n = -dt     X = dx      Y = dy     Z = dz     =#
            # T components -T_α = (1, 0, 0, 0)
            inve[1,1] = 1.0
            # X components X_α = (0, 1, 0, 0)
            inve[2,2] = 1.0
            # Y components Y_α = (0, 0, 1, 0)
            inve[3,3] = 1.0
            # Z components Z_α = (0, 0, 0, 1)
            inve[4,4] = 1.0

            return nothing

        end
        @inline function CoordinateFluxSpaceAIntegrand!(xyzt::MVector{4,T},A::MVector{4,T},::Minkowski,::Cartesian,::StationaryObserverTetrad) where T 
            #= A=A_a=e_a^{~0}χ =#
            A[1] = 1.0

        end
        @inline function CoordinateFluxSpaceBIntegrand!(yztx::MVector{4,T},B::MVector{4,T},::Minkowski,::Cartesian,::StationaryObserverTetrad) where T 
            #= B=B_a=e_a^{~1}χ =#
            B[2] = 1.0

        end
        @inline function CoordinateFluxSpaceCIntegrand!(ztxy::MVector{4,T},C::MVector{4,T},::Minkowski,::Cartesian,::StationaryObserverTetrad) where T 
            #= C=C_a=e_a^{~2}χ =#
            C[3] = 1.0

        end
        @inline function CoordinateFluxSpaceDIntegrand!(txyz::MVector{4,T},D::MVector{4,T},::Minkowski,::Cartesian,::StationaryObserverTetrad) where T 
            #= D=D_a=e_a^{~3}χ =#
            D[4] = 1.0

        end

    # ============ Minkowski Cylindrical ========== #
        @inline function TetradComponents!(pos::MVector{4,T},e::MMatrix{4,4,T,16},::Minkowski,::Cylindrical,::StationaryObserverTetrad) where T

            ρ = pos[2]
            #= n = -dt     X = dρ     Y = ρdϕ     Z = dz =#
            # T components T^α = g^αβT_β = (1, 0, 0, 0)
            e[1,1] = 1.0
            # X components X^α = g^αβX_β = (0, 1, 0, 0)
            e[2,2] = 1.0
            # Y components Y^α = g^αβY_β = (0, 0, 1/ρ, 0)
            e[3,3] = 1.0/ρ
            # Z components Z^α = g^αβZ_β = (0, 0, 0, 1)
            e[4,4] = 1.0

            return nothing

        end
        @inline function InverseTetradComponents!(pos::MVector{4,T},inve::MMatrix{4,4,T,16},::Minkowski,::Cylindrical,::StationaryObserverTetrad) where T

            ρ = pos[2]
            #=     n = -dt     X = dρ       Y = ρdϕ      Z = dz    =#
            # T components -T_α = (1, 0, 0, 0)
            inve[1,1] = 1.0
            # X components X_α = (0, 1, 0, 0)
            inve[2,2] = 1.0
            # Y components Y_α = (0, 0, ρ, 0)
            inve[3,3] = ρ
            # Z components Z_α = (0, 0, 0, 1)
            inve[4,4] = 1.0

            return nothing

        end
        @inline function CoordinateFluxSpaceAIntegrand!(xyzt::MVector{4,T},A::MVector{4,T},::Minkowski,::Cylindrical,::StationaryObserverTetrad) where T 
            #= A=A_a=e_a^{~0}χ =#
            ρ = xyzt[2]
            A[1] = ρ

        end
        @inline function CoordinateFluxSpaceBIntegrand!(yztx::MVector{4,T},B::MVector{4,T},::Minkowski,::Cylindrical,::StationaryObserverTetrad) where T 
            #= B=B_a=e_a^{~1}χ =#
            ρ = yztx[1]
            B[2] = ρ

        end
        @inline function CoordinateFluxSpaceCIntegrand!(ztxy::MVector{4,T},C::MVector{4,T},::Minkowski,::Cylindrical,::StationaryObserverTetrad) where T 
            #= C=C_a=e_a^{~2}χ =#
            C[3] = 1.0 # 1/ρ * ρ

        end
        @inline function CoordinateFluxSpaceDIntegrand!(txyz::MVector{4,T},D::MVector{4,T},::Minkowski,::Cylindrical,::StationaryObserverTetrad) where T 
            #= D=D_a=e_a^{~3}χ =#
            ρ = txyz[3]
            D[4] = ρ

        end

    # ============ Minkowski Spherical ============ #
        @inline function TetradComponents!(pos::MVector{4,T},e::MMatrix{4,4,T,16},::Minkowski,::Spherical,::StationaryObserverTetrad) where T
            # align momentum with radial coordinate
            t = pos[1]
            r = pos[2]
            θ = pos[3]
            ϕ = pos[4]
            #= n = -dt     X = rdθ     Y = r*sin(θ)*dϕ     Z = dr =#
            # T components T^α = g^αβT_β = (1, 0, 0, 0)
            e[1,1] = 1.0
            # X components X^α = g^αβX_β = (0, 0, 1/r, 0)
            e[2,3] = 1.0/r
            # Y components Y^α = g^αβY_β = (0, 0, 0, 1/rsinθ)
            e[3,4] = 1.0/(r*sin(θ))
            # Z components Z^α = g^αβZ_β = (0, 1, 0, 0)
            e[4,2] = 1.0

            return nothing

        end
        @inline function InverseTetradComponents!(pos::MVector{4,T},inve::MMatrix{4,4,T,16},::Minkowski,::Spherical,::StationaryObserverTetrad) where T

            t = pos[1]
            r = pos[2]
            θ = pos[3]
            ϕ = pos[4]
            #= n = -dt     X = rdθ     Y = r*sin(θ)*dϕ     Z = dr =#
            # T components -T_α = (1, 0, 0, 0)
            inve[1,1] = 1.0
            # X components X_α = (0, 0, r, 0)
            inve[3,2] = r
            # Y components Y_α = (0, 0, 0, r*sin(θ))
            inve[4,3] = r * sin(θ)
            # Z components Z_α = (0, 1, 0, 0)
            inve[2,4] = 1.0

            return nothing

        end
        @inline function CoordinateFluxSpaceAIntegrand!(xyzt::MVector{4,T},A::MVector{4,T},::Minkowski,::Spherical,::StationaryObserverTetrad) where T 
            #= A=A_a=e_a^{~0}χ =#
            r = xyzt[1]
            θ = xyzt[2]
            A[1] = r^2 * sin(θ)

        end
        @inline function CoordinateFluxSpaceBIntegrand!(yztx::MVector{4,T},B::MVector{4,T},::Minkowski,::Spherical,::StationaryObserverTetrad) where T 
            #= B=B_a=e_a^{~1}χ =#
            r = yztx[4]
            θ = yztx[1]
            B[4] = r^2 * sin(θ)

        end
        @inline function CoordinateFluxSpaceCIntegrand!(ztxy::MVector{4,T},C::MVector{4,T},::Minkowski,::Spherical,::StationaryObserverTetrad) where T 
            #= C=C_a=e_a^{~2}χ =#
            r = ztxy[3]
            θ = ztxy[4]
            C[2] = r * sin(θ)

        end
        @inline function CoordinateFluxSpaceDIntegrand!(txyz::MVector{4,T},D::MVector{4,T},::Minkowski,::Spherical,::StationaryObserverTetrad) where T 
            #= D=D_a=e_a^{~3}χ =#
            r = txyz[2]
            D[3] = r

        end
  
    # ============ Minkowski Paraboloidal ========= #
        @inline function TetradComponents!(pos::MVector{4,T},e::MMatrix{4,4,T,16},::Minkowski,::Paraboloidal,::StationaryObserverTetrad) where T

            t = pos[1]
            ϕ = pos[2]
            u = pos[3]
            v = pos[4]
            #= n = -dt     X = uvdϕ     Y = sqrt(u^2+v^2)dv     Z = sqrt(u^2+v^2)dv =#
            # T components T^α = g^αβT_β = (1, 0, 0, 0)
            e[1,1] = 1.0
            # X components X^α = g^αβX_β = (0, 1/uv, 0, 0)
            e[2,2] = 1.0/(u*v)
            # Y components Y^α = g^αβY_β = (0, 0, 1/sqrt(u^2+v^2), 0)
            e[3,3] = 1.0/sqrt(u^2+v^2)
            # Z components Z^α = g^αβZ_β = (0, 0, 0, 1/sqrt(u^2+v^2))
            e[4,4] = 1.0/sqrt(u^2+v^2)

            return nothing

        end
        @inline function InverseTetradComponents!(pos::MVector{4,T},inve::MMatrix{4,4,T,16},::Minkowski,::Paraboloidal,::StationaryObserverTetrad) where T

            t = pos[1]
            ϕ = pos[2]
            u = pos[3]
            v = pos[4]
            #= n = -dt     X = uvdϕ     Y = sqrt(u^2+v^2)dv     Z = sqrt(u^2+v^2)dv =#
            # T components -T_α = (1, 0, 0, 0)
            inve[1,1] = 1.0
            # X components X_α = (0, uv, 0, 0)
            inve[2,2] = u*v
            # Y components Y_α = (0, 0, sqrt(u^2+v^2), 0)
            inve[3,3] = sqrt(u^2+v^2)
            # Z components Z_α = (0, 0, 0, sqrt(u^2+v^2))
            inve[4,4] = sqrt(u^2+v^2)

            return nothing

        end
        @inline function CoordinateFluxSpaceAIntegrand!(xyzt::MVector{4,T},A::MVector{4,T},::Minkowski,::Paraboloidal,::StationaryObserverTetrad) where T 
            #= A=A_a=e_a^{~0}χ =#
            u = xyzt[2]
            v = xyzt[3]
            A[1] = u*v*(u^2+v^2)

        end
        @inline function CoordinateFluxSpaceBIntegrand!(yztx::MVector{4,T},B::MVector{4,T},::Minkowski,::Paraboloidal,::StationaryObserverTetrad) where T 
            #= B=B_a=e_a^{~1}χ =#
            u = yztx[1]
            v = yztx[2]
            B[2] = u^2+ v^2

        end
        @inline function CoordinateFluxSpaceCIntegrand!(ztxy::MVector{4,T},C::MVector{4,T},::Minkowski,::Paraboloidal,::StationaryObserverTetrad) where T 
            #= C=C_a=e_a^{~2}χ =#
            u = ztxy[4]
            v = ztxy[1]
            C[3] = u*v*sqrt(u^2+v^2)

        end
        @inline function CoordinateFluxSpaceDIntegrand!(txyz::MVector{4,T},D::MVector{4,T},::Minkowski,::Paraboloidal,::StationaryObserverTetrad) where T 
            #= D=D_a=e_a^{~3}χ =#
            u = txyz[3]
            v = txyz[4]
            D[4] = u*v*sqrt(u^2+v^2)

        end


#============= Uniform ElectromagneticField =============#
#========================================================#
    # ============ Minkowski Cartesian ============ #
        @inline function TetradComponents!(pos::MVector{4,T},e::MMatrix{4,4,T,16},::Minkowski,::Cartesian,tetrad::UniformElectromagneticFieldTetrad) where T

            # B field in Cartesian z direction, E field in Cartesian y direction
            B = tetrad.B0
            E = tetrad.E0
            #=     Y = dy      Z = dz      n = -dt      U_perp = (E/B)dx     T = γ(n-U_perp) = γ(-dt - (E/B)dx)     X = *(T∧Y∧Z) = γ(E/B)dt + γdx    =#
            γ = sqrt(B^2/(B^2 - E^2))
            v = E/B
            # T components T^α = (γ, -γE/B, 0, 0)
            e[1,1] = γ
            e[1,2] = -γ*v
            # X components X^α = (-γE/B, γ, 0, 0)
            e[2,1] = -γ*v
            e[2,2] = γ
            # Y components Y^α = (0, 0, 1, 0)
            e[3,3] = 1.0
            # Z components Z^α = (0, 0, 0, 1)
            e[4,4] = 1.0

            return nothing

        end
        @inline function InverseTetradComponents!(pos::MVector{4,T},inve::MMatrix{4,4,T,16},::Minkowski,::Cartesian,tetrad::UniformElectromagneticFieldTetrad) where T
            
            # B field in z direction, E field in y direction
            B = tetrad.B0
            E = tetrad.E0
            #=     Y = dy      Z = dz      n = -dt      U_perp = (E/B)dx     T = γ(n-U_perp) = γ(-dt - (E/B)dx)     X = *(T∧Y∧Z) = γ(E/B)dt + γdx    =#
            γ = sqrt(B^2/(B^2 - E^2))
            v = E/B
            # T components -T_α = (γ, γE/B, 0, 0)
            inve[1,1] = γ
            inve[2,1] = γ*v
            # X components X_α = (-γE/B, γ, 0, 0)
            inve[1,2] = -γ*v
            inve[2,2] = γ
            # Y components Y_α = (0, 0, 1, 0)
            inve[3,3] =  1.0
            # Z components Z_α = (0, 0, 0, 1)
            inve[4,4] = 1.0

            return nothing

        end
        @inline function CoordinateFluxSpaceAIntegrand!(xyzt::MVector{4,T},A::MVector{4,T},::Minkowski,::Cartesian,tetrad::UniformElectromagneticFieldTetrad) where T 
            # B field in Cartesian z direction, E field in Cartesian y direction
            B = tetrad.B0
            E = tetrad.E0
            #=     Y = dy      Z = dz      n = -dt      U_perp = (E/B)dx     T = γ(n-U_perp) = γ(-dt - (E/B)dx)     X = *(T∧Y∧Z) = γ(E/B)dt + γdx    =#
            γ = sqrt(B^2/(B^2 - E^2))
            v = E/B

            A[1] = γ
            A[2] = -γ*v

        end
        @inline function CoordinateFluxSpaceBIntegrand!(yztx::MVector{4,T},B::MVector{4,T},::Minkowski,::Cartesian,tetrad::UniformElectromagneticFieldTetrad) where T 
            # B field in Cartesian z direction, E field in Cartesian y direction
            B = tetrad.B0
            E = tetrad.E0
            #=     Y = dy      Z = dz      n = -dt      U_perp = (E/B)dx     T = γ(n-U_perp) = γ(-dt - (E/B)dx)     X = *(T∧Y∧Z) = γ(E/B)dt + γdx    =#
            γ = sqrt(B^2/(B^2 - E^2))
            v = E/B

            B[1] = -γ*v
            B[2] = γ 

        end
        @inline function CoordinateFluxSpaceCIntegrand!(ztxy::MVector{4,T},C::MVector{4,T},::Minkowski,::Cartesian,tetrad::UniformElectromagneticFieldTetrad) where T 

            C[3] = 1.0

        end
        @inline function CoordinateFluxSpaceDIntegrand!(txyz::MVector{4,T},D::MVector{4,T},::Minkowski,::Cartesian,tetrad::UniformElectromagneticFieldTetrad) where T 

            D[4] = 1.0

        end

    # ============ Minkowski Cylindrical ========== #
        @inline function TetradComponents!(pos::MVector{4,T},e::MMatrix{4,4,T,16},::Minkowski,::Cylindrical,tetrad::UniformElectromagneticFieldTetrad) where T

        ρ = pos[2]
        ϕ = pos[3]
        # B field in z direction, E field in y direction
        B = tetrad.B0
        E = tetrad.E0
        #=     Y = dy = sinϕ dρ + ρcosϕ dϕ      Z = dz      n = -dt      U_perp = (E/B)dx = (E/B)(cosϕ dρ - ρsinϕ dϕ)     T = γ(n-U_perp) = γ(-dt - (E/B)dx) = γ(-dt - (E/B)(cosϕ dρ - ρsinϕ dϕ))     X = *(T∧Y∧Z) = γ(E/B)dt + γdx = γ(E/B)dt + γ(cosϕ dρ - ρsinϕ dϕ)    =#
        γ = sqrt(B^2/(B^2 - E^2))
        v = E/B
        sϕ,cϕ = sincos(ϕ)
        # T components T^α = g^αβT_β = (γ, -γE/Bcosϕ, +γE/Bsinϕ/ρ, 0)
        e[1,1] = γ
        e[1,2] = -γ*v*cϕ
        e[1,3] = γ*v*sϕ/ρ
        # X components X^α = g^αβX_β (-γE/B, γcosϕ, -γsinϕ/ρ, 0)
        e[2,1] = -γ*v
        e[2,2] = γ*cϕ
        e[2,3] = -γ*sϕ/ρ
        # Y components Y^α = g^αβY_β (0, sinϕ, cosϕ/ρ, 0)
        e[3,2] = sϕ
        e[3,3] = cϕ/ρ
        # Z components Z^α = g^αβZ_β = (0, 0, 0, 1)
        e[4,4] = 1.0

        return nothing

        end
        @inline function InverseTetradComponents!(pos::MVector{4,T},inve::MMatrix{4,4,T,16},::Minkowski,::Cylindrical,tetrad::UniformElectromagneticFieldTetrad) where T

            ρ = pos[2]
            ϕ = pos[3]
            # B field in z direction, E field in y direction
            B = tetrad.B0
            E = tetrad.E0
            #=     Y = dy = sinϕ dρ + ρcosϕ dϕ      Z = dz      n = -dt      U_perp = (E/B)dx = (E/B)(cosϕ dρ - ρsinϕ dϕ)     T = γ(n-U_perp) = γ(-dt - (E/B)dx) = γ(-dt - (E/B)(cosϕ dρ - ρsinϕ dϕ))     X = *(T∧Y∧Z) = γ(E/B)dt + γdx = γ(E/B)dt + γ(cosϕ dρ - ρsinϕ dϕ)    =#
            γ = sqrt(B^2/(B^2 - E^2))
            v = E/B
            sϕ,cϕ = sincos(ϕ)
            # T components -T_α = (γ, γE/Bcosϕ, -γρE/Bsinϕ, 0)
            inve[1,1] = γ
            inve[2,1] = γ*v*cϕ
            inve[3,1] = -γ*v*ρ*sϕ
            # X components X_α = (γE/B, γcosϕ, -γρsinϕ, 0)
            inve[1,2] = γ*v
            inve[2,2] = γ*cϕ
            inve[3,2] = -γ*ρ*sϕ
            # Y components Y_α = (0, sinϕ, ρcosϕ, 0)
            inve[2,3] = sϕ
            inve[3,3] = ρ*cϕ
            # Z components Z_α = (0, 0, 0, 1)
            inve[4,4] = 1.0

            return nothing

        end
        @inline function CoordinateFluxSpaceAIntegrand!(xyzt::MVector{4,T},A::MVector{4,T},::Minkowski,::Cylindrical,tetrad::UniformElectromagneticFieldTetrad) where T 
            ρ = pos[2]
            B = tetrad.B0
            E = tetrad.E0
            γ = sqrt(B^2/(B^2 - E^2))
            v = E/B
            A[1] = γ*ρ
            A[2] = -γ*v*ρ


        end
        @inline function CoordinateFluxSpaceBIntegrand!(yztx::MVector{4,T},B::MVector{4,T},::Minkowski,::Cylindrical,tetrad::UniformElectromagneticFieldTetrad) where T 
            ρ = pos[2]
            ϕ = pos[3]
            B = tetrad.B0
            E = tetrad.E0
            γ = sqrt(B^2/(B^2 - E^2))
            v = E/B
            sϕ,cϕ = sincos(ϕ)
            B[1] = -γ*v*cϕ*ρ
            B[2] = γ*cϕ*ρ
            B[3] = sϕ*ρ

        end
        @inline function CoordinateFluxSpaceCIntegrand!(ztxy::MVector{4,T},C::MVector{4,T},::Minkowski,::Cylindrical,tetrad::UniformElectromagneticFieldTetrad) where T 
            ϕ = pos[3]
            B = tetrad.B0
            E = tetrad.E0
            γ = sqrt(B^2/(B^2 - E^2))
            v = E/B
            sϕ,cϕ = sincos(ϕ)
            C[1] = γ*v*sϕ
            C[2] = -γ*sϕ
            C[3] = cϕ
        end
        @inline function CoordinateFluxSpaceDIntegrand!(txyz::MVector{4,T},D::MVector{4,T},::Minkowski,::Cylindrical,tetrad::UniformElectromagneticFieldTetrad) where T 
            ρ = pos[2]
            D[4] = ρ
        end

#=============== ParabolicForceFreeField ================#
#========================================================#
    # ============ Minkowski Paraboloidal ========= #
        @inline function TetradComponents!(pos::MVector{4,T},e::MMatrix{4,4,T,16},::Minkowski,::Paraboloidal,tetrad::ParabolicForceFreeFieldTetrad) where T

            fill!(e, zero(T)) # needed as components that are non-zero change with sign z
            t = pos[1]
            ϕ = pos[2]
            u = pos[3]
            v = pos[4]
            #= See Mathematica notebook for details of how these are derived. T corresponds to the drifting observer with T = γ(n-U⟂) then Y is in E field direction, Z in the B field direction and X orthogonal to T,X,Y
            These directions also depend on the sign of the z=(u^2-v^2)/2 cylindrical coordinate =#
            signz = sign((u^2-v^2))
            @inline Ω = tetrad.Ω
            if signz == 1
                # T components T^α = g^αβT_β
                e[1,1] = sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2)/sqrt(1 + v^4*Ω(v)^2)
                e[1,2] = Ω(v)/(sqrt(1 + v^4*Ω(v)^2) * sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2))
                e[1,3] = u*v^2*Ω(v)^2/(sqrt(1 + v^4*Ω(v)^2) * sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2))
                # X components X^α = g^αβX_β
                e[2,1] = -u*v*Ω(v)/sqrt(1 + v^4*Ω(v)^2)
                e[2,2] = -1/(u*v*sqrt(1 + v^4*Ω(v)^2))
                e[2,3] = -v*Ω(v)/sqrt(1 + v^4*Ω(v)^2)
                # Y components Y^α = g^αβY_β
                e[3,4] = -1/sqrt(u^2 + v^2)
                # Z components Z^α = g^αβZ_β
                e[4,2] = -Ω(v)*sqrt(u^2 + v^2)/(u*sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2))
                e[4,3] = 1/(sqrt(u^2 + v^2) * sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2))
            else
                # T components T^α = g^αβT_β
                e[1,1] = sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2)/sqrt(1 + u^4*Ω(u)^2)
                e[1,2] = Ω(u)/(sqrt(1 + u^4*Ω(u)^2) * sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2))
                e[1,4] = -u^2*v*Ω(u)^2/(sqrt(1 + u^4*Ω(u)^2) * sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2))
                # X components X^α = g^αβX_β
                e[2,1] = -u*v*Ω(u)/sqrt(1 + u^4*Ω(u)^2)
                e[2,2] = -1/(u*v*sqrt(1 + u^4*Ω(u)^2))
                e[2,4] = u*Ω(u)/sqrt(1 + u^4*Ω(u)^2)
                # Y components Y^α = g^αβY_β
                e[3,4] = -1/sqrt(u^2 + v^2)
                # Z components Z^α = g^αβZ_β
                e[4,2] = -Ω(u)*sqrt(u^2 + v^2)/(v*sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2))
                e[4,4] = -1/(sqrt(u^2 + v^2) * sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2))
            end

            return nothing

        end
        @inline function InverseTetradComponents!(pos::MVector{4,T},inve::MMatrix{4,4,T,16},::Minkowski,::Paraboloidal,tetrad::ParabolicForceFreeFieldTetrad) where T

            fill!(inve, zero(T)) # needed as components that are non-zero change with sign z
            t = pos[1]
            ϕ = pos[2]
            u = pos[3]
            v = pos[4]
            #= See Mathematica notebook for details of how these are derived. T corresponds to the drifting observer with T = γ(n-U⟂) then Y is in E field direction, Z in the B field direction and X orthogonal to T,X,Y
            These directions also depend on the sign of the z=(u^2-v^2)/2 cylindrical coordinate =#
            signz = sign((u^2-v^2))
            @inline Ω = tetrad.Ω
            if signz == 1
                # T components -T_α = (1, 0, 0, 0)
                inve[1,1] = sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2)/sqrt(1 + v^4*Ω(v)^2)
                inve[2,1] = -u^2*v^2*Ω(v)/(sqrt(1 + v^4*Ω(v)^2)*sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2))
                inve[3,1] = -u*v^2*(u^2 + v^2)*Ω(v)^2/(sqrt(1 + v^4*Ω(v)^2) * sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2))
                # X components X_α = (0, uv, 0, 0)
                inve[1,2] = u*v*Ω(v)/sqrt(1 + v^4*Ω(v)^2)
                inve[2,2] = -u*v/sqrt(1 + v^4*Ω(v)^2)
                inve[3,2] = -v*(u^2 + v^2)*Ω(v)/sqrt(1 + v^4*Ω(v)^2)            
                # Y components Y_α = (0, 0, sqrt(u^2+v^2), 0)
                inve[4,3] = -sqrt(u^2 + v^2)
                # Z components Z_α = (0, 0, 0, sqrt(u^2+v^2))
                inve[2,4] = -u*v^2*Ω(v)*sqrt(u^2 + v^2)/sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2)
                inve[3,4] = sqrt(u^2 + v^2)/sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2)
            else
                # T components -T_α = (1, 0, 0, 0)
                inve[1,1] = sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2)/sqrt(1 + u^4*Ω(u)^2)
                inve[2,1] = -u^2*v^2*Ω(u)/(sqrt(1 + u^4*Ω(u)^2)*sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2))
                inve[4,1] = u^2*v*(u^2 + v^2)*Ω(u)^2/(sqrt(1 + u^4*Ω(u)^2) * sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2))
                # X components X_α = (0, uv, 0, 0)
                inve[1,2] = u*v*Ω(u)/sqrt(1 + u^4*Ω(u)^2)
                inve[2,2] = -u*v/sqrt(1 + u^4*Ω(u)^2)
                inve[4,2] = u*(u^2 + v^2)*Ω(u)/sqrt(1 + u^4*Ω(u)^2)            
                # Y components Y_α = (0, 0, sqrt(u^2+v^2), 0)
                inve[3,3] = -sqrt(u^2 + v^2)
                # Z components Z_α = (0, 0, 0, sqrt(u^2+v^2))
                inve[2,4] = -u^2*v*Ω(u)*sqrt(u^2 + v^2)/sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2)
                inve[4,4] = -sqrt(u^2 + v^2)/sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2)
            end

            return nothing

        end
        @inline function CoordinateFluxSpaceAIntegrand!(xyzt::MVector{4,T},A::MVector{4,T},::Minkowski,::Paraboloidal,tetrad::ParabolicForceFreeFieldTetrad) where T 
            #= A=A_a=e_a^{~0}χ =#
            u = xyzt[2]
            v = xyzt[3]
            Ω = tetrad.Ω
            signz = sign((u^2-v^2))
            if signz == 1
                A[1] = u*v*(u^2 + v^2)*sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2)/sqrt(1 + v^4*Ω(v)^2)
                A[2] = -u^2*v^2*(u^2 + v^2)*Ω(v)/sqrt(1 + v^4*Ω(v)^2)
            else
                A[1] = u*v*(u^2 + v^2)*sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2)/sqrt(1 + u^4*Ω(u)^2)
                A[2] = -u^2*v^2*(u^2 + v^2)*Ω(u)/sqrt(1 + u^4*Ω(u)^2)
            end
        end
        @inline function CoordinateFluxSpaceBIntegrand!(yztx::MVector{4,T},B::MVector{4,T},::Minkowski,::Paraboloidal,tetrad::ParabolicForceFreeFieldTetrad) where T 
            #= B=B_a=e_a^{~1}χ =#
            u = yztx[1]
            v = yztx[2]
            Ω = tetrad.Ω
            signz = sign((u^2-v^2))
            if signz == 1
                B[1] = u*v*(u^2 + v^2)*Ω(v)/(sqrt(1 + v^4*Ω(v)^2) * sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2))
                B[2] = -(u^2 + v^2)/sqrt(1 + v^4*Ω(v)^2)
                B[4] = -v*(u^2 + v^2)^(3/2)*Ω(v)/sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2)
            else
                B[1] = u*v*(u^2 + v^2)*Ω(u)/(sqrt(1 + u^4*Ω(u)^2) * sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2))
                B[2] = -(u^2 + v^2)/sqrt(1 + u^4*Ω(u)^2)
                B[4] = -u*(u^2 + v^2)^(3/2)*Ω(u)/sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2)
            end
        end
        @inline function CoordinateFluxSpaceCIntegrand!(ztxy::MVector{4,T},C::MVector{4,T},::Minkowski,::Paraboloidal,tetrad::ParabolicForceFreeFieldTetrad) where T 
            #= C=C_a=e_a^{~2}χ =#
            fill!(C, zero(T)) # needed as components that are non-zero change with sign z
            u = ztxy[4]
            v = ztxy[1]
            Ω = tetrad.Ω
            signz = sign((u^2-v^2))
            if signz == 1
                C[1] = u^2*v^3*(u^2 + v^2)*Ω(v)^2/(sqrt(1 + v^4*Ω(v)^2) * sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2))
                C[2] = -u*v^2*(u^2 + v^2)*Ω(v)/sqrt(1 + v^4*Ω(v)^2)
                C[4] = u*v*sqrt(u^2 + v^2)/sqrt(1 + v^2*(u^2 + v^2)*Ω(v)^2)
            else
                C[3] = -u*v*sqrt(u^2 + v^2)
            end

        end
        @inline function CoordinateFluxSpaceDIntegrand!(txyz::MVector{4,T},D::MVector{4,T},::Minkowski,::Paraboloidal,tetrad::ParabolicForceFreeFieldTetrad) where T 
            #= D=D_a=e_a^{~3}χ =#
            fill!(D, zero(T)) # needed as components that are non-zero change with sign z
            u = txyz[3]
            v = txyz[4]
            Ω = tetrad.Ω
            signz = sign((u^2-v^2))
            if signz == 1
                D[3] = -u*v*sqrt(u^2 + v^2)
            else
                D[1] = -u^3*v^2*(u^2 + v^2)*Ω(u)^2/(sqrt(1 + u^4*Ω(u)^2) * sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2))
                D[2] = u^2*v*(u^2 + v^2)*Ω(u)/sqrt(1 + u^4*Ω(u)^2)
                D[4] = -u*v*sqrt(u^2 + v^2)/sqrt(1 + u^2*(u^2 + v^2)*Ω(u)^2)
            end
        end
