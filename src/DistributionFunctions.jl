"""
    MaxwellJuttner_Distribution(PhaseSpace,species,T;n=1.0)

Generates a Maxwell-Juttner distribution in momentum space: 
```math
f(\\vec{p})\\mathrm{d}^3\\vec{p} = f(p,u,ϕ)\\mathrm{d}p\\mathrm{d}u\\mathrm{d}ϕ = \\frac{p^2}{4π m^3c^3θK_2(1/θ)}e^{-γ(p)/θ}\\mathrm{d}p\\mathrm{d}u\\mathrm{d}ϕ
```
where ``γ(p) = \\sqrt{1+(p/mc)^2}`` and ``1/θ = m c^2/(k_B T)``. 

As f(p,u,ϕ)ΔpΔuΔϕ is the normalised distribution function in a single momentum-space bin. The value of the distribution function in each bin is given by the mean value of the Maxwell-Juttner distribution at the mean momentum of the bin. 

"""
function MaxwellJuttner_Distribution(PhaseSpace::PhaseSpaceStruct,species::String,T::Float64;n::Float64=1e0)
    # Generates the height of the MJ distribution at positions of the mean momentum per bin

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    species_index = findfirst(x->x==species,name_list)

    dp = Grids.dpx_list[species_index]
    #du = Grids.dpy_list[species_index] # angle averaged therefore not used and dudh = 4π
    #dh = Grids.dpz_list[species_index]
    meanp = Grids.mpx_list[species_index]
    mass = Grids.mass_list[species_index]

    mEle = 9.11e-31
    c = 3e8
    kb = 1.38e-23

    MJPoints = zeros(Float64,size(meanp))
    # besselk doesn't do well with small T besselkx = besselk * e^x does better.
    invtheta = mEle*c^2/(kb*T)
    @. MJPoints = (meanp^2)*(invtheta/(2*mass^2)) * (1/besselkx(2,mass*invtheta)) * exp((mass-sqrt(mass^2+meanp^2))*invtheta) * dp

    num = sum(MJPoints)
    MJPoints = MJPoints ./ num .*n # scale to number density

    return MJPoints

end

function BlackBody_Distribution(PhaseSpace::PhaseSpaceStruct,species::String,T::Float64;n::Float64=1e0)
    # Generates the height of the MJ distribution at positions of the mean momentum per bin

    name_list = PhaseSpace.name_list
    Momentum = PhaseSpace.Momentum
    Grids = PhaseSpace.Grids

    species_index = findfirst(x->x==species,name_list)

    dp = Grids.dpx_list[species_index]
    #du = Grids.dpy_list[species_index] # angle averaged therefore not used and dudh = 4π
    #dh = Grids.dpz_list[species_index]
    meanp = Grids.mpx_list[species_index]
    mass = Grids.mass_list[species_index]

    mEle = 9.11e-31
    c = 3e8
    kb = 1.38e-23
    h = 6.62607015e-34

    BBPoints = zeros(Float64,size(meanp))
    invtheta = mEle*c^2/(kb*T)
    @. BBPoints = (meanp^2) / (exp(meanp*invtheta)-1) * dp

    num = sum(MJPoints)
    BBPoints = BBPoints ./ num .*n # scale to number density

    return BBPoints

end