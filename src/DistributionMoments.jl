"""
    FourFlow(u,dp,du,dE,du2)

Returns the four-flow vector Ua 'vector{Float64}' from the flattened axisymmetric distribution function f1D.
"""
function FourFlow(f1D::Vector{T},p_num,u_num,h_num,pr,ur,hr,m) where T<:AbstractFloat

    f3D = zeros(Float64,p_num,u_num,h_num)
    f3D = Float64.(reshape(f1D,(p_num,u_num,h_num)))

    Na = zeros(Float64,4)
    
    for p in 1:p_num, u in 1:u_num, h in 1:h_num

        Na[1] += f3D[p,u,h]
        Na[2] += ((sqrt(m^2 + pr[p]^2) - sqrt(m^2 + pr[p+1]^2)) * (ur[u]*sqrt(1 - ur[u]^2) - ur[u+1]*sqrt(1 - ur[u+1]^2) - 2acot_mod(ur[u]) + 2acot_mod(ur[u+1])) * (sin(hr[h]) - sin(hr[h+1])))/(2(pr[p] - pr[p+1])*(ur[u] - ur[u+1])*(hr[h] - hr[h+1]))
        Na[3] += ((sqrt(m^2 + pr[p]^2) - sqrt(m^2 + pr[p+1]^2)) * (ur[u]*sqrt(1 - ur[u]^2) - ur[u+1]*sqrt(1 - ur[u+1]^2) - 2acot_mod(ur[u]) + 2acot_mod(ur[u+1])) * (cos(hr[h]) - cos(hr[h+1])))/(2(pr[p] - pr[p+1])*(ur[u] - ur[u+1])*(hr[h] - hr[h+1]))
        Na[4] += (sqrt(pr[p+1]^2+m^2) - sqrt(pr[p]^2+m^2)) * (ur[u+1]+ur[u]) / (pr[p+1]-pr[p]) / 2

    end

    return Na

end

function HydroFourVelocity(Na::Vector{Float64})

    Ua = zeros(Float64,4)

    metric = zeros(Float64,4,4)
    metric[1,1] = -1
    metric[2,2] = 1
    metric[3,3] = 1
    metric[4,4] = 1

    sqrtNa2 = sqrt(abs(Na'*metric*Na))
    Ua .= Na/sqrtNa2

    return Ua

end

function ProjectionTensor(Uₐ::Vector{Float64})

    Δab = zeros(Float64,4,4)

    metric = zeros(Float64,4,4)
    metric[1,1] = -1
    metric[2,2] = 1
    metric[3,3] = 1
    metric[4,4] = 1

    for i in 1:4, j in 1:4
        Δab[i,j] = metric[i,j] + Uₐ[i]*Uₐ[j]
    end

    return Δab

end

function StressEnergyTensor(f1D::Vector{T},p_num,u_num,h_num,pr,ur,hr,m) where T<:AbstractFloat

    f3D = zeros(Float64,p_num,u_num,h_num)
    f3D = Float64.(reshape(f1D,(p_num,u_num,h_num)))

    #dp = zeros(Float64,p_num)
    #du = zeros(Float64,u_num)
    #dp2 = zeros(Float64,p_num)
    #du2 = zeros(Float64,u_num)
    #du3 = zeros(Float64,u_num)
    #duplusu3 = zeros(Float64,u_num)
    #dpfunc = zeros(Float64,p_num)
    #dE = DC.deltaEVector(pr,m)

    pr64 = Float64.(pr)

    Tab = zeros(Float64,4,4)

    #=for i in 1:u_num
        du[i] = ur[i+1]-ur[i]
        du2[i] = (ur[i+1]^2-ur[i]^2)/2
        du3[i] = (ur[i+1]^3-ur[i]^3)/3
        duplusu3[i] = ur[i+1]-ur[i] - (ur[i+1]^3-ur[i]^3)/3
    end
    for i in 1:p_num 
        dp[i] = (pr[i+1]-pr[i])
        dp2[i] = (pr[i+1]^2-pr[i]^2)/2
        if m == 0e0
            dpfunc[i] = pr64[i+1]^2 - pr64[i]^2
        else
            if pr[i+1]/m < 1e-3
                dpfunc[i] = 2/3 * pr[i+1]^3/m - 1/5 * pr[i+1]^5/m^3
                dpfunc[i] -= 2/3 * pr[i]^3/m - 1/5 * pr[i]^5/m^3
            else
                dpfunc[i] = pr64[i+1]*sqrt(pr64[i+1]^2+m^2) - m^2*atanh(pr64[i+1]/sqrt(pr64[i+1]^2+m^2))
                dpfunc[i] -= pr64[i]*sqrt(pr64[i]^2+m^2) - m^2*atanh(pr64[i]/sqrt(pr64[i]^2+m^2))
            end
        end

    end
    =#

    for p in 1:p_num, u in 1:u_num, h in 1:h_num

        Tab[1,1] += (pr[p]sqrt(m^2 + pr[p]^2) - pr[p+1]sqrt(m^2 + pr[p+1]^2) + m^2(asinh(pr[p]/m) - asinh(pr[p+1]/m)))/(2(pr[p] - pr[p+1]))

        Tab[1,2] += ((pr[p] + pr[p+1])*(ur[u]sqrt(1 - ur[u]^2) - ur[u+1]sqrt(1 - ur[u+1]^2) - 2acot_mod(ur[u]) + 2acot_mod(ur[u+1]))*(sin(hr[h]) - sin(hr[h+1])))/(4(ur[u] - ur[u+1])*(hr[h] - hr[h+1]))
        Tab[1,3] += -(((pr[p] + pr[p+1])*(ur[u]sqrt(1 - ur[u]^2) - ur[u+1]sqrt(1 - ur[u+1]^2) - 2acot_mod(ur[u]) + 2acot_mod(ur[u+1]))*(cos(hr[h]) - cos(hr[h+1])))/(4(ur[u] - ur[u+1])*(hr[h] - hr[h+1])))
        Tab[1,4] += 1/4 * (pr[p] + pr[p+1])*(ur[u] + ur[u+1])

        Tab[2,2] += -(((-3 + ur[u]^2 + ur[u]ur[u+1] + ur[u+1]^2)*(pr[p]sqrt(m^2 + pr[p]^2) - pr[p+1]sqrt(m^2 + pr[p+1]^2) + m^2 * (-asinh(pr[p]/m) + asinh(pr[p+1]/m)))*(hr[h] - hr[h+1] + cos(hr[h])sin(hr[h]) - cos(hr[h+1])sin(hr[h+1])))/(12(pr[p] - pr[p+1])*(hr[h] - hr[h+1])))

        Tab[2,3] += ((-3 + ur[u]^2 + ur[u]ur[u+1] + ur[u+1]^2)*(pr[p]sqrt(m^2 + pr[p]^2) - pr[p+1]sqrt(m^2 + pr[p+1]^2) + m^2 * (-asinh(pr[p]/m) + asinh(pr[p+1]/m)))*(cos(2hr[h]) - cos(2hr[h+1])))/(24(pr[p] - pr[p+1])*(hr[h] - hr[h+1]))

        Tab[3,3] += -(((-3 + ur[u]^2 + ur[u]ur[u+1] + ur[u+1]^2)*(pr[p]sqrt(m^2 + pr[p]^2) - pr[p+1]sqrt(m^2 + pr[p+1]^2) + m^2 * (-asinh(pr[p]/m) + asinh(pr[p+1]/m)))*(hr[h] - hr[h+1] - cos(hr[h])sin(hr[h]) + cos(hr[h+1])sin(hr[h+1])))/(12(pr[p] - pr[p+1])*(hr[h] - hr[h+1])))

        Tab[3,4] += -(((-sqrt(1 - ur[u]^2) + ur[u]^2 * sqrt(1 - ur[u]^2) + sqrt(1 - ur[u+1]^2) - ur[u+1]^2 * sqrt(1 - ur[u+1]^2))*(pr[p]sqrt(m^2 + pr[p]^2) - pr[p+1]sqrt(m^2 + pr[p+1]^2) + m^2 * (-asinh(pr[p]/m) + asinh(pr[p+1]/m)))*(cos(hr[h]) - cos(hr[h+1])))/(6(pr[p] - pr[p+1])*(ur[u] - ur[u+1])*(hr[h] - hr[h+1])))

        Tab[4,4] += (1/(6(pr[p] - pr[p+1])))*(ur[u]^2 + ur[u]ur[u+1] + ur[u+1]^2)*(pr[p]sqrt(m^2 + pr[p]^2) - pr[p+1]sqrt(m^2 + pr[p+1]^2) + m^2 * (-asinh(pr[p]/m) + asinh(pr[p+1]/m)))

    end

    Tab[2,1] = Tab[1,2]
    Tab[3,1] = Tab[1,3]
    Tab[4,1] = Tab[1,4]

    Tab[3,2] = Tab[2,3]
    Tab[4,2] = Tab[2,4]

    Tab[4,3] = Tab[3,4]

    # unscale by dp*du
    #for i in axes(f2D,1), j in axes(f2D,2)
    #    f2D[i,j] /= dp[i] * du[j]
    #end

    #Tab[1,1] = dE' * f2D * du
    #Tab[2,2] = 1/4 * dpfunc' * f2D * duplusu3
    #Tab[3,3] = Tab[2,2]
    #Tab[4,4] = 1/2 * dpfunc' * f2D * du3 
    #Tab[1,4] = dp2' * f2D * du2
    #Tab[4,1] = Tab[1,4]

    # unscale by dp*du
    #Tab[1,1] = (dE./dp)' * f2D * (du./du)
    #Tab[2,2] = 1/4 * (dpfunc./dp)' * f2D * (duplusu3./du)
    #Tab[3,3] = Tab[2,2]
    #Tab[4,4] = 1/2 * (dpfunc./dp)' * f2D * (du3./du) 
    #Tab[1,4] = (dp2./dp)' * f2D * (du2./du)
    #Tab[4,1] = Tab[1,4]
    
    return Tab
end

function ScalarNumberDensity(Nᵃ::Vector{Float64},Uₐ::Vector{Float64})

    n::Float64 = 0.0

    if Nᵃ[1] == 0
        n = 0.0
    else
        n = - dot(Nᵃ,Uₐ)
    end

    return n

end

function ScalarEnergyDensity(Tᵃᵇ::Matrix{Float64},Uₐ::Vector{Float64},n::Float64)

    e::Float64 = 0.0
    en::Float64 = Uₐ' * Tᵃᵇ * Uₐ

    if n == 0.0
        e = 0.0
    else
        e = en/n
    end

    return en # not sure if this is correct 

end

function ScalarPressure(Tᵃᵇ::Matrix{Float64},Δab::Matrix{Float64})

    p::Float64 = 0.0

    metric = zeros(Float64,4,4)
    metric[1,1] = -1
    metric[2,2] = 1
    metric[3,3] = 1
    metric[4,4] = 1

    #p = (1/3) * sum(Tab .* (metric * Δab * metric))
    p = (1/3) * sum(Tᵃᵇ .* Δab)

    return p

end

function ScalarTemperature(p::Float64,n::Float64)

    T::Float64 = 0.0

    kb = 1.38f-23
    c = 3f8
    mEle = 9.11e-31
    T = p/(n*kb) * mEle * c^2

    return T

end


