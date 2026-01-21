module RootFinding

using Optim
using LinearAlgebra

using ..Parameters
using ..Laser

export vL, find_roots, find_roots_log10

function δc(vz, vd; Δ, params::SystemParams)
    ħ = params.ħ
    α = params.α

    Δ -α/ħ * sum(Ec(0, yi, zi; Δ=Δ, params=params)^2 for(yi, zi) in zip(vd, vz))
end

function α_c_eq(vz, vd, vϕ; Δ, κ, params::SystemParams)
    ħ = params.ħ

    α = params.α

    δc_ = δc(vz, vd; Δ=Δ, params=params)
    α/ħ * sum(Ec(0, yi, zi; Δ=Δ, params=params)*Et(0, 0, zi; params=params)*exp(im*ϕi) for(yi, zi, ϕi) in zip(vd, vz, vϕ))/(δc_ - im*κ/2)
end

function vL(vz, vd, vϕ; Δ, κ, params::SystemParams)
    c = params.c
    k0 = params.k0
    ω0 = params.ω0
    kc = let
        ωc = Δ + ω0
        ωc/c
    end
    Wc = params.Wc
    zR = params.zR

    fz(z) = (im*(k0 - zR/(z^2 + zR^2)) - z/(z^2 + zR^2))
    gz(y, z) = -2*z/Wc^2 
    α_c = α_c_eq(vz, vd, vϕ; Δ=Δ, κ=κ, params=params)
    cα_c = conj(α_c)

    [begin
          t1 = (
                cα_c * Ec(0, yi, zi; Δ=Δ, params=params) 
                + conj(Et(0, 0, zi; params=params))*exp(-im*ϕi)
               )
          t2 = (
                α_c * gz(yi, zi)*Ec(0,yi, zi; Δ=Δ, params=params)
                + fz(zi)*Et(0, 0, zi; params=params) * exp(im * ϕi)
          )
          real(t1*t2)
     end
    for (yi, zi, ϕi) in zip(vd, vz, vϕ)]
end


function find_roots(zguess, vd, vϕ; Δ, κ, params::SystemParams, method=NelderMead())

    vL_ = zg-> norm(vL(zg ./ 1e6, vd, vϕ; Δ=Δ, κ=κ, params=params))
    root = optimize(vL_, collect(zguess), method)

    minimum(root)
end

end # module RootFinding
