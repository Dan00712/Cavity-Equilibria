module Laser

using ..Parameters

export Ec, Et

function E0c(;Δ, params::SystemParams)
    ħ = params.ħ
    ϵ0 = params.ϵ0

    ω0 = params.ω0
    Vc = params.Vc

    # Δω = ωc - ω0
    ωc = Δ + ω0

    sqrt(ħ*ωc/2/ϵ0/Vc)
end


function Ec(x, y, z; Δ, params::SystemParams)
    c = params.c

    ω0 = params.ω0
    Wc = params.Wc

    ζ = params.ζ

    kc = let
        ωc = Δ + ω0
        ωc/c
    end
    E0c(;Δ=Δ, params=params) * exp(-(x^2+z^2)/Wc^2) * ζ * cos(kc*y)
end

function ϕt(x, y, z; params::SystemParams)
    zR = params.zR
    k0 = params.k0

    -atan(z/zR) + k0*z/2 * (x^2+y^2)/(z^2+zR^2)
end

function W(z; params::SystemParams)
    Wt = params.Wt
    zR = params.zR

    Wt * sqrt(1+(z/zR)^2)
end

expi(z) = exp(im*z)

function Et(x, y, z; params)
    E0 = params.E0
    k0 = params.k0
    Wt = params.Wt
    Ax, Ay = params.Ax, params.Ay

    (E0/2 
       * expi(k0*z+ϕt(x, y, z; params=params)) 
       * Wt/W(z; params=params) 
       * exp(-(y/Ay/W(z; params=params))^2) * exp(-(x/Ax/W(z; params=params))^2)
    )
end

end # module Laser
