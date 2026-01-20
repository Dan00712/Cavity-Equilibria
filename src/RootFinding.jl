module RootFinding

using LinearAlgebra
using JuMP, Ipopt

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


function find_roots(zguess, vd, vϕ; Δ, κ, params::SystemParams, silent::Bool=true)
    model = Model(Ipopt.Optimizer)

    if silent
        set_silent(model)
    end

    nvars = length(zguess)
    @variable(model, vz[1:nvars])
    set_start_value.(vz, zguess)

    function objective(vz...)
        vz_ = collect(vz)
        F = vL(vz_, vd, vϕ; Δ=Δ, κ=κ, params=params)
        norm(F)
    end
    register(model, :obj_func, nvars, objective, autodiff=true)
    @NLobjective(model, Min, obj_func(vz...))

    optimize!(model)

    solution = value.(vz)
    residual_norm = sqrt(objective_value(model))

    solution
end

end # module RootFinding
